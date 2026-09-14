"""새 근육 k개(다리당)로 걷기 — 사람 근육 80개를 전부 떼고, 시너지당 하나씩 '토크형 새 근육'을 단다.

새 근육 i = 한 개의 신호 u_i(t) ∈ [0, u_max] 가 고관절 굴곡·내전·무릎·발목에 **고정 비율** V_i 의
토크를 내는 액추에이터. `--v-from SUBJ TRIAL` 로 V 를 다른 걸음에서(TRIAL=ALL 이면 그 피험자의
모든 걸음을 합쳐) 만들어 하드웨어 고정 교차 검증을 한다. `--reserve-weight 100` 을 권장 — 1 이면
u 와 보조(Nm)가 같은 비용이라 새 근육이 낼 수 있어도 보조에 떠넘긴다(NOTES 3절). 구현은 좌표당 CoordinateActuator(optimal_force = 1 Nm) 4개를
SynergyController 로 묶고 시너지 벡터 = V_i (부호 있음). 운동학·지면반력은 측정값 그대로.
답은 보조 액추에이터(reserve)가 떠맡는 몫과 각 새 근육의 최대 신호.

    python report/moco_new_muscles.py --k 8
    python report/moco_new_muscles.py --k 5
"""
import argparse, csv, glob, os, sys, time
import numpy as np
import opensim as osim

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import addb_io as io                                   # noqa: E402
from moco_synergy_inverse import strip_locked, complete_kinematics   # noqa: E402
import lumped_actuator_test as lat                     # noqa: E402

osim.Logger.setLevelString("Error")
COORDS = lat.COORDS


def build_model(model_path, ext_loads, k, sub, trial, out_dir, name, umax, v_from=None):
    """v_from=(sub, trial): 새 근육 방향 V 를 다른 trial 에서 만든다 (하드웨어 고정 → 교차 검증)."""
    pre = os.path.join(out_dir, f"{name}_pre.osim")
    welds, removed = strip_locked(model_path, pre)
    proc = osim.ModelProcessor(pre)
    sv = osim.StdVectorString()
    for j in welds: sv.append(j)
    proc.append(osim.ModOpReplaceJointsWithWelds(sv))
    proc.append(osim.ModOpAddExternalLoads(os.path.abspath(ext_loads)))
    proc.append(osim.ModOpRemoveMuscles())                 # 사람 근육 80개 전부 제거
    proc.append(osim.ModOpAddReserves(1.0))                # 모든 좌표에 보조 — 새 근육이 못 내는 몫을 잰다
    model = proc.process(); model.initSystem()

    V_all = {}
    for side in ("r", "l"):
        vs, vt = v_from if v_from else (sub, trial)
        if v_from and vt.endswith(".npz"):                      # 설계 파일: V_nm (8×4, 최대 신호 u=1 일 때 Nm), 좌우 같은 설계
            V = np.load(vt)["V_nm"]
        elif vt == "ALL":                                     # 그 피험자의 모든 걸음을 합쳐서 V
            trials = sorted(os.path.basename(p) for p in glob.glob(os.path.join(ROOT, "out/real", vs, "walking*")))
            A, F, R, tau = lat.load_pooled(vs, trials, side); V = lat.directions(A, F, R, k)
        else:
            A, F, R, tau = lat.load(vs, vt, side); V = lat.directions(A, F, R, k)   # (k, 4) Nm, 단위 신호당
        V_all[side] = V
        acts = []
        for c in COORDS:
            coord = f"{c}_{side}"
            a = osim.CoordinateActuator(coord); a.setName(f"newmus_{coord}")
            a.setOptimalForce(1.0); a.setMinControl(-1e4); a.setMaxControl(1e4)
            model.addForce(a); acts.append(a)
        ctrl = osim.SynergyController(); ctrl.setName(f"newmuscles_{side}")
        model.finalizeConnections()
        for a in acts: ctrl.addActuator(a)
        for i in range(k):
            w = V[i].copy(); w[np.abs(w) < 1e-8] = 0.0
            ctrl.addSynergyVector(osim.Vector(w.tolist()))
        model.addController(ctrl)
    model.finalizeConnections(); model.initSystem()
    out = os.path.join(out_dir, f"{name}_model.osim"); model.printToXML(out)
    print(f"[{name}] 근육 제거, 새 근육 {k}개×2  (weld {welds}, 뗀 팔 액추에이터 {len(removed)}개)", flush=True)
    return out, V_all


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subject", default="subject3"); ap.add_argument("--trial", default="walking1_segment_0")
    ap.add_argument("--k", type=int, default=8); ap.add_argument("--umax", type=float, default=3.0)
    ap.add_argument("--mesh", type=float, default=0.02); ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--tag", default="")
    ap.add_argument("--reserve-weight", type=float, default=1.0,
                    help="보조 액추에이터 비용 가중치. 1 이면 u 와 보조(Nm)가 같은 값이라 u 가 크면 보조에 떠넘긴다; 100 이면 정말 못 내는 몫만 보조로 남는다")
    ap.add_argument("--v-from", nargs=2, metavar=("SUBJECT", "TRIAL"), default=None,
                    help="새 근육 방향 V 를 이 (subject, trial) 로 만든다. TRIAL=ALL 이면 그 피험자 전체, .npz 면 설계 파일(V_nm). 기본은 자기 자신")
    args = ap.parse_args()
    sub, trial, k = args.subject, args.trial, args.k
    name = f"newmus_{sub}_{trial.replace('_segment_0','')}_k{k}{args.tag}"
    out_dir = os.path.join(ROOT, "report", "moco"); os.makedirs(out_dir, exist_ok=True)
    base = os.path.join(ROOT, "out/real", sub, trial)
    model_path = os.path.join(ROOT, "data", sub, f"{sub}_scaled.osim")
    ik_raw = os.path.join(ROOT, "data", sub, f"{trial}_ik.mot")
    ext = os.path.join(base, f"{trial}_external_loads.xml")
    idm = io.read_mot(os.path.join(base, f"{trial}_id.sto")); t0, t1 = float(idm.t0), float(idm.t1)
    ik_path, _ = complete_kinematics(model_path, ik_raw, os.path.join(out_dir, f"{name}_ik_full.mot"))

    model_osim, V_all = build_model(model_path, ext, k, sub, trial, out_dir, name, args.umax, v_from=args.v_from)
    if args.v_from: print(f"[{name}] V 출처: {args.v_from[0]}/{args.v_from[1]}  (교차 검증)", flush=True)
    inv = osim.MocoInverse(); inv.setName(name)
    inv.setModel(osim.ModelProcessor(model_osim))
    tp = osim.TableProcessor(ik_path); tp.append(osim.TabOpLowPassFilter(6.0)); inv.setKinematics(tp)
    inv.set_kinematics_allow_extra_columns(True)
    inv.set_initial_time(t0); inv.set_final_time(t1); inv.set_mesh_interval(args.mesh)
    inv.set_convergence_tolerance(args.tol); inv.set_constraint_tolerance(args.tol); inv.set_max_iterations(2000)
    study = inv.initialize()
    study.updProblem().setInputControlInfoPattern(".*synergy_excitation.*", osim.MocoBounds(0.0, args.umax))
    # 새 근육이 내는 토크(Nm)가 보조 액추에이터와 같은 제곱 비용으로 들어가면 둘을 반반으로 나눈다.
    # 비용은 신호 u 와 보조에만 매긴다 — 새 근육 출력의 가중치는 0.
    goal = osim.MocoControlGoal.safeDownCast(study.updProblem().updGoal("excitation_effort"))
    goal.setWeightForControlPattern("/forceset/newmus_.*", 0.0)
    goal.setWeightForControlPattern("/forceset/reserve_.*", args.reserve_weight)
    tic = time.time(); msol = study.solve(); dt = time.time() - tic
    ok = msol.success(); msol.unseal()
    print(f"[{name}] 성공={ok}  {dt/60:.1f} 분  목적함수={msol.getObjective():.3f}  반복 {msol.getNumIterations()}", flush=True)
    msol.write(os.path.join(out_dir, f"{name}_solution.sto"))

    names = list(msol.getControlNames()); rows = []
    print(f"[{name}] 좌표            reserve 피크(Nm)   ID 피크(Nm)   reserve/ID")
    for c in COORDS:
        for side in ("r", "l"):
            coord = f"{c}_{side}"; cands = [n for n in names if "reserve" in n and n.endswith("_" + coord)]
            if not cands: continue
            res = float(np.abs(msol.getControlMat(cands[0])).max()); idpk = float(np.abs(idm.column(f"{coord}_moment")).max())
            rows.append((coord, res, idpk, res / idpk * 100)); print(f"  {coord:16s} {res:12.2f} {idpk:12.1f} {res/idpk*100:10.1f} %")
    worst = max(r[3] for r in rows)
    # 새 근육 신호 크기 — 하드웨어 용량 감각
    try:
        inames = list(msol.getInputControlNames())
        U = np.column_stack([np.asarray(msol.getInputControlMat(n)).ravel() for n in inames]) if inames else None
    except Exception:
        try:
            U = msol.getInputControlsTrajectoryMat(); inames = list(msol.getInputControlNames())
        except Exception as exc:
            U = None; print("  (입력 신호 조회 불가:", exc, ")")
    if U is not None:
        umax_each = U.max(axis=0)
        print(f"  새 근육 최대 신호 (u_max 한도 {args.umax}): " + "  ".join(f"{n.split('/')[-2][-1]}{n.split('_')[-1]}={v:.2f}" for n, v in zip(inames, umax_each)))
    with open(os.path.join(out_dir, "results_newmus.csv"), "a", newline="") as fh:
        csv.writer(fh).writerow([name, k, ok, f"{dt/60:.1f}", f"{msol.getObjective():.3f}", f"{worst:.1f}"] + [f"{r[3]:.1f}" for r in rows])
    print(f"[{name}] reserve 최악 {worst:.1f} %  → {'새 근육만으로 걸음 재현' if worst < 10 else '보조 개입 큼'}", flush=True)


if __name__ == "__main__":
    main()
