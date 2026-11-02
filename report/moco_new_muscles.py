"""새 근육 k개(다리당)로 걷기 — 사람 근육 80개를 전부 떼고, 시너지당 하나씩 '토크형 새 근육'을 단다.

새 근육 i = 한 개의 신호 u_i(t) ∈ [0, u_max] 가 고관절 굴곡·내전·무릎·발목에 **고정 비율** V_i 의
토크를 내는 액추에이터. `--v-from SUBJ TRIAL` 로 V 를 다른 걸음에서(TRIAL=ALL 이면 그 피험자의
모든 걸음을 합쳐) 만들어 하드웨어 고정 교차 검증을 한다. `--reserve-weight 100` 을 권장 — 1 이면
u 와 보조(Nm)가 같은 비용이라 새 근육이 낼 수 있어도 보조에 떠넘긴다(NOTES 3절). 구현은 좌표당 CoordinateActuator(optimal_force = 1 Nm) 4개를
SynergyController 로 묶고 시너지 벡터 = V_i (부호 있음). 운동학·지면반력은 측정값 그대로.
답은 보조 액추에이터(reserve)가 떠맡는 몫과 각 새 근육의 최대 신호.

    python report/moco_new_muscles.py --k 8
    python report/moco_new_muscles.py --k 5
    python report/moco_new_muscles.py --space arm_r --v-from design report/newmus_V_design_arm_r_k6.npz --ext none \
        --model data/moore2015/subject3/subject3_scaled.osim --ik data/moore2015/subject3/T048_segment_26_ik.mot \
        --id out/moore2015/subject3/T048_segment_26/T048_segment_26_id.sto --v-scale 0.96 --reserve-weight 100 --name arm_s3_T048_26
"""
import argparse, csv, glob, os, sys, time
import numpy as np
import opensim as osim

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import addb_io as io                                   # noqa: E402
from moco_synergy_inverse import strip_locked, complete_kinematics   # noqa: E402
import lumped_actuator_test as lat                     # noqa: E402
import design_space as ds                              # noqa: E402
import path_design as pd                               # noqa: E402

osim.Logger.setLevelString("Error")
SP = ds.get("lower"); COORDS = SP["base"]
EDGE = 0.05   # s. 보조 % 를 잴 때 구간 양 끝에서 뺄 폭 (MocoInverse 는 마지막 몇 노드에서 토크가 튄다)


def build_model(model_path, ext_loads, k, sub, trial, out_dir, name, umax, v_from=None, v_scale=1.0):
    """v_from=(sub, trial): 새 근육 방향 V 를 다른 trial 에서 만든다 (하드웨어 고정 → 교차 검증)."""
    pre = os.path.join(out_dir, f"{name}_pre.osim")
    welds, removed = strip_locked(model_path, pre)
    proc = osim.ModelProcessor(pre)
    sv = osim.StdVectorString()
    for j in welds: sv.append(j)
    proc.append(osim.ModOpReplaceJointsWithWelds(sv))
    if ext_loads: proc.append(osim.ModOpAddExternalLoads(os.path.abspath(ext_loads)))   # 상체 전용이면 외력 없음
    proc.append(osim.ModOpRemoveMuscles())                 # 사람 근육 80개 전부 제거
    proc.append(osim.ModOpAddReserves(1.0))                # 모든 좌표에 보조 — 새 근육이 못 내는 몫을 잰다
    model = proc.process(); model.initSystem()
    pd.remove_design_actuators(model)                      # 설계 좌표의 원래 토크 모터(shoulder_flex_r 등) 제거

    V_all = {}
    for side in SP["sides"]:
        vs, vt = v_from if v_from else (sub, trial)
        if v_from and vt.endswith(".npz"):                      # 설계 파일: V_nm (8×4, 최대 신호 u=1 일 때 Nm), 좌우 같은 설계
            V = np.load(vt if os.path.isabs(vt) else os.path.join(ROOT, vt))["V_nm"] * v_scale   # 몸 크기 비율(질량×키 / 피험자3)로 용량을 맞춘다
        elif vt == "ALL" and SP["name"] == "lower":                                     # 그 피험자의 모든 걸음을 합쳐서 V
            trials = sorted(os.path.basename(p) for p in glob.glob(os.path.join(ROOT, "out/real", vs, "walking*")))
            A, F, R, tau = lat.load_pooled(vs, trials, side); V = lat.directions(A, F, R, k)
        else:
            A, F, R, tau = lat.load(vs, vt, side); V = lat.directions(A, F, R, k)   # (k, 4) Nm, 단위 신호당
        V_all[side] = V
        acts = []
        for c in COORDS:
            coord = ds.full(c, side)
            a = osim.CoordinateActuator(coord); a.setName(f"newmus_{coord}")
            a.setOptimalForce(1.0); a.setMinControl(-1e4); a.setMaxControl(1e4)
            model.addForce(a); acts.append(a)
        ctrl = osim.SynergyController(); ctrl.setName(f"newmuscles_{side or 'c'}")
        model.finalizeConnections()
        for a in acts: ctrl.addActuator(a)
        for i in range(k):
            w = V[i].copy(); w[np.abs(w) < 1e-8] = 0.0
            ctrl.addSynergyVector(osim.Vector(w.tolist()))
        model.addController(ctrl)
    model.finalizeConnections(); model.initSystem()
    out = os.path.join(out_dir, f"{name}_model.osim"); model.printToXML(out)
    print(f"[{name}] 근육 제거, 새 근육 {k}개×{len(SP['sides'])}  (weld {welds}, 뗀 잠긴 액추에이터 {len(removed)}개)", flush=True)
    return out, V_all


def msol_names(study):
    """문제에 들어간 제어 이름 목록 (보조 가중치를 좌표별로 매기려고)."""
    rep = study.getProblem().createRep(); m = rep.getModelBase()
    return [f"/forceset/{m.getForceSet().get(i).getName()}" for i in range(m.getForceSet().getSize())]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subject", default="subject3"); ap.add_argument("--trial", default="walking1_segment_0")
    ap.add_argument("--k", type=int, default=8); ap.add_argument("--umax", type=float, default=3.0)
    ap.add_argument("--mesh", type=float, default=0.02); ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--tag", default="")
    ap.add_argument("--reserve-weight", type=float, default=1.0,
                    help="보조 액추에이터 비용 가중치. 1 이면 u 와 보조(Nm)가 같은 값이라 u 가 크면 보조에 떠넘긴다; 100 이면 정말 못 내는 몫만 보조로 남는다")
    ap.add_argument("--prefilter", action="store_true", help="운동학을 전체 길이에서 6 Hz 로 미리 필터해 넘긴다 (Moco 내부 필터는 창으로 자른 뒤 걸어 양 끝이 튄다)")
    ap.add_argument("--v-scale", type=float, default=None, help="설계 V(Nm) 에 곱할 비율. 기본: 피험자 (질량×키)/(63.5×1.69) 자동")
    ap.add_argument("--v-from", nargs=2, metavar=("SUBJECT", "TRIAL"), default=None,
                    help="새 근육 방향 V 를 이 (subject, trial) 로 만든다. TRIAL=ALL 이면 그 피험자 전체, .npz 면 설계 파일(V_nm). 기본은 자기 자신")
    ap.add_argument("--space", default="lower", help=f"설계 공간 {list(ds.SPACES)}")
    ap.add_argument("--model", default=None); ap.add_argument("--ik", default=None); ap.add_argument("--id", default=None)
    ap.add_argument("--ext", default=None, help="외력 xml. 'none' 이면 외력 없이 (상체 전용)"); ap.add_argument("--name", default=None)
    ap.add_argument("--t0", type=float, default=None); ap.add_argument("--t1", type=float, default=None)
    args = ap.parse_args()
    global SP, COORDS
    SP = ds.get(args.space); COORDS = SP["base"]; pd.configure(args.space)
    sub, trial, k = args.subject, args.trial, args.k
    if args.v_from and args.v_from[1].endswith(".npz"): k = len(np.load(args.v_from[1] if os.path.isabs(args.v_from[1]) else os.path.join(ROOT, args.v_from[1]))["V_nm"])
    name = args.name or f"newmus_{sub}_{trial.replace('_segment_0','')}_k{k}{args.tag}"
    out_dir = os.path.join(ROOT, "report", "moco"); os.makedirs(out_dir, exist_ok=True)
    base = os.path.join(ROOT, "out/real", sub, trial)
    model_path = args.model or os.path.join(ROOT, "data", sub, f"{sub}_scaled.osim")
    ik_raw = args.ik or os.path.join(ROOT, "data", sub, f"{trial}_ik.mot")
    ext = None if (args.ext or "").lower() == "none" else (args.ext or os.path.join(base, f"{trial}_external_loads.xml"))
    idm = io.read_mot(args.id or os.path.join(base, f"{trial}_id.sto")); t0, t1 = float(idm.t0), float(idm.t1)
    if args.t0 is not None: t0 = args.t0
    if args.t1 is not None: t1 = args.t1
    ik_path, _ = complete_kinematics(model_path, ik_raw, os.path.join(out_dir, f"{name}_ik_full.mot"))

    v_scale = args.v_scale
    if v_scale is None:
        meta = os.path.join(os.path.dirname(ROOT), "LabValidation_withoutVideos", sub, "sessionMetadata.yaml")
        if os.path.exists(meta):
            kv = dict(l.strip().split(":") for l in open(meta) if l.startswith(("mass_kg", "height_m")))
            v_scale = float(kv["mass_kg"]) * float(kv["height_m"]) / (63.5 * 1.69)
        else: v_scale = 1.0
    model_osim, V_all = build_model(model_path, ext, k, sub, trial, out_dir, name, args.umax, v_from=args.v_from, v_scale=v_scale)
    if args.v_from and args.v_from[1].endswith(".npz"): print(f"[{name}] 설계 V 사용, 몸 크기 비율 {v_scale:.2f}", flush=True)
    if args.v_from: print(f"[{name}] V 출처: {args.v_from[0]}/{args.v_from[1]}  (교차 검증)", flush=True)
    inv = osim.MocoInverse(); inv.setName(name)
    inv.setModel(osim.ModelProcessor(model_osim))
    if args.prefilter:
        from scipy.signal import butter, filtfilt
        mot = io.read_mot(ik_path); fs = 1.0 / np.median(np.diff(mot.time)); b, a_ = butter(4, 6.0 / (fs / 2))
        D = np.column_stack([filtfilt(b, a_, mot.column(n)) for n in mot.names])
        ik_path = ik_path.replace(".mot", "_f6.mot"); io.write_mot(ik_path, mot.names, mot.time, D)
        tp = osim.TableProcessor(ik_path)
    else:
        tp = osim.TableProcessor(ik_path); tp.append(osim.TabOpLowPassFilter(6.0))
    inv.setKinematics(tp)
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
    # 골반 6개 보조는 ID 의 잔차와 같은 역할 — 힘판 접촉 시작 순간의 가짜 모멘트를 여기서 흡수해야지 다리로 밀리면 안 된다.
    goal.setWeightForControlPattern("/forceset/reserve_jointset_ground_pelvis_.*", 0.01)
    if ext is None:                                        # 외력 없음: 다리·골반·비설계 좌표의 보조는 싸게 — 설계 좌표만 잰다
        design = {ds.full(c, s) for c in COORDS for s in SP["sides"]}
        for i in range(study.getProblem().createRep().getModelBase().getForceSet().getSize()):
            pass
        for cn in list(msol_names(study)):
            if "reserve" in cn and cn.split("_")[-1] not in {d.split("_")[-1] for d in design} and not any(cn.endswith("_" + d) for d in design):
                goal.setWeightForControlPattern(cn, 0.01)
    tic = time.time(); msol = study.solve(); dt = time.time() - tic
    ok = msol.success(); msol.unseal()
    print(f"[{name}] 성공={ok}  {dt/60:.1f} 분  목적함수={msol.getObjective():.3f}  반복 {msol.getNumIterations()}", flush=True)
    msol.write(os.path.join(out_dir, f"{name}_solution.sto"))

    names = list(msol.getControlNames()); rows = []
    print(f"[{name}] 좌표            reserve 피크(Nm)   ID 피크(Nm)   reserve/ID")
    for c in COORDS:
        for side in SP["sides"]:
            coord = ds.full(c, side); cands = [n for n in names if "reserve" in n and n.endswith("_" + coord)]
            if not cands: continue
            tt = np.asarray(msol.getTimeMat()).ravel(); inner = (tt >= t0 + EDGE) & (tt <= t1 - EDGE)   # 구간 양 끝은 뺀다 — 최적화 경계에서 토크가 튄다
            res = float(np.abs(np.asarray(msol.getControlMat(cands[0])).ravel()[inner]).max()); idpk = float(np.abs(idm.column(f"{coord}_moment")).max())
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
