"""시너지 신호 k개로 근육 80개를 구동해 측정된 걸음을 재현할 수 있나 — MocoInverse + SynergyController.

운동학(IK)과 지면반력은 측정값 그대로 주고, 근육 흥분을 **k개 시너지 신호의 조합으로만**
허용한다(OpenSim 4.6 내장 SynergyController, Moco가 자동 감지). 근육 동역학
(활성화 지연, 힘-길이-속도)이 들어간 진짜 시뮬레이션이다. 답은 보조 액추에이터(reserve)가
얼마나 떠맡느냐 — 0에 가까우면 k개 신호로 그 걸음을 낼 수 있다는 뜻.

    python report/moco_synergy_inverse.py --k 0   # 기준선: 근육 80개 자유 제어
    python report/moco_synergy_inverse.py --k 5   # 다리당 시너지 5개
    python report/moco_synergy_inverse.py --k 8
"""
import argparse, csv, os, sys, time
import numpy as np
import opensim as osim

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import synergy as sy, addb_io as io  # noqa: E402

osim.Logger.setLevelString("Error")
COLCAP = True
MIN_PEAK = sy.MIN_PEAK
ACT_FLOOR0 = False
REPORT_COORDS = ("hip_flexion", "hip_adduction", "knee_angle", "ankle_angle")


def synergy_matrix(sub, trial, side, k, model_muscles, colcap=True):
    """물리 단위 시너지 행렬. 흥분 e_j(t) = Σ_i u_i(t)·W_ij, u_i ∈ [0,1] 이면 SO 활성도를 재현한다."""
    act = os.path.join(ROOT, "out/real", sub, trial, "so", f"{trial}_StaticOptimization_activation.sto")
    mot = io.read_mot(act)
    mus = [n for n in model_muscles if n.endswith(f"_{side}")]
    A_full = np.column_stack([np.clip(mot.column(n), 0, None) for n in mus])
    peaks = A_full.max(axis=0)
    keep = peaks >= MIN_PEAK
    A = A_full[:, keep] / peaks[keep]
    C, W = sy.nmf(A, k, seed=0); C, W, _ = sy.sort_synergies(C, W)
    cmax = np.maximum(C.max(axis=0), 1e-9)
    Wphys = np.zeros((k, len(mus)))
    Wphys[:, keep] = W * peaks[keep][None, :] * cmax[:, None]
    Wphys[np.abs(Wphys) < 1e-8] = 0.0            # 4.9e-324 같은 비정규 값은 XML 왕복에서 속성을 통째로 비운다
    if colcap:                                   # u ∈ [0,1] 일 때 흥분 Σ_i u_i W_ij ≤ 1 이 되도록 열 합을 1로 제한
        colsum = Wphys.sum(axis=0)
        scale = np.where(colsum > 1.0, 1.0 / np.maximum(colsum, 1e-12), 1.0)
        Wphys = Wphys * scale[None, :]
    return mus, Wphys, C / cmax, mot.time


def strip_locked(model_path, out_path):
    """잠긴 좌표를 Moco 가 못 다루므로 (1) 그 좌표를 쓰는 CoordinateActuator 를 떼고
    (2) 그 좌표가 속한 관절 이름을 돌려준다 — 이어서 WeldJoint 로 바꾼다."""
    m = osim.Model(os.path.abspath(model_path)); m.initSystem()
    cs = m.getCoordinateSet(); locked, joints = set(), set()
    for i in range(cs.getSize()):
        c = cs.get(i)
        if c.get_locked():
            locked.add(c.getName()); joints.add(c.getJoint().getName())
    fs = m.updForceSet(); removed = []
    for i in reversed(range(fs.getSize())):
        ca = osim.CoordinateActuator.safeDownCast(fs.get(i))
        if ca is not None and ca.get_coordinate() in locked:
            removed.append(ca.getName()); fs.remove(i)
    m.finalizeConnections(); m.printToXML(out_path)
    return sorted(joints), removed


def complete_kinematics(model_path, ik_path, out_path):
    """IK 에 없는 종속 좌표(예: knee_angle_r_beta)를 커플러 구속 함수로 계산해 붙인다.
    nimblephysics 가 내보낸 IK 는 독립 DOF 37개뿐이라 모델 좌표 39개보다 적다."""
    m = osim.Model(os.path.abspath(model_path)); m.initSystem()
    ik = io.read_mot(ik_path)
    names, cols = list(ik.names), [ik.data[:, j] for j in range(len(ik.names))]
    added = []
    cs = m.getConstraintSet()
    for i in range(cs.getSize()):
        cc = osim.CoordinateCouplerConstraint.safeDownCast(cs.get(i))
        if cc is None: continue
        dep = cc.getDependentCoordinateName()
        if dep in names: continue
        ind = cc.getIndependentCoordinateNames(); ind = [ind.get(j) for j in range(ind.getSize())]
        f = cc.getFunction()
        vals = np.array([f.calcValue(osim.Vector([float(ik.column(n)[t]) for n in ind])) for t in range(len(ik.time))])
        names.append(dep); cols.append(vals); added.append(dep)
    io.write_mot(out_path, names, ik.time, np.column_stack(cols), name="ik_full", in_degrees=ik.in_degrees)
    return out_path, added


def build_model(model_path, ext_loads, k, sub, trial, out_dir, name):
    pre = os.path.join(out_dir, f"{name}_pre.osim")
    welds, removed = strip_locked(model_path, pre)
    print(f"[{name}] 잠긴 좌표의 관절 → weld: {welds}   뗀 액추에이터: {removed}", flush=True)
    proc = osim.ModelProcessor(pre)
    if welds:
        sv = osim.StdVectorString()
        for j in welds: sv.append(j)
        proc.append(osim.ModOpReplaceJointsWithWelds(sv))
    proc.append(osim.ModOpAddExternalLoads(os.path.abspath(ext_loads)))
    proc.append(osim.ModOpIgnoreTendonCompliance())
    proc.append(osim.ModOpReplaceMusclesWithDeGrooteFregly2016())
    proc.append(osim.ModOpIgnorePassiveFiberForcesDGF())
    proc.append(osim.ModOpScaleActiveFiberForceCurveWidthDGF(1.5))
    proc.append(osim.ModOpAddReserves(1.0))
    model = proc.process()
    model.initSystem()
    if ACT_FLOOR0:
        # Muscle.min_control 기본값 0.01 이 흥분·활성도 하한으로 걸린다. 시너지 가중치가 0인
        # 근육은 흥분이 0 이라 이 하한을 영원히 못 지켜 문제 전체가 '불가능'이 된다 → 0 으로.
        for i in range(model.getMuscles().getSize()):
            model.getMuscles().get(i).setMinControl(0.0)
        model.finalizeConnections(); model.initSystem()
    muscles = [model.getMuscles().get(i).getName() for i in range(model.getMuscles().getSize())]
    W_used = {}
    if k > 0:
        for side in ("r", "l"):
            mus, Wphys, C, t = synergy_matrix(sub, trial, side, k, muscles, colcap=COLCAP)
            ctrl = osim.SynergyController(); ctrl.setName(f"synergy_{side}")
            for m in mus:
                ctrl.addActuator(model.getMuscles().get(m))
            for i in range(k):
                # addSynergyVector 는 SimTK::Vector 를 직접 받는다 (SynergyVector 객체가 아니다)
                ctrl.addSynergyVector(osim.Vector(Wphys[i].tolist()))
            model.addController(ctrl)
            W_used[side] = (mus, Wphys, C, t)
    model.finalizeConnections()
    model.initSystem()
    out_osim = os.path.join(out_dir, f"{name}_model.osim")
    model.printToXML(out_osim)
    return out_osim, W_used


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subject", default="subject3"); ap.add_argument("--trial", default="walking1_segment_0")
    ap.add_argument("--k", type=int, default=5, help="0 = 시너지 없이 근육 80개 자유 제어 (기준선)")
    ap.add_argument("--mesh", type=float, default=0.02); ap.add_argument("--tol", type=float, default=1e-3)
    ap.add_argument("--max-iter", type=int, default=2000)
    ap.add_argument("--umax", type=float, default=1.0, help="시너지 입력 신호 상한")
    ap.add_argument("--no-colcap", action="store_true", help="열 합 ≤ 1 스케일링 끄기")
    ap.add_argument("--tag", default="")
    ap.add_argument("--min-peak", type=float, default=sy.MIN_PEAK, help="0 이면 모든 근육 포함")
    ap.add_argument("--act-floor0", action="store_true", help="DGF 최소 활성도를 0 으로")
    args = ap.parse_args()
    sub, trial, k = args.subject, args.trial, args.k
    name = f"{sub}_{trial.replace('_segment_0','')}_k{k}{args.tag}"
    global COLCAP, MIN_PEAK, ACT_FLOOR0; COLCAP = not args.no_colcap; MIN_PEAK = args.min_peak; ACT_FLOOR0 = args.act_floor0
    out_dir = os.path.join(ROOT, "report", "moco"); os.makedirs(out_dir, exist_ok=True)
    base = os.path.join(ROOT, "out/real", sub, trial)
    model_path = os.path.join(ROOT, "data", sub, f"{sub}_scaled.osim")
    ik_path = os.path.join(ROOT, "data", sub, f"{trial}_ik.mot")
    ext = os.path.join(base, f"{trial}_external_loads.xml")
    idm = io.read_mot(os.path.join(base, f"{trial}_id.sto")); t0, t1 = float(idm.t0), float(idm.t1)
    print(f"[{name}] 구간 {t0:.3f}~{t1:.3f} s, mesh {args.mesh} s", flush=True)

    ik_full, added = complete_kinematics(model_path, ik_path, os.path.join(out_dir, f"{name}_ik_full.mot"))
    print(f"[{name}] IK 에 종속 좌표 보강: {added}", flush=True)
    ik_path = ik_full
    model_osim, W_used = build_model(model_path, ext, k, sub, trial, out_dir, name)
    # 모델은 파일 경로로, 테이블 연산은 append 로 — 파이썬 임시객체를 C++ 에 넘기면 segfault 가 났다.
    inv = osim.MocoInverse(); inv.setName(name)
    inv.setModel(osim.ModelProcessor(model_osim))
    tp = osim.TableProcessor(os.path.abspath(ik_path)); tp.append(osim.TabOpLowPassFilter(6.0))
    inv.setKinematics(tp)
    inv.set_kinematics_allow_extra_columns(True)
    inv.set_initial_time(t0); inv.set_final_time(t1); inv.set_mesh_interval(args.mesh)
    inv.set_convergence_tolerance(args.tol); inv.set_constraint_tolerance(args.tol)
    inv.set_max_iterations(args.max_iter)

    study = inv.initialize()
    if k > 0:                                     # 시너지 입력 신호 범위 [0,1] 명시
        try:
            study.updProblem().setInputControlInfoPattern(".*synergy_excitation.*", osim.MocoBounds(0.0, args.umax))
        except Exception as exc:
            print(f"[{name}] 입력 범위 설정 불가 ({exc}) — 기본값 사용", flush=True)
    tic = time.time(); msol = study.solve(); dt = time.time() - tic
    ok = msol.success(); msol.unseal()       # 실패(반복 한도 등)해도 내용은 본다 — unseal 이 먼저
    print(f"[{name}] 성공={ok}  {dt/60:.1f} 분  목적함수={msol.getObjective():.4f}  반복 {msol.getNumIterations()}", flush=True)
    msol.write(os.path.join(out_dir, f"{name}_solution.sto"))

    # 보조 액추에이터가 떠맡은 몫 — 좌표별 피크 |reserve| / 피크 |ID 모멘트|
    names = list(msol.getControlNames()); T = msol.getTimeMat()
    rows = []
    for c in REPORT_COORDS:
        for side in ("r", "l"):
            coord = f"{c}_{side}"
            # ModOpAddReserves 는 reserve_<좌표 경로를 _ 로 이은 이름> 으로 짓는다 → 접미로 찾는다
            cands = [n for n in names if "reserve" in n and n.endswith("_" + coord)]
            if not cands: continue
            res = np.abs(msol.getControlMat(cands[0])).max()                # optimal_force = 1 → Nm
            idpk = np.abs(idm.column(f"{coord}_moment")).max()
            rows.append((coord, res, idpk, res / idpk * 100))
    print(f"\n[{name}] 좌표            reserve 피크(Nm)   ID 피크(Nm)   reserve/ID")
    for coord, res, idpk, pct in rows:
        print(f"  {coord:16s} {res:12.2f} {idpk:12.1f} {pct:10.1f} %")
    worst = max(r[3] for r in rows) if rows else float("nan")
    # 시너지 입력 신호
    try:
        inames = list(msol.getInputControlNames())
        if inames:
            print(f"  시너지 입력 {len(inames)}개: " + ", ".join(n.split('/')[-2] + '/' + n.split('/')[-1] for n in inames[:4]) + " …")
            U = np.column_stack([msol.getInputControlMat(n) for n in inames]); print(f"  입력 범위 {U.min():.2f}~{U.max():.2f}")
    except Exception as exc:
        print("  (입력 신호 조회 불가:", exc, ")")
    with open(os.path.join(out_dir, "results.csv"), "a", newline="") as fh:
        w = csv.writer(fh); w.writerow([name, k, ok, f"{dt/60:.1f}", f"{msol.getObjective():.4f}", f"{worst:.1f}"]
                                       + [f"{r[3]:.1f}" for r in rows])
    print(f"[{name}] reserve 최악 {worst:.1f} %  → {'근육만으로 재현' if worst < 10 else '보조 개입 큼'}", flush=True)


if __name__ == "__main__":
    main()
