"""케이블 경로 근육 8개/다리 (PathActuator) 로 걸음 재현 — MocoInverse.

    python report/moco_path_muscles.py --subject subject3 --trial walking1_segment_0 [--umax 1.0] [--tag _x]
지표: (1) 보조 % (창 양 끝 0.05 s 제외), (2) 오른다리 새 근육 토크 Σ F_i u_i r_i(q) vs 우리 ID (Moco 내부 토크의 튐과 무관한 지표).
"""
import argparse, csv, json, os, sys, time
import numpy as np
import opensim as osim
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import addb_io as io                                                   # noqa: E402
from moco_synergy_inverse import strip_locked, complete_kinematics    # noqa: E402
import path_design as pd                                               # noqa: E402
import design_space as ds                                              # noqa: E402
COORDS = ("hip_flexion", "hip_adduction", "knee_angle", "ankle_angle"); EDGE = 0.05


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subject", default="subject3"); ap.add_argument("--trial", default="walking1_segment_0")
    ap.add_argument("--model", default=os.path.join(ROOT, "report/moco/newmus_paths_subject3_model.osim"))
    ap.add_argument("--umax", type=float, default=1.0); ap.add_argument("--mesh", type=float, default=0.02); ap.add_argument("--tag", default="")
    ap.add_argument("--ik", default=None, help="운동학 .mot (기본: data/<subject>/<trial>_ik.mot). 도 단위면 자동 변환")
    ap.add_argument("--ext", default=None, help="외력 설정 xml (기본: out/real 의 것)"); ap.add_argument("--id", default=None, help="기준 ID .sto (기본: out/real 의 것)")
    ap.add_argument("--t0", type=float, default=None); ap.add_argument("--t1", type=float, default=None); ap.add_argument("--name", default=None)
    ap.add_argument("--space", default="lower"); ap.add_argument("--paths", default=None, help="경로 json (기본 report/newmus_paths_k8.json)")
    ap.add_argument("--prefilter", action="store_true", help="운동학을 전체 길이에서 6 Hz 로 미리 거른다 (Moco 내부 필터는 창으로 자른 뒤 걸어 양 끝·급변 구간이 튄다)")
    a = ap.parse_args(); sub, trial = a.subject, a.trial
    global COORDS
    SP = pd.configure(a.space); COORDS = SP["base"]; dside = SP["design_side"]
    paths_json = a.paths or os.path.join(ROOT, "report/newmus_paths_k8.json")
    if a.space != "lower" and a.model == os.path.join(ROOT, "report/moco/newmus_paths_subject3_model.osim"):
        a.model = os.path.join(ROOT, f"report/moco/newmus_paths_{a.space}_k{len(json.load(open(paths_json))['muscles'])}_model.osim")
    name = a.name or f"paths_{sub}_{trial.replace('_segment_0','')}{a.tag}"; out_dir = os.path.join(ROOT, "report/moco")
    base = os.path.join(ROOT, "out/real", sub, trial); ext = None if (a.ext or "").lower() == "none" else (a.ext or os.path.join(base, f"{trial}_external_loads.xml"))
    idm = io.read_mot(a.id or os.path.join(base, f"{trial}_id.sto")); t0 = a.t0 if a.t0 is not None else float(idm.t0); t1 = a.t1 if a.t1 is not None else float(idm.t1)
    ik_raw = a.ik or os.path.join(ROOT, "data", sub, f"{trial}_ik.mot")
    if a.ik and "inDegrees=yes" in open(ik_raw).read(2000):       # 논문 IK 는 도 단위 → 라디안으로
        mot = io.read_mot(ik_raw); D = np.column_stack([mot.column(n) if n in ("pelvis_tx", "pelvis_ty", "pelvis_tz") else np.deg2rad(mot.column(n)) for n in mot.names])   # pelvis_tilt 는 회전이다 (2026-11-03 버그)
        ik_raw = os.path.join(out_dir, f"{name}_ik_rad.mot"); io.write_mot(ik_raw, mot.names, mot.time, D)
    ik_path, _ = complete_kinematics(a.model, ik_raw, os.path.join(out_dir, f"{name}_ik_full.mot"))
    pre = os.path.join(out_dir, f"{name}_pre.osim"); welds, removed = strip_locked(a.model, pre)
    proc = osim.ModelProcessor(pre); sv = osim.StdVectorString()
    for j in welds: sv.append(j)
    proc.append(osim.ModOpReplaceJointsWithWelds(sv))
    if ext: proc.append(osim.ModOpAddExternalLoads(os.path.abspath(ext)))
    proc.append(osim.ModOpAddReserves(1.0))
    model = proc.process(); model.initSystem(); mp = os.path.join(out_dir, f"{name}_model.osim"); model.printToXML(mp)
    inv = osim.MocoInverse(); inv.setName(name); inv.setModel(osim.ModelProcessor(mp))
    if a.prefilter:
        from scipy.signal import butter, filtfilt
        mot = io.read_mot(ik_path); fs = 1.0 / np.median(np.diff(mot.time)); b_, a_ = butter(4, 6.0 / (fs / 2))
        D = np.column_stack([filtfilt(b_, a_, mot.column(n)) for n in mot.names])
        ik_path = ik_path.replace(".mot", "_f6.mot"); io.write_mot(ik_path, mot.names, mot.time, D); tp = osim.TableProcessor(ik_path)
    else:
        tp = osim.TableProcessor(ik_path); tp.append(osim.TabOpLowPassFilter(6.0))
    inv.setKinematics(tp); inv.set_kinematics_allow_extra_columns(True)
    inv.set_initial_time(t0); inv.set_final_time(t1); inv.set_mesh_interval(a.mesh); inv.set_convergence_tolerance(1e-3); inv.set_constraint_tolerance(1e-3); inv.set_max_iterations(2000)
    study = inv.initialize(); prob = study.updProblem()
    prob.setControlInfoPattern("/forceset/M[0-9]+(_[rl])?", osim.MocoBounds(0.0, a.umax))
    goal = osim.MocoControlGoal.safeDownCast(prob.updGoal("excitation_effort"))
    goal.setWeightForControlPattern("/forceset/reserve_.*", 100.0); goal.setWeightForControlPattern("/forceset/reserve_jointset_ground_pelvis_.*", 0.01)
    if ext is None:                                        # 외력 없음: 설계 좌표 외의 보조는 싸게 (다리는 어차피 틀린다)
        design = {ds.full(c, s) for c in COORDS for s in SP["sides"]}
        fs = model.getForceSet()
        for i in range(fs.getSize()):
            nm = fs.get(i).getName()
            if nm.startswith("reserve_") and not any(nm.endswith("_" + d) for d in design): goal.setWeightForControlPattern(f"/forceset/{nm}", 0.01)
    tic = time.time(); msol = study.solve(); ok = msol.success(); msol.unseal(); msol.write(os.path.join(out_dir, f"{name}_solution.sto"))
    print(f"[{name}] 성공={ok}  {(time.time()-tic)/60:.1f} 분  목적함수={msol.getObjective():.1f}  반복 {msol.getNumIterations()}", flush=True)
    names = list(msol.getControlNames()); tt = np.asarray(msol.getTimeMat()).ravel(); inner = (tt >= t0 + EDGE) & (tt <= t1 - EDGE)
    print(f"[{name}] 좌표            보조 피크(Nm)  ID 피크  보조/ID    | 설계 쪽: 새 근육 토크 vs ID 최악")
    # 설계 쪽 토크 = Σ F_i u_i r_i(q)
    paths = json.load(open(paths_json)); ctrl_name = (lambda nm: f"/forceset/{nm}_{dside}" if dside else f"/forceset/{nm}")
    ik = io.read_mot(ik_path); Q = np.column_stack([np.interp(tt, ik.time, ik.column(c)) for c in list(pd.COORDS) + list(pd.EXTRA)])
    fr = pd.Frames(a.model, Q)
    tau_new = np.zeros((len(tt), len(COORDS)))
    for m in paths["muscles"]:
        u = np.asarray(msol.getControlMat(ctrl_name(m["name"]))).ravel(); chain = [pd.BODIES.index(b) for b in m["chain"]]
        r, _ = pd.moment_arms(fr, chain, np.array(m["points"])); tau_new += (m["force_N"] * u)[:, None] * r
    rows = []; worst = 0; worst_err = 0; worst_e95 = 0; worst_r95 = 0
    for j, c in enumerate(COORDS):
        for side in SP["sides"]:
            coord = ds.full(c, side); cands = [n for n in names if "reserve" in n and n.endswith("_" + coord)]
            r_t = np.abs(np.asarray(msol.getControlMat(cands[0])).ravel()[inner]); res = float(r_t.max()); idc = idm.column(f"{coord}_moment"); idpk = float(np.abs(idc).max())
            pct = res / idpk * 100; worst = max(worst, pct); r95 = float(np.percentile(r_t, 95)) / idpk * 100; worst_r95 = max(worst_r95, r95); extra = f"  (95퍼센타일 {r95:5.1f} %)"
            if side == dside:
                e_t = np.abs(tau_new[:, j] - np.interp(tt, idm.time, idc))[inner]; err = e_t.max() / idpk * 100; e95 = np.percentile(e_t, 95) / idpk * 100
                worst_err = max(worst_err, err); worst_e95 = max(worst_e95, e95); extra = f"   | 최대 {err:5.1f} %  95퍼센타일 {e95:5.1f} %"
            rows.append((coord, pct)); print(f"  {coord:16s} {res:10.2f} {idpk:8.1f} {pct:8.1f} %{extra}")
    upk = {m["name"]: float(np.asarray(msol.getControlMat(ctrl_name(m["name"]))).ravel().max()) for m in paths["muscles"]}
    print(f"  설계 쪽 신호 피크 (u_max {a.umax}): " + "  ".join(f"{k}={v:.2f}" for k, v in upk.items()))
    # 최대값은 Moco 가 창을 잘라 6 Hz 로 거른 운동학과 ID(전체 길이 필터)가 어긋나는 한두 순간에 튄다 — 95퍼센타일로 판정 (2026-11-03)
    print(f"[{name}] 보조 최악 {worst:.1f} % (95퍼센타일 {worst_r95:.1f} %)   새 근육 토크 vs ID 최대 {worst_err:.1f} % / 95퍼센타일 {worst_e95:.1f} %  → {'재현' if worst_r95 < 10 else '미달'}", flush=True)
    with open(os.path.join(out_dir, "results_paths.csv"), "a", newline="") as fh:
        csv.writer(fh).writerow([name, ok, f"{worst:.1f}", f"{worst_r95:.1f}", f"{worst_err:.1f}", f"{worst_e95:.1f}"] + [f"{r[1]:.1f}" for r in rows] + [f"{v:.2f}" for v in upk.values()])


if __name__ == "__main__":
    main()
