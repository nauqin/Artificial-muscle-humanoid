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
COORDS = ("hip_flexion", "hip_adduction", "knee_angle", "ankle_angle"); EDGE = 0.05


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subject", default="subject3"); ap.add_argument("--trial", default="walking1_segment_0")
    ap.add_argument("--model", default=os.path.join(ROOT, "report/moco/newmus_paths_subject3_model.osim"))
    ap.add_argument("--umax", type=float, default=1.0); ap.add_argument("--mesh", type=float, default=0.02); ap.add_argument("--tag", default="")
    a = ap.parse_args(); sub, trial = a.subject, a.trial
    name = f"paths_{sub}_{trial.replace('_segment_0','')}{a.tag}"; out_dir = os.path.join(ROOT, "report/moco")
    base = os.path.join(ROOT, "out/real", sub, trial); ext = os.path.join(base, f"{trial}_external_loads.xml")
    idm = io.read_mot(os.path.join(base, f"{trial}_id.sto")); t0, t1 = float(idm.t0), float(idm.t1)
    ik_path, _ = complete_kinematics(a.model, os.path.join(ROOT, "data", sub, f"{trial}_ik.mot"), os.path.join(out_dir, f"{name}_ik_full.mot"))
    pre = os.path.join(out_dir, f"{name}_pre.osim"); welds, removed = strip_locked(a.model, pre)
    proc = osim.ModelProcessor(pre); sv = osim.StdVectorString()
    for j in welds: sv.append(j)
    proc.append(osim.ModOpReplaceJointsWithWelds(sv)); proc.append(osim.ModOpAddExternalLoads(os.path.abspath(ext))); proc.append(osim.ModOpAddReserves(1.0))
    model = proc.process(); model.initSystem(); mp = os.path.join(out_dir, f"{name}_model.osim"); model.printToXML(mp)
    inv = osim.MocoInverse(); inv.setName(name); inv.setModel(osim.ModelProcessor(mp))
    tp = osim.TableProcessor(ik_path); tp.append(osim.TabOpLowPassFilter(6.0)); inv.setKinematics(tp); inv.set_kinematics_allow_extra_columns(True)
    inv.set_initial_time(t0); inv.set_final_time(t1); inv.set_mesh_interval(a.mesh); inv.set_convergence_tolerance(1e-3); inv.set_constraint_tolerance(1e-3); inv.set_max_iterations(2000)
    study = inv.initialize(); prob = study.updProblem()
    prob.setControlInfoPattern("/forceset/M[0-9]_.*", osim.MocoBounds(0.0, a.umax))
    goal = osim.MocoControlGoal.safeDownCast(prob.updGoal("excitation_effort"))
    goal.setWeightForControlPattern("/forceset/reserve_.*", 100.0); goal.setWeightForControlPattern("/forceset/reserve_jointset_ground_pelvis_.*", 0.01)
    tic = time.time(); msol = study.solve(); ok = msol.success(); msol.unseal(); msol.write(os.path.join(out_dir, f"{name}_solution.sto"))
    print(f"[{name}] 성공={ok}  {(time.time()-tic)/60:.1f} 분  목적함수={msol.getObjective():.1f}  반복 {msol.getNumIterations()}", flush=True)
    names = list(msol.getControlNames()); tt = np.asarray(msol.getTimeMat()).ravel(); inner = (tt >= t0 + EDGE) & (tt <= t1 - EDGE)
    print(f"[{name}] 좌표            보조 피크(Nm)  ID 피크  보조/ID    | 오른다리: 새 근육 토크 vs ID 최악")
    # 오른다리 토크 = Σ F_i u_i r_i(q)
    paths = json.load(open(os.path.join(ROOT, "report/newmus_paths_k8.json")))
    ik = io.read_mot(ik_path); Q = np.column_stack([np.interp(tt, ik.time, ik.column(c)) for c in list(pd.COORDS) + list(pd.EXTRA)])
    fr = pd.Frames(a.model, Q)
    tau_new = np.zeros((len(tt), 4))
    for m in paths["muscles"]:
        u = np.asarray(msol.getControlMat(f"/forceset/{m['name']}_r")).ravel(); chain = [pd.BODIES.index(b) for b in m["chain"]]
        r, _ = pd.moment_arms(fr, chain, np.array(m["points"])); tau_new += (m["force_N"] * u)[:, None] * r
    rows = []; worst = 0; worst_err = 0
    for j, c in enumerate(COORDS):
        for side in "rl":
            coord = f"{c}_{side}"; cands = [n for n in names if "reserve" in n and n.endswith("_" + coord)]
            res = float(np.abs(np.asarray(msol.getControlMat(cands[0])).ravel()[inner]).max()); idc = idm.column(f"{coord}_moment"); idpk = float(np.abs(idc).max())
            pct = res / idpk * 100; worst = max(worst, pct); extra = ""
            if side == "r":
                err = np.abs(tau_new[:, j] - np.interp(tt, idm.time, idc))[inner].max() / idpk * 100; worst_err = max(worst_err, err); extra = f"   | {err:5.1f} %"
            rows.append((coord, pct)); print(f"  {coord:16s} {res:10.2f} {idpk:8.1f} {pct:8.1f} %{extra}")
    upk = {m["name"]: float(np.asarray(msol.getControlMat(f"/forceset/{m['name']}_r")).ravel().max()) for m in paths["muscles"]}
    print(f"  오른다리 신호 피크 (u_max {a.umax}): " + "  ".join(f"{k}={v:.2f}" for k, v in upk.items()))
    print(f"[{name}] 보조 최악 {worst:.1f} %   새 근육 토크 vs ID 최악 {worst_err:.1f} %  → {'재현' if worst_err < 10 else '미달'}", flush=True)
    with open(os.path.join(out_dir, "results_paths.csv"), "a", newline="") as fh:
        csv.writer(fh).writerow([name, ok, f"{worst:.1f}", f"{worst_err:.1f}"] + [f"{r[1]:.1f}" for r in rows] + [f"{v:.2f}" for v in upk.values()])


if __name__ == "__main__":
    main()
