"""OpenSim GUI 에서 케이블 근육을 활성도 색으로 보기 위한 파일.

GUI 는 Muscle 만 활성도로 색을 바꾸므로 PathActuator 를 같은 경로의 Millard 근육으로 바꾼 모델을 만들고,
모션 파일에 관절 각도 + 각 근육의 activation(=Moco 신호 u) 을 넣는다.
    python report/gui_export.py paths_subject3_squats1 paths_subject3_walking1
→ report/moco/gui_muscles_model.osim, report/moco/gui_<name>.mot
"""
import json, os, sys
import numpy as np, opensim as osim
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import addb_io as io                                   # noqa: E402
COLORS = [(0.75, 0.22, 0.17), (0.16, 0.5, 0.73), (0.15, 0.68, 0.38), (0.9, 0.49, 0.13), (0.56, 0.27, 0.68), (0.09, 0.63, 0.52), (0.83, 0.33, 0.0), (0.5, 0.55, 0.55)]


def build():
    paths = json.load(open(os.path.join(ROOT, "report/newmus_paths_k8.json")))
    proc = osim.ModelProcessor(os.path.join(ROOT, "data/subject3/subject3_scaled.osim")); proc.append(osim.ModOpRemoveMuscles()); model = proc.process()
    for side, sgn in (("r", 1.0), ("l", -1.0)):
        for i, m in enumerate(paths["muscles"]):
            mus = osim.Millard2012EquilibriumMuscle(); mus.setName(f"{m['name']}_{side}")
            mus.setMaxIsometricForce(m["force_N"]); L = m.get("L_mean", 0.8)
            mus.setOptimalFiberLength(0.6 * L); mus.setTendonSlackLength(0.4 * L); mus.setPennationAngleAtOptimalFiberLength(0.0)
            mus.set_ignore_tendon_compliance(True); mus.set_ignore_activation_dynamics(True)
            for k, (bn, pt) in enumerate(zip(m["chain"], m["points"])):
                body = bn if bn == "pelvis" else bn[:-1] + side
                mus.addNewPathPoint(f"{m['name']}_{side}_p{k}", model.getBodySet().get(body), osim.Vec3(pt[0], pt[1], sgn * pt[2]))
            mus.updGeometryPath().setDefaultColor(osim.Vec3(*COLORS[i]))
            model.addForce(mus)
    model.finalizeConnections(); model.initSystem()
    p = os.path.join(ROOT, "report/moco/gui_muscles_model.osim"); model.printToXML(p); print("모델", p); return paths


def motion(name, paths):
    mdir = os.path.join(ROOT, "report/moco"); sol = io.read_mot(os.path.join(mdir, f"{name}_solution.sto")); ik = io.read_mot(os.path.join(mdir, f"{name}_ik_full.mot"))
    t = sol.time; cols = [n for n in ik.names]; D = [np.interp(t, ik.time, ik.column(n)) for n in cols]
    names = list(cols)
    for side in "rl":
        for m in paths["muscles"]:
            names.append(f"/forceset/{m['name']}_{side}/activation"); D.append(np.clip(sol.column(f"/forceset/{m['name']}_{side}"), 0, 1))
    p = os.path.join(mdir, f"gui_{name}.mot"); io.write_mot(p, names, t, np.column_stack(D)); print("모션", p, f"({t[0]:.2f}~{t[-1]:.2f} s)")


if __name__ == "__main__":
    paths = build()
    for n in sys.argv[1:]: motion(n, paths)
