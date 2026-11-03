"""상체 케이블 근육(몸통 + 어깨 + 팔꿈치)을 OpenSim GUI 에서 활성도 색으로 보기 위한 파일.

`gui_export.py`(하체)와 같은 방식: PathActuator 를 같은 경로의 Millard 근육으로 바꾼 모델 + 관절 각도와 activation(=Moco 신호 u)
열을 넣은 모션. 세 설계 공간의 경로 json 과 Moco 해를 한 모델·한 모션으로 합친다.

    python report/gui_export_upper.py --base data/carrt/subject1/subject1_load3.79.osim \\
        --paths trunk=report/newmus_paths_trunk_k4.json shoulder_r=report/newmus_paths_shoulder_r_k5.json elbow_r=report/newmus_paths_elbow_r_k2.json \\
        --sol trunk=paths_trunk_k4_lift_s1_PFWJ2 shoulder_r=paths_shoulder_k5_lift_s1_PFWJ2 elbow_r=paths_elbow_k2_lift_s1_PFWJ2 \\
        --name lift_s1_PFWJ2
→ report/moco/gui_upper_<name>_model.osim, report/moco/gui_upper_<name>.mot
GUI: File → Open Model 로 .osim, File → Load Motion 으로 .mot. 뼈 메시는 report/moco/Geometry.
"""
import argparse, json, os, sys
import numpy as np, opensim as osim
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import addb_io as io                                   # noqa: E402
import design_space as ds                              # noqa: E402
import path_design as pd                               # noqa: E402
COLORS = {"trunk": [(0.75, 0.22, 0.17), (0.9, 0.49, 0.13), (0.56, 0.27, 0.68), (0.09, 0.63, 0.52)],
          "shoulder_r": [(0.16, 0.5, 0.73), (0.15, 0.68, 0.38), (0.83, 0.33, 0.0), (0.5, 0.55, 0.55), (0.2, 0.2, 0.6), (0.6, 0.6, 0.1)],
          "elbow_r": [(0.9, 0.1, 0.5), (0.1, 0.4, 0.1)]}


def kv(items):
    return dict(s.split("=", 1) for s in items)


def build(base, paths_by_space, out_path):
    proc = osim.ModelProcessor(base); proc.append(osim.ModOpRemoveMuscles()); model = proc.process()
    for space, pj in paths_by_space.items():
        sp = pd.configure(space); pd.remove_design_actuators(model)
        paths = json.load(open(pj)); wrap = sp.get("wrap")
        sides = (("r", 1.0), ("l", -1.0)) if sp["mirror"] else (("", 1.0),)
        for side, sgn in sides:
            wrap_obj = None
            if wrap:
                wbody = model.getBodySet().get(ds.mirror_body(sp, wrap["body"], side)); c = wrap["center"]
                wrap_obj = osim.WrapSphere(); wrap_obj.setName(f"wrap_{space}_{side}"); wrap_obj.set_radius(float(wrap["radius"]))
                wrap_obj.set_translation(osim.Vec3(c[0], c[1], sgn * c[2])); wbody.addWrapObject(wrap_obj)
            for i, m in enumerate(paths["muscles"]):
                nm = f"{space.replace('_r', '')}_{m['name']}" + (f"_{side}" if side else "")
                mus = osim.Millard2012EquilibriumMuscle(); mus.setName(nm)
                mus.setMaxIsometricForce(m["force_N"]); L = m.get("L_mean", 0.3)
                mus.setOptimalFiberLength(0.6 * L); mus.setTendonSlackLength(0.4 * L); mus.setPennationAngleAtOptimalFiberLength(0.0)
                mus.set_ignore_tendon_compliance(True); mus.set_ignore_activation_dynamics(True)
                for k, (bn, pt) in enumerate(zip(m["chain"], m["points"])):
                    mus.addNewPathPoint(f"{nm}_p{k}", model.getBodySet().get(ds.mirror_body(sp, bn, side)), osim.Vec3(pt[0], pt[1], sgn * pt[2]))
                if wrap_obj is not None and pd.wrap_segment([sp["bodies"].index(b) for b in m["chain"]]) is not None:
                    mus.updGeometryPath().addPathWrap(wrap_obj)
                mus.updGeometryPath().setDefaultColor(osim.Vec3(*COLORS[space][i % len(COLORS[space])]))
                model.addForce(mus)
    model.finalizeConnections(); model.initSystem(); model.printToXML(out_path); print("모델", out_path)


def motion(sols_by_space, paths_by_space, out_path):
    mdir = os.path.join(ROOT, "report/moco"); first = next(iter(sols_by_space.values()))
    ik = io.read_mot(os.path.join(mdir, f"{first}_ik_full.mot")); sol0 = io.read_mot(os.path.join(mdir, f"{first}_solution.sto"))
    t = sol0.time; names = list(ik.names); D = [np.interp(t, ik.time, ik.column(n)) for n in names]
    for space, sname in sols_by_space.items():
        sp = ds.get(space); sol = io.read_mot(os.path.join(mdir, f"{sname}_solution.sto")); paths = json.load(open(paths_by_space[space]))
        for side in (("r", "l") if sp["mirror"] else ("",)):
            for m in paths["muscles"]:
                col = f"/forceset/{m['name']}" + (f"_{side}" if side else ""); nm = f"{space.replace('_r', '')}_{m['name']}" + (f"_{side}" if side else "")
                names.append(f"/forceset/{nm}/activation"); D.append(np.clip(np.interp(t, sol.time, sol.column(col)), 0, 1))
    io.write_mot(out_path, names, t, np.column_stack(D), in_degrees=ik.in_degrees); print("모션", out_path, f"({t[0]:.2f}~{t[-1]:.2f} s)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True); ap.add_argument("--paths", nargs="+", required=True, help="space=json …")
    ap.add_argument("--sol", nargs="+", required=True, help="space=Moco 이름 (paths_… ) …"); ap.add_argument("--name", required=True)
    a = ap.parse_args(); paths = {k: (v if os.path.isabs(v) else os.path.join(ROOT, v)) for k, v in kv(a.paths).items()}; sols = kv(a.sol)
    build(a.base if os.path.isabs(a.base) else os.path.join(ROOT, a.base), paths, os.path.join(ROOT, f"report/moco/gui_upper_{a.name}_model.osim"))
    motion(sols, paths, os.path.join(ROOT, f"report/moco/gui_upper_{a.name}.mot"))
