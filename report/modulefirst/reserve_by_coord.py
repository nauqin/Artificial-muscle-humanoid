"""Moco 해에서 좌표별 보조 95퍼센타일(% of ID 피크)을 다시 잰다 — moco_path_muscles 와 같은 창(양 끝 0.05 s 제외)·같은 정규화.
results_paths.csv 의 넷째 열은 좌표 중 최악값만 남기므로, 어느 좌표가 정했는지 보려고 (2026-11-05).
    python report/modulefirst/reserve_by_coord.py <space> <walk|squat|lift> <Moco 이름> [...]
"""
import os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import opensim as osim, addb_io as io, design_space as ds   # noqa: E402
IDS = {("lower", "walk"): ("out/lower40_h160/moore/T047_segment_12_id.sto", 1, 4), ("lower", "squat"): ("out/lower40_h160/subject3/squats1_segment_0_id.sto", 2.6, 4.4),
       ("upper", "walk"): ("out/upper40_h160/moore2015/subject3/T047_segment_12/T047_segment_12_id.sto", 1, 4), ("upper", "lift"): ("out/upper40_h160/carrt/subject1/PFWJ2/PFWJ2_id.sto", 0.3, 5.5)}
space, task = sys.argv[1], sys.argv[2]; SP = ds.get(space)
idp, t0, t1 = IDS[("lower" if space == "lower" else "upper", task)]; idm = io.read_mot(os.path.join(ROOT, idp))
coords = [ds.full(c, s) for c in SP["base"] for s in SP["sides"]]
print("이름".ljust(28) + "".join(c.rjust(17) for c in coords))
for name in sys.argv[3:]:
    tr = osim.MocoTrajectory(os.path.join(ROOT, "report/moco", f"{name}_solution.sto")); tt = np.asarray(tr.getTimeMat()).ravel(); inner = (tt >= t0 + 0.05) & (tt <= t1 - 0.05)
    names = list(tr.getControlNames()); row = []
    for c in coords:
        rn = [n for n in names if "reserve" in n and n.endswith("_" + c)][0]; r = np.abs(np.asarray(tr.getControlMat(rn)).ravel()[inner])
        pk = float(np.abs(idm.column(f"{c}_moment")).max()); row.append(f"{np.percentile(r, 95) / pk * 100:6.1f} % ({pk:4.1f})")
    print(name.ljust(28) + "".join(x.rjust(17) for x in row))
