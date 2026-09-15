"""사용: python report/artifact_metric.py <moco 이름> <paths json> <trial> [t0 t1]
케이블 Moco 해에서 '새 근육 토크 vs ID' 를 Moco 튐 노드(새 근육+보조 = Moco 내부 토크가 ID 와 5 Nm 넘게 다른 노드) 제외로 다시 계산."""
import sys, os, json, numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/report")
import addb_io as io, path_design as pd, opensim as osim
name, jsonfile, trial = sys.argv[1], sys.argv[2], sys.argv[3]
t_lo, t_hi = (float(sys.argv[4]), float(sys.argv[5])) if len(sys.argv) > 5 else (None, None)
sol = osim.MocoTrajectory(f"{ROOT}/report/moco/{name}_solution.sto"); tt = np.asarray(sol.getTimeMat()).ravel(); names = list(sol.getControlNames())
idm = io.read_mot(f"{ROOT}/out/real/subject3/{trial}_segment_0/{trial}_segment_0_id.sto")
lo = max(tt[0] + 0.05, idm.time[0], t_lo if t_lo else -1); hi = min(tt[-1] - 0.05, idm.time[-1], t_hi if t_hi else 1e9); inner = (tt >= lo) & (tt <= hi)
ik = io.read_mot(f"{ROOT}/report/moco/{name}_ik_full.mot"); Q = np.column_stack([np.interp(tt, ik.time, ik.column(c)) for c in list(pd.COORDS) + list(pd.EXTRA)])
fr = pd.Frames(f"{ROOT}/data/subject3/subject3_scaled.osim", Q)
paths = json.load(open(jsonfile)); tau = np.zeros((len(tt), 4))
for m in paths["muscles"]:
    u = np.asarray(sol.getControlMat(f"/forceset/{m['name']}_r")).ravel(); chain = [pd.BODIES.index(b) for b in m["chain"]]
    r, _ = pd.moment_arms(fr, chain, np.array(m["points"])); tau += (m["force_N"] * u)[:, None] * r
print(f"[{name}] 창 {lo:.2f}~{hi:.2f} s, 노드 {inner.sum()}")
worst_all = worst_ex = 0; n_art = 0
for j, c in enumerate(pd.COORDS):
    res = np.asarray(sol.getControlMat([n for n in names if "reserve" in n and n.endswith("_" + c)][0])).ravel()
    idc = np.interp(tt, idm.time, idm.column(c + "_moment")); pk = np.abs(idm.column(c + "_moment")).max()
    art = np.abs(tau[:, j] + res - idc) > 5.0                    # Moco 내부 토크 ≠ ID
    err = np.abs(tau[:, j] - idc) / pk * 100
    e_all = err[inner].max(); e_ex = err[inner & ~art].max(); n_art = max(n_art, (art & inner).sum())
    worst_all = max(worst_all, e_all); worst_ex = max(worst_ex, e_ex)
    print(f"  {c:16s} ID피크 {pk:6.1f}  오차 최악 {e_all:6.1f} %  튐 노드 {(art & inner).sum():3d} 개 제외하면 {e_ex:5.1f} %")
print(f"  → 전체 {worst_all:.1f} %, 튐 노드 제외 {worst_ex:.1f} %")
