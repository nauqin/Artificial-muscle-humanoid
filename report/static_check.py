"""케이블 경로의 실제 모멘트암 r_i(q_t)로, ID 토크를 u∈[0,1]로 정적으로 재현할 수 있나 (Moco 없이, 프레임별 box-lsq).
현재 경로(report/newmus_paths_k8.json)와 v1(archive_v1)을 스쿼트·걷기 전 프레임에서 비교한다. 사용: python report/static_check.py"""
import sys, os, json, numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/report")
import addb_io as io, path_design as pd
from scipy.optimize import lsq_linear
MODEL = ROOT + "/data/subject3/subject3_scaled.osim"
def run(jsonfile, npzfile, trial, t0, t1, step=4):
    paths = json.load(open(jsonfile)); d = np.load(npzfile)
    ik = io.read_mot(f"{ROOT}/data/subject3/{trial}_segment_0_ik.mot"); idm = io.read_mot(f"{ROOT}/out/real/subject3/{trial}_segment_0/{trial}_segment_0_id.sto")
    sel = (ik.time >= t0) & (ik.time <= t1); tt = ik.time[sel][::step]
    Q = np.column_stack([np.interp(tt, ik.time, ik.column(c)) for c in list(pd.COORDS) + list(pd.EXTRA)])
    tau = np.column_stack([np.interp(tt, idm.time, idm.column(c + "_moment")) for c in pd.COORDS]); pk = np.abs(tau).max(axis=0)
    fr = pd.Frames(MODEL, Q)
    A = np.zeros((len(tt), 4, 8))
    for i, m in enumerate(paths["muscles"]):
        chain = [pd.BODIES.index(b) for b in m["chain"]]; r, _ = pd.moment_arms(fr, chain, np.array(m["points"])); A[:, :, i] = m["force_N"] * r
    Vd = d["V_nm"]  # (8,4) 설계 토크방향×용량
    res_c = np.zeros((len(tt), 4)); res_d = np.zeros((len(tt), 4)); umax_c = np.zeros(8)
    for t in range(len(tt)):
        u = lsq_linear(A[t], tau[t], bounds=(0, 1)).x; res_c[t] = A[t] @ u - tau[t]; umax_c = np.maximum(umax_c, u)
        u2 = lsq_linear(Vd.T, tau[t], bounds=(0, 1)).x; res_d[t] = Vd.T @ u2 - tau[t]
    print(f"{os.path.basename(os.path.dirname(jsonfile))}/{trial} {t0}-{t1}s  프레임 {len(tt)}")
    for j, c in enumerate(pd.COORDS):
        print(f"  {c:16s} ID피크 {pk[j]:6.1f}  케이블 잔차 최악 {np.abs(res_c[:, j]).max()/pk[j]*100:5.1f} %   설계방향(토크형) 잔차 {np.abs(res_d[:, j]).max()/pk[j]*100:5.1f} %")
    print("  u 최대(케이블):", np.round(umax_c, 2))
    worst = np.abs(res_c).max(axis=1) / pk.max(); k = worst.argmax(); print(f"  최악 프레임 t={tt[k]:.2f}s q(deg)={np.round(np.rad2deg(Q[k,:4]),0)} 잔차 {np.round(res_c[k],1)} Nm  ID {np.round(tau[k],1)}")
for tag, jf, nf in [("v2", ROOT + "/report/newmus_paths_k8.json", ROOT + "/report/newmus_V_design_k8.npz"), ("v1", ROOT + "/report/archive_v1/newmus_paths_k8.json", ROOT + "/report/archive_v1/newmus_V_design_k8.npz")]:
    print("=====", tag)
    run(jf, nf, "squats1", 2.6, 4.4); run(jf, nf, "walking1", 0.03, 1.27)
