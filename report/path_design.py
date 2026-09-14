"""설계 방향 8개 → 피험자3 뼈대 위의 케이블 경로 (부착점 최적화).

케이블 i 의 모멘트암 벡터 r(q) = −∂L/∂q (tendon excursion, 중앙차분) 가 걷기·스쿼트·앉았다서기의
모든 자세에서 설계 방향 v_i 와 나란하도록 부착점(시작·경유·끝)을 고른다. 길이 계산은 뼈 변환을
미리 저장해 두고 numpy 로만 하므로(마이크로초) 전역 방법(differential evolution)을 쓸 수 있다.

    python report/path_design.py            → report/newmus_paths_k8.json, report/moco/newmus_paths_subject3_model.osim,
                                              report/fig_path_design.png
"""
import glob, json, os, sys, time
import numpy as np
from scipy.optimize import differential_evolution, minimize
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import addb_io as io                                   # noqa: E402
LAB = os.path.join(os.path.dirname(ROOT), "LabValidation_withoutVideos")
COORDS = ("hip_flexion_r", "hip_adduction_r", "knee_angle_r", "ankle_angle_r")
EXTRA = ("hip_rotation_r",)                           # 자세에는 반영하되 설계 방향에는 없음
BODIES = ("pelvis", "femur_r", "tibia_r", "calcn_r")
DELTA = 0.01                                          # rad, 중앙차분
STROKE = 0.40                                         # 모듈 스트로크 (변형률)
F_MODULE = 5000.0                                     # N, 초대형 모듈 상한
MARGIN = 1.1                                          # 용량 여유 10 %
#: 부착 허용 상자 (바디 좌표계, m). 오른다리 기준 (z+ 가 바깥쪽).
BOX = {"pelvis": ((-0.12, 0.08), (-0.12, 0.05), (0.02, 0.16)),
       "femur_r": ((-0.07, 0.07), (-0.42, -0.02), (-0.06, 0.07)),
       "tibia_r": ((-0.07, 0.06), (-0.40, -0.02), (-0.05, 0.05)),
       "calcn_r": ((-0.03, 0.18), (-0.02, 0.06), (-0.05, 0.05))}
for f in font_manager.findSystemFonts():
    if "AppleSDGothicNeo" in f: plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
plt.rcParams["axes.unicode_minus"] = False


def poses(n_walk=10, n_other=30, seed=0):
    """자세 표본 (N × 5): 4 설계 좌표 + hip_rotation, rad."""
    rng = np.random.default_rng(seed); rows = []
    cols = list(COORDS) + list(EXTRA)
    for f in sorted(glob.glob(os.path.join(ROOT, "data/subject3/walking*_ik.mot"))):
        m = io.read_mot(f); idx = rng.choice(len(m.time), n_walk, replace=False)
        rows.append(np.column_stack([m.column(c) for c in cols])[idx])
    for task in ("squats1", "STS1"):
        m = io.read_mot(os.path.join(LAB, "subject3/OpenSimData/Mocap/IK", f"{task}.mot")); idx = rng.choice(len(m.time), n_other, replace=False)
        rows.append(np.deg2rad(np.column_stack([m.column(c) for c in cols])[idx]))
    return np.vstack(rows)


class Frames:
    """자세 × (기준, 좌표 j ±δ) × 바디 의 지면 변환을 미리 계산해 둔다."""
    def __init__(self, model_path, Q):
        import opensim as osim
        osim.Logger.setLevelString("Error")
        proc = osim.ModelProcessor(model_path); proc.append(osim.ModOpRemoveMuscles()); model = proc.process(); state = model.initSystem()
        cs = model.getCoordinateSet(); names = list(COORDS) + list(EXTRA)
        nP = 1 + 2 * len(COORDS)
        self.R = np.zeros((len(Q), nP, len(BODIES), 3, 3)); self.p = np.zeros((len(Q), nP, len(BODIES), 3))
        e = [osim.Vec3(1, 0, 0), osim.Vec3(0, 1, 0), osim.Vec3(0, 0, 1)]; o = osim.Vec3(0, 0, 0)
        for t, q in enumerate(Q):
            for k in range(nP):
                qq = q.copy()
                if k > 0: j = (k - 1) // 2; qq[j] += DELTA if k % 2 == 1 else -DELTA
                for name, v in zip(names, qq): cs.get(name).setValue(state, float(v), False)
                model.realizePosition(state)
                for b, bn in enumerate(BODIES):
                    body = model.getBodySet().get(bn)
                    p0 = np.array([body.findStationLocationInGround(state, o).get(i) for i in range(3)])
                    self.p[t, k, b] = p0
                    for a in range(3):
                        self.R[t, k, b, :, a] = np.array([body.findStationLocationInGround(state, e[a]).get(i) for i in range(3)]) - p0
        self.model = model


def path_lengths(fr, chain, P):
    """chain: 바디 인덱스 목록, P: (len(chain), 3) 바디 좌표. 반환 L (poses, nP)."""
    pts = np.stack([np.einsum("tkij,j->tki", fr.R[:, :, b], P[n]) + fr.p[:, :, b] for n, b in enumerate(chain)], axis=2)  # (T, nP, n, 3)
    return np.linalg.norm(np.diff(pts, axis=2), axis=3).sum(axis=2)


def moment_arms(fr, chain, P):
    L = path_lengths(fr, chain, P)
    r = np.stack([-(L[:, 1 + 2 * j] - L[:, 2 + 2 * j]) / (2 * DELTA) for j in range(len(COORDS))], axis=1)  # (T, 4)
    return r, L[:, 0]


def chain_for(v, thr=0.12):
    hip, knee, ankle = max(abs(v[0]), abs(v[1])) > thr, abs(v[2]) > thr, abs(v[3]) > thr
    first = 0 if hip else (1 if knee else 2); last = 3 if ankle else (2 if knee else 1)
    return list(range(first, last + 1))


def cost(x, fr, chain, v, cap, detail=False):
    P = x.reshape(len(chain), 3); r, L = moment_arms(fr, chain, P)
    nr = np.linalg.norm(r, axis=1); ok = nr > 1e-4
    if ok.sum() < len(r) * 0.5: return 1e3 if not detail else None
    dir_err = float((np.linalg.norm(r[ok] / nr[ok, None] - v, axis=1) ** 2).mean())
    arm_eff = (r @ v)                                          # 설계 방향으로의 유효 모멘트암 (m)
    arm_min = float(arm_eff.min()); arm_mean = float(arm_eff.mean())
    stroke = float((L.max() - L.min()) / L.mean())
    force = MARGIN * cap / max(arm_min, 1e-4)                  # 설계 용량을 가장 불리한 자세에서 내는 데 필요한 힘
    pen = 20 * max(0.0, stroke - STROKE) ** 2 + 5 * max(0.0, force / F_MODULE - 1) ** 2 + 1e4 * max(0.0, 0.015 - arm_min) ** 2
    if detail: return dict(dir_err=dir_err, arm_min=arm_min, arm_mean=arm_mean, stroke=stroke, force=force, L_mean=float(L.mean()), pen=pen)
    return dir_err + pen


def main():
    d = np.load(os.path.join(ROOT, "report/newmus_V_design_k8.npz")); V, c = d["V"], d["c"]
    Q = poses(); print(f"자세 표본 {len(Q)}개", flush=True)
    tic = time.time(); fr = Frames(os.path.join(ROOT, "data/subject3/subject3_scaled.osim"), Q); print(f"뼈 변환 저장 {time.time()-tic:.0f}s", flush=True)
    out = {"coords": COORDS, "bodies": BODIES, "muscles": []}
    for i in range(8):
        v = V[i]; chain = chain_for(v)
        bounds = [b for n in chain for b in BOX[BODIES[n]]]
        tic = time.time()
        res = differential_evolution(cost, bounds, args=(fr, chain, v, c[i]), seed=i, maxiter=400, popsize=25, tol=1e-7, polish=False)
        res2 = minimize(cost, res.x, args=(fr, chain, v, c[i]), method="Powell", bounds=bounds)
        x = res2.x if res2.fun < res.fun else res.x
        det = cost(x, fr, chain, v, c[i], detail=True)
        r, _ = moment_arms(fr, chain, x.reshape(len(chain), 3)); rm = r.mean(axis=0)
        ang = np.degrees(np.arccos(np.clip(rm @ v / np.linalg.norm(rm), -1, 1)))
        print(f"M{i}: 바디 {[BODIES[n] for n in chain]}  방향오차 {det['dir_err']:.3f} (평균 모멘트암과 설계 방향 각도 {ang:.1f}°)  유효 모멘트암 {det['arm_min']*100:.1f}~{det['arm_mean']*100:.1f} cm  스트로크 {det['stroke']:.2f}  필요 힘 {det['force']:.0f} N  길이 {det['L_mean']*100:.0f} cm  벌점 {det['pen']:.3f}  {time.time()-tic:.0f}s", flush=True)
        print(f"     설계 v = {np.round(v, 2)}   평균 r/|r| = {np.round(rm / np.linalg.norm(rm), 2)}   r 범위(cm) = {np.round(r.min(axis=0)*100, 1)} ~ {np.round(r.max(axis=0)*100, 1)}", flush=True)
        out["muscles"].append(dict(name=f"M{i}", v=v.tolist(), capacity_Nm=float(c[i]), chain=[BODIES[n] for n in chain],
                                   points=x.reshape(len(chain), 3).tolist(), force_N=float(det["force"]), arm_min=det["arm_min"], arm_mean=det["arm_mean"],
                                   stroke=det["stroke"], dir_err=det["dir_err"], angle_deg=float(ang), r_mean=rm.tolist()))
    json.dump(out, open(os.path.join(ROOT, "report/newmus_paths_k8.json"), "w"), indent=1, ensure_ascii=False)
    build_model(out); plot(out, fr, V)


def build_model(out):
    """피험자3 모델(근육 제거)에 PathActuator 8개 × 좌우. 왼쪽은 z 부호를 뒤집는다."""
    import opensim as osim
    proc = osim.ModelProcessor(os.path.join(ROOT, "data/subject3/subject3_scaled.osim")); proc.append(osim.ModOpRemoveMuscles()); model = proc.process()
    colors = [(0.9, 0.2, 0.2), (0.2, 0.6, 0.9), (0.2, 0.8, 0.3), (0.9, 0.6, 0.1), (0.6, 0.3, 0.8), (0.1, 0.7, 0.7), (0.8, 0.3, 0.6), (0.5, 0.5, 0.1)]
    for side, sgn in (("r", 1.0), ("l", -1.0)):
        for i, m in enumerate(out["muscles"]):
            pa = osim.PathActuator(); pa.setName(f"{m['name']}_{side}"); pa.setOptimalForce(m["force_N"]); pa.setMinControl(0.0); pa.setMaxControl(1.0)
            for k, (bn, pt) in enumerate(zip(m["chain"], m["points"])):
                body = bn if bn == "pelvis" else bn[:-1] + side
                pa.addNewPathPoint(f"{m['name']}_{side}_p{k}", model.getBodySet().get(body), osim.Vec3(pt[0], pt[1], sgn * pt[2]))
            pa.updGeometryPath().setDefaultColor(osim.Vec3(*colors[i]))
            model.addForce(pa)
    model.finalizeConnections(); model.initSystem()
    p = os.path.join(ROOT, "report/moco/newmus_paths_subject3_model.osim"); model.printToXML(p); print("모델 저장", p)


def plot(out, fr, V):
    fig, axes = plt.subplots(2, 4, figsize=(16, 7))
    labels = ["고관절 굴곡", "고관절 내전", "무릎", "발목"]
    for i, (ax, m) in enumerate(zip(axes.ravel(), out["muscles"])):
        chain = [BODIES.index(b) for b in m["chain"]]; r, _ = moment_arms(fr, chain, np.array(m["points"]))
        x = np.arange(4); w = 0.38
        ax.bar(x - w / 2, np.array(m["v"]) * m["arm_mean"] * 100, w, color="0.6", label="설계 방향 × 평균 모멘트암")
        ax.bar(x + w / 2, r.mean(axis=0) * 100, w, color="#1e8449", label="경로의 평균 모멘트암")
        ax.errorbar(x + w / 2, r.mean(axis=0) * 100, yerr=[r.mean(axis=0) * 100 - r.min(axis=0) * 100, r.max(axis=0) * 100 - r.mean(axis=0) * 100], fmt="none", ecolor="k", lw=0.8)
        ax.axhline(0, color="k", lw=0.6); ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
        ax.set_title(f"{m['name']}: {'→'.join(b.replace('_r','') for b in m['chain'])}  각도차 {m['angle_deg']:.0f}°  {m['force_N']:.0f} N  스트로크 {m['stroke']:.2f}", fontsize=9)
        if i == 0: ax.legend(fontsize=7); ax.set_ylabel("모멘트암 (cm)")
    fig.suptitle("케이블 경로 설계 — 설계 방향 vs 경로가 실제로 내는 모멘트암 (막대 = 자세 평균, 선 = 자세에 따른 범위)", fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(ROOT, "report/fig_path_design.png"), dpi=150); print("그림 저장")


if __name__ == "__main__":
    main()
