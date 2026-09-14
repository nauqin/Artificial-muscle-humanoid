"""새 근육 8개의 토크 방향을 ID 토크에서 직접 설계한다 — SO 를 거치지 않는다.

문제: 단위 방향 v_i (8×4, 각 행 노름 1) 와 용량 c_i 를 골라, 모든 프레임의 관절 토크 τ_t 를
    τ_t = Σ_i v_i u_it,   0 ≤ u_it ≤ c_i
로 내면서 Σ c_i (총 액추에이터 용량, Nm) 를 최소로. 안쪽(V 고정)은 LP, 바깥(V)은 Powell 다중 시작.
데이터: 7명 걷기 42걸음(우리 ID, 지지 구간) + 논문 Mocap ID 의 squats1·STS1, 좌우 합침,
토크는 (질량×키) 로 정규화해 피험자3 몸(63.5 kg, 1.69 m) 크기로 맞춘다.

    python report/design_directions.py [--frames 600] [--maxfev 300] [--with-dj]
→ report/newmus_V_design_k8.npz (V 단위방향, c 용량, V_nm = V*c), design_directions.csv, fig_design_directions.png
"""
import argparse, csv, glob, os, sys, time
import numpy as np
from scipy.optimize import linprog, minimize
from scipy.sparse import coo_matrix
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import addb_io as io                                   # noqa: E402
LAB = os.path.join(os.path.dirname(ROOT), "LabValidation_withoutVideos")
KEEP = ("subject3", "subject4", "subject5", "subject7", "subject8", "subject10", "subject11")
COORDS = ("hip_flexion", "hip_adduction", "knee_angle", "ankle_angle")
REF = 63.5 * 1.69                                      # 피험자3 질량×키
for f in font_manager.findSystemFonts():
    if "AppleSDGothicNeo" in f: plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
plt.rcParams["axes.unicode_minus"] = False


def body(sub):
    m = h = None
    for line in open(os.path.join(LAB, sub, "sessionMetadata.yaml")):
        if line.startswith("mass_kg"): m = float(line.split(":")[1])
        if line.startswith("height_m"): h = float(line.split(":")[1])
    return m * h


def torque_rows(mot, scale):
    out = []
    for side in "rl":
        cols = [f"{c}_{side}_moment" for c in COORDS]
        if all(c in mot.names for c in cols):
            out.append(np.column_stack([mot.column(c) for c in cols]) * scale)
    return out


def dataset(with_dj=False):
    parts = {}
    for sub in KEEP:
        scale = REF / body(sub)
        for f in sorted(glob.glob(os.path.join(ROOT, "out/real", sub, "walking*", "*_id.sto"))):
            parts.setdefault("walking", []).extend(torque_rows(io.read_mot(f), scale))
        for task in ("squats1", "STS1") + (("DJ1", "DJ2", "DJ3", "DJ4") if with_dj else ()):
            f = os.path.join(LAB, sub, "OpenSimData/Mocap/ID", f"{task}.sto")
            if os.path.exists(f):
                parts.setdefault(task.rstrip("1234"), []).extend(torque_rows(io.read_mot(f), scale))
    return {k: np.vstack(v) for k, v in parts.items()}


def capacity_lp(V, tau, slack_w=50.0):
    """V 고정. min Σc + w·Σ|slack|  s.t.  Vᵀu_t + slack_t = τ_t, 0 ≤ u_t ≤ c. 반환 (Σc, c, slack 최대)."""
    T, n = tau.shape; k = V.shape[0]
    nu, ns = T * k, T * n
    # 변수: u (T*k), s+ (T*n), s- (T*n), c (k)
    nv = nu + 2 * ns + k
    cost = np.concatenate([np.zeros(nu), np.full(2 * ns, slack_w), np.ones(k)])
    # 등식: 각 (t, j): Σ_i V[i,j] u_ti + s+_tj − s−_tj = τ_tj
    ri = np.repeat(np.arange(T * n), k); ci = (np.arange(T)[:, None, None] * k + np.arange(k)[None, None, :]).repeat(n, axis=1).ravel()
    vals = np.tile(V.T.ravel(), T)   # V.T[j,i] 순서 = (j 바깥, i 안쪽)
    rows = [ri, np.arange(ns), np.arange(ns)]; cols = [ci, nu + np.arange(ns), nu + ns + np.arange(ns)]; data = [vals, np.ones(ns), -np.ones(ns)]
    Aeq = coo_matrix((np.concatenate(data), (np.concatenate(rows), np.concatenate(cols))), shape=(ns, nv)).tocsr()
    beq = tau.ravel()
    # 부등식: u_ti − c_i ≤ 0
    ri2 = np.arange(nu); ci2 = np.arange(nu); ci3 = nu + 2 * ns + np.tile(np.arange(k), T)
    Aub = coo_matrix((np.concatenate([np.ones(nu), -np.ones(nu)]), (np.concatenate([ri2, ri2]), np.concatenate([ci2, ci3]))), shape=(nu, nv)).tocsr()
    res = linprog(cost, A_ub=Aub, b_ub=np.zeros(nu), A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs")
    if not res.success: return np.inf, None, np.inf
    c = res.x[nu + 2 * ns:]; s = res.x[nu:nu + 2 * ns]
    return float(c.sum() + slack_w * s.sum()), c, float(s.max())


def unit(V):
    return V / np.maximum(np.linalg.norm(V, axis=1, keepdims=True), 1e-9)


def axis_aligned():
    return np.vstack([np.eye(4), -np.eye(4)])


def subsample(tau, n, seed=0):
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(tau), min(n, len(tau)), replace=False)
    ext = np.concatenate([tau.argmax(axis=0), tau.argmin(axis=0)])   # 좌표별 극값은 항상 포함
    return tau[np.unique(np.concatenate([idx, ext]))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", type=int, default=600); ap.add_argument("--maxfev", type=int, default=300)
    ap.add_argument("--starts", type=int, default=3); ap.add_argument("--with-dj", action="store_true")
    a = ap.parse_args()
    data = dataset(a.with_dj)
    tau_all = np.vstack(list(data.values()))
    print("데이터:", {k: v.shape[0] for k, v in data.items()}, "합계", tau_all.shape[0], "프레임 (좌우 합침, 피험자3 크기)")
    print("좌표별 피크 (Nm):", dict(zip(COORDS, np.round(np.abs(tau_all).max(axis=0), 1))))
    tau_fit = subsample(tau_all, a.frames)
    print(f"최적화 프레임 {len(tau_fit)}개")

    def evaluate(V, tau=tau_fit):
        return capacity_lp(unit(V), tau)

    results = {}
    # 기준 1: 관절마다 굴근·신근
    cA = evaluate(axis_aligned(), tau_all); results["축 정렬 8개 (관절마다 굴근·신근)"] = (axis_aligned(), cA)
    print(f"축 정렬 8개: 총용량 {cA[0]:.0f} Nm, slack {cA[2]:.2f}")
    # 기준 2: NMF 에서 나온 8개 (subject3 6걸음 합침, 오른다리)
    npz = os.path.join(ROOT, "report/newmus_V_subject3_pooled_k8.npz")
    if os.path.exists(npz):
        Vn = np.load(npz)["V_r"]; cN = evaluate(Vn, tau_all); results["NMF 8개 (subject3 합침)"] = (unit(Vn), cN)
        print(f"NMF 8개: 총용량 {cN[0]:.0f} Nm, slack {cN[2]:.2f}")

    # 최적화: 다중 시작
    starts = [axis_aligned()] + ([unit(Vn)] if os.path.exists(npz) else [])
    rng = np.random.default_rng(1)
    while len(starts) < a.starts + 2: starts.append(unit(rng.normal(size=(8, 4))))
    best = (np.inf, None)
    for si, V0 in enumerate(starts):
        tic = time.time()
        f = lambda x: evaluate(x.reshape(8, 4))[0]
        r = minimize(f, unit(V0).ravel(), method="Powell", options={"maxfev": a.maxfev, "xtol": 1e-2, "ftol": 1e-3})
        V = unit(r.x.reshape(8, 4)); full = evaluate(V, tau_all)
        print(f"  시작 {si}: {r.nfev}회 {time.time()-tic:.0f}s  부분집합 {r.fun:.0f} → 전체 {full[0]:.0f} Nm (slack {full[2]:.2f})", flush=True)
        if full[0] < best[0]: best = (full[0], V, full)
    V, full = best[1], best[2]; results["최적화 8개"] = (V, full)

    # 정리
    print("\n설계안 (행 = 근육, 열 =", COORDS, ") 단위방향 × 용량 = Nm")
    c = full[1]; Vnm = V * c[:, None]
    for i in range(8): print(f"  M{i}: 용량 {c[i]:6.1f} Nm   방향 {np.round(V[i], 2)}   Nm {np.round(Vnm[i], 1)}")
    # 과제별 필요 용량 (설계 V 고정)
    print("\n과제별 (설계 V 고정) 필요 총용량:")
    per_task = {}
    for k, v in data.items():
        r = capacity_lp(V, v); per_task[k] = r[0]; print(f"  {k:10s} {r[0]:7.0f} Nm  slack {r[2]:.2f}")
    np.savez(os.path.join(ROOT, "report/newmus_V_design_k8.npz"), V=V, c=c, V_nm=Vnm, coords=np.array(COORDS),
             axis_c=cA[1], nmf_c=(cN[1] if os.path.exists(npz) else np.zeros(8)))
    with open(os.path.join(ROOT, "report/design_directions.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["design", "total_capacity_Nm", "slack_max"] + [f"c{i}" for i in range(8)])
        for name, (Vx, rx) in results.items(): w.writerow([name, f"{rx[0]:.1f}", f"{rx[2]:.3f}"] + [f"{x:.1f}" for x in rx[1]])
        for k, v in per_task.items(): w.writerow([f"design on {k}", f"{v:.1f}", ""])

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw={"width_ratios": [1, 1.3, 1.3]})
    ax = axes[0]; names = list(results); tot = [results[n][1][0] for n in names]
    ax.bar(range(len(names)), tot, color=["#95a5a6", "#e59866", "#1e8449"][:len(names)])
    for i, v in enumerate(tot): ax.text(i, v, f"{v:.0f}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(range(len(names))); ax.set_xticklabels([n.split(" (")[0] for n in names], fontsize=9); ax.set_ylabel("총 액추에이터 용량 Σc (Nm)")
    ax.set_title("같은 토크를 내는 데 필요한 용량", fontsize=10)
    ax = axes[1]; im = ax.imshow(Vnm, cmap="RdBu_r", vmin=-np.abs(Vnm).max(), vmax=np.abs(Vnm).max(), aspect="auto")
    ax.set_xticks(range(4)); ax.set_xticklabels(["고관절 굴곡", "고관절 내전", "무릎", "발목"], fontsize=9); ax.set_yticks(range(8)); ax.set_yticklabels([f"M{i} ({c[i]:.0f} Nm)" for i in range(8)], fontsize=9)
    for i in range(8):
        for j in range(4): ax.text(j, i, f"{Vnm[i, j]:.0f}", ha="center", va="center", fontsize=8)
    ax.set_title("설계안: 근육별 관절 토크 (Nm, 최대 신호일 때)", fontsize=10); plt.colorbar(im, ax=ax, fraction=0.04)
    ax = axes[2]; ks = list(per_task); ax.bar(range(len(ks)), [per_task[k] for k in ks], color="#1e8449")
    ax.set_xticks(range(len(ks))); ax.set_xticklabels(ks, fontsize=9); ax.set_ylabel("Nm"); ax.set_title("설계 V 고정, 과제별 필요 총용량", fontsize=10)
    fig.suptitle("새 근육 8개 방향 — ID 토크에서 직접 설계 (7명, 걷기 42 + 스쿼트 + 앉았다서기, 피험자3 크기)", fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(ROOT, "report/fig_design_directions.png"), dpi=150); print("그림 저장")


if __name__ == "__main__":
    main()
