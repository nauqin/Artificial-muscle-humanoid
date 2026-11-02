"""새 근육 8개의 토크 방향을 ID 토크에서 직접 설계한다 — SO 를 거치지 않는다.

문제: 단위 방향 v_i (8×4, 각 행 노름 1) 와 용량 c_i 를 골라, 모든 프레임의 관절 토크 τ_t 를
    τ_t = Σ_i v_i u_it,   0 ≤ u_it ≤ c_i
로 내면서 Σ c_i (총 액추에이터 용량, Nm) 를 최소로. 안쪽(V 고정)은 LP, 바깥(V)은 Powell 다중 시작.
데이터: 7명 걷기 42걸음(우리 ID, 지지 구간) + 논문 Mocap ID 의 squats1·STS1, 좌우 합침,
토크는 (질량×키) 로 정규화해 피험자3 몸(63.5 kg, 1.69 m) 크기로 맞춘다.

    python report/design_directions.py [--frames 600] [--maxfev 300] [--with-dj]            # 하체 (기본, 9월 설계 재현)
    python report/design_directions.py --space trunk --tasks report/upper/tasks.csv --k 4,5,6  # 상체
→ report/newmus_V_design_k8.npz (V 단위방향, c 용량, V_nm = V*c), design_directions.csv, fig_design_directions.png
  (--space 가 lower 가 아니면 newmus_V_design_<space>_k<k>.npz 등으로 저장)

`--tasks` CSV 열: dataset, subject, trial, ik, id, mass_kg, height_m, task, hand_load_kg — 한 줄이 trial 하나.
id 는 `_id.sto` 경로(ROOT 기준 상대 가능). 토크는 REF/(mass×height) 로 정규화해 같은 몸 크기로 맞춘다.
"""
import argparse, csv, glob, os, sys, time
import numpy as np
from scipy.optimize import linprog, minimize
from scipy.sparse import coo_matrix
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import addb_io as io                                   # noqa: E402
import design_space as ds                              # noqa: E402
LAB = os.path.join(os.path.dirname(ROOT), "LabValidation_withoutVideos")
KEEP = ("subject3", "subject4", "subject5", "subject7", "subject8", "subject10", "subject11")
SP = ds.get("lower"); COORDS = SP["base"]              # main() 에서 --space 로 바뀐다
REF = 63.5 * 1.69                                      # 피험자3 질량×키
for f in font_manager.findSystemFonts():
    if "AppleSDGothicNeo" in f: plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
plt.rcParams["axes.unicode_minus"] = False


def mass_kg(sub):
    for line in open(os.path.join(LAB, sub, "sessionMetadata.yaml")):
        if line.startswith("mass_kg"): return float(line.split(":")[1])


def body(sub):
    m = h = None
    for line in open(os.path.join(LAB, sub, "sessionMetadata.yaml")):
        if line.startswith("mass_kg"): m = float(line.split(":")[1])
        if line.startswith("height_m"): h = float(line.split(":")[1])
    return m * h


def torque_rows(mot, scale):
    out = []
    for side in SP["sides"]:
        cols = [ds.full(c, side) + "_moment" for c in COORDS]
        if all(c in mot.names for c in cols):
            out.append(np.column_stack([mot.column(c) for c in cols]) * scale)
    return out


def dataset(with_dj=False):
    parts = {}
    for sub in KEEP:
        scale = REF / body(sub)
        for f in sorted(glob.glob(os.path.join(ROOT, "out/real", sub, "walking*", "*_id.sto"))):
            parts.setdefault("walking", []).extend(torque_rows(io.read_mot(f), scale))
        for task in ("squats1",) + (("DJ1", "DJ2", "DJ3", "DJ4") if with_dj else ()):   # STS 는 뺐다(9/15 결정): 의자 힘이 힘판에 없어 앉은 구간 ID 가 틀리고, 남는 구간은 스쿼트와 같다
            f = os.path.join(LAB, sub, "OpenSimData/Mocap/ID", f"{task}.sto")
            if os.path.exists(f):
                mot = io.read_mot(f); rows = torque_rows(mot, scale)
                if task.startswith("STS"):
                    # 의자는 힘판 위에 없었다 (힘판 3 = 0). 앉아 있는 동안 발 힘판은 체중의 20~35 % 뿐이라 ID 가 틀리다
                    # → 발 힘판 합이 체중의 85 % 이상인 프레임(의자에서 떨어진 상태)만 쓴다.
                    fz = io.read_mot(os.path.join(LAB, sub, "ForceData", f"{task}_forces.mot"))
                    tot = np.interp(mot.time, fz.time, fz.column("R_ground_force_vy") + fz.column("L_ground_force_vy"))
                    keep = tot >= 0.85 * mass_kg(sub) * 9.80665; rows = [r[keep] for r in rows]
                parts.setdefault(task.rstrip("1234"), []).extend(rows)
    return {k: np.vstack(v) for k, v in parts.items()}


def dataset_from_tasks(path):
    """tasks.csv 의 trial 들을 과제별로 모은다. 좌우(sides)를 같이 쌓고 몸 크기를 REF 로 맞춘다."""
    parts, n_rows = {}, 0
    for row in csv.DictReader(open(path, encoding="utf-8")):
        if not row.get("id") or row["id"].startswith("#"): continue
        f = row["id"] if os.path.isabs(row["id"]) else os.path.join(ROOT, row["id"])
        scale = REF / (float(row["mass_kg"]) * float(row["height_m"]))
        rows = torque_rows(io.read_mot(f), scale)
        if not rows: print(f"  [건너뜀] {row['trial']}: 열 {COORDS} 없음"); continue
        parts.setdefault(row["task"], []).extend(rows); n_rows += 1
    print(f"tasks.csv: trial {n_rows}개, 과제 {list(parts)}")
    return {k: np.vstack(v) for k, v in parts.items()}


def capacity_lp(V, tau, slack_w=50.0, time_limit=None):
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
    res = linprog(cost, A_ub=Aub, b_ub=np.zeros(nu), A_eq=Aeq, b_eq=beq, bounds=(0, None), method="highs", options=({"time_limit": time_limit} if time_limit else {}))
    if not res.success: return np.inf, None, np.inf
    c = res.x[nu + 2 * ns:]; s = res.x[nu:nu + 2 * ns]
    return float(c.sum() + slack_w * s.sum()), c, float(s.max())


def unit(V):
    return V / np.maximum(np.linalg.norm(V, axis=1, keepdims=True), 1e-9)


def axis_aligned(n=None):
    n = n or len(COORDS)
    return np.vstack([np.eye(n), -np.eye(n)])


def subsample(tau, n, seed=0):
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(tau), min(n, len(tau)), replace=False)
    ext = np.concatenate([tau.argmax(axis=0), tau.argmin(axis=0)])   # 좌표별 극값은 항상 포함
    return tau[np.unique(np.concatenate([idx, ext]))]


def main():
    global SP, COORDS
    ap = argparse.ArgumentParser()
    ap.add_argument("--frames", type=int, default=600); ap.add_argument("--maxfev", type=int, default=300)
    ap.add_argument("--starts", type=int, default=3); ap.add_argument("--with-dj", action="store_true")
    ap.add_argument("--space", default="lower", help=f"설계 공간 {list(ds.SPACES)}")
    ap.add_argument("--tasks", default=None, help="trial 목록 CSV (없으면 하체 기본 데이터)")
    ap.add_argument("--k", default=None, help="근육 수 (쉼표로 여러 개: 4,5,6). 기본 = 공간의 k")
    a = ap.parse_args()
    SP = ds.get(a.space); COORDS = SP["base"]
    ks = [int(x) for x in a.k.split(",")] if a.k else [SP["k"]]
    data = dataset_from_tasks(a.tasks if os.path.isabs(a.tasks) else os.path.join(ROOT, a.tasks)) if a.tasks else dataset(a.with_dj)
    summary = []
    for k in ks:
        summary.append(design(a, data, k))
    if len(ks) > 1:
        print("\n근육 수별 요약 (slack 0 이면 모든 과제 토크를 낼 수 있다):")
        for k, tot, slack, worst in summary: print(f"  k={k}: 총용량 {tot:.0f} Nm  slack {slack:.2f}  전 프레임 최대 부족 {worst:.2f} Nm")


def design(a, data, k):
    n = len(COORDS); legacy = (SP["name"] == "lower" and k == 8 and not a.tasks)
    tag = "" if legacy else f"_{SP['name']}_k{k}"
    tau_all = np.vstack(list(data.values()))
    print(f"\n=== {SP['name']} 근육 {k}개 ===  데이터:", {kk: v.shape[0] for kk, v in data.items()}, "합계", tau_all.shape[0], "프레임 (side 합침, 피험자3 크기)")
    print("좌표별 피크 (Nm):", dict(zip(COORDS, np.round(np.abs(tau_all).max(axis=0), 1))))
    tau_fit = subsample(tau_all, a.frames)
    tau_eval = subsample(tau_all, 3000, seed=7)
    print(f"최적화 프레임 {len(tau_fit)}개, 평가 프레임 {len(tau_eval)}개 (전체 LP 는 조건 나쁜 V 에서 한없이 길어져 표본으로 잰다)")

    def evaluate(V, tau=tau_fit):
        if tau is tau_all: tau = tau_eval
        return capacity_lp(unit(V), tau, time_limit=(5.0 if tau is tau_fit else 120.0))

    def verify_all(V, c):
        """용량 c 로 전 프레임을 정말 낼 수 있나 — 프레임별 0 ≤ u ≤ c 최소제곱, 남는 토크 최대 (Nm)."""
        from scipy.optimize import lsq_linear
        Vu = unit(V); worst = 0.0
        for row in tau_all[::3]:
            r = lsq_linear(Vu.T, row, bounds=(np.zeros(len(c)), np.maximum(c, 1e-9)), lsmr_tol="auto")
            worst = max(worst, float(np.abs(Vu.T @ r.x - row).max()))
        return worst

    results = {}
    # 기준 1: 관절마다 굴근·신근
    cA = evaluate(axis_aligned(), tau_all); results[f"축 정렬 {2*n}개 (관절마다 굴근·신근)"] = (axis_aligned(), cA)
    print(f"축 정렬 {2*n}개: 총용량 {cA[0]:.0f} Nm, slack {cA[2]:.2f}")
    # 기준 2: NMF 에서 나온 8개 (subject3 6걸음 합침, 오른다리) — 하체 k=8 에서만
    npz = os.path.join(ROOT, "report/newmus_V_subject3_pooled_k8.npz") if (SP["name"] == "lower" and k == 8) else "/nonexistent"
    if os.path.exists(npz):
        Vn = np.load(npz)["V_r"]; cN = evaluate(Vn, tau_all); results["NMF 8개 (subject3 합침)"] = (unit(Vn), cN)
        print(f"NMF 8개: 총용량 {cN[0]:.0f} Nm, slack {cN[2]:.2f}")

    # 최적화: 다중 시작
    rng = np.random.default_rng(1)
    ax = axis_aligned()                                      # 2n 행. k 가 다르면 앞 k 행을 쓰거나 난수로 채운다
    V0a = ax[:k] if k <= 2 * n else np.vstack([ax, unit(rng.normal(size=(k - 2 * n, n)))])
    starts = [V0a] + ([unit(Vn)] if os.path.exists(npz) else [])
    while len(starts) < a.starts + 2: starts.append(unit(rng.normal(size=(k, n))))
    class Stop(Exception): pass
    best = (np.inf, None)
    for si, V0 in enumerate(starts):
        tic = time.time(); state = {"n": 0, "best": (np.inf, None), "t_lp": 0.0}
        def f(x):
            t0 = time.time(); v = evaluate(x.reshape(k, n))[0]; state["t_lp"] += time.time() - t0
            state["n"] += 1
            if v < state["best"][0]: state["best"] = (v, x.copy())
            if state["n"] % 100 == 0:
                print(f"    시작 {si} 평가 {state['n']}회  현재 최선 {state['best'][0]:.0f} Nm  LP 평균 {state['t_lp']/state['n']*1000:.0f} ms", flush=True)
            if state["n"] >= a.maxfev: raise Stop
            return v
        try:
            minimize(f, unit(V0).ravel(), method="Powell", options={"maxfev": a.maxfev, "xtol": 1e-2, "ftol": 1e-3})
        except Stop:
            pass
        V = unit(state["best"][1].reshape(k, n)); full = evaluate(V, tau_all)
        print(f"  시작 {si}: {state['n']}회 {time.time()-tic:.0f}s  부분집합 {state['best'][0]:.0f} → 전체 {full[0]:.0f} Nm (slack {full[2]:.2f})", flush=True)
        if full[0] < best[0]: best = (full[0], V, full)
    V, full = best[1], best[2]; results[f"최적화 {k}개"] = (V, full)
    worst_all = verify_all(V, full[1])
    print(f"전 프레임 검산 (1/3 표본): 설계 용량으로 못 내는 토크 최대 {worst_all:.2f} Nm  (축 정렬: {verify_all(axis_aligned(), cA[1]):.2f} Nm)", flush=True)

    # 정리
    print("\n설계안 (행 = 근육, 열 =", COORDS, ") 단위방향 × 용량 = Nm")
    c = full[1]; Vnm = V * c[:, None]
    for i in range(k): print(f"  M{i}: 용량 {c[i]:6.1f} Nm   방향 {np.round(V[i], 2)}   Nm {np.round(Vnm[i], 1)}")
    # 과제별 필요 용량 (설계 V 고정)
    print("\n과제별 (설계 V 고정) 필요 총용량:")
    per_task = {}
    for kk, v in data.items():
        r = capacity_lp(V, subsample(v, 3000, seed=3), time_limit=120.0); per_task[kk] = r[0]; print(f"  {kk:10s} {r[0]:7.0f} Nm  slack {r[2]:.2f}")
    np.savez(os.path.join(ROOT, f"report/newmus_V_design{tag or '_k8'}.npz"), V=V, c=c, V_nm=Vnm, coords=np.array(COORDS), space=SP["name"],
             axis_c=cA[1], nmf_c=(cN[1] if os.path.exists(npz) else np.zeros(k)))
    with open(os.path.join(ROOT, f"report/design_directions{tag}.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["design", "total_capacity_Nm", "slack_max"] + [f"c{i}" for i in range(k)])
        for name, (Vx, rx) in results.items(): w.writerow([name, f"{rx[0]:.1f}", f"{rx[2]:.3f}"] + [f"{x:.1f}" for x in rx[1]])
        for kk, v in per_task.items(): w.writerow([f"design on {kk}", f"{v:.1f}", ""])

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4), gridspec_kw={"width_ratios": [1, 1.3, 1.3]})
    ax = axes[0]; names = list(results); tot = [results[n][1][0] for n in names]
    ax.bar(range(len(names)), tot, color=["#95a5a6", "#e59866", "#1e8449"][:len(names)])
    for i, v in enumerate(tot): ax.text(i, v, f"{v:.0f}", ha="center", va="bottom", fontsize=9)
    ax.set_xticks(range(len(names))); ax.set_xticklabels([n.split(" (")[0] for n in names], fontsize=9); ax.set_ylabel("총 액추에이터 용량 Σc (Nm)")
    ax.set_title("같은 토크를 내는 데 필요한 용량", fontsize=10)
    ax = axes[1]; im = ax.imshow(Vnm, cmap="RdBu_r", vmin=-np.abs(Vnm).max(), vmax=np.abs(Vnm).max(), aspect="auto")
    ax.set_xticks(range(n)); ax.set_xticklabels(SP["labels"], fontsize=9); ax.set_yticks(range(k)); ax.set_yticklabels([f"M{i} ({c[i]:.0f} Nm)" for i in range(k)], fontsize=9)
    for i in range(k):
        for j in range(n): ax.text(j, i, f"{Vnm[i, j]:.0f}", ha="center", va="center", fontsize=8)
    ax.set_title("설계안: 근육별 관절 토크 (Nm, 최대 신호일 때)", fontsize=10); plt.colorbar(im, ax=ax, fraction=0.04)
    ax = axes[2]; ks = list(per_task); ax.bar(range(len(ks)), [per_task[k] for k in ks], color="#1e8449")
    ax.set_xticks(range(len(ks))); ax.set_xticklabels(ks, fontsize=9); ax.set_ylabel("Nm"); ax.set_title("설계 V 고정, 과제별 필요 총용량", fontsize=10)
    fig.suptitle(f"새 근육 {k}개 방향 ({SP['name']}) — ID 토크에서 직접 설계, 피험자3 크기", fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(ROOT, f"report/fig_design_directions{tag}.png"), dpi=150); print("그림 저장")
    return k, full[0], full[2], worst_all


if __name__ == "__main__":
    main()
