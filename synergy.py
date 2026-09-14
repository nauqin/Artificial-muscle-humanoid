"""NMF 시너지 분해, VAF로 개수 선택, 표·그림.

NMF는 Lee & Seung 곱셈 업데이트(Frobenius norm)를 직접 구현했다. scikit-learn
기본 solver와 같은 알고리즘이고, 의존성을 늘리지 않으려고 직접 짰다. 초기값에
따라 국소해에 빠지므로 여러 번 재시작해 최선을 고른다.

개수는 VAF로 고르되 **두 기준을 나란히 보고한다.**

  global  전체 VAF ≥ 90 %
  dual    전체 VAF ≥ 90 % 그리고 모든 근육의 개별 VAF ≥ 75 %

전체 VAF만 보면 활성도가 큰 근육 몇 개에 가려 조용한 근육이 전혀 설명되지
않아도 통과한다. 반대로 `dual`은 SO가 모멘트암 겹치는 근육 사이에 부하를
임의로 나눈 것까지 재현하라고 요구해서 실제 데이터에서 5~12로 흔들린다.
그래서 **본 결과는 `global`, `dual`은 "이 근육은 개별 해석 금지" 진단 지표**다.
"""

from __future__ import annotations

import _env  # noqa: F401  환경 확인 — 다른 import 보다 먼저

import argparse
import os

import numpy as np

import addb_io as io

VAF_TOTAL = 0.90
VAF_MUSCLE = 0.75
K_MAX = 10
#: 이보다 조용한 근육은 분석에서 제외 (정규화 전 활성도 피크)
MIN_PEAK = 0.05
#: 표에 표시할 가중치 하한
WEIGHT_SHOW = 0.30
N_RESTARTS = 20
MAX_ITER = 3000
TOL = 1e-8
#: 수렴 검사 주기. 25로 두면 k가 클 때 덜 수렴한 채로 멈춘다.
CHECK_EVERY = 10


# ─────────────────────────────────────────────────────────────────────────────
# NMF
# ─────────────────────────────────────────────────────────────────────────────


def nmf(A, k, *, seed=0, n_restarts=N_RESTARTS, max_iter=MAX_ITER, tol=TOL):
    """A(T×M) ≈ C(T×k) @ W(k×M). 최선의 재시작 결과를 돌려준다."""
    rng = np.random.default_rng(seed)
    best = None
    scale = np.sqrt(A.mean() / k) if A.mean() > 0 else 1.0
    for _ in range(n_restarts):
        C = rng.random((A.shape[0], k)) * scale + 1e-6
        W = rng.random((k, A.shape[1])) * scale + 1e-6
        prev = np.inf
        for it in range(max_iter):
            W *= (C.T @ A) / (C.T @ C @ W + 1e-10)
            C *= (A @ W.T) / (C @ W @ W.T + 1e-10)
            if it % CHECK_EVERY == 0:
                err = np.linalg.norm(A - C @ W)
                if prev - err < tol * max(prev, 1e-12):
                    break
                prev = err
        err = np.linalg.norm(A - C @ W)
        if best is None or err < best[0]:
            best = (err, C.copy(), W.copy())
    return best[1], best[2]


def vaf(A, recon):
    """전체 VAF와 근육별 VAF. 평균을 빼지 않는 정의(비음수 데이터 관례)."""
    resid = A - recon
    total = 1.0 - float(np.sum(resid ** 2) / max(np.sum(A ** 2), 1e-12))
    denom = np.sum(A ** 2, axis=0)
    per = 1.0 - np.sum(resid ** 2, axis=0) / np.maximum(denom, 1e-12)
    return total, per


def sort_synergies(C, W):
    """기여도(활성 총합 × 가중치 합) 내림차순으로."""
    contrib = C.sum(axis=0) * W.sum(axis=1)
    order = np.argsort(-contrib)
    return C[:, order], W[order, :], contrib[order] / max(contrib.sum(), 1e-12)


# ─────────────────────────────────────────────────────────────────────────────
# 개수 선택
# ─────────────────────────────────────────────────────────────────────────────


def scan_k(A, kmax=K_MAX, *, seed=0):
    """k=1..kmax 를 전부 돌려 (k, 전체VAF, 최저근육VAF, 미달 근육 수)."""
    rows = []
    fits = {}
    for k in range(1, kmax + 1):
        C, W = nmf(A, k, seed=seed)
        total, per = vaf(A, C @ W)
        rows.append((k, total, float(per.min()), int(np.sum(per < VAF_MUSCLE))))
        fits[k] = (C, W, total, per)
    return rows, fits


def choose_k(rows, criterion="dual"):
    """기준을 만족하는 최소 k. 못 찾으면 None.

    `--kmax`까지 기준을 못 채우면 **조용히 아무 k나 고르지 않는다.** 예전에
    kmax=6에 두었더니 기준을 만족하는 k가 없는데도 말없이 6을 골라 그럴듯한
    표를 그린 적이 있다. 실제로는 k=7에서 충족됐다.
    """
    for k, total, worst_muscle, _n_bad in rows:
        if total < VAF_TOTAL:
            continue
        if criterion == "global" or worst_muscle >= VAF_MUSCLE:
            return k
    return None


# ─────────────────────────────────────────────────────────────────────────────
# 데이터 준비
# ─────────────────────────────────────────────────────────────────────────────


def muscle_names(model_path):
    """모델의 근육 이름 집합.

    활성도 파일에는 근육 말고 팔 토크 액추에이터(`shoulder_flex_r`,
    `elbow_flex_r` …)도 함께 들어있다. 그것들을 시너지에 넣으면 근육이 아닌
    성분이 만들어져 개수와 구성이 통째로 어긋난다.
    """
    import opensim as osim

    model = osim.Model(os.path.abspath(model_path))
    model.initSystem()
    muscles = model.getMuscles()
    return {muscles.get(i).getName() for i in range(muscles.getSize())}


def load_activations(paths, side, *, min_peak=MIN_PEAK, muscles=None, log=print):
    """활성도 파일(들)을 이어붙여 (A, 근육이름, 시간, 제외된 근육)."""
    if isinstance(paths, str):
        paths = [paths]

    mats, times, names = [], [], None
    offset = 0.0
    for path in paths:
        mot = io.read_mot(path)
        cols = [n for n in mot.names
                if not n.startswith("reserve_") and _matches_side(n, side)
                and (muscles is None or n in muscles)]
        if names is None:
            names = cols
        else:
            cols = [n for n in names if n in mot.names]
            if len(cols) != len(names):
                raise ValueError(f"근육 집합이 다르다: {path}")
        mats.append(np.column_stack([mot.column(n) for n in cols]))
        times.append(mot.time - mot.t0 + offset)
        offset += mot.duration + (mot.time[1] - mot.time[0] if len(mot.time) > 1 else 0)

    A = np.vstack(mats)
    time = np.concatenate(times)
    A = np.clip(A, 0.0, None)

    peaks = A.max(axis=0)
    keep = peaks >= min_peak
    dropped = [n for n, k in zip(names, keep) if not k]
    if dropped:
        log(f"  조용한 근육 {len(dropped)}개 제외 (피크 < {min_peak}): "
            + ", ".join(dropped[:8]) + ("…" if len(dropped) > 8 else ""))
    A, names = A[:, keep], [n for n, k in zip(names, keep) if k]

    # 근육마다 피크로 나눈다 — 안 그러면 큰 근육이 VAF를 독점한다.
    A = A / A.max(axis=0)
    return A, names, time, dropped


def _matches_side(name, side):
    if side == "both":
        return True
    return name.endswith(f"_{side}")


# ─────────────────────────────────────────────────────────────────────────────
# 보고
# ─────────────────────────────────────────────────────────────────────────────


def render_table(rows, chosen, criterion) -> str:
    lines = [f"  {'k':>3s} {'전체VAF':>8s} {'최저근육':>9s} {'미달':>5s}"]
    for k, total, worst_muscle, n_bad in rows:
        mark = ""
        if k == chosen.get("global"):
            mark += " ←global"
        if k == chosen.get("dual"):
            mark += " ←dual"
        lines.append(f"  {k:3d} {total*100:7.1f}% {worst_muscle*100:8.1f}% {n_bad:5d}{mark}")
    return "\n".join(lines)


def describe(W, names, contrib, weight_show=WEIGHT_SHOW) -> str:
    lines = []
    for i in range(W.shape[0]):
        w = W[i] / max(W[i].max(), 1e-12)
        members = [(names[j], w[j]) for j in np.argsort(-w) if w[j] >= weight_show]
        body = "  ".join(f"{n} {v:.2f}" for n, v in members) or "(표시 기준 미만)"
        flag = "   ← 근육 1개짜리 성분" if len(members) == 1 else ""
        lines.append(f"  시너지 {i+1}  기여 {contrib[i]*100:4.1f}%  근육 {len(members)}개   {body}{flag}")
    return "\n".join(lines)


def plot(path, C, W, names, time, contrib, title=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    k = W.shape[0]
    fig, axes = plt.subplots(k, 2, figsize=(13, 1.9 * k + 1),
                             gridspec_kw={"width_ratios": [2, 1]}, squeeze=False)
    for i in range(k):
        ax = axes[i][0]
        w = W[i] / max(W[i].max(), 1e-12)
        ax.bar(range(len(names)), w, color="tab:blue")
        ax.axhline(WEIGHT_SHOW, color="0.6", lw=0.6, ls="--")
        ax.set_ylim(0, 1.05)
        ax.set_ylabel(f"W{i+1}\n{contrib[i]*100:.0f}%")
        if i == k - 1:
            ax.set_xticks(range(len(names)))
            ax.set_xticklabels(names, rotation=90, fontsize=6)
        else:
            ax.set_xticks([])
        axes[i][1].plot(time, C[:, i], color="tab:red")
        axes[i][1].set_ylabel(f"C{i+1}")
        if i == k - 1:
            axes[i][1].set_xlabel("time [s]")
    axes[0][0].set_title(title or os.path.basename(path))
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    fig.savefig(path, dpi=110)
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# 한 판
# ─────────────────────────────────────────────────────────────────────────────


def analyze(activation_paths, out_dir, name, side, *, kmax=K_MAX, criterion="dual",
            seed=0, min_peak=MIN_PEAK, model_path=None, muscles=None,
            plot_png=True, log=print):
    if muscles is None and model_path:
        muscles = muscle_names(model_path)
    A, names, time, dropped = load_activations(activation_paths, side,
                                               min_peak=min_peak, muscles=muscles,
                                               log=log)
    if A.shape[1] < 2:
        raise ValueError(f"분석할 근육이 부족하다 ({A.shape[1]}개)")
    log(f"  {side}: 근육 {A.shape[1]}개, 표본 {A.shape[0]}개")

    rows, fits = scan_k(A, kmax, seed=seed)
    chosen = {c: choose_k(rows, c) for c in ("global", "dual")}
    k = chosen.get(criterion)

    verdict = "PASS"
    if k is None:
        verdict = "WARN"
        fallback = chosen.get("global") or rows[-1][0]
        log(f"  [WARN] '{criterion}' 기준을 k={kmax}까지 채우지 못했다 — "
            f"k={fallback}로 그리되 PASS를 주지 않는다.")
        k = fallback

    C, W, total, per = fits[k]
    C, W, contrib = sort_synergies(C, W)

    log(render_table(rows, chosen, criterion))
    log(f"  선택 : global k={chosen['global']}  dual k={chosen['dual']}  "
        f"→ 기준 '{criterion}' 사용, k={k}, 전체 VAF {total*100:.1f}%")
    log(describe(W, names, contrib))
    log(f"  판정 : {verdict}")

    os.makedirs(out_dir, exist_ok=True)
    npz = os.path.join(out_dir, f"{name}_synergy_{side}.npz")
    np.savez(npz, W=W, C=C, time=time, muscles=np.array(names),
             vaf_total=total, vaf_muscle=per, contrib=contrib, k=k,
             k_global=chosen["global"] if chosen["global"] else -1,
             k_dual=chosen["dual"] if chosen["dual"] else -1,
             scan=np.array(rows), dropped=np.array(dropped))
    png = None
    if plot_png:
        png = os.path.join(out_dir, f"{name}_synergy_{side}.png")
        plot(png, C, W, names, time, contrib,
             title=f"{name} [{side}]  k={k} ({criterion})  VAF {total*100:.1f}%")

    return {
        "k": k, "k_global": chosen["global"], "k_dual": chosen["dual"],
        "vaf": round(total, 4), "vaf_muscle_min": round(float(per.min()), 4),
        "n_muscles": A.shape[1], "n_samples": A.shape[0],
        "verdict": verdict, "npz": npz, "png": png, "dropped": dropped,
    }


def main():
    ap = argparse.ArgumentParser(description="NMF 근육 시너지")
    ap.add_argument("--activation", nargs="+", required=True,
                    help="*_StaticOptimization_activation.sto (여러 개면 이어붙인다)")
    ap.add_argument("--model", default=None,
                    help="근육 목록을 뽑을 .osim — 팔 토크 액추에이터를 걸러내려면 필요하다")
    ap.add_argument("--out", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--side", default="both", choices=("r", "l", "both"))
    ap.add_argument("--kmax", type=int, default=K_MAX)
    ap.add_argument("--criterion", default="dual", choices=("dual", "global"))
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--min-peak", type=float, default=MIN_PEAK)
    args = ap.parse_args()

    if args.model is None:
        print("  [WARN] --model 이 없다 — 팔 토크 액추에이터가 근육으로 섞일 수 있다")
    sides = ("r", "l") if args.side == "both" else (args.side,)
    for side in sides:
        print(f"[{args.name}] side={side}")
        analyze(args.activation, args.out, args.name, side, kmax=args.kmax,
                criterion=args.criterion, seed=args.seed, min_peak=args.min_peak,
                model_path=args.model)


if __name__ == "__main__":
    main()
