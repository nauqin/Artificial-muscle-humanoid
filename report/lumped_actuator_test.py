"""새 근육 k개(시너지당 하나, 토크 비율 고정)로 측정된 관절 모멘트를 재현할 수 있나 — 준정적 검사.

새 근육 i = 고정 토크 방향벡터 V_i ∈ R^4 (hip_flex, hip_add, knee, ankle) × 흥분 u_i(t) ≥ 0.
프레임마다  min ||Σ u_i V_i − τ_ID(t)||²  s.t. u ≥ 0  (NNLS) 를 풀고 남는 토크를 잰다.
남는 토크 = 보조 액추에이터가 떠맡아야 할 몫 = 하드웨어가 못 내는 부분.

    python report/lumped_actuator_test.py            # 7명 42 trial × 좌우, k=5 vs 8, in-sample + cross-trial
"""
import csv, glob, os, sys
import numpy as np
from scipy.optimize import nnls

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import synergy as sy, addb_io as io, run_muscle_analysis as rma  # noqa: E402

COORDS = ("hip_flexion", "hip_adduction", "knee_angle", "ankle_angle")
KEEP = ("subject3", "subject4", "subject5", "subject7", "subject8", "subject10", "subject11")
KS = (5, 8)
_mus_cache = {}


def load(sub, trial, side, min_peak=None):
    base = os.path.join(ROOT, "out/real", sub, trial)
    frc = io.read_mot(f"{base}/so/{trial}_StaticOptimization_force.sto")
    idm = io.read_mot(f"{base}/{trial}_id.sto")
    if sub not in _mus_cache:
        _mus_cache[sub] = sy.muscle_names(os.path.join(ROOT, "data", sub, f"{sub}_scaled.osim"))
    cs = [f"{c}_{side}" for c in COORDS]
    ma = rma.load_moment_arms(f"{base}/ma", trial, cs)
    kw = {} if min_peak is None else {"min_peak": min_peak}
    A, names, _, _ = sy.load_activations(f"{base}/so/{trial}_StaticOptimization_activation.sto",
                                         side, muscles=_mus_cache[sub], log=lambda *a: None, **kw)
    t = frc.time
    F = np.column_stack([frc.column(n) for n in names])
    R = np.stack([np.column_stack([np.interp(t, ma[c].time, ma[c].column(n)) for n in names]) for c in cs], axis=2)  # (T,M,4)
    tau = np.column_stack([np.interp(t, idm.time, idm.column(c + "_moment")) for c in cs])                          # (T,4)
    return A, F, R, tau


def load_pooled(sub, trials, side):
    """여러 trial 을 시간축으로 이어 붙인다 — 방향 V 를 걸음 하나가 아니라 피험자 전체에서 뽑을 때."""
    parts = [load(sub, t, side, min_peak=0.0) for t in trials]      # 근육 집합을 통일하려고 문턱 0
    return tuple(np.concatenate([p[i] for p in parts], axis=0) for i in range(4))


def directions(A, F, R, k):
    """시너지 k개의 고정 토크 방향벡터 V (k×4) — 근육 힘을 시너지에 귀속시킨 토크의 C-가중 평균."""
    C, W = sy.nmf(A, k, seed=0); C, W, _ = sy.sort_synergies(C, W)
    share = C[:, :, None] * W[None, :, :]; share /= np.maximum(share.sum(axis=1, keepdims=True), 1e-12)
    T = np.einsum("tkm,tm,tmc->ktc", share, F, R)                     # (k,T,4)
    V = np.stack([(T[i] * C[:, i:i + 1]).sum(axis=0) / max(C[:, i].sum(), 1e-9) for i in range(k)])
    return V


def track(V, tau):
    """프레임별 NNLS. 반환: 남는 토크 RMS/피크 (%) 좌표별, 최악 프레임 (%)."""
    resid = np.empty_like(tau)
    for i, row in enumerate(tau):
        u, _ = nnls(V.T, row); resid[i] = V.T @ u - row
    peak = np.abs(tau).max(axis=0)
    rms_pct = np.sqrt((resid ** 2).mean(axis=0)) / peak * 100
    worst_pct = np.abs(resid).max(axis=0) / peak * 100
    return rms_pct, worst_pct


def main():
    rows = []
    for sub in KEEP:
        trials = sorted(os.path.basename(p) for p in glob.glob(os.path.join(ROOT, "out/real", sub, "walking*")))
        cache = {}
        for trial in trials:
            for side in ("r", "l"):
                try:
                    A, F, R, tau = load(sub, trial, side)
                except Exception as exc:
                    print(f"  [skip] {sub}/{trial}/{side}: {exc}"); continue
                cache[(trial, side)] = (A, F, R, tau)
        first = {s: min(tr for (tr, sd) in cache if sd == s) for s in ("r", "l")}
        for (trial, side), (A, F, R, tau) in cache.items():
            for k in KS:
                V_in = directions(A, F, R, k)
                rms_in, worst_in = track(V_in, tau)
                A0, F0, R0, _ = cache[(first[side], side)]
                V_x = directions(A0, F0, R0, k)                        # 그 subject 첫 walking trial 로 만든 근육
                rms_x, worst_x = track(V_x, tau)
                rows.append(dict(subject=sub, trial=trial, side=side, k=k,
                                 **{f"in_rms_{c}": r for c, r in zip(COORDS, rms_in)},
                                 **{f"in_worst_{c}": r for c, r in zip(COORDS, worst_in)},
                                 **{f"x_rms_{c}": r for c, r in zip(COORDS, rms_x)},
                                 **{f"x_worst_{c}": r for c, r in zip(COORDS, worst_x)},
                                 cross_from=first[side]))
        print(f"  {sub}: {len(cache)} 분해 완료")

    out = os.path.join(ROOT, "report", "lumped_actuator_test.csv")
    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

    def summarize(prefix, label):
        print(f"\n=== {label}  (남는 토크 RMS / 피크 ID 모멘트, %. 42 trial × 좌우 평균) ===")
        print(f"  {'k':>2s}  " + "  ".join(f"{c:>13s}" for c in COORDS) + "   최악프레임 평균   >10% 인 분해 비율")
        for k in KS:
            sub_rows = [r for r in rows if r["k"] == k]
            rms = [np.mean([r[f"{prefix}_rms_{c}"] for r in sub_rows]) for c in COORDS]
            worst = np.mean([max(r[f"{prefix}_worst_{c}"] for c in COORDS) for r in sub_rows])
            bad = np.mean([max(r[f"{prefix}_rms_{c}"] for c in COORDS) > 10 for r in sub_rows]) * 100
            print(f"  {k:2d}  " + "  ".join(f"{v:12.1f}%" for v in rms) + f"   {worst:12.1f}%   {bad:8.0f}%")
    summarize("in", "in-sample: 그 trial 자신의 시너지로 만든 근육")
    summarize("x", "cross-trial: 첫 walking trial 로 만든 근육을 같은 사람의 다른 trial(walkingTS 포함)에")
    print(f"\n저장: {out}")


if __name__ == "__main__":
    main()
