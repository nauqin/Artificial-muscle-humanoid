"""실측 EMG(OpenCap LabValidation) 로 SO 활성도와 시너지를 검증한다.

    python report/emg_check.py            → report/emg_check.csv, fig_emg_subject3_walking1.png, fig_emg_summary.png

EMG 는 100 Hz 포락선(정규화 단위, 약간의 음수 = 기저선 뺀 것), 시간축은 AddB trial 과 같다(둘 다 0 부터).
비교는 SO 창(ID 창) 안에서, 둘 다 창 내 피크로 나눈 뒤: 지연 0 의 Pearson r, EMG 를 0~120 ms 앞당겼을 때
최선 r 과 그 지연(전기역학 지연이면 30~100 ms 가 정상). 시너지: 같은 근육 부분집합으로 EMG 와 SO 각각
NMF 해서 전체 VAF 90 % 에 드는 k 를 비교하고, 시간 곡선 C 를 최대 상관으로 짝지어 본다.
"""
import csv, glob, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
from scipy.optimize import linear_sum_assignment

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import addb_io as io, synergy as sy                      # noqa: E402
LAB = os.path.join(os.path.dirname(ROOT), "LabValidation_withoutVideos")
KEEP = ("subject3", "subject4", "subject5", "subject7", "subject8", "subject10", "subject11")
LAGS = np.arange(0.0, 0.121, 0.01)
for f in font_manager.findSystemFonts():
    if "AppleSDGothicNeo" in f: plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
plt.rcParams["axes.unicode_minus"] = False


def load_pair(sub, trial):
    """(t, {muscle: (emg, so)}) — SO 창 위에서, 둘 다 피크 정규화."""
    so = io.read_mot(os.path.join(ROOT, "out/real", sub, f"{trial}_segment_0/so/{trial}_segment_0_StaticOptimization_activation.sto"))
    emg = io.read_mot(os.path.join(LAB, sub, "EMGData", f"{trial}_EMG.sto"))
    t = so.time
    pairs = {}
    for n in emg.names:
        m = n.replace("_activation", "")
        if m not in so.names: continue
        e = np.clip(np.interp(t, emg.time, emg.column(n)), 0, None); s = so.column(m)
        if e.max() < 1e-6 or s.max() < 1e-6: continue
        pairs[m] = (e / e.max(), s / s.max(), emg.time, np.clip(emg.column(n), 0, None))
    return t, pairs


def corr(a, b):
    if a.std() < 1e-9 or b.std() < 1e-9: return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def best_lag(t, emg_t, emg_full, s):
    """EMG 를 lag 만큼 앞당겨(emg(t-lag) 를 so(t) 와) 최선 r."""
    best = (-2, 0.0)
    for lag in LAGS:
        e = np.interp(t - lag, emg_t, emg_full); e = e / max(e.max(), 1e-9)
        r = corr(e, s)
        if r > best[0]: best = (r, lag)
    return best


def k_at_vaf(A, thr=0.90, kmax=6):
    for k in range(1, kmax + 1):
        C, W = sy.nmf(A, k, seed=0, n_restarts=5, max_iter=1000)
        tot, _ = sy.vaf(A, C @ W)
        if tot >= thr: return k, C, W
    return kmax, C, W


def match_C(C1, C2):
    R = np.array([[corr(C1[:, i], C2[:, j]) for j in range(C2.shape[1])] for i in range(C1.shape[1])])
    R = np.nan_to_num(R, nan=-1)
    ri, ci = linear_sum_assignment(-R)
    return R[ri, ci]


def main():
    rows, syn_rows = [], []
    per_muscle = {}
    for sub in KEEP:
        trials = sorted(os.path.basename(p).replace("_segment_0", "") for p in glob.glob(os.path.join(ROOT, "out/real", sub, "walking*")))
        for trial in trials:
            try: t, pairs = load_pair(sub, trial)
            except Exception as exc: print(f"[skip] {sub}/{trial}: {exc}"); continue
            for m, (e, s, et, ef) in pairs.items():
                r0 = corr(e, s); rb, lag = best_lag(t, et, ef, s)
                rows.append((sub, trial, m, r0, rb, lag)); per_muscle.setdefault(m[:-2], []).append((r0, rb, lag))
            # 시너지: 같은 근육 부분집합, 좌우 따로
            for side in "rl":
                ms = [m for m in pairs if m.endswith("_" + side)]
                if len(ms) < 5: continue
                E = np.column_stack([pairs[m][0] for m in ms]); S = np.column_stack([pairs[m][1] for m in ms])
                ke, Ce, We = k_at_vaf(E); ks, Cs, Ws = k_at_vaf(S)
                # SO 전체 근육(파이프라인 k=5) 의 C 와 EMG 시너지 C 짝짓기
                A_full, names, _, _ = sy.load_activations(os.path.join(ROOT, "out/real", sub, f"{trial}_segment_0/so/{trial}_segment_0_StaticOptimization_activation.sto"), side, log=lambda *a: None)
                C5, W5 = sy.nmf(A_full, 5, seed=0, n_restarts=5, max_iter=1000)
                mr = match_C(Ce, C5)
                syn_rows.append((sub, trial, side, len(ms), ke, ks, float(np.mean(mr)), float(np.min(mr))))
            print(f"{sub} {trial}: 근육 {len(pairs)}개, r0 평균 {np.mean([corr(v[0], v[1]) for v in pairs.values()]):.2f}", flush=True)

    with open(os.path.join(ROOT, "report/emg_check.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["subject", "trial", "muscle", "r_lag0", "r_best", "lag_s"]); w.writerows(rows)
    with open(os.path.join(ROOT, "report/emg_check_synergy.csv"), "w", newline="") as fh:
        w = csv.writer(fh); w.writerow(["subject", "trial", "side", "n_muscles", "k_emg_90", "k_so_90_same_muscles", "C_match_mean_r", "C_match_min_r"]); w.writerows(syn_rows)

    print("\n근육별 (전 피험자·걸음 평균)   r(지연0)   r(최선)   지연 ms   n")
    order = sorted(per_muscle, key=lambda m: -np.nanmean([v[1] for v in per_muscle[m]]))
    for m in order:
        v = np.array(per_muscle[m]); print(f"  {m:8s} {np.nanmean(v[:,0]):8.2f} {np.nanmean(v[:,1]):8.2f} {np.nanmean(v[:,2])*1000:8.0f} {len(v):4d}")
    sr = np.array([(r[4], r[5], r[6], r[7]) for r in syn_rows])
    print(f"\n시너지 k (VAF 90 %, 같은 근육 부분집합): EMG 중앙값 {np.median(sr[:,0]):.0f} (범위 {sr[:,0].min():.0f}~{sr[:,0].max():.0f}), SO 중앙값 {np.median(sr[:,1]):.0f} (범위 {sr[:,1].min():.0f}~{sr[:,1].max():.0f}), n={len(sr)}")
    print(f"EMG 시너지 C ↔ SO(k=5) 시너지 C 짝지은 상관: 평균 {sr[:,2].mean():.2f}, 최저의 평균 {sr[:,3].mean():.2f}")

    # 그림 1: subject3 walking1 겹쳐 그리기
    t, pairs = load_pair("subject3", "walking1")
    ms = sorted(pairs, key=lambda m: (m[-1], m)); n = len(ms)
    fig, axes = plt.subplots(2, (n + 1) // 2, figsize=(2.6 * ((n + 1) // 2), 5.2), sharex=True, sharey=True)
    for ax, m in zip(axes.ravel(), ms):
        e, s, et, ef = pairs[m]; rb, lag = best_lag(t, et, ef, s)
        ax.plot(t, e, color="#1f77b4", lw=1.8, label="EMG"); ax.plot(t, s, color="#d62728", lw=1.4, label="SO")
        ax.set_title(f"{m}  r={corr(e, s):.2f} (지연 {lag*1000:.0f} ms: {rb:.2f})", fontsize=8)
    for ax in axes.ravel()[n:]: ax.axis("off")
    axes[0, 0].legend(fontsize=8); fig.suptitle("subject3 walking1 — 실측 EMG(파랑) vs SO 활성도(빨강), 각각 피크 정규화", fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(ROOT, "report/fig_emg_subject3_walking1.png"), dpi=150)

    # 그림 2: 근육별 r 분포
    fig, ax = plt.subplots(figsize=(10, 4))
    data = [[v[1] for v in per_muscle[m]] for m in order]
    ax.boxplot(data, tick_labels=order, showfliers=False); ax.axhline(0.5, color="0.5", ls="--", lw=0.8)
    ax.set_ylabel("EMG–SO 상관 r (최선 지연)"); ax.set_ylim(-0.6, 1.05)
    ax.set_title(f"7명 × 걸음 6개 — 근육별 EMG 대 SO 활성도 상관 (n={len(rows)})", fontsize=11)
    fig.tight_layout(); fig.savefig(os.path.join(ROOT, "report/fig_emg_summary.png"), dpi=150)
    print("그림 저장")


if __name__ == "__main__":
    main()
