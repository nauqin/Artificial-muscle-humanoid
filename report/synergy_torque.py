"""관절 토크 시너지 분석 → 관절 배분비 → 직접 최적화 설계와 비교 (2026-11-06).

상체 모델에는 팔·몸통 근육이 없어(뼈대의 근육 80개는 전부 다리) 근육 시너지를 뽑을 수 없다. 그래서 설계 입력과 같은
관절 토크에 NMF 를 건다. 토크는 부호가 있으므로 좌표마다 +쪽·−쪽 두 채널로 나눠 비음수로 만들고(τ⁺ = max(τ,0), τ⁻ = max(−τ,0)),
채널 피크로 정규화한 뒤(synergy.py 의 근육 정규화와 같은 이유) A ≈ C·W 로 분해한다.
시너지 i 의 관절 배분비 = (W_i⁺ − W_i⁻) × 채널 피크 를 단위벡터로 만든 것.

개수는 synergy.py 와 같은 VAF 기준으로 고른다: global = 전체 VAF ≥ 90 %, dual = 전체 ≥ 90 % 이고 모든 채널 ≥ 75 %.
비교는 design_directions.py 와 같은 함수·같은 평가 표본으로 한다 (총 토크 용량 LP, 전체 프레임 1/3 검산).
검증: --check-npz 로 준 직접 최적화 V 의 총 토크 용량이 그 설계의 csv 값과 같게 나와야 한다.

    python report/synergy_torque.py --space trunk --tasks report/upper/tasks_40kg_h160.csv --ref-mass 40 --ref-height 1.60 \\
        --check-npz report/newmus_V_design_trunk_k4_40kg_h160.npz --kmax 8
→ report/synergy_torque_<space>_<tag>.csv, report/newmus_V_design_<space>_k<k><tag>_syn.npz (path_design --design 입력과 같은 형식)
"""
import argparse, csv, os, sys, time
import numpy as np
from scipy.optimize import lsq_linear

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import design_directions as dd                         # noqa: E402
import design_space as ds                              # noqa: E402
import synergy as sy                                   # noqa: E402

N_FIT = 6000                                           # NMF 에 쓸 프레임 수 (곱셈 업데이트가 프레임 수에 비례해 느리다). 좌표별 극값은 항상 포함


def channels(tau):
    """(T, n) 부호 있는 토크 → (T, 2n) 비음수 채널과 채널 피크. 열 순서 = + 채널 n 개, − 채널 n 개."""
    A = np.hstack([np.maximum(tau, 0.0), np.maximum(-tau, 0.0)])
    peak = np.maximum(A.max(axis=0), 1e-9)
    return A / peak, peak


def signed_allocation(W, peak):
    """시너지 가중치 W (k, 2n) → 단위 관절 배분비 V (k, n). 채널 정규화를 되돌린 뒤 +쪽 − −쪽."""
    n = W.shape[1] // 2; Wn = W * peak
    return dd.unit(Wn[:, :n] - Wn[:, n:])


def worst_shortfall(V, c, tau_all):
    """토크 용량 c 로 전체 프레임의 1/3 을 낼 때 못 내는 토크 최대 (N·m). design_directions.verify_all 과 같은 계산."""
    worst = 0.0
    for row in tau_all[::3]:
        r = lsq_linear(V.T, row, bounds=(np.zeros(len(c)), np.maximum(c, 1e-9)), lsmr_tol="auto")
        worst = max(worst, float(np.abs(V.T @ r.x - row).max()))
    return worst


def evaluate(V, tau_eval, tau_all):
    """같은 평가 표본에서 총 토크 용량·slack·검산. LP 가 시간 안에 못 풀면 c=None."""
    total, c, slack = dd.capacity_lp(dd.unit(V), tau_eval, time_limit=120.0)
    if c is None: return dict(total=np.nan, cap=np.nan, slack=np.nan, worst=np.nan, c=None)
    return dict(total=total, cap=float(c.sum()), slack=slack, worst=worst_shortfall(dd.unit(V), c, tau_all), c=c)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--space", required=True); ap.add_argument("--tasks", default=None)
    ap.add_argument("--ref-mass", type=float, default=40.0); ap.add_argument("--ref-height", type=float, default=1.60)
    ap.add_argument("--kmax", type=int, default=8); ap.add_argument("--tag", default="_40kg_h160")
    ap.add_argument("--check-npz", default=None, help="직접 최적화 설계 npz (같은 평가로 다시 재서 비교)")
    ap.add_argument("--restarts", type=int, default=20)
    a = ap.parse_args()
    dd.REF = a.ref_mass * a.ref_height; dd.SP = ds.get(a.space); dd.COORDS = dd.SP["base"]; n = len(dd.COORDS)
    data = dd.dataset_from_tasks(os.path.join(ROOT, a.tasks)) if a.tasks else dd.dataset()
    tau_all = np.vstack(list(data.values())); tau_eval = dd.subsample(tau_all, 3000, seed=7); tau_fit = dd.subsample(tau_all, N_FIT, seed=11)
    print(f"[{a.space}] 좌표 {n}개 {list(dd.COORDS)}  프레임 {len(tau_all)} (과제 {[(k, len(v)) for k, v in data.items()]})  NMF 표본 {len(tau_fit)}  평가 표본 {len(tau_eval)}", flush=True)
    print("좌표별 피크 N·m:", np.round(np.abs(tau_all).max(axis=0), 1), flush=True)
    A, peak = channels(tau_fit)

    ax = evaluate(dd.axis_aligned(n), tau_eval, tau_all); print(f"축 정렬 {2*n}개: 총 토크 용량 {ax['total']:.1f}  slack {ax['slack']:.2f}  검산 {ax['worst']:.2f}", flush=True)
    direct = None
    if a.check_npz:
        Vd = np.load(os.path.join(ROOT, a.check_npz))["V"]; direct = evaluate(Vd, tau_eval, tau_all); direct["k"] = len(Vd)
        print(f"직접 최적화 {len(Vd)}개: 총 토크 용량 {direct['total']:.1f}  slack {direct['slack']:.2f}  검산 {direct['worst']:.2f}", flush=True)

    rows = []; print(f"\n{'k':>2s} {'전체VAF':>8s} {'최저채널VAF':>11s} {'토크용량합':>9s} {'LP목적':>9s} {'slack':>7s} {'못내는토크':>9s} {'s':>5s}", flush=True)
    for k in range(1, a.kmax + 1):
        tic = time.time(); C, W = sy.nmf(A, k, seed=k, n_restarts=a.restarts); total, per = sy.vaf(A, C @ W)
        V = signed_allocation(W, peak); ev = evaluate(V, tau_eval, tau_all)
        rows.append(dict(k=k, vaf=total, vaf_min=float(per.min()), cap=ev["cap"], lp=ev["total"], slack=ev["slack"], worst=ev["worst"]))
        print(f"{k:2d} {total*100:7.1f}% {per.min()*100:10.1f}% {ev['cap']:9.1f} {ev['total']:9.1f} {ev['slack']:7.2f} {ev['worst']:9.2f} {time.time()-tic:5.0f}", flush=True)
        if ev["c"] is not None:
            np.savez(os.path.join(ROOT, f"report/newmus_V_design_{a.space}_k{k}{a.tag}_syn.npz"), V=V, c=ev["c"], V_nm=V * ev["c"][:, None], W=W, peak=peak, coords=np.array(dd.COORDS), space=a.space, vaf=total)
    glob_k = next((r["k"] for r in rows if r["vaf"] >= sy.VAF_TOTAL), None); dual_k = next((r["k"] for r in rows if r["vaf"] >= sy.VAF_TOTAL and r["vaf_min"] >= sy.VAF_MUSCLE), None)
    print(f"\nVAF 로 고른 개수: global(전체 ≥ 90 %) {glob_k}, dual(+ 모든 채널 ≥ 75 %) {dual_k}.  당기기만 하는 케이블의 최소 개수 n+1 = {n+1}", flush=True)
    out = os.path.join(ROOT, f"report/synergy_torque_{a.space}{a.tag}.csv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(["design", "k", "vaf_total", "vaf_channel_min", "capacity_sum_Nm", "lp_objective", "slack_max_Nm", "shortfall_max_Nm"])
        w.writerow(["축 정렬", 2 * n, "", "", f"{ax['cap']:.1f}", f"{ax['total']:.1f}", f"{ax['slack']:.3f}", f"{ax['worst']:.3f}"])
        if direct: w.writerow(["직접 최적화", direct["k"], "", "", f"{direct['cap']:.1f}", f"{direct['total']:.1f}", f"{direct['slack']:.3f}", f"{direct['worst']:.3f}"])
        for r in rows: w.writerow(["토크 시너지", r["k"], f"{r['vaf']:.4f}", f"{r['vaf_min']:.4f}", f"{r['cap']:.1f}", f"{r['lp']:.1f}", f"{r['slack']:.3f}", f"{r['worst']:.3f}"])
    print("→", out, flush=True)


if __name__ == "__main__":
    main()
