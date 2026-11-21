"""정해진 걷기 동작에서 SMA 근육이 필요한 토크를 제때 내는가 — 근육만 빠르게 계산 (2026-11-21). conda env `fit` (numba).

입력: report/mujoco160/walk_table.npz (tools/mj_walk_sma.py --dump-table: 한 주기의 자세별 케이블 길이·모멘트암·ID 토크),
      report/mujoco160/sma_fit.json (tools/sma_fit.py 가 기계연 원자료에 맞춘 근육 값).
근육 모델: tools/sma_fast.py (가닥 모델). 등급 환산·묶음·다는 방식·가열 제어는 tools/mj_walk_sma.py 와 같다
  (등급 g = M3 의 힘·열용량·가열·냉각 × 규격 힘 비율, 길이 × 규격 수축 비율; 케이블이 가장 길 때 변형률 = 최대 수축; 듀티 = clip(KU × 힘 오차 ÷ 모듈 최대 힘, 0, 1)).

    FIT=/c/Users/USER/miniforge3/envs/fit/python.exe
    $FIT tools/walk_kin_fast.py                 # 모델(원래/1단계/2단계) × 걸음 주기(1·2·4·8배) 표와 그림
"""
import json, math, os, sys
os.environ.setdefault("MKL_THREADING_LAYER", "SEQUENTIAL")
import numpy as np                                           # noqa: E402
from scipy.optimize import lsq_linear                        # noqa: E402
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "tools"))
import sma_fast as SF                                        # noqa: E402
from sma_fast import njit, step_strands                      # noqa: E402
from actuator_model import BY_NAME                           # noqa: E402
OUT = os.path.join(ROOT, "report", "mujoco160"); W_ROWS = np.array([1.0, 1.0, 1.0, 2.0])


@njit(cache=True)
def kin_sim(pars, offs, nb, Pq, hq, honq, fmax1, Tmaxc, eps_tab, FD, lead, ku, dt, ncyc, ki=0.0):
    """케이블 nc 개를 ncyc 주기 돌린다. 마지막 주기의 힘 F, 평균 온도 T, 듀티 U (각 (n1, nc))."""
    n1, nc = eps_tab.shape; ns = offs.shape[1]; T = np.empty((nc, ns)); xi = np.ones((nc, ns)); ht = np.ones((nc, ns)); Ts = np.empty((nc, ns)); xs = np.ones((nc, ns))
    for i in range(nc):
        for j in range(ns): T[i, j] = pars[i, 11]; Ts[i, j] = pars[i, 11]
    nsub = int(math.ceil(dt / 0.001))                                                             # 열 식은 1 ms 이하로 잘라 푼다 (냉각이 센 모델은 시간 상수가 0.03 s)
    I = np.zeros(nc); F = np.zeros(nc); Fo = np.zeros((n1, nc)); To = np.zeros((n1, nc)); Uo = np.zeros((n1, nc))
    for k in range(ncyc * n1):
        kc = k % n1; kl = (kc + lead) % n1; last = k >= (ncyc - 1) * n1
        for i in range(nc):
            h = dt / nsub; um = 0.0; Tm = 0.0
            for _s in range(nsub):                                                                  # 제어와 열 식을 1 ms 이하 간격으로 (걸음을 늦춰도 제어 간격은 그대로)
                e = (FD[kl, i] / nb[i] - F[i] / nb[i]) / fmax1[i]; I[i] += ki * e * h                 # 적분 항 (ki = 0 이면 비례 제어만), 0~1 로 묶어 쌓임 방지
                if I[i] < 0.0: I[i] = 0.0
                if I[i] > 1.0: I[i] = 1.0
                u = ku * e + I[i]
                if u < 0.0: u = 0.0
                if u > 1.0: u = 1.0
                Tm = 0.0
                for j in range(ns): Tm += T[i, j]
                Tm /= ns; cap = (Tmaxc - Tm) / 10.0
                if cap < 0.0: cap = 0.0
                if cap < 1.0: u *= cap
                F[i] = nb[i] * step_strands(T[i], xi[i], ht[i], Ts[i], xs[i], offs[i], pars[i], u * Pq[i], h, eps_tab[kc, i], honq[i] if u > 0.0 else hq[i]); um += u / nsub
            if last: Fo[kc, i] = F[i]; To[kc, i] = Tm; Uo[kc, i] = um
    return Fo, To, Uo


def build(model, tab, speed=1.0, cool=1.0):
    """케이블별 값 묶음 + 원하는 힘 표 (배분 상한 = 그 자세에서 다 데웠을 때의 힘)."""
    nc = tab["TL"].shape[1]; ns = max(int(model.get("n", 1)), 1); pars = np.zeros((nc, 12)); offs = np.zeros((nc, ns)); nb = np.zeros(nc); Pq = np.zeros(nc); hq = np.zeros(nc); honq = np.zeros(nc); fmax1 = np.zeros(nc); L0 = np.zeros(nc)
    for i in range(nc):
        G = BY_NAME[str(tab["grade"][i])]; r = G.force_n[1] / model["FMAX"]; d = dict(model); d["C"] = model["C"] * r; d["FMAX"] = G.force_n[1]
        pars[i] = SF.par_array(d); offs[i] = SF.offsets(ns, model.get("sigma", 0.0)); nb[i] = float(tab["bundles"][i]); Pq[i] = model["P"] * r * speed; hq[i] = model["h"] * r * speed * cool; honq[i] = model["h_on"] * r * speed
        fmax1[i] = G.force_n[1]; L0[i] = model["L0"] * G.stroke_mm[1] / BY_NAME["M3"].stroke_mm[1]
    epsL = model["epsL"]; eps = np.clip(epsL - (tab["Lmax"][None, :] - tab["TL"]) / L0[None, :], 0.0, None); cap = fmax1 * nb; fhot = cap[None, :] * np.minimum(eps / epsL, 1.0)
    n1 = len(eps); FD = np.zeros((n1, nc)); idx = (tab["act_r"], tab["act_l"])
    for k in range(n1):
        for si in range(2):
            R = tab["RR"][k, si]; tau = tab["TD"][k, 4 * si:4 * si + 4]; Fm = np.maximum(fhot[k, idx[si]], 1e-3)
            A = np.vstack([W_ROWS[:, None] * R, np.sqrt(1e-4) * np.diag(1.0 / Fm)]); b = np.concatenate([W_ROWS * tau, np.zeros(8)]); FD[k, idx[si]] = lsq_linear(A, b, bounds=(np.zeros(8), Fm), method="bvls", tol=1e-8).x
    return dict(pars=pars, offs=offs, nb=nb, Pq=Pq, hq=hq, honq=honq, fmax1=fmax1, eps=eps, FD=FD, idx=idx)


def run(B, tab, slow=1.0, lead_s=0.0, ku=8.0, ncyc=4, Tmax=120.0, ki=0.0):
    dt = float(tab["dt"]) * slow; lead = int(round(lead_s / dt)); Fo, To, Uo = kin_sim(B["pars"], B["offs"], B["nb"], B["Pq"], B["hq"], B["honq"], B["fmax1"], Tmax, B["eps"], B["FD"], lead, ku, dt, ncyc, ki)
    ta = np.concatenate([np.einsum("kij,kj->ki", tab["RR"][:, si], Fo[:, B["idx"][si]]) for si in range(2)], axis=1); td = tab["TD"]; pk = np.abs(td).max(axis=0); e = ta - td
    tc = np.concatenate([np.einsum("kij,kj->ki", tab["RR"][:, si], B["FD"][:, B["idx"][si]]) for si in range(2)], axis=1)
    return dict(slow=slow, period_s=float(tab["T"]) * slow, lead_s=lead_s, err_rms_pct=(100 * np.sqrt(np.mean(e ** 2, axis=0)) / pk), err_p95_pct=(100 * np.percentile(np.abs(e), 95, axis=0) / pk),
                static_p95_pct=(100 * np.percentile(np.abs(tc - td), 95, axis=0) / pk), Tmin=To.min(axis=0), Tmax=To.max(axis=0), heat_W=float((Uo * B["Pq"][None, :] * B["nb"][None, :]).mean(axis=0).sum()), tau_act=ta, tau_des=td, T=To)


def best_lead(B, tab, slow, leads=(0.0, 0.05, 0.1, 0.15, 0.2, 0.3), Tmax=120.0):
    best = None
    for ld in leads:
        if ld > 0.3 * float(tab["T"]) * slow: continue
        for ku, ki in ((8.0, 0.0), (4.0, 5.0), (20.0, 5.0)):                                            # 제어 값 세 벌 중 가장 좋은 것 (비례만 / 비례 + 적분 두 가지)
            r = run(B, tab, slow, ld, ku=ku, ki=ki, Tmax=Tmax); r["ku"] = ku; r["ki"] = ki
            if best is None or r["err_rms_pct"].mean() < best["err_rms_pct"].mean(): best = r
    return best


def main():
    tab = dict(np.load(os.path.join(OUT, "walk_table.npz"), allow_pickle=True)); J = json.load(open(os.path.join(OUT, "sma_fit.json"), encoding="utf-8")) if os.path.exists(os.path.join(OUT, "sma_fit.json")) else {}
    models = [("원래 모델", dict(SF.BASE))] + [(lab, J[k]["model"]) for k, lab in (("fit1", "1단계 맞춤 (근육 하나 자료)"), ("fit2", "2단계 맞춤 (다섯 자료)"), ("fit3", "3단계 맞춤 (범위 넓힘)"), ("fit4", "4단계 맞춤 (As ≥ 40 °C)"), ("fit5", "5단계 맞춤 (실온에서 모두 풀림)")) if k in J]
    res = {}; keep = {}
    for name, model in models:
        B = build(model, tab); res[name] = []
        for slow in (1.0, 2.0, 4.0, 8.0):
            r0 = run(B, tab, slow, 0.0); r = best_lead(B, tab, slow); keep[(name, slow)] = r
            res[name].append(dict(period_s=r["period_s"], err_rms_nolead=float(r0["err_rms_pct"].mean()), lead_s=r["lead_s"], err_rms=float(r["err_rms_pct"].mean()), err_rms_min=float(r["err_rms_pct"].min()), err_rms_max=float(r["err_rms_pct"].max()),
                                  err_p95=float(r["err_p95_pct"].mean()), static_p95=float(r["static_p95_pct"].max()), Tmin=[float(r["Tmin"].min()), float(r["Tmin"].max())], Tmax=float(r["Tmax"].max()), heat_W=r["heat_W"]))
            x = res[name][-1]; print(f"{name:28s} 주기 {x['period_s']:4.1f} s | 토크 오차 {x['err_rms']:4.0f} % ({x['err_rms_min']:.0f}~{x['err_rms_max']:.0f}), 앞당김 없으면 {x['err_rms_nolead']:.0f} %, 앞당김 {x['lead_s']:.2f} s | 정적 한계 {x['static_p95']:.0f} % | 가장 식은 온도 {x['Tmin'][0]:.0f}~{x['Tmin'][1]:.0f} °C, 최고 {x['Tmax']:.0f} °C | 가열 {x['heat_W'] / 1000:.1f} kW", flush=True)
    json.dump(res, open(os.path.join(OUT, "walk_kin_fast.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for fnt in ("Malgun Gothic", "NanumGothic"):
        if any(x.name == fnt for x in font_manager.fontManager.ttflist): plt.rcParams["font.family"] = fnt; break
    plt.rcParams["axes.unicode_minus"] = False
    cols = ["#2980b9", "#e67e22", "#c0392b", "#27ae60", "#8e44ad", "#16a085"]; JN = [(0, "오른 고관절 굽힘"), (2, "오른 무릎"), (3, "오른 발목")]; slows = (1.0, 4.0); fig, ax = plt.subplots(3, len(slows) + 1, figsize=(7 * (len(slows) + 1), 11))
    for c, slow in enumerate(slows):
        for r_, (j, jn) in enumerate(JN):
            a = ax[r_, c]; x = np.linspace(0, 100, len(tab["TD"])); a.plot(x, tab["TD"][:, j], "k--", lw=1.5, label="걷기에 필요한 토크")
            for (name, _), col in zip(models, cols): rr = keep[(name, slow)]; a.plot(x, rr["tau_act"][:, j], color=col, lw=1.7, label=f"{name} (이 관절 오차 {rr['err_rms_pct'][j]:.0f} %)")
            a.grid(alpha=0.3); a.set_ylabel(f"{jn} [N·m]"); a.legend(fontsize=7.5)
            if r_ == 0: a.set_title(f"걸음 주기 {float(tab['T']) * slow:.1f} s", weight="bold")
            if r_ == 2: a.set_xlabel("걸음 주기 [%]")
    a = ax[0, -1]
    for (name, _), col in zip(models, cols): a.plot([x["period_s"] for x in res[name]], [x["err_rms"] for x in res[name]], "o-", color=col, lw=2, label=name)
    a.set_xscale("log"); a.set_xlabel("걸음 주기 [s]"); a.set_ylabel("토크 오차, 관절 8개 평균 [%]"); a.set_title("걸음이 느릴수록 근육이 따라간다", weight="bold"); a.grid(alpha=0.3, which="both"); a.legend(fontsize=8); a.set_ylim(0, None)
    a = ax[1, -1]
    for (name, _), col in zip(models, cols): a.plot([x["period_s"] for x in res[name]], [x["heat_W"] / 1000 for x in res[name]], "o-", color=col, lw=2, label=name)
    a.set_xscale("log"); a.set_xlabel("걸음 주기 [s]"); a.set_ylabel("평균 가열 전력, 케이블 16개 합 [kW]"); a.grid(alpha=0.3, which="both"); a.legend(fontsize=8); a.set_ylim(0, None); ax[2, -1].axis("off")
    fig.suptitle("160 cm 하체 설계 — 정해진 걷기 동작에서 SMA 근육이 낸 토크 (근육 모델별)", fontsize=13, weight="bold"); fig.tight_layout(rect=(0, 0, 1, 0.96)); fig.savefig(os.path.join(OUT, "fig_walk_kin_models.png"), dpi=110); print("그림 → report/mujoco160/fig_walk_kin_models.png")


if __name__ == "__main__":
    main()
