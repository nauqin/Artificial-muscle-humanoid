"""서서히 시작하는 수축을 넣은 SMA 모델 — 여러 가닥(분포) 모델 (2026-11-21).

김은서 님 모델(sma_arm_kimm/sma_model.py, 근육 전체 = 온도 하나)은 55 °C 에 닿기 전에는 전혀 줄지 않다가 한꺼번에 준다.
기계연 원자료(근육 하나 5 kg)는 전류를 넣고 0.07 s 만에 움직이기 시작해 처음 0.25 s 동안 22 % 를 천천히 줄고 그 뒤에 빨라진다.
→ 근육을 가닥 N 개로 나누고 가닥마다 변태 온도를 조금씩 다르게 둔다 (정규분포, 평균 이동 mu, 표준편차 sigma °C).
  가닥 하나는 김은서 님 SMA 클래스 그대로 (열용량·가열·냉각·힘을 1/N 씩). 변형률은 모든 가닥이 같고 힘은 더한다.
  sigma = 0, mu = 0 이면 원래 모델과 같다.

맞추는 자료: 근육 하나 + 5 kg, 0.5 Hz, 켜는 0.5 s, 36 A 의 변위 파형 (8 kHz, 주기 평균). 팔 실험 4개는 맞추는 데 쓰지 않고 검증에만 쓴다.
맞추는 값: 가열 P [W], 강제 냉각 h [W/K], sigma, mu, 추 감쇠 c [N·s/m]. 켜는 동안 냉각 0.72 W/K 는 그대로 [가정].

    MJ=/c/Users/USER/miniforge3/envs/mj/python.exe
    $MJ tools/sma_dist.py fit            # 맞추기 → report/mujoco160/sma_dist_fit.json
    $MJ tools/sma_dist.py check          # 근육 하나 파형 + 팔 실험 4개 비교 그림·표 (원래 모델 vs 새 모델)
"""
import json, math, os, sys
os.environ.setdefault("MKL_THREADING_LAYER", "SEQUENTIAL")
import numpy as np                                           # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kimm_data_check as KD                                 # noqa: E402
ROOT, OUT = KD.ROOT, KD.OUT; FIT = os.path.join(OUT, "sma_dist_fit.json"); N_STRANDS = 7
_SM = None


def sm_mod():
    global _SM
    if _SM is None: _SM = KD.load_model()[0]
    return _SM


def _norm_quantiles(n):
    """정규분포를 같은 확률 n 칸으로 나눈 가운데 값 (표준편차 1)."""
    from statistics import NormalDist
    return np.array([NormalDist().inv_cdf((i + 0.5) / n) for i in range(n)])


class SMAStrands:
    """가닥 N 개 묶음. 김은서 님 SMA 와 같은 쓰임 (step → 힘, .T, .xi, u_cap, power_duty, .p)."""

    def __init__(self, p=None, n=N_STRANDS, sigma=0.0, mu=0.0, **kw):
        sm = sm_mod(); p = dict(sm.M3 if p is None else p); self.p = p; self.n = n; self.strands = []
        for dq in (_norm_quantiles(n) if n > 1 else np.zeros(1)):
            q = dict(p); d = mu + sigma * dq
            q.update(m=p["m"] / n, hS=p["hS"] / n, Qon=p["Qon"] / n, FMAX=p["FMAX"] / n, As=p["As"] + d, Af=p["Af"] + d, Ms=p["Ms"] + d, Mf=p["Mf"] + d)
            q["A_wire"] = q["FMAX"] / 400e6; q["K_mod"] = q["FMAX"] / q["epsL"]; self.strands.append(sm.SMA(q, **kw))
        self.hS = p["hS"]

    def step(self, P, dt, eps=None, hS=None):
        h = (self.hS if hS is None else hS) / self.n; F = 0.0
        for s in self.strands:
            f = s.step(P / self.n, dt, eps=eps, hS=h); F += f if f is not None else 0.0
        return F if eps is not None else None

    def force(self, eps): return sum(s.force(eps) for s in self.strands)
    def power_duty(self, u): return u * self.p["Qon"]
    def u_cap(self, margin=10.0): return min(1.0, max(0.0, (self.p["Tmax"] - self.T) / margin))
    T = property(lambda s: float(np.mean([x.T for x in s.strands])))
    xi = property(lambda s: float(np.mean([x.xi for x in s.strands])))


def make(par=None, p=None, n=N_STRANDS):
    """맞춘 값(dict: P, h, sigma, mu)으로 근육 하나. par 없으면 저장된 맞춤, 그것도 없으면 원래 모델."""
    sm = sm_mod(); base = dict(sm.M3 if p is None else p)
    if par is None and os.path.exists(FIT): par = json.load(open(FIT, encoding="utf-8"))["par"]
    if par is None: return sm.SMA(base)
    r = base["FMAX"] / sm.M3["FMAX"]; base["Qon"] = par["P"] * r; base["hS"] = par["h"] * r
    return SMAStrands(base, n=n, sigma=par["sigma"], mu=par["mu"])


# ───────────────────────── 근육 하나 + 추 (기계연 시험) ─────────────────────────
def sim_single(mus, c=30.0, mass=5.0, f=0.5, on=0.5, cycles=6, dt=0.001, h_on=0.72, P=None, h=None):
    """근육 하나에 추를 달고 주기 1/f 중 on 초 동안 데움. 반환 (t, 수축 변위 mm). 0 = 다 식은 평형."""
    p = mus.p; L0 = p["L0"]; g = 9.81; P = p["Qon"] if P is None else P; h = p["hS"] if h is None else h
    KM = p["K_mod"] * p["EM"] / p["EA"]; eps0 = p["epsL"] + mass * g / KM; y = v = 0.0; n = int(round(cycles / f / dt)); Y = np.zeros(n); per = 1.0 / f
    for k in range(n):
        t = k * dt; heat = (t % per) < on; eps = max(eps0 - y / L0, 0.0)
        F = mus.step(P if heat else 0.0, dt, eps=eps, hS=h_on if heat else h)
        v += (F - mass * g - c * v) / mass * dt; y += v * dt; Y[k] = y * 1000.0
    return np.arange(n) * dt, Y


def wave_err(wm, wo):
    """두 주기 파형(같은 칸 수)의 RMS 차이 — 시작 시점을 모르므로 가장 잘 맞는 자리로 돌려서."""
    best = (1e9, 0)
    for s in range(len(wm)):
        e = float(np.sqrt(np.mean((np.roll(wo, s) - wm) ** 2)))
        if e < best[0]: best = (e, s)
    return best


_MEAS = {}


def meas_single(nb=200):
    if "s" not in _MEAS:
        d = KD.read_xlsx(os.path.join(KD.DD, "5kg_반복_Displacement_Time_Data.xlsx"), 2); w = KD.cycle_avg(d[:, 0], d[:, 1], 0.5, skip=3, nb=nb); _MEAS["s"] = w - w.min()
    return _MEAS["s"]


def model_wave_single(par, nb=200, n=N_STRANDS):
    mus = make(par, n=n) if par else sm_mod().SMA()
    t, Y = sim_single(mus, c=(par or {}).get("c", 30.0)); m = t >= 4 * 2.0; w = KD.cycle_avg(t[m], Y[m], 0.5, skip=0, nb=nb)
    return w - w.min()


def fit():
    from scipy.optimize import minimize
    wm = meas_single(); names = ("P", "h", "sigma", "mu", "c"); x0 = np.array([267.0, 4.62, 6.0, -3.0, 30.0]); sc = np.array([267.0, 4.62, 6.0, 6.0, 30.0]); hist = []

    def cost(z):
        par = dict(P=abs(z[0]) * sc[0], h=abs(z[1]) * sc[1], sigma=abs(z[2]) * sc[2], mu=z[3] * sc[3], c=abs(z[4]) * sc[4])
        e, s = wave_err(wm, model_wave_single(par)); hist.append(e)
        print(f"  {len(hist):3d}: P {par['P']:6.1f} W  h {par['h']:5.2f} W/K  sigma {par['sigma']:5.2f}  mu {par['mu']:+6.2f} °C  c {par['c']:5.1f}  → 파형 오차 {e:.2f} mm", flush=True)
        return e
    e0, _ = wave_err(wm, model_wave_single(None)); print(f"원래 모델 (가닥 1개): 파형 오차 {e0:.2f} mm (실측 변위 {wm.max():.1f} mm)", flush=True)
    r = minimize(cost, x0 / sc, method="Nelder-Mead", options=dict(maxfev=int(os.environ.get("MAXFEV", 140)), xatol=5e-3, fatol=5e-3, initial_simplex=None))
    z = r.x; par = dict(P=abs(z[0]) * sc[0], h=abs(z[1]) * sc[1], sigma=abs(z[2]) * sc[2], mu=z[3] * sc[3], c=abs(z[4]) * sc[4])
    json.dump(dict(par=par, n_strands=N_STRANDS, wave_err_mm=float(r.fun), wave_err_mm_original=e0, measured_pp_mm=float(wm.max()), h_on=0.72,
                   note="근육 하나 5 kg·0.5 Hz·36 A 변위 파형에 맞춤. 팔 실험은 쓰지 않음."), open(FIT, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print("맞춤 →", os.path.relpath(FIT, ROOT), json.dumps(par), f"오차 {r.fun:.2f} mm (원래 {e0:.2f})")


# ───────────────────────── 팔 실험 (검증) — 김은서 kimm_arm_new.py 와 같은 설정, 근육만 바꿈 ─────────────────────────
def sim_arm(f, par, CX=0.6, MF=0.2, DAMP=0.02, dt=0.002, cycles=None):
    sm = sm_mod(); R, LOAD, ARM, g = 0.05, 0.5, 0.32, 9.81; L0, EPSL = sm.M3["L0"], sm.M3["epsL"]; QMAX = L0 * EPSL / R
    CUR = {0.125: 25.0, 0.25: 29.0, 0.5: 29.0, 1.0: 29.5}; P36 = (par or {}).get("P", sm.KIMM_HEAT["P"]); HF = (par or {}).get("h", sm.KIMM_HEAT["h_forced"]) * CX; H_ON, H_HALF = 0.72, HF / 2
    I = MF * ARM ** 2 / 3 + LOAD * ARM ** 2; Gm = (MF * ARM / 2 + LOAD * ARM) * g; P = P36 * (CUR[f] / 36.0) ** 2; per = 1 / f; n = cycles or max(6, int(round(12 * f)))
    mus = [make(par) if par else sm.SMA(), make(par) if par else sm.SMA()]; q = w = 0.0; N = int(round(n * per / dt)); Q = np.zeros(N)
    for k in range(N):
        t = k * dt; x = (t % per) / per
        b = (True, H_ON) if x < .25 else (False, H_HALF) if x < .5 else (False, H_HALF if x >= .75 else HF)
        tr = (False, H_HALF if .25 <= x < .5 else HF) if x < .5 else (True, H_ON) if x < .75 else (False, H_HALF)
        eps = (max(EPSL - R * q / L0, 0.0), EPSL + R * q / L0); F = [mus[j].step(P if c[0] else 0.0, dt, eps=eps[j], hS=c[1]) for j, c in enumerate((b, tr))]
        w += (R * (F[0] - F[1]) - Gm * math.sin(q) - DAMP * w) / I * dt; q += w * dt
        if q < 0: q, w = 0.0, max(w, 0.0)
        if q > QMAX: q, w = QMAX, min(w, 0.0)
        Q[k] = math.degrees(q)
    return np.arange(N) * dt, Q


def check():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for fnt in ("Malgun Gothic", "NanumGothic"):
        if any(x.name == fnt for x in font_manager.fontManager.ttflist): plt.rcParams["font.family"] = fnt; break
    plt.rcParams["axes.unicode_minus"] = False
    J = json.load(open(FIT, encoding="utf-8")); par = J["par"]; res = dict(par=par); fig, ax = plt.subplots(1, 5, figsize=(25, 4.8))
    wm = meas_single(); w0 = model_wave_single(None); w1 = model_wave_single(par); e0, s0 = wave_err(wm, w0); e1, s1 = wave_err(wm, w1); x = np.linspace(0, 2, len(wm), endpoint=False)
    ax[0].plot(x, wm, "k-", lw=2.5, label="기계연 실측"); ax[0].plot(x, np.roll(w0, s0), color="#2980b9", lw=1.6, label=f"원래 모델 (오차 {e0:.1f} mm)"); ax[0].plot(x, np.roll(w1, s1), color="#c0392b", lw=1.8, label=f"가닥 모델 (오차 {e1:.1f} mm)")
    ax[0].set_title("근육 하나 + 5 kg, 0.5 Hz (맞춘 자료)", weight="bold"); ax[0].set_xlabel("시간 [s]"); ax[0].set_ylabel("수축 변위 [mm]"); ax[0].grid(alpha=0.3); ax[0].legend(fontsize=8)
    res["single"] = dict(err_mm_original=e0, err_mm_new=e1, pp_meas=float(wm.max()), pp_original=float(w0.max()), pp_new=float(w1.max()))
    print(f"근육 하나: 실측 {wm.max():.1f} mm | 원래 {w0.max():.1f} mm (파형 오차 {e0:.2f}) | 가닥 {w1.max():.1f} mm (파형 오차 {e1:.2f})", flush=True); rows = {}
    for i, (f, fn) in enumerate(KD.FILES.items()):
        d = KD.read_xlsx(os.path.join(KD.DD, fn)); t, y = d[:, 0], d[:, -1]; on = np.argmax(y > 1.0); m = t >= t[on]; nb = int(min(100, 10 / f))
        wmm = KD.cycle_avg(t[m], y[m], f, nb=nb); pm = float(KD.pp_per_cycle(t[m], y[m], f).mean()); out = {}
        for key, pr, col in (("original", None, "#2980b9"), ("new", par, "#c0392b")):
            tt, Q = sim_arm(f, pr); mm = tt >= 2 / f; wo = KD.cycle_avg(tt[mm], Q[mm], f, skip=0, nb=nb); e, s = wave_err(wmm - wmm.min(), wo - wo.min()); pp = float(KD.pp_per_cycle(tt[mm], Q[mm], f, skip=0).mean())
            out[key] = dict(bend=pp, wave_err_deg=e); ax[i + 1].plot(np.linspace(0, 1 / f, nb, endpoint=False), np.roll(wo - wo.min(), s), color=col, lw=1.7, label=f"{'원래' if key == 'original' else '가닥'} 모델 {pp:.1f}° (파형 오차 {e:.1f}°)")
        ax[i + 1].plot(np.linspace(0, 1 / f, nb, endpoint=False), wmm - wmm.min(), "k-", lw=2.5, label=f"기계연 실측 {pm:.1f}°"); ax[i + 1].set_title(f"팔 실험 {f:g} Hz (검증 — 맞추는 데 안 씀)", weight="bold"); ax[i + 1].set_xlabel("시간 [s]"); ax[i + 1].set_ylabel("굽힘 [°]"); ax[i + 1].grid(alpha=0.3); ax[i + 1].legend(fontsize=8)
        rows[f] = dict(meas=pm, **out); print(f"팔 {f:5.3f} Hz: 실측 {pm:5.1f}° | 원래 {out['original']['bend']:5.1f}° (파형 오차 {out['original']['wave_err_deg']:.1f}°) | 가닥 {out['new']['bend']:5.1f}° (파형 오차 {out['new']['wave_err_deg']:.1f}°)", flush=True)
    res["arm"] = {str(k): v for k, v in rows.items()}
    for key in ("original", "new"):
        res[f"arm_bend_rmse_{key}"] = float(np.sqrt(np.mean([(rows[f][key]["bend"] - rows[f]["meas"]) ** 2 for f in rows]))); res[f"arm_wave_err_mean_{key}"] = float(np.mean([rows[f][key]["wave_err_deg"] for f in rows]))
    print(f"팔 굽힘 오차 RMS: 원래 {res['arm_bend_rmse_original']:.1f}° → 가닥 {res['arm_bend_rmse_new']:.1f}° | 파형 오차 평균: 원래 {res['arm_wave_err_mean_original']:.1f}° → 가닥 {res['arm_wave_err_mean_new']:.1f}°")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_sma_dist_check.png"), dpi=110); json.dump(res, open(os.path.join(OUT, "sma_dist_check.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)


if __name__ == "__main__":
    {"fit": fit, "check": check}[sys.argv[1] if len(sys.argv) > 1 else "check"]()
