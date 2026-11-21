"""가닥(분포) SMA 모델을 기계연 원자료에 맞추기 (2026-11-21). conda env `fit` (numba) 로 돌린다.

자료 (저장소에 올리지 않는 '기계연구원 실험 파일 공유 건/'):
  S  근육 하나 + 5 kg, 0.5 Hz, 켜는 0.5 s, 36 A 의 수축 변위 (8 kHz → 주기 평균 200 칸)
  A  팔 실험 0.125 / 0.25 / 0.5 / 1 Hz 의 굽힘각 (10 Hz → 주기 평균)
전류를 켠 시각이 기록에 없으므로, 파형은 가장 잘 맞는 자리로 돌려서(주기 안 이동) 비교한다. 오차 = 파형 RMS 차이 ÷ 실측 굽힘(변위) 크기.

    FIT=/c/Users/USER/miniforge3/envs/fit/python.exe
    $FIT tools/sma_fit.py verify            # 빠른 구현(가닥 1개) = 김은서 님 모델인지
    $FIT tools/sma_fit.py fit1              # 1단계: S 만으로 근육 값 맞춤 → 팔 실험은 검증
    $FIT tools/sma_fit.py fit2              # 2단계: S + A 다섯 개를 한꺼번에 (팔 팬 배율·감쇠 포함)
    $FIT tools/sma_fit.py report            # 그림·표 (원래 / 1단계 / 2단계)
결과: report/mujoco160/sma_fit.json, fig_sma_fit.png
"""
import json, os, sys, time
os.environ.setdefault("MKL_THREADING_LAYER", "SEQUENTIAL")
import numpy as np                                           # noqa: E402
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import kimm_data_check as KD                                 # noqa: E402
import sma_fast as SF                                        # noqa: E402
OUT = KD.OUT; JS = os.path.join(OUT, "sma_fit.json"); NS = 9; CUR = {0.125: 25.0, 0.25: 29.0, 0.5: 29.0, 1.0: 29.5}; FS = (0.125, 0.25, 0.5, 1.0)
ARM0 = dict(CX=0.6, DAMP=0.02, MF=0.2)                       # 김은서 님 새 설정 [가설]
# 맞추는 값: 이름, 아래, 위, 시작 (원래 모델 값)
P1 = [("P", 80, 900, 267), ("h", 1.0, 20, 4.62), ("As", 26, 68, 55), ("dA", 3, 45, 17), ("dH", 1, 45, 20), ("dM", 3, 45, 14), ("sigma", 0, 16, 4), ("c", 5, 300, 30), ("h_on", 0.0, 4, 0.72), ("rEM", 0.2, 0.9, 28 / 65)]
P2 = P1 + [("CX", 0.25, 1.5, 0.6), ("DAMP", 0.003, 1.0, 0.02), ("MF", 0.1, 0.5, 0.2)]
# 3단계: 1·2단계에서 범위 끝에 붙은 값(h, sigma, c, rEM, h_on, CX)의 범위를 넓힘. 4단계: 3단계 + 변태 시작 온도 As ≥ 40 °C (실온 근처 변태를 막은 것)
P3 = [("P", 80, 900, 267), ("h", 1.0, 80, 4.62), ("As", 26, 68, 55), ("dA", 3, 45, 17), ("dH", 1, 45, 20), ("dM", 3, 45, 14), ("sigma", 0, 30, 4), ("c", 5, 600, 30), ("h_on", 0.0, 8, 0.72), ("rEM", 0.1, 0.9, 28 / 65),
      ("CX", 0.08, 1.5, 0.6), ("DAMP", 0.003, 1.0, 0.02), ("MF", 0.1, 0.5, 0.2)]
QMAX = 1.6                                                    # 가닥 9개에서 가장 낮은(높은) 가닥의 온도 차이 = 1.59 sigma
P4 = [(n if n[0] != "As" else ("As", 40, 68, 55)) for n in P3]


def to_model(v):
    d = dict(SF.BASE); d.update(P=v["P"], h=v["h"], h_on=v["h_on"], As=v["As"], Af=v["As"] + v["dA"], Ms=v["As"] + v["dA"] - v["dH"], Mf=v["As"] + v["dA"] - v["dH"] - v["dM"], sigma=v["sigma"], rEM=v["rEM"], n=NS if v["sigma"] > 0 else 1)
    return d


_M = {}


def data():
    if not _M:
        d = KD.read_xlsx(os.path.join(KD.DD, "5kg_반복_Displacement_Time_Data.xlsx"), 2); w = KD.cycle_avg(d[:, 0], d[:, 1], 0.5, skip=3, nb=200); _M["S"] = w - w.min()
        for f in FS:
            a = KD.read_xlsx(os.path.join(KD.DD, KD.FILES[f])); t, y = a[:, 0], a[:, -1]; on = np.argmax(y > 1.0); m = t >= t[on]; nb = int(round(10 / f))
            w = KD.cycle_avg(t[m], y[m], f, nb=nb); _M[f] = (w - w.min(), float(KD.pp_per_cycle(t[m], y[m], f).mean()))
    return _M


def fold(Y, dt, per, ncyc_last, nb):
    """마지막 ncyc_last 주기를 nb 칸 평균 파형으로."""
    n1 = int(round(per / dt)); seg = Y[-ncyc_last * n1:].reshape(ncyc_last, n1).mean(axis=0); k = n1 // nb
    return seg[:nb * k].reshape(nb, k).mean(axis=1)


def best_shift(wm, wo):
    e = np.array([np.sqrt(np.mean((np.roll(wo, s) - wm) ** 2)) for s in range(len(wm))]); s = int(np.argmin(e)); return float(e[s]), s


def run_S(d, c):
    Y = SF.sim_single(SF.par_array(d), SF.offsets(d["n"], d["sigma"]), d["P"], d["h"], d["h_on"], c, 5.0, 0.5, 0.5, 6.0, 0.001, d["L0"]); w = fold(Y, 0.001, 2.0, 2, 200); return w - w.min()


def run_A(d, f, arm):
    P = d["P"] * (CUR[f] / 36.0) ** 2; cyc = max(5.0, round(8 * f)); nb = int(round(10 / f))
    Q = SF.sim_arm(SF.par_array(d), SF.offsets(d["n"], d["sigma"]), P, d["h"] * arm["CX"], d["h_on"], f, cyc, 0.002, d["L0"], 0.05, 0.5, 0.32, arm["MF"], arm["DAMP"])
    w = fold(Q, 0.002, 1.0 / f, 2, nb); n1 = int(round(1 / f / 0.002)); last = Q[-2 * n1:].reshape(2, n1); return w - w.min(), float((last.max(axis=1) - last.min(axis=1)).mean())


def evaluate(v, arm=None, full=False):
    M = data(); d = to_model(v); arm = dict(ARM0 if arm is None else arm); out = {}
    wS = run_S(d, v["c"]); eS, sS = best_shift(M["S"], wS); out["S"] = dict(err=eS / M["S"].max(), pp=float(wS.max()), pp_meas=float(M["S"].max()), wave=np.roll(wS, sS) if full else None)
    for f in FS:
        w, pp = run_A(d, f, arm); wm, pm = M[f]; e, s = best_shift(wm, w); out[f] = dict(err=e / pm, pp=pp, pp_meas=pm, wave=np.roll(w, s) if full else None)
    return out


def vec(names, z): return {n[0]: float(np.clip(x, n[1], n[2])) for n, x in zip(names, z)}


def fit(stage):
    from scipy.optimize import differential_evolution, minimize
    names = {1: P1, 2: P2, 3: P3, 4: P4, 5: P4}[stage]; lo = np.array([n[1] for n in names], float); hi = np.array([n[2] for n in names], float); best = [1e9, None]; t0 = time.time(); cnt = [0]
    prev = json.load(open(JS, encoding="utf-8")) if os.path.exists(JS) else {}

    def cost(z):
        v = vec(names, z)
        if v["dH"] > v["dA"] + 30 or v["As"] + v["dA"] - v["dH"] - v["dM"] < 20: return 5.0          # 말이 안 되는 온도 배치
        if stage == 5 and v["As"] + v["dA"] - v["dH"] - v["dM"] - QMAX * v["sigma"] < SF.BASE["T_inf"] + 2.0:    # 5단계: 가장 낮은 가닥도 실온(+2 °C)에서 다 풀려야 한다 (Mf 가 실온보다 높다)
            return 5.0 + 0.05 * (SF.BASE["T_inf"] + 2.0 - (v["As"] + v["dA"] - v["dH"] - v["dM"] - QMAX * v["sigma"]))
        try:
            if stage == 1:
                d = to_model(v); wS = run_S(d, v["c"]); M = data(); e = best_shift(M["S"], wS)[0] / M["S"].max()
            else:
                o = evaluate(v, arm=dict(CX=v["CX"], DAMP=v["DAMP"], MF=v["MF"])); errs = [o["S"]["err"]] + [o[f]["err"] for f in FS]; amp = [abs(o[k]["pp"] - o[k]["pp_meas"]) / o[k]["pp_meas"] for k in ("S",) + FS]
                e = float(np.mean(errs) + 0.5 * np.mean(amp))                                       # 파형 오차 평균 + 굽힘 크기 오차 평균의 절반
        except Exception: return 5.0
        cnt[0] += 1
        if e < best[0]:
            best[0], best[1] = e, v
            print(f"  [{cnt[0]:5d} · {time.time() - t0:5.0f} s] {e:.4f}  " + " ".join(f"{k} {x:.3g}" for k, x in v.items()), flush=True)
        return e
    x0 = np.array([n[3] for n in names], float)
    if stage == 2 and "fit1" in prev: x0 = np.array([prev["fit1"]["v"].get(n[0], n[3]) for n in names], float)
    cost(x0); rng = np.random.default_rng(0); init = lo + (hi - lo) * rng.random((int(os.environ.get("POP", 12)) * len(names), len(names))); init[0] = x0
    if stage >= 3:                                                                                  # 앞 단계 결과들을 처음 무리에 넣는다
        for i, k in enumerate(("fit1", "fit2", "fit3", "fit4")):
            if k in prev: init[i + 1] = np.clip([prev[k]["v"].get(n[0], n[3]) for n in names], lo, hi)
    if stage == 5:                                                                                  # 앞 결과에서 sigma 만 조건에 맞게 줄인 것도 넣는다
        js = [n[0] for n in names].index("sigma")
        for i in range(1, 5):
            z = vec(names, init[i]); init[i + 5] = init[i]; init[i + 5, js] = max(0.0, (z["As"] + z["dA"] - z["dH"] - z["dM"] - SF.BASE["T_inf"] - 3.0) / QMAX)
    r = differential_evolution(cost, list(zip(lo, hi)), init=init, maxiter=int(os.environ.get("GEN", 60)), tol=1e-4, seed=1, polish=False, mutation=(0.4, 0.9), recombination=0.8)
    r2 = minimize(cost, r.x, method="Nelder-Mead", options=dict(maxfev=1500, xatol=1e-3, fatol=1e-5)); v = best[1]
    arm = dict(CX=v.get("CX", ARM0["CX"]), DAMP=v.get("DAMP", ARM0["DAMP"]), MF=v.get("MF", ARM0["MF"])); o = evaluate(v, arm=arm)
    prev = json.load(open(JS, encoding="utf-8")) if os.path.exists(JS) else {}
    prev[f"fit{stage}"] = dict(v=v, model=to_model(v), arm=arm, cost=best[0], n_eval=cnt[0], n_strands=NS, result={str(k): dict(err=x["err"], pp=x["pp"], pp_meas=x["pp_meas"]) for k, x in o.items()})
    json.dump(prev, open(JS, "w", encoding="utf-8"), indent=1, ensure_ascii=False); show(f"{stage}단계 맞춤", o); print("저장 →", os.path.relpath(JS, KD.ROOT))


def show(label, o):
    print(f"{label}: 근육 하나 {o['S']['pp']:.1f} mm (실측 {o['S']['pp_meas']:.1f}, 파형 오차 {100 * o['S']['err']:.0f} %) | 팔 굽힘 " + " / ".join(f"{o[f]['pp']:.1f}" for f in FS) + "° (실측 " + " / ".join(f"{o[f]['pp_meas']:.1f}" for f in FS)
          + "°) 파형 오차 " + " / ".join(f"{100 * o[f]['err']:.0f}" for f in FS) + " %", flush=True)


def original_v():
    return {n[0]: float(n[3]) for n in P1} | dict(sigma=0.0)


def verify():
    import subprocess, types
    sm = types.ModuleType("sma_model"); exec(compile(subprocess.run(["git", "-C", KD.ROOT, "show", "origin/main:sma_arm_kimm/sma_model.py"], capture_output=True, check=True).stdout.decode("utf-8"), "sma_model.py", "exec"), sm.__dict__); mus = sm.SMA(); p = mus.p; KM = p["K_mod"] * p["EM"] / p["EA"]; eps0 = p["epsL"] + 49.05 / KM; y = v = 0.0; Y = []
    for k in range(6000):
        t = k * 0.001; ht = (t % 2) < 0.5; eps = max(eps0 - y / p["L0"], 0); F = mus.step(267.0 if ht else 0.0, 0.001, eps=eps, hS=0.72 if ht else 4.62); v += (F - 49.05 - 30 * v) / 5 * 0.001; y += v * 0.001; Y.append(y * 1000)
    d = to_model(original_v()); Z = SF.sim_single(SF.par_array(d), SF.offsets(1, 0.0), 267.0, 4.62, 0.72, 30.0, 5.0, 0.5, 0.5, 3.0, 0.001, 0.18); Y = np.array(Y)
    print(f"빠른 구현 vs 김은서 님 클래스 (근육 하나 시험 3주기): 변위 차이 최대 {np.abs(Y - Z).max():.4f} mm (변위 크기 {Y.max():.1f} mm), numba = {SF.HAVE_NUMBA}")
    show("원래 모델 (가닥 1개, 팔 새 설정)", evaluate(original_v()))


def report():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for fnt in ("Malgun Gothic", "NanumGothic"):
        if any(x.name == fnt for x in font_manager.fontManager.ttflist): plt.rcParams["font.family"] = fnt; break
    plt.rcParams["axes.unicode_minus"] = False
    J = json.load(open(JS, encoding="utf-8")); M = data(); cases = [("원래 모델", original_v(), ARM0, "#2980b9")]
    if "fit1" in J: cases.append(("1단계 (근육 하나 자료만으로 맞춤)", J["fit1"]["v"], ARM0, "#e67e22"))
    if "fit2" in J: cases.append(("2단계 (다섯 자료 모두로 맞춤)", J["fit2"]["v"], J["fit2"]["arm"], "#c0392b"))
    if "fit3" in J: cases.append(("3단계 (다섯 자료, 범위 넓힘)", J["fit3"]["v"], J["fit3"]["arm"], "#27ae60"))
    if "fit4" in J: cases.append(("4단계 (3단계 + As ≥ 40 °C)", J["fit4"]["v"], J["fit4"]["arm"], "#8e44ad"))
    if "fit5" in J: cases.append(("5단계 (4단계 + 실온에서 모두 풀림)", J["fit5"]["v"], J["fit5"]["arm"], "#16a085"))
    fig, ax = plt.subplots(1, 5, figsize=(26, 5)); summ = {}
    ax[0].plot(np.linspace(0, 2, 200, endpoint=False), M["S"], "k-", lw=2.8, label=f"기계연 실측 {M['S'].max():.1f} mm")
    for i, f in enumerate(FS): ax[i + 1].plot(np.linspace(0, 1 / f, len(M[f][0]), endpoint=False), M[f][0], "k-", lw=2.8, label=f"기계연 실측 {M[f][1]:.1f}°")
    for label, v, arm, col in cases:
        o = evaluate(v, arm=arm, full=True); show(label, o); summ[label] = {str(k): dict(err_pct=100 * x["err"], pp=x["pp"], pp_meas=x["pp_meas"]) for k, x in o.items()}
        ax[0].plot(np.linspace(0, 2, 200, endpoint=False), o["S"]["wave"], color=col, lw=1.7, label=f"{label} {o['S']['pp']:.1f} mm (오차 {100 * o['S']['err']:.0f} %)")
        for i, f in enumerate(FS): ax[i + 1].plot(np.linspace(0, 1 / f, len(M[f][0]), endpoint=False), o[f]["wave"], color=col, lw=1.7, label=f"{label.split(' (')[0]} {o[f]['pp']:.1f}° (오차 {100 * o[f]['err']:.0f} %)")
    ax[0].set_title("근육 하나 + 5 kg, 0.5 Hz, 36 A", weight="bold"); ax[0].set_ylabel("수축 변위 [mm]")
    for i, f in enumerate(FS): ax[i + 1].set_title(f"팔 실험 {f:g} Hz ({CUR[f]:g} A)", weight="bold"); ax[i + 1].set_ylabel("굽힘 [°]")
    for a in ax: a.set_xlabel("시간 [s] (한 주기)"); a.grid(alpha=0.3); a.legend(fontsize=7.5, loc="upper right")
    fig.suptitle("SMA 근육 모델과 기계연 원자료 — 오차 = 파형 RMS 차이 ÷ 실측 크기", fontsize=13, weight="bold"); fig.tight_layout(rect=(0, 0, 1, 0.95)); fig.savefig(os.path.join(OUT, "fig_sma_fit.png"), dpi=110)
    J["summary"] = summ; json.dump(J, open(JS, "w", encoding="utf-8"), indent=1, ensure_ascii=False); print("그림 → report/mujoco160/fig_sma_fit.png")


if __name__ == "__main__":
    a = sys.argv[1] if len(sys.argv) > 1 else "report"
    {"verify": verify, "fit1": lambda: fit(1), "fit2": lambda: fit(2), "fit3": lambda: fit(3), "fit4": lambda: fit(4), "fit5": lambda: fit(5), "report": report}[a]()
