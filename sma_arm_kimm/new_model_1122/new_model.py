"""SMA 근육 새 모델 (11/22, 김은서) — arm.py · sma_model.py 위에 새 모델 값만 바꿔 끼운 파일

같은 폴더에 arm.py · sma_model.py 만 있으면 돌아감 (기계연 원본 자료 필요 없음).
필요: pip install mujoco numpy

전 모델 (arm.py 그대로) 에서 바꾼 것 — 값의 근거는 노션 "SMA 근육 모델 — 새 모델 (11/22)"
  1. 변태 온도 구간 2.5배: 데울 때 55→72 → 42.25→84.75 °C, 식힐 때 52→38 → 62.5→27.5 °C (유효 구간, 11쪽 모양에 맞춤)
  2. 근육 최대 힘 450 → 89.3 N [가정: 5 · 6 kg 추세를 늘린 짐작]. 단면적 = 힘 / 400 MPa → 무게에 따라 줄어드는 온도가 오르는 정도
  3. 무거우면 덜 줄어듦: 최대 수축 0.42, 힘이 45.95 N 넘으면 1 N 당 0.0097 씩 줄어듦 (5 kg 39.0 % · 6 kg 29.5 % 를 지나는 선)
     0.42 = 3 · 4 kg 하한 (41.2 · 41.6 % 이상) 과 안 어긋나는 가장 작은 값
  4. 가열 · 강제 냉각 (근육 하나 11쪽에 맞춤): 267 W · 4.62 W/K → 344.5 W (@ 36 A) · 3.005 W/K
  5. 팔: 차가운 삼두가 버티기 시작하는 각도 40° · 팔 팬 냉각 ×1.4 (팔 굽힘으로 고름, 한 빠르기씩 빼도 40°)
  모두 맞춘 값 · 가정값이라 "맞춘 모델" 이지 검증된 모델은 아님

돌리기:  python new_model.py          → 팔 네 빠르기 + 무게별 + 근육 하나 (11쪽) 요약 (2코어 약 5분)
         python new_model.py --arm    → 팔만
쓰기:    import new_model as N;  r = N.NewArm().run(0.5);  r["bend"], r["log"]
"""
import sys, math, numpy as np, mujoco
import arm as A, sma_model as sm

# ── 새 모델 값 ──
FMAX = 48.0 + 0.40 / (0.095 / 9.81)          # 89.3 N [가정]
SLOPE = 0.095 / 9.81                          # 무거우면 덜 줄어드는 기울기 [1/N] (5 kg 39.0 % · 6 kg 29.5 %)
EPS = 0.42                                    # 최대 수축
F1 = 49.05 - (EPS - 0.390) / SLOPE            # 45.95 N 넘으면 덜 줄어듦
WIDTH = 2.5                                   # 변태 구간 배율
P36, HF = 344.5, 3.005                        # 가열 @ 36 A [W] · 강제 냉각 [W/K] (11쪽에 맞춤)
QT, CX = 40.0, 1.4                            # 삼두가 버티기 시작하는 각도 [°] · 팔 팬 냉각 배율

P = dict(sm.M3)
P["FMAX"] = FMAX; P["A_wire"] = FMAX / 400e6; P["K_mod"] = FMAX / P["epsL"]          # 강성은 0.40 기준 (만든 순서 그대로)
P["As"], P["Af"] = 63.5 - 8.5 * WIDTH, 63.5 + 8.5 * WIDTH; P["Ms"], P["Mf"] = 45 + 7 * WIDTH, 45 - 7 * WIDTH
P["epsL"] = EPS; A.EPSL = EPS                                                         # 팔에 근육 다는 기준도 0.42

def eff(F):
    """근육 힘 F [N] → 줄어들 수 있는 최대 비율"""
    return EPS if F <= F1 else max(EPS - SLOPE * (F - F1), 0.0)

class NewArm(A.Arm):
    """기계연 1자유도 팔 + 새 모델 근육 두 개"""
    def __init__(s, qt=QT, cx=CX, p36=P36, hf=HF):
        super().__init__(mode="ours", qt=qt, cx=cx); s.P36 = p36; s.HF = hf * cx; s.H_HALF = s.HF / 2

    def run(s, f, n=None):
        n = n or max(6, int(round(30 * f))); per = 1 / f; Pw = s.P36 * (A.CUR[f] / 36.0) ** 2
        mus = [sm.SMA(dict(P)), sm.SMA(dict(P))]; Fl = [0.0, 0.0]; d = mujoco.MjData(s.m); mujoco.mj_forward(s.m, d)
        t = 0.0; L = []; DT = s.m.opt.timestep; Fm = 0.0
        while t < n * per:
            cmd = s.phase((t % per) / per); eps = s.strains(d.qpos[0])
            for j in range(2):
                on, h = cmd[j]; mus[j].p["epsL"] = eff(Fl[j])
                F = mus[j].step(Pw if on else 0.0, DT, eps=eps[j], hS=h); Fl[j] = F; Fm = max(Fm, F); d.ctrl[j] = min(max(F, 0.0) / 1000.0, 1.0)
            mujoco.mj_step(s.m, d); t += DT; L.append((t, math.degrees(d.qpos[0]), mus[0].T, mus[1].T))
        L = np.array(L); r = []
        for k in range(n):
            sl = (L[:, 0] >= k * per) & (L[:, 0] < (k + 1) * per); r.append((L[sl, 1].max() - L[sl, 1].min(), L[sl, 1].min()))
        r = np.array(r)                                      # 처음 두 주기는 빼고 평균
        return dict(bend=r[2:, 0].mean(), lo=r[2:, 1].mean(), Tmax=L[:, 2:].max(), Fmax=Fm, log=L)

# ── 근육 하나 (기계연 11쪽: 추 5 kg, 0.5 Hz, 0.5 s 켜기, 25회) · 무게별 시험 ──
PEAK11 = np.array([39.4, 43.1, 46.4, 47.3, 46.7, 47.2, 47.2, 46.8, 46.2, 46.2, 46.0, 45.6, 45.2,
                   45.7, 45.2, 45.4, 44.5, 45.1, 45.3, 45.3, 45.2, 44.8, 44.6, 44.3, 44.3])   # 실측 회마다 최고 길이 [mm]
DT11, H_ON = 0.002, 0.72

def muscle(Pw, hf, W=5 * 9.81, f=0.5, t_on=0.5, n=25, cont=False, dur=None):
    """근육 하나에 추 W [N]. cont=True 면 계속 데움 (무게별 시험: 120 °C 닿으면 유지). 반환: (시간, 줄어든 길이 mm, 온도)"""
    st = sm.SMA(dict(P), latent=True); st.sig = W / P["A_wire"]; per = 1 / f; t = 0.0; out = []; hold = False
    for _ in range(int((dur if dur else n * per) / DT11)):
        if cont:
            hold = hold or st.T >= 120.0; st.step(min(Pw, H_ON * (120.5 - P["T_inf"])) if hold else Pw, DT11, hS=H_ON)
        else:
            on = (t % per) < t_on; st.step(Pw if on else 0.0, DT11, hS=H_ON if on else hf)
        t += DT11; out.append((t, eff(W) * 0.18 * 1000 * (1 - st.xi), st.T))
    return np.array(out)

def measure11(L, f=0.5, n=25):
    """회마다 최고 길이 · 끈 뒤 절반 · 80 % 풀리는 시간 평균"""
    t, d = L[:, 0], L[:, 1]; per = 1 / f; pk, r50, r20 = [], [], []
    for k in range(n):
        s = (t >= k * per) & (t < (k + 1) * per); ts, ds = t[s], d[s]; i = ds.argmax(); p = ds[i]; pk.append(p)
        for frac, lst in ((0.5, r50), (0.2, r20)):
            j = np.nonzero(ds[i:] <= frac * p)[0]; lst.append(ts[i + j[0]] - ts[i] if len(j) else per)
    return np.array(pk), float(np.mean(r50)), float(np.mean(r20))

def weight_test(kg):
    """무게별 시험: 계속 데워 6 s 뒤 줄어든 비율 [%]"""
    return float(muscle(P36, 1.0, W=kg * 9.81, cont=True, dur=6.0)[-500:, 1].mean() / 1.8)

def _arm(f):
    r = NewArm().run(f); return f, r["bend"], r["lo"], r["Tmax"], r["Fmax"]

if __name__ == "__main__":
    from multiprocessing import Pool
    print(f"새 모델: 최대 힘 {FMAX:.1f} N · 최대 수축 {EPS} ({F1:.2f} N 넘으면 덜 줄어듦) · 변태 구간 ×{WIDTH} "
          f"({P['As']:.2f}→{P['Af']:.2f} / {P['Ms']:.1f}→{P['Mf']:.1f} °C) · {P36} W · {HF} W/K · 삼두 {QT:.0f}° · 팔 냉각 ×{CX}", flush=True)
    with Pool(2) as pl: R = pl.map(_arm, list(A.MEAS))
    print("\n팔 (기계연 15쪽) | 굽힘 · 실측 | 쉴 때 | 최고 온도 | 근육 최대 힘")
    for f, b, lo, T, F in R: print(f"  {f:5.3f} Hz | {b:5.1f}° · {A.MEAS[f]}° | {lo:5.1f}° | {T:4.0f} °C | {F:3.0f} N", flush=True)
    print(f"  평균 오차 {A.err([x[1] for x in R]):.1f}°")
    if "--arm" in sys.argv: sys.exit()
    print("\n무게별 (실측 3 kg 41.2 % 이상 · 4 kg 41.6 % 이상 · 5 kg 39.0 % · 6 kg 29.5 %): "
          + " / ".join(f"{kg} kg {weight_test(kg):.1f} %" for kg in (3, 4, 5, 6)), flush=True)
    L = muscle(P36, HF); pk, r50, r20 = measure11(L)
    print(f"근육 하나 (11쪽): 1회 {pk[0]:.1f} · 4회 {pk[3]:.1f} · 25회 {pk[24]:.1f} mm (실측 39.4 · 47.3 · 44.3) | "
          f"끈 뒤 절반 풀림 {r50:.2f} s (0.38) · 80 % 풀림 {r20:.2f} s (0.71) | 최고 길이 오차 {np.sqrt(np.mean((pk - PEAK11) ** 2)):.2f} mm")
