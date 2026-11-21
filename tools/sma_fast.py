"""SMA 근육 열 모델 — 가닥(분포) 모델의 빠른 구현 (2026-11-21).

식은 김은서 님 sma_arm_kimm/sma_model.py (v3) 와 같다:
  열  C dT/dt = P − hS (T − T∞) + 잠열,  상태 ξ = Brinson cos 식 (방향이 바뀐 순간의 ξ 에서 이어감),
  문턱 이동 As′ = As + σ/C_sig,  힘 F = K(ξ)(ε − ε_L ξ)  (K(ξ) = K_mod · (ξ·E_M + (1−ξ)·E_A)/E_A),  온도와 ξ 는 이분법으로 같이 푼다.
달라진 것: 근육을 가닥 n 개로 나누고 가닥마다 변태 온도 네 개를 off_i [°C] 만큼 옮긴다 (정규분포 분위수 × sigma).
  가닥은 열용량·가열·냉각·힘을 1/n 씩 나눠 가지므로 가닥 하나의 식은 근육 전체 값으로 쓴 식과 같다 (off 만 다르다).
  n = 1 (sigma = 0) 이면 원래 모델과 같다 — tools/sma_fit.py verify 로 확인.
numba 가 있으면 컴파일해서 쓰고 (맞춤용 conda env `fit`), 없으면 같은 코드를 파이썬으로 돈다 (mujoco 가 있는 env `mj`).

par (길이 12): [C 열용량 J/K, k_lat 잠열/비열 K, As, Af, Ms, Mf, C_sig Pa/K, A 단면적 m², K_mod N, rEM = E_M/E_A, epsL, T_inf]
"""
import math
import numpy as np
try:
    from numba import njit
    HAVE_NUMBA = True
except ImportError:                                           # 파이썬으로 그대로
    HAVE_NUMBA = False
    def njit(*a, **k):
        return a[0] if a and callable(a[0]) else (lambda f: f)

BASE = dict(C=1.8, k_lat=24180.0 / 600.0, As=55.0, Af=72.0, Ms=52.0, Mf=38.0, C_sig=8.0e6, FMAX=450.0, rEM=28.0 / 65.0, epsL=0.40, T_inf=25.0,
            P=267.0, h=4.62, h_on=0.72, L0=0.18, Tmax=120.0, sigma=0.0, n=1)          # 김은서 님 M3 값 (기계연 실측 기반 + 교수님 T3.1)


def par_array(d):
    return np.array([d["C"], d["k_lat"], d["As"], d["Af"], d["Ms"], d["Mf"], d["C_sig"], d["FMAX"] / 400e6, d["FMAX"] / d["epsL"], d["rEM"], d["epsL"], d["T_inf"]], dtype=np.float64)


def offsets(n, sigma):
    from statistics import NormalDist
    return np.array([sigma * NormalDist().inv_cdf((i + 0.5) / n) for i in range(n)], dtype=np.float64) if n > 1 else np.zeros(1)


@njit(cache=True)
def _force1(xi, eps, par):
    f = par[8] * (xi * par[9] + (1.0 - xi)) * (eps - par[10] * xi)
    return f if f > 0.0 else 0.0


@njit(cache=True)
def _xi_of(TT, heating, x, xi_s, As, Af, Ms, Mf):
    if heating:
        f = 1.0 if TT <= As else (0.0 if TT >= Af else 0.5 * (math.cos(math.pi * (TT - As) / (Af - As)) + 1.0))
        xb = xi_s * f
        return x if xb > x else xb
    f = 0.0 if TT >= Ms else (1.0 if TT <= Mf else 0.5 * math.cos(math.pi * (TT - Mf) / (Ms - Mf)) + 0.5)
    xb = xi_s + (1.0 - xi_s) * f
    return x if xb < x else xb


@njit(cache=True)
def step_strands(T, xi, heat, Ts, xis, off, par, P, dt, eps, hS):
    """가닥 상태(T, xi, heat, Ts, xis: 길이 n 배열)를 dt 만큼 진행하고 근육 전체 힘 [N] 을 돌려준다. P·hS 는 근육 전체 값."""
    n = len(T); C = par[0]; k = par[1]; Ftot = 0.0
    for i in range(n):
        x = xi[i]; Ti = T[i]; sh = _force1(x, eps, par) / par[7] / par[6]
        As = par[2] + off[i] + sh; Af = par[3] + off[i] + sh; Ms = par[4] + off[i] + sh; Mf = par[5] + off[i] + sh
        q = P - hS * (Ti - par[11]); hn = q > 0.0
        if hn != (heat[i] > 0.5):
            heat[i] = 1.0 if hn else 0.0; Ts[i] = Ti; xis[i] = x
        Tp = Ti + dt * q / C; lo = Tp - k - 1.0; hi = Tp + k + 1.0
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if mid - Tp - k * (_xi_of(mid, hn, x, xis[i], As, Af, Ms, Mf) - x) > 0.0: hi = mid
            else: lo = mid
        Tn = 0.5 * (lo + hi); xn = _xi_of(Tn, hn, x, xis[i], As, Af, Ms, Mf); T[i] = Tn; xi[i] = xn
        Ftot += _force1(xn, eps, par)
    return Ftot / n


@njit(cache=True)
def sim_single(par, off, P, h, h_on, c, mass, f, on, cycles, dt, L0):
    """근육 하나 + 추 (기계연 근육 하나 시험). 반환: 수축 변위 [mm] 배열."""
    n = len(off); T = np.full(n, par[11]); xi = np.ones(n); heat = np.ones(n); Ts = T.copy(); xis = xi.copy()
    g = 9.81; KM = par[8] * par[9]; eps0 = par[10] + mass * g / KM; y = 0.0; v = 0.0; N = int(round(cycles / f / dt)); Y = np.zeros(N); per = 1.0 / f
    for k in range(N):
        t = k * dt; ph = t - math.floor(t / per) * per; ht = ph < on; eps = eps0 - y / L0
        if eps < 0.0: eps = 0.0
        F = step_strands(T, xi, heat, Ts, xis, off, par, P if ht else 0.0, dt, eps, h_on if ht else h)
        v += (F - mass * g - c * v) / mass * dt; y += v * dt; Y[k] = y * 1000.0
    return Y


@njit(cache=True)
def sim_arm(par, off, P, HF, h_on, f, cycles, dt, L0, R, LOAD, ARM, MF, DAMP):
    """기계연 1자유도 팔 (이두·삼두, 블레이드 냉각 시간표) — 김은서 님 kimm_arm_new.py 와 같은 설정. 반환: 굽힘각 [°] 배열."""
    n = len(off); g = 9.81; EPSL = par[10]; QMAX = L0 * EPSL / R; I = MF * ARM * ARM / 3.0 + LOAD * ARM * ARM; Gm = (MF * ARM / 2.0 + LOAD * ARM) * g; HH = HF / 2.0
    Tb = np.full(n, par[11]); xb = np.ones(n); hb = np.ones(n); Tsb = Tb.copy(); xsb = xb.copy(); Tt = Tb.copy(); xt = xb.copy(); ht = hb.copy(); Tst = Tb.copy(); xst = xb.copy()
    per = 1.0 / f; N = int(round(cycles * per / dt)); Q = np.zeros(N); q = 0.0; w = 0.0
    for k in range(N):
        t = k * dt; x = (t - math.floor(t / per) * per) / per
        if x < 0.25: Pb = P; cb = h_on
        elif x < 0.5: Pb = 0.0; cb = HH
        elif x < 0.75: Pb = 0.0; cb = HF
        else: Pb = 0.0; cb = HH
        if x < 0.25: Pt = 0.0; ct = HF
        elif x < 0.5: Pt = 0.0; ct = HH
        elif x < 0.75: Pt = P; ct = h_on
        else: Pt = 0.0; ct = HH
        eb = EPSL - R * q / L0
        if eb < 0.0: eb = 0.0
        et = EPSL + R * q / L0
        Fb = step_strands(Tb, xb, hb, Tsb, xsb, off, par, Pb, dt, eb, cb); Ft = step_strands(Tt, xt, ht, Tst, xst, off, par, Pt, dt, et, ct)
        w += (R * (Fb - Ft) - Gm * math.sin(q) - DAMP * w) / I * dt; q += w * dt
        if q < 0.0:
            q = 0.0
            if w < 0.0: w = 0.0
        if q > QMAX:
            q = QMAX
            if w > 0.0: w = 0.0
        Q[k] = q * 180.0 / math.pi
    return Q


class FastSMA:
    """근육 하나 (가닥 n 개). 김은서 님 SMA 와 같은 쓰임: step(P, dt, eps, hS) → 힘, .T, .xi, .p, power_duty, u_cap.
    d: BASE 와 같은 키의 dict (근육 전체 값). scale: 힘·열용량·가열·냉각을 이 배수로 (다른 등급 환산), L0 는 따로 준다."""

    def __init__(self, d=None, scale=1.0, L0=None):
        d = dict(BASE if d is None else d); self.d = d; n = int(d.get("n", 1)); self.n = n
        dd = dict(d); dd["C"] = d["C"] * scale; dd["FMAX"] = d["FMAX"] * scale; self.par = par_array(dd); self.off = offsets(n, d.get("sigma", 0.0))
        self.p = dict(FMAX=dd["FMAX"], Qon=d["P"] * scale, hS=d["h"] * scale, h_on=d["h_on"] * scale, epsL=d["epsL"], L0=d["L0"] if L0 is None else L0, Tmax=d.get("Tmax", 120.0), K_mod=dd["FMAX"] / d["epsL"])
        self._T = np.full(n, d["T_inf"]); self._xi = np.ones(n); self._h = np.ones(n); self._Ts = self._T.copy(); self._xs = self._xi.copy()

    def step(self, P, dt, eps=None, hS=None):
        return step_strands(self._T, self._xi, self._h, self._Ts, self._xs, self.off, self.par, float(P), float(dt), float(self.p["epsL"] if eps is None else eps), float(self.p["hS"] if hS is None else hS))

    def force(self, eps): return float(np.mean([_force1(x, eps, self.par) for x in self._xi]))
    def power_duty(self, u): return u * self.p["Qon"]
    def u_cap(self, margin=10.0): return min(1.0, max(0.0, (self.p["Tmax"] - self.T) / margin))
    T = property(lambda s: float(s._T.mean()))
    xi = property(lambda s: float(s._xi.mean()))
