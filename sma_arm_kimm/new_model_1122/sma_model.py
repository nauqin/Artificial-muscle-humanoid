"""
SMA 근육 열 모델 (M3 모듈, 기계연 실제 값) — 한 칸씩 (step) 계산, MuJoCo 에 끼우기용
김은서 · v3 (11/20): 기계연 실제 값만 남김 (예전 값 · 스프링 힘 식 뺌). v2 (11/18): 온도 튐 수정

식
  ① 열    : m·c_p·dT/dt = P − hS·(T − T∞) + 잠열   (온도와 상태 ξ 를 같이 풂)
  ② 상태  : Brinson cos 식, 방향 바뀐 순간의 ξ 에서 이어감 (부분 루프)
  ③ 하중  : 응력만큼 변태 문턱이 올라감  As′ = As + σ/C  (교수님 T3.1)
  ④ 힘    : F = K(ξ)·(ε − ε_L·ξ)  (교수님 T3.1)

쓰는 법 (근육마다 객체 하나, 매 칸 호출)
  sma = SMA()                                         # = SMA(M3), 기계연 실제 값
  F = sma.step(열[W], dt, eps=ε, hS=냉각세기)          # → 힘 [N]
  열: 267 × (전류 / 36)² [W]  또는  sma.power_duty(u)  (u 0~1, 1 = 36 A = 267 W, 상한 sma.u_cap())
  냉각세기 hS [W/K]: 강제 4.62 (기본값) · 반냉각 2.31 · 켜는 동안 0.72 [가정]
  상태: sma.T (온도 °C), sma.xi (1 = 안 줄어듦, 0 = 다 줄어듦)

값 출처
  기계연 회신 · 원본 데이터: 근육 유효 길이 0.18 m, 최대 수축 0.40 (하중 곡선 3~5 kg 39~42 %)
  11쪽 (M3 하나 + 5 kg, 0.5 Hz, 켜는 0.5 s, 36 A) 에 맞춘 열 값: 가열 267 W, 강제 냉각 4.62 W/K
     (변위 오차 1.0 mm, 50 % · 20 % 풀림 0.40 · 0.57 s, 실측 0.38 · 0.71) — 15쪽 (팔) 은 맞추는 데 안 씀
  교수님 T3.1: 열용량 1.8 J/K, 변태 온도 As/Af/Ms/Mf 55/72/52/38 °C, 문턱 이동 8 MPa/K, 최고 온도 120 °C, 힘 450 N
  문헌 (Rodinò 2023): 잠열 24 180 J/kg, 비열 600 J/(kg·K)
  ※ 가열 W 는 모델 안에서 근육에 들어간 열. 기계연 전기 입력은 22.2 V × 36 A ≈ 800 W (효율 약 1/3, 어림)
  ※ 36 A 가 장치 최대 전류인지는 확인 안 됨 (기계연: 5 · 6 kg 하중 시험에서 "거의 최대 전류")
"""
import math

# ════════════════════ 값 (기계연 실제 값) ════════════════════
_Ca, _cp = 1.8, 600.0                 # 열용량 [J/K] (교수님 T3.1), 비열 [J/(kg·K)] (Rodinò) → 질량 3 g
M3 = dict(
    name="M3 (기계연 실제 값)",
    m=_Ca / _cp, cpA=_cp, cpM=_cp, lam=24180.0,
    As=55.0, Af=72.0, Ms=52.0, Mf=38.0, T_inf=25.0,
    C_sig=8.0e6,                      # 응력 1 Pa 당 문턱 이동 [Pa/K]
    hS=4.62,                          # 기본 냉각 세기 = 강제 냉각 [W/K] (11쪽에 맞춤)
    Qon=267.0,                        # 듀티 1 일 때 열 [W] = 36 A (11쪽에 맞춤)
    Tmax=120.0,                       # 최고 온도 (u_cap 상한, 교수님)
    FMAX=450.0, epsL=0.40, EA=65e9, EM=28e9, L0=0.18,
)
M3["A_wire"] = M3["FMAX"] / 400e6     # 응력 환산 단면적 (교수님 가정)
M3["K_mod"] = M3["FMAX"] / M3["epsL"] # 모듈 강성 (T3.1 힘 식)
KIMM_HEAT = dict(P=M3["Qon"], h_forced=M3["hS"], h_half=M3["hS"] / 2, h_on=0.72)   # 열 [W] @ 36 A · 강제 · 반냉각 · 켜는 동안 [가정] [W/K]
M3_KIMM = M3                          # 예전 이름 (v2 코드 호환)


# ════════════════════ 모델 ════════════════════
class SMA:
    def __init__(s, p=None, latent=True, T0=None, xi0=1.0, hS=None, solver="joint"):
        """p: 값 묶음 (기본 M3 = 기계연 실제 값) · latent: 잠열 넣기 · T0/xi0: 처음 온도 · 상태 · hS: 기본 냉각 세기"""
        p = M3 if p is None else p
        s.p = p
        s.lam = p["lam"] if latent else 0.0
        s.T = p["T_inf"] if T0 is None else T0
        s.xi = xi0
        s.heating = True; s.T_s = s.T; s.xi_s = xi0      # 방향 바뀐 순간 값 (Brinson ξ_M0 / ξ_A0)
        s.hS = p.get("hS") if hS is None else hS          # 냉각 세기 [W/K]
        s.sig = 0.0                                        # 지금 응력 [Pa] (문턱 이동용)
        s.solver = solver                                  # "joint" (온도·ξ 같이 풂, 11/18) / "old" (예전 방식)

    # ── 가열 입력 ──
    def power_duty(s, u):
        """듀티 0~1 → 근육에 들어가는 열 [W]. 1 = 36 A (267 W). 전류로 주려면 열 = 267 × (I / 36)²"""
        return u * s.p["Qon"]

    # ── 문턱 (하중만큼 이동) ──
    def _thr(s):
        p = s.p
        sh = 0.0 if p["C_sig"] is None else s.sig / p["C_sig"]
        return p["As"] + sh, p["Af"] + sh, p["Ms"] + sh, p["Mf"] + sh

    # ── ② 온도 → ξ  ──
    def _xi_branch(s, T):
        As, Af, Ms, Mf = s._thr()
        if s.heating:
            f = 1.0 if T <= As else 0.0 if T >= Af else 0.5 * (math.cos(math.pi * (T - As) / (Af - As)) + 1)
            return s.xi_s * f
        f = 0.0 if T >= Ms else 1.0 if T <= Mf else 0.5 * math.cos(math.pi * (T - Mf) / (Ms - Mf)) + 0.5
        return s.xi_s + (1 - s.xi_s) * f

    def _dxi_dT(s, T):
        As, Af, Ms, Mf = s._thr()
        if s.heating:
            if T <= As or T >= Af: return 0.0
            a = math.pi / (Af - As)
            return -(s.xi_s / 2) * a * math.sin(a * (T - As))
        if T >= Ms or T <= Mf: return 0.0
        a = math.pi / (Ms - Mf)
        return -((1 - s.xi_s) / 2) * a * math.sin(a * (T - Mf))

    # ── ① 한 칸 계산 ──
    def step(s, P, dt, eps=None, hS=None):
        """P: 넣는 열 [W], dt: 시간 [s], eps: 지금 변형률 (힘 계산·문턱 이동용, 없으면 무시)
        hS: 이번 칸 냉각 세기 (팬 켜고 끄기 등, 없으면 기본값). 반환: 힘 [N] (eps 없으면 None)"""
        p = s.p
        if eps is not None and p["C_sig"] is not None:
            s.sig = s.force(eps) / p["A_wire"]
        hS = s.hS if hS is None else hS
        x, T = s.xi, s.T
        q = P - hS * (T - p["T_inf"])                     # 줄열 − 냉각

        h_now = q > 0                                     # 데우는 중 / 식는 중 
        if h_now != s.heating:
            s.heating, s.T_s, s.xi_s = h_now, T, x

        cp = p["cpA"] * (1 - x) + p["cpM"] * x
        if s.solver == "old":                             # 예전 방식 (비교용)
            dxi_dT = s._dxi_dT(T)
            c_eff = cp - s.lam * min(dxi_dT, 0.0)
            T_new = T + dt * q / (p["m"] * c_eff)
            xi_new = s._xi_branch(T_new)
            xi_new = min(x, xi_new) if s.heating else max(x, xi_new)
            T_new += s.lam * ((xi_new - x) - dxi_dT * (T_new - T)) / cp
        else:                                             # 온도와 ξ 를 같이 풂 (에너지 맞춤)
            T_new, xi_new = s._solve(T + dt * q / (p["m"] * cp), x, cp)
        s.T, s.xi = T_new, xi_new
        return s.force(eps) if eps is not None else None

    # ── ④ 힘 ──
    def force(s, eps):
        """교수님 T3.1: F = K(ξ)·(ε − ε_L·ξ)"""
        p = s.p
        Kx = p["K_mod"] * (s.xi * p["EM"] + (1 - s.xi) * p["EA"]) / p["EA"]
        return max(0.0, Kx * (eps - p["epsL"] * s.xi))

    def _solve(s, T_pred, x, cp):
        """잠열을 넣은 뒤의 온도 T 와 ξ 를 함께 구함 (11/18 수정)
        조건: T = T_pred + (λ/c_p)·(ξ(T) − x),  ξ(T) = 지금 방향의 cos 식 (데울 땐 x 이하, 식힐 땐 x 이상)
        왼쪽 − 오른쪽은 T 에 대해 늘어나기만 하므로 답이 하나 → 이분법으로 찾음.
        예전 방식은 문턱이 한 칸에 크게 움직이면 ξ 가 0↔1 로 통째로 바뀌고 그 잠열 (최대 ±40 °C) 이
        한 번에 들어가 온도가 문턱 밖까지 튀었음. 이 방식은 잠열을 넣은 뒤의 온도에서 ξ 를 다시 맞추므로
        온도가 문턱 안에 머물고, ξ 는 그 온도에 맞는 만큼만 바뀜."""
        k = s.lam / cp
        def xi_of(TT):
            xb = s._xi_branch(TT)
            return min(x, xb) if s.heating else max(x, xb)
        if k == 0.0:
            return T_pred, xi_of(T_pred)
        g = lambda TT: TT - T_pred - k * (xi_of(TT) - x)
        lo, hi = T_pred - k - 1.0, T_pred + k + 1.0          # ξ 변화는 ±1 안 → 답은 이 범위 안
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            if g(mid) > 0: hi = mid
            else: lo = mid
        T_new = 0.5 * (lo + hi)
        return T_new, xi_of(T_new)

    def u_cap(s, margin=10.0):
        """열 상한 (교수님과 같음)"""
        return min(1.0, max(0.0, (s.p["Tmax"] - s.T) / margin))


def heat_from_current(I):
    """전류 [A] → 근육에 들어가는 열 [W]  (36 A = 267 W, 전류 제곱 비례)"""
    return M3["Qon"] * (I / 36.0) ** 2


def hS_from_tau(p, tau):
    """냉각 시간상수 τ → 냉각 세기 h·S """
    return p["m"] * p["cpA"] / tau
