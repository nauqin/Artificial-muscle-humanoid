"""시간 영역 힘 설계 공통 함수 (2026-11-19).

지금까지(path_design 기본 `--force-mode amin`)는 케이블 힘을 `1.1 × 토크 용량 ÷ 최소 모멘트암` 으로 정했다. 토크가 가장 큰 순간과
모멘트암이 가장 작은 순간을 한데 묶는 최악 조합이라, 모멘트암이 크게 변하는 팔꿈치·어깨에서 힘이 3배쯤 부풀었다 (README 12절).
여기서는 프레임마다
    T_i(t) = 그 순간 케이블 i 가 내야 하는 장력 (그 순간의 모멘트암으로 나눔),
    s_i(t) = 그 순간 모듈이 줄어들어 있어야 하는 길이 = max_t L_i − L_i(t)
를 구하고, 같은 순간의 (s, T) 쌍이 모듈 등급의 힘–수축 곡선 안에 드는지로 묶음 수를 정한다.

모듈 크기 = 묶음 수 × 등급 최대 힘 × 등급 최대 스트로크 (J, "일 용량"). SMA 선의 부피는 힘(단면적) × 스트로크(길이) 에 비례하므로
질량도 이것에 비례한다. 힘만으로 질량을 매기면(F/ρ) 최적화가 스트로크가 가장 긴 등급으로만 도망간다.
보고용으로 예전 방식 질량(묶음 × 최대 힘 ÷ 10.8 N/g, park2019)도 같이 낸다.

힘–수축 곡선 (report/module_curves.json, **임시**): 등급 g 한 묶음이 하중 f 에서 낼 수 있는 수축 = S_g · (1 − β f / F_g).
즉 무부하면 최대 스트로크 S_g, 최대 힘 F_g 에서는 (1−β)·S_g. 거꾸로, 수축 s 에서 낼 수 있는 힘 = F_g · min(1, (1 − s/S_g)/β).
β = 0 이면 지금까지처럼 힘·스트로크를 동시에 다 쓴다고 보는 것. SMA 는 하중이 커질수록 수축률이 줄어 최대 힘과 최대 변위를 동시에 쓰기 어렵다.
"""
import json, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
from actuator_model import BY_NAME                       # noqa: E402  과제 계획서 모듈표 (힘 N, 스트로크 mm)

RHO_PARK = 130.0 / 12.0                                  # N/g, park2019 SMA 스프링 묶음 (발표 질량 추정과 같은 기준)
CURVES_DEFAULT = os.path.join(ROOT, "report", "module_curves.json")


# ─────────────────────────────── 힘–수축 곡선 ───────────────────────────────
def load_curves(path=None, beta=None):
    """등급 이름 → dict(F 최대 힘 N, S 최대 스트로크 m, beta). beta 를 주면 모든 등급에 그 값을 쓴다 (민감도)."""
    js = json.load(open(path or CURVES_DEFAULT, encoding="utf-8"))
    out = {}
    for name, m in BY_NAME.items():
        g = js.get("grades", {}).get(name, {})
        b = float(beta if beta is not None else g.get("beta", js.get("beta", 0.0)))
        # 등급별 덮어쓰기 (2026-11-20, 실측 곡선을 넣을 때): grades.<등급>.F_N / S_mm / cap_mm (cap = 이 수축을 넘으면 힘 0 — 실측 범위 밖은 쓰지 않는다)
        out[name] = dict(F=float(g.get("F_N", m.force_n[1])), S=float(g.get("S_mm", m.stroke_mm[1])) / 1000.0, beta=b)
        if "cap_mm" in g: out[name]["cap"] = float(g["cap_mm"]) / 1000.0
    return out


def f_avail(cv, s):
    """한 묶음이 수축 s (m, 배열) 에서 낼 수 있는 힘 (N). s > S 이면 0."""
    x = np.asarray(s, float) / cv["S"]
    f = np.where(x <= 1.0, cv["F"], 0.0) if cv["beta"] <= 0 else cv["F"] * np.clip((1.0 - x) / cv["beta"], 0.0, 1.0)
    return np.where(np.asarray(s, float) > cv["cap"], 0.0, f) if "cap" in cv else f


MAX_BUNDLES = 4                                          # 한 케이블에 묶는 모듈 수 상한 (예전 가정 "M3 450 N × 4 묶음 = 1,800 N" 과 같은 4). 일 용량만 보면 M4 를 10개 넘게 묶는 답이 나온다


def size_cable(T, s, curves, grades, max_bundles=None):
    """한 케이블의 프레임별 (장력 T N, 수축 s m) → 허용 등급 중 일 용량(묶음 올림)이 가장 작은 등급 (묶음 수 ≤ max_bundles).
    반환 dict(grade, n, n_cont, work_J, work_cont_J, mass_park_g, T_peak, s_max, obj, feasible, per_grade).
    obj 는 최적화용 연속값(J): 가능한 등급 중 묶음 수를 올림하지 않은 일 용량의 최소. 가능한 등급이 없으면
    스트로크가 가장 긴 허용 등급으로 계산하고 초과량 벌점을 곱한다 (최적화가 그쪽으로 돌아오도록)."""
    T = np.maximum(np.asarray(T, float), 0.0); s = np.maximum(np.asarray(s, float), 0.0)
    Tpk = float(T.max()) if len(T) else 0.0; smax = float(s.max()) if len(s) else 0.0
    nb = MAX_BUNDLES if max_bundles is None else max_bundles
    per = {}
    for g in grades:
        cv = curves[g]
        if smax > min(cv["S"], cv.get("cap", cv["S"])): per[g] = None; continue
        fa = f_avail(cv, s); need = T > 1e-6
        if np.any(need & (fa <= 0)): per[g] = None; continue
        n_cont = float((T[need] / fa[need]).max()) if need.any() else 0.0
        n = max(1, int(np.ceil(n_cont - 1e-9)))
        if n > nb: per[g] = None; continue
        # 최적화용 연속값: 묶음은 최소 1 개라 max(1, n_cont) — 그냥 n_cont 면 큰 등급(M1)을 일부만 쓰는 해가 싸 보여 최적화가 M1·최대 스트로크로 간다
        # (175 cm 하체 M0: 767 N 에 M1 600 J 를 골랐다, 2026-11-19). + 0.1 n_cont 는 같은 등급 안에서 장력이 작은 쪽을 고르게 하는 기울기
        per[g] = dict(n_cont=n_cont, n=n, work_cont=(max(1.0, n_cont) + 0.1 * n_cont) * cv["F"] * cv["S"], work=n * cv["F"] * cv["S"], mass_park=n * cv["F"] / RHO_PARK)
    ok = {g: v for g, v in per.items() if v is not None}
    if ok:
        g = min(ok, key=lambda k: (ok[k]["work"], ok[k]["n"], curves[k]["F"]))
        obj = min(v["work_cont"] for v in ok.values())
        b = ok[g]
        return dict(grade=g, n=b["n"], n_cont=b["n_cont"], work_J=b["work"], work_cont_J=b["work_cont"], mass_park_g=b["mass_park"],
                    T_peak=Tpk, s_max=smax, obj=obj, feasible=True, per_grade=per)
    gL = max(grades, key=lambda k: (curves[k]["S"], curves[k]["F"])); cv = curves[gL]   # 어느 등급도 안 됨 → 스트로크 가장 긴 등급 + 벌점 (스트로크 초과, 묶음 상한 초과)
    sc = np.minimum(s, 0.999 * cv["S"]); fa = f_avail(cv, sc); need = T > 1e-6
    n_cont = float((T[need] / np.maximum(fa[need], 1e-6 * cv["F"])).max()) if need.any() else 0.0
    # 배정 불가는 어떤 가능한 해보다 비싸야 한다: 허용 등급 중 가장 큰 모듈 묶음 상한 × 1.5 를 바닥으로 깔고, 넘친 정도만큼 더한다
    # (바닥 없이 약한 벌점만 주었더니 2단계가 방향오차를 줄이려고 배정 불가로 넘어갔다 — 175 cm 팔꿈치 브래킷 20·50 mm, 2026-11-19)
    floor = 1.5 * nb * max(curves[g]["F"] * curves[g]["S"] for g in grades)
    obj = floor * (1.0 + 10.0 * max(0.0, smax / cv["S"] - 1.0) + max(0.0, n_cont / nb - 1.0)) + 0.1 * n_cont * cv["F"] * cv["S"]
    return dict(grade=None, n=None, n_cont=n_cont, work_J=None, work_cont_J=None, mass_park_g=None, T_peak=Tpk, s_max=smax, obj=obj, feasible=False, per_grade=per)


# ─────────────────────────────── 프레임별 몫·장력 ───────────────────────────────
def shares(V, c, tau, lam=1e-4, iters=600):
    """토크 공간 몫 (1단계용, 경로가 아직 없을 때): 프레임마다 min ‖Σ_i v_i u_i − τ_t‖² + λ Σ (u_i/c_i)², 0 ≤ u_i ≤ c_i.
    V (k, n), c (k,), tau (N, n) → u (N, k) N·m. design_directions 의 LP 와 같은 제약에 정칙화를 더해 프레임마다 유일한 몫."""
    A = (np.asarray(V, float) * np.asarray(c, float)[:, None]).T                 # (n, k), a = u/c ∈ [0, 1]
    AtA = A.T @ A; sig = float(np.linalg.eigvalsh(AtA).max()); lamA = lam * sig
    H = AtA + lamA * np.eye(A.shape[1]); step = 1.0 / (sig + lamA); Atb = tau @ A   # (N, k)
    a = np.zeros((len(tau), A.shape[1])); y = a.copy(); tk = 1.0
    for _ in range(iters):
        g = y @ H - Atb; an = np.clip(y - step * g, 0.0, 1.0)
        tn = (1 + np.sqrt(1 + 4 * tk * tk)) / 2; y = an + ((tk - 1) / tn) * (an - a); a, tk = an, tn
    return a * np.asarray(c, float)[None, :]


def allocate(R, tau, C, lam=1e-4, iters=400, a0=None):
    """실제 모멘트암으로 장력 분배 (2단계·평가용): 프레임마다 min ‖R_t T_t − τ_t‖² + λ_t Σ (T_i/C_i)², T ≥ 0.
    R (N, n, k) m, tau (N, n) N·m, C (k,) N (정규화 크기, 예: 1단계 T_peak). → T (N, k) N, 잔차 (N, n) N·m, a (N, k) (다음 호출 warm start).
    Moco 의 Σu² (u = T/F) 와 같은 꼴이라 설계와 검증의 분배가 맞는다."""
    C = np.maximum(np.asarray(C, float), 1.0); A = R * C[None, None, :]
    AtA = np.einsum("tij,tik->tjk", A, A); Atb = np.einsum("tij,ti->tj", A, tau)
    sig = np.linalg.eigvalsh(AtA)[:, -1] + 1e-12; lamt = lam * sig
    H = AtA + lamt[:, None, None] * np.eye(A.shape[2])[None]; step = 1.0 / (sig + lamt)
    a = np.zeros(Atb.shape) if a0 is None else np.maximum(a0, 0.0).copy(); y = a.copy(); tk = 1.0
    for _ in range(iters):
        g = np.einsum("tjk,tk->tj", H, y) - Atb
        an = np.maximum(y - step[:, None] * g, 0.0)
        tn = (1 + np.sqrt(1 + 4 * tk * tk)) / 2; y = an + ((tk - 1) / tn) * (an - a); a, tk = an, tn
    T = a * C[None, :]; res = np.einsum("tij,tj->ti", R, T) - tau
    return T, res, a


def contraction(L):
    """케이블 길이 L (N,) 또는 (N, k) → 모듈 수축 s = max L − L (가장 늘어났을 때 0)."""
    L = np.asarray(L, float); return L.max(axis=0) - L


def fmt_size(z):
    """size_cable 결과 한 줄 요약."""
    if not z["feasible"]: return f"배정 불가 (수축 {z['s_max']*1000:.0f} mm, 장력 최대 {z['T_peak']:.0f} N)"
    return f"{z['grade']} ×{z['n']} (연속 {z['n_cont']:.2f}, 일 용량 {z['work_J']:.0f} J, park {z['mass_park_g']:.0f} g)  장력 최대 {z['T_peak']:.0f} N, 수축 최대 {z['s_max']*1000:.0f} mm"
