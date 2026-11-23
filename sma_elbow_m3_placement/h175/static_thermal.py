"""내 정적 최적화 2 — '열까지 따진 정적' (정어진, 11/23)

MuJoCo(팔 운동) 없이 계산식만으로 자리를 고르되, 유빈식 · 내 첫 정적처럼 '힘'만 보지 않고
SMA 가 그 힘을 내려면 몇 도까지, 얼마나 빨리 데우고 식혀야 하는지까지 계산해서 '못 하는 만큼'을 줄인다.

들기 동작의 매 순간 t, 케이블 하나 (M3 × n 묶음, 은서 v3 값):
  ① 필요한 장력      T(t) = τ(t) / r(t)                                  (τ: 역동역학 팔꿈치 토크, r: 모멘트암)
  ② 근육 변형률      ε(t) = εL − s(t) / L0                               (s: 가장 길 때보다 줄어든 길이)
  ③ 필요한 변태      n·K(ξ)·(ε − εL·ξ) = T  를 ξ 에 대해 풂 → ξ_req(t)  (은서 힘 식 T3.1, K(ξ) = K·(EM ξ + EA (1−ξ))/EA)
  ④ 필요한 온도      데울 때: T_h = As' + (Af' − As')·acos(2ξ−1)/π,  식힐 때: T_c = Mf' + (Ms' − Mf')·acos(2ξ−1)/π
                     (' = 응력만큼 문턱 이동 σ/C_sig, σ = T / (n·A_wire))
  ⑤ 필요한 열        P = m c · dT/dt(+) + h_on (T − 25) + λ m · (−dξ/dt)(+)   → 267 W 넘는 만큼이 '가열 부족'
  ⑥ 식힐 수 있는 빠르기  −dT/dt ≤ h_forced (T − 25) / (m c)                  → 넘는 만큼이 '냉각 부족'
  목적 = 평균(가열 부족/267 + 냉각 부족/(100 K/s) + 최고 온도 초과/10) + 0.1·평균(1 − ξ_req)
  (힘을 못 내는 순간 = 큰 벌점). 조건(스트로크 72 mm, 근육 0.18 m, 모멘트암, 뼈 거리)은 다른 정적과 같음.
돌리기:  NM3=4 python static_thermal.py   → out/static_thermal_<n>.json
"""
import sys, os, json, math, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); os.chdir(D); sys.path.insert(0, D); sys.argv = sys.argv[:1]
g = {'__file__': os.path.join(D, 'place_elbow.py'), '__name__': 'pe'}
exec(open('place_elbow.py', encoding='utf-8').read().split("if __name__ == '__main__':")[0], g)
NM3 = int(os.environ.get('NM3', '4'))
P = dict(m=1.8 / 600, c=600.0, lam=24180.0, As=55.0, Af=72.0, Ms=52.0, Mf=38.0, Csig=8.0e6, Q=267.0, Tmax=120.0,
         K=450.0 / 0.40, EA=65e9, EM=28e9, epsL=0.40, L0=0.18, A=450.0 / 400e6, h_on=0.72, h_f=4.62)
S = P['epsL'] * P['L0']


def thermal_cost(G, X, cable, n=NM3):
    sign = g['SIGN'][cable]; L, r = g['lengths'](G, X); rr = sign * r; nt = G['nt']
    s = (L.max(1, keepdims=True) - L)[:, :nt]; eps = np.clip(P['epsL'] - s / P['L0'], 0, None)
    tau = np.clip(sign * G['tau'], 0, None)[None]; T = np.where(tau > 1e-6, tau / np.maximum(rr[:, :nt], 1e-3), 0.0)
    # ③ n K (1 − b ξ)(ε − εL ξ) = T  →  b εL ξ² − (εL + b ε) ξ + ε − T/(nK) = 0 , 작은 근
    b = 1 - P['EM'] / P['EA']; nK = n * P['K']
    qa = b * P['epsL']; qb = -(P['epsL'] + b * eps); qc = eps - T / nK
    disc = np.clip(qb ** 2 - 4 * qa * qc, 0, None); xi = (-qb - np.sqrt(disc)) / (2 * qa)
    short = np.clip(T / nK - eps, 0, None) / 0.05                      # ξ = 0 이어도 힘이 모자람 (스트로크 · 힘 부족)
    xi = np.clip(xi, 0.0, 1.0)
    sh = T / (n * P['A']) / P['Csig']
    a = np.arccos(np.clip(2 * xi - 1, -1, 1)) / math.pi
    Th = P['As'] + sh + (P['Af'] - P['As']) * a; Tc = P['Mf'] + sh + (P['Ms'] - P['Mf']) * a
    dt = np.gradient(G['t'])[None]; dxi = np.gradient(xi, axis=1)
    Tn = np.where(dxi > 1e-6, Tc, Th); Tn = np.where(xi >= 1 - 1e-6, P['As'] + sh - 1.0, Tn)   # 다 늘어나 있으면 As 아래면 됨
    dT = np.gradient(Tn, axis=1) / dt
    Pneed = P['m'] * P['c'] * np.clip(dT, 0, None) + P['h_on'] * (Tn - 25) + P['lam'] * P['m'] * np.clip(-dxi / dt, 0, None)
    heat = np.clip(Pneed - P['Q'], 0, None) / P['Q']
    cool = np.clip(-dT - P['h_f'] * (Tn - 25) / (P['m'] * P['c']), 0, None) / 100.0
    hot = np.clip(Tn - P['Tmax'], 0, None) / 10.0
    J = np.minimum(heat + cool + hot, 5.0).mean(1) + 0.1 * (1 - xi).mean(1) + 10 * short.mean(1)
    return J, dict(heat=heat.mean(1), cool=cool.mean(1), hot=hot.mean(1), xi_min=xi.min(1), short=short.mean(1))


def objective(G, X, cable):
    _, F, use, rmin, Lmin = g['evaluate'](G, X, cable, S, True)            # 기하 조건 (다른 정적과 같음)
    v = np.where(np.isfinite(F), 0.0, 1.0)
    J, _ = thermal_cost(G, X, cable)
    return np.where(v > 0, 1e3 + g['evaluate'](G, X, cable, S, True)[0] / 1e2, J)


if __name__ == '__main__':
    from scipy.optimize import differential_evolution
    G, yubin = g['task']('path'); out = {'NM3': NM3, 'opt': {}, 'yubin': {}}
    for c in g['CABLES']:
        best = None
        for sd in (0, 1):
            r = differential_evolution(lambda X: objective(G, X.T, c), list(zip(g['LO'], g['HI'])), seed=sd, maxiter=300, popsize=30,
                                       tol=1e-10, vectorized=True, updating='deferred', polish=False, init='sobol')
            if best is None or r.fun < best.fun: best = r
        J, d = thermal_cost(G, best.x[None], c); Jy, dy = thermal_cost(G, yubin[c][None], c)
        out['opt'][c] = dict(x=best.x.tolist(), J=float(J[0]), **{k: float(v[0]) for k, v in d.items()})
        out['yubin'][c] = dict(J=float(Jy[0]), **{k: float(v[0]) for k, v in dy.items()})
        print(c, '열 정적 최적', {k: round(v, 3) for k, v in out['opt'][c].items() if k != 'x'}, '| 유빈 자리', {k: round(v, 3) for k, v in out['yubin'][c].items()}, flush=True)
    json.dump(out, open(f'out/static_thermal_{NM3}.json', 'w'), indent=1)
