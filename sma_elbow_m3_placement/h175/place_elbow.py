# -*- coding: utf-8 -*-
"""③ SMA 성질을 넣은 팔꿈치 케이블 배치 최적화 (논문 핵심 실험의 첫 단계)

질문: 유빈 최종 A안 부착점은 "이상적 근육"(OpenSim/Moco) 기준으로 고른 것이다.
      SMA 근육의 세 가지 성질을 넣고 부착점을 다시 고르면 필요한 근육이 얼마나 작아지나?
        ① 스트로크 한계  : 줄어드는 길이 ≤ 스트로크 S = εL(0.33) × 근육 길이 L0
        ② 줄수록 힘 빠짐 : 낼 수 있는 힘 = 최대 힘 × (1 − 줄어든 길이 / S)   ← 3상태 힘 식 (뜨거울 때)
        ③ 근육이 들어갈 자리 : 케이블 길이(가장 짧을 때) ≥ 근육 길이 L0 = S / 0.33   [가정: 근육이 케이블 경로에 일자로 놓임]
      + 모멘트암 ≥ 5 mm (방향 맞게, 유빈 path_design 과 같은 하한)
      + 뼈와 거리 ≥ 20 mm : 케이블이 위팔·아래팔 뼈 중심선에서 20 mm 이상 떨어짐 (뼈 뚫고 지나가기 방지. 유빈 배치는 24~37 mm)

방법 (정적 계산 — 자세마다 따로 계산해서 빠름)
  필요한 팔꿈치 토크 τ(t) = 유빈 Moco 들기 결과의 케이블 활성도 × 설계 힘 × 유빈 배치 모멘트암 (두 케이블 합)
  부착점 후보 = 유빈 design_space.py 의 elbow_r 상자 (위팔 점 3개 좌표 + 아래팔 점 3개 좌표 = 6개 변수)
  후보마다: 필요한 최대 힘 F_need = max_t  τ(t) / r(t) / (1 − s(t)/S)
  → F_need 가 가장 작은 부착점을 차분진화(DE, scipy)로 찾음. 굽힘(M0)·폄(M1) 따로, 스트로크 S 마다 따로.

실행 (04_코드\\elbow_finalA 에서, conda activate osim) — 자세한 순서는 place_elbow_실행방법.txt
  python place_elbow.py search [path]  → ① 정적 최적화 (기본: 들기 경로 + 버티기 조건 / path: 들기 경로만)  약 5분
  python place_elbow.py lift [path]    → ② 유빈 배치 vs SMA 반영 배치 물병 들기 (MuJoCo)  약 1분
  python place_elbow.py dyn            → ③ 2단계: 기하 조건으로 후보 추리기 → MuJoCo 로 직접 돌려 고르기  약 3~5분
  python place_elbow.py view dyn        → MuJoCo 창: dyn 으로 고른 배치 (2단계, 실제 속도 재생)
  python place_elbow.py view opt        → MuJoCo 창: 3단계 최적화 배치 (opt_dyn.py → validate_opt.py 결과)
  3단계 최적화:  pip install cma  →  python opt_dyn.py 25 0.03  (약 10~20분)  →  python validate_opt.py  (약 2~3분)
  python place_elbow.py view 1 [path]  → MuJoCo 창 (1 유빈 · 2 SMA 반영 같은 M1 · 3 정적 최적 배치 · 필요 힘)
필요: elbow_finalA.py, lift_finalA.py, sma_model.py, 유빈 osim · mot, finalA_elbow_r.xml (elbow_finalA.py build 로 만든 것)
"""
import sys, os, math, json, copy, re, time, logging
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
logging.getLogger('matplotlib').setLevel(logging.ERROR)

NAME = 'gui_upper_robot40h160_finalA_lift_s6_PFWJ_2'
CABLES = ['elbow_M0_r', 'elbow_M1_r']
SIGN = {'elbow_M0_r': 1, 'elbow_M1_r': -1}                         # 굽힘 = 토크 + 담당, 폄 = − 담당
KOR = {'elbow_M0_r': '굽힘(M0)', 'elbow_M1_r': '폄(M1)'}
DESIGN_F = {'elbow_M0_r': 2327.05, 'elbow_M1_r': 1834.86}
BOX = {'humerus_r': ((-0.07, 0.05), (-0.27, -0.02), (-0.05, 0.05)),   # 유빈 report/design_space.py elbow_r (m)
       'ulna_r': ((-0.08, 0.05), (-0.20, 0.04), (-0.04, 0.04))}
EPSL = 0.40
STROKES = [0.050, 0.060, 0.070, 0.086]                             # 50 = M3(보정, L0 150 mm) · 60 = M3 최대 · 86 = M1 크기(L0 260 mm)
S_MAIN = 0.086
F_M3 = 450.0
R_MIN = 0.005
MARGIN = 1.2                                                       # 버티기 토크 여유 (가속·제어 오차 몫)
CLEAR = 0.020                                                      # 뼈 중심선과 최소 거리
WRIST = np.array([-0.0164, -0.2463, 0.0399])   # 175 cm: radioulnar + radius_hand (osim)                           # 아래팔(척골) 좌표의 손목 근처 (osim radioulnar + radius_hand 위치)
LO = np.array([BOX['humerus_r'][k][0] for k in range(3)] + [BOX['ulna_r'][k][0] for k in range(3)])
HI = np.array([BOX['humerus_r'][k][1] for k in range(3)] + [BOX['ulna_r'][k][1] for k in range(3)])
OUT = os.path.join(HERE, 'out')


def find(fname):
    for c in [os.path.join(HERE, '..', '..', 'report', 'moco', fname), os.path.join(HERE, fname), os.path.join('report', 'moco', fname), os.path.join(r'C:\work\addb_opensim\report\moco', fname)]:
        if os.path.exists(c): return c
    sys.exit(f'{fname} 을 못 찾았어요 → 이 폴더로 복사:  copy C:\\work\\addb_opensim\\report\\moco\\{fname} .')


def font():
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['Malgun Gothic', 'NanumGothic', 'DejaVu Sans']; plt.rcParams['axes.unicode_minus'] = False
    return plt


# ───────────────────────── 들기 동작 · 필요한 토크 ─────────────────────────
def rot_set(J, q):
    c, s = np.cos(q), np.sin(q); z, o = np.zeros_like(q), np.ones_like(q)
    Rz = np.stack([np.stack([c, -s, z], -1), np.stack([s, c, z], -1), np.stack([z, z, o], -1)], -2)
    dRz = np.stack([np.stack([-s, -c, z], -1), np.stack([c, -s, z], -1), np.stack([z, z, z], -1)], -2)
    R = J['R_p'] @ Rz @ J['R_c'].T; dR = J['R_p'] @ dRz @ J['R_c'].T
    return dict(R=R, dR=dR, T0=J['t_p'] - R @ J['t_c'], dT0=-dR @ J['t_c'])


def hold_torque(q_grid):
    """각 팔꿈치 각도에서 물병을 들고 '버티는' 데 필요한 토크 — 들기 동작 중 나온 모든 위팔 방향(중력 방향) 중 최악값.
    MuJoCo 로 계산: 속도 0 에서 qfrc_bias = 중력 토크 (+ = 굽힘 토크가 필요)"""
    import mujoco
    src = open(os.path.join(HERE, 'lift_finalA.py'), encoding='utf-8').read()
    ns = {'__file__': os.path.join(HERE, 'lift_finalA.py'), '__name__': 'lift'}; exec(src[:src.index('mode = sys.argv')], ns)
    REF = ns['prepare']()
    m = mujoco.MjModel.from_xml_string(open(os.path.join(HERE, 'finalA_elbow_r.xml'), encoding='utf-8').read()); d = mujoco.MjData(m)
    qa = m.jnt_qposadr[m.joint('elbow_flex_r').id]; va = m.jnt_dofadr[m.joint('elbow_flex_r').id]
    B = np.zeros((len(REF['grav'][::5]), len(q_grid)))
    for i, g in enumerate(REF['grav'][::5]):
        m.opt.gravity[:] = g
        for j, qq in enumerate(q_grid):
            d.qpos[qa] = qq; d.qvel[:] = 0; mujoco.mj_forward(m, d); B[i, j] = d.qfrc_bias[va]
    return B.max(0), -B.min(0)


def arm175():
    """유빈 175 cm 팔꿈치 설계 osim (레포 report/moco) → 케이블 이름만 elbow_M0_r / elbow_M1_r 로 바꾼 사본 (out/)"""
    src = os.path.join(HERE, '..', '..', 'report', 'moco', 'newmus_paths_elbow_r_k2_40kg_h175sk_td_br50_rb_load3.79_model.osim')
    dst = os.path.join(HERE, 'out', 'arm175.osim')
    if not os.path.exists(dst):
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        open(dst, 'w', encoding='utf-8').write(re.sub(r'name="M([01])_r(_p[0-9]+)?"', r'name="elbow_M\1_r\2"', open(src, encoding='utf-8').read()))
    return dst


def task(req='range'):
    """req='range' : 들기 경로의 토크 + 어느 각도에서든 물병을 버틸 토크 × 여유 MARGIN (기본, 튼튼한 설계)
                    굽힘 근육은 0°~143°(목표 최대), 폄 근육은 6°(목표 최소)~150° 에서 버텨야 함
                    = 목표 범위 안에서는 어느 쪽으로든 버티고, 범위를 넘어가면 되돌릴 수 있어야 함
       req='path'  : 들기 경로(Moco)의 토크만 (1차 시도 — 동적 시뮬레이션에서 실패함)"""
    src = open(os.path.join(HERE, 'elbow_finalA.py'), encoding='utf-8').read()
    ns = {'__file__': os.path.join(HERE, 'elbow_finalA.py'), '__name__': 'elbow_geom'}
    exec(src[:src.index('# ───────────────────────── MuJoCo xml 만들기')], ns)          # osim 읽기 함수만 가져옴
    bodies, joints, cables = ns['read_osim'](arm175())   # 유빈 175 cm 팔꿈치 설계 (M3×3, td_br50_rb)
    L = open(find(NAME + '.mot'), encoding='utf-8').read().splitlines(); i = [k for k, l in enumerate(L) if l.strip() == 'endheader'][0]
    hdr = L[i + 1].split('\t'); D = np.array([[float(x) for x in l.split('\t')] for l in L[i + 2:] if l.strip()])
    deg = any('inDegrees=yes' in l for l in L[:i])
    qt = D[:, hdr.index('elbow_flex_r')] * (math.pi / 180 if deg else 1.0); t = D[:, 0] - D[0, 0]
    J = joints['elbow_r']; nt = len(qt)
    yubin = {n: np.concatenate([cables[n]['points'][0][1], cables[n]['points'][1][1]]) for n in CABLES}
    Ls = open(os.path.join(HERE, '..', '..', 'out', 'upper40_h175sk', 'carrt', 'subject6', 'PFWJ_2', 'PFWJ_2_id.sto'), encoding='utf-8').read().splitlines(); k = [i for i, l in enumerate(Ls) if l.strip() == 'endheader'][0]
    hs = Ls[k + 1].split('\t'); DS = np.array([[float(x) for x in l.split()] for l in Ls[k + 2:] if l.strip()])
    tau = np.interp(D[:, 0], DS[:, 0], DS[:, hs.index('elbow_flex_r_moment')])     # 175 cm 역동역학 팔꿈치 토크 (물병 3.79 kg, 유빈 out/upper40_h175sk)
    q_grid = np.radians(np.arange(0, 151, 5.0))
    q = np.concatenate([qt, q_grid])
    G = dict(rot_set(J, q), q=q, t=t, joint=J, nt=nt, tau=tau, req=req)
    if req == 'range':
        hf, he = hold_torque(q_grid)
        hf = np.where(q_grid <= qt.max() + 1e-6, hf, 0.0)        # 굽힘 근육: 0° ~ 목표 최대(143°) 에서 버티기 (그 위는 되돌릴 일 없음)
        he = np.where(q_grid >= qt.min() - 1e-6, he, 0.0)        # 폄 근육: 목표 최소(6°) ~ 150° 에서 버티기
        G['hold'] = (hf, he)
        G['tauF'] = np.concatenate([np.clip(tau, 0, None), MARGIN * hf]); G['tauE'] = np.concatenate([np.clip(-tau, 0, None), MARGIN * he])
    else:
        G['tauF'] = np.concatenate([np.clip(tau, 0, None), 0 * q_grid]); G['tauE'] = np.concatenate([np.clip(-tau, 0, None), 0 * q_grid])
    return G, yubin


def lengths(G, X):
    """X: 후보 N×6 (위팔 점 xyz, 아래팔 점 xyz, m) → 케이블 길이 N×T, 모멘트암 N×T (= −dL/dq, + 면 굽힘)"""
    P0, P1 = X[:, :3], X[:, 3:]
    A = np.einsum('tij,nj->nti', G['R'], P1) + G['T0'][None]
    dA = np.einsum('tij,nj->nti', G['dR'], P1) + G['dT0'][None]
    d = A - P0[:, None, :]; L = np.linalg.norm(d, axis=-1)
    return L, -np.einsum('nti,nti->nt', d, dA) / L


def arm(G, X): return lengths(G, X)[1]


def seg_dist(P0, P1, Q0, Q1):
    """선분 P0-P1 과 선분 Q0-Q1 사이 최소 거리 (배열 한꺼번에)"""
    d1, d2, r = P1 - P0, Q1 - Q0, P0 - Q0
    a, e, f, c, b = (d1 * d1).sum(-1), (d2 * d2).sum(-1), (d2 * r).sum(-1), (d1 * r).sum(-1), (d1 * d2).sum(-1)
    den = a * e - b * b
    s = np.where(den > 1e-12, np.clip((b * f - c * e) / np.maximum(den, 1e-12), 0, 1), 0.0)
    t = (b * s + f) / e
    s = np.where(t < 0, np.clip(-c / a, 0, 1), np.where(t > 1, np.clip((b - c) / a, 0, 1), s)); t = np.clip(t, 0, 1)
    return np.linalg.norm(P0 + d1 * s[..., None] - (Q0 + d2 * t[..., None]), axis=-1)


def clearance(G, X):
    """케이블과 뼈 중심선(위팔: 어깨~팔꿈치, 아래팔: 팔꿈치~손목) 사이 최소 거리, 후보마다 (m)"""
    P0 = np.broadcast_to(X[:, None, :3], (len(X), len(G['q']), 3))
    A = np.einsum('tij,nj->nti', G['R'], X[:, 3:]) + G['T0'][None]
    dh = seg_dist(P0, A, np.zeros(3), G['joint']['t_p'])
    du = seg_dist(P0, A, G['T0'][None], (G['T0'] + G['R'] @ WRIST)[None])
    return np.minimum(dh, du).min(1)


def evaluate(G, X, cable, S, fit=True):
    """반환: 목적값(제약 위반이면 큰 값), 필요한 최대 힘 N, 쓰는 스트로크 m, 최소 모멘트암 m, 최소 케이블 길이 m"""
    sign = SIGN[cable]; L, r = lengths(G, X); rr = sign * r
    s = L.max(1, keepdims=True) - L; use = s.max(1)                     # 가장 늘어난 자세 기준 줄어든 길이
    tau_c = (G['tauF'] if sign > 0 else G['tauE'])[None]
    BETA = float(os.environ.get('BETA', '1'))                          # 힘–수축 곡선: 낼 수 있는 힘 = F × min(1, (1 − s/S)/β). β 1 = 내 정적 (은서 모델처럼 직선), β 0.25 = 유빈 시간 영역 설계 (임시 곡선)
    F = np.where(tau_c > 1e-6, tau_c / np.maximum(rr, 1e-3) / np.maximum(np.minimum(1.0, (1 - s / S) / BETA), 1e-3), 0.0).max(1)
    v = np.maximum(0, use - S) / S + np.maximum(0, R_MIN - rr.min(1)) / R_MIN + np.maximum(0, CLEAR - clearance(G, X)) / CLEAR
    if fit: v = v + np.maximum(0, S / EPSL - L.min(1)) / (S / EPSL)
    return np.where(v > 0, 1e5 * (1 + 10 * v), F), np.where(v > 0, np.inf, F), use, rr.min(1), L.min(1)


def optimize(G, cable, S, fit=True, seeds=(0, 1)):
    from scipy.optimize import differential_evolution
    best = None
    for sd in seeds:
        r = differential_evolution(lambda X: evaluate(G, X.T, cable, S, fit)[0], list(zip(LO, HI)), seed=sd, maxiter=300, popsize=30,
                                   tol=1e-10, vectorized=True, updating='deferred', polish=False, init='sobol')
        if best is None or r.fun < best.fun: best = r
    _, F, use, rmin, Lmin = evaluate(G, best.x[None], cable, S, fit)
    return dict(x=best.x.tolist(), F_need=float(F[0]), stroke_used=float(use[0]), r_min=float(rmin[0]), L_min=float(Lmin[0]))


def fmt(F): return '불가능' if not np.isfinite(F) else f'{F:,.0f} N'


# ───────────────────────── search ─────────────────────────
def json_path(req): return os.path.join(OUT, f'place_elbow_{req}.json')


def search(req):
    G, yubin = task(req); os.makedirs(OUT, exist_ok=True); t0 = time.time(); nt = G['nt']
    print(f"[{'경로 + 전 범위 버티기 × ' + str(MARGIN) if req == 'range' else '들기 경로만'}] 들기 {G['t'][-1]:.1f} 초, 팔꿈치 {math.degrees(G['q'][:nt].min()):.0f}~{math.degrees(G['q'][:nt].max()):.0f}° (스트로크는 0~150° 전체로 계산)")
    print(f"필요한 토크: 굽힘 쪽 최대 {G['tauF'].max():.1f} N·m, 폄 쪽 최대 {G['tauE'].max():.1f} N·m"
          + (f"  (들기 경로만이면 {G['tau'].max():.1f} / {-G['tau'].min():.1f})" if req == 'range' else '') + '\n')
    res = dict(req=req, tauF=float(G['tauF'].max()), tauE=float(G['tauE'].max()), yubin={}, opt={})
    for n in CABLES:
        Ly, ry = lengths(G, yubin[n][None])
        res['yubin'][n] = dict(x=yubin[n].tolist(), stroke_used=float(Ly.max() - Ly.min()), L_min=float(Ly.min()), r_min=float((SIGN[n] * ry).min()), F={})
        for S in STROKES:
            res['yubin'][n]['F'][f'{S * 1000:.0f}'] = [float(evaluate(G, yubin[n][None], n, S, fit)[1][0]) for fit in (False, True)]
            print(f'  {KOR[n]} 스트로크 {S * 1000:.0f} mm 최적화 중 ...      ', end='\r', flush=True)
            res['opt'][f'{n}|{S * 1000:.0f}'] = optimize(G, n, S, True)
    print(' ' * 60, end='\r')
    best = min(STROKES, key=lambda S: res['opt'][f'elbow_M0_r|{S * 1000:.0f}']['F_need'] + res['opt'][f'elbow_M1_r|{S * 1000:.0f}']['F_need'])
    res['best_S'] = best
    json.dump(res, open(json_path(req), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)

    print('필요한 최대 근육 힘 — 작을수록 좋음. "불가능" = 스트로크가 모자라거나 근육이 들어갈 자리가 없음')
    print(f"{'케이블':8s} {'스트로크 S':>9s} {'근육 길이':>8s} | {'유빈 배치':>10s} {'(자리 조건 빼면)':>14s} | {'SMA 반영 배치':>12s}  M3(450 N) 몇 개?")
    for n in CABLES:
        y = res['yubin'][n]
        for S in STROKES:
            k = f'{S * 1000:.0f}'; o = res['opt'][f'{n}|{k}']
            nm3 = f"{math.ceil(o['F_need'] / F_M3)}개" if np.isfinite(o['F_need']) and S <= 0.06 else ''
            print(f"{KOR[n]:8s} {k:>7s} mm {S / EPSL * 1000:6.0f} mm | {fmt(y['F'][k][1]):>10s} {'(' + fmt(y['F'][k][0]) + ')':>14s} | {fmt(o['F_need']):>12s}  {nm3}")
        print(f"{'':8s} 유빈 배치: 케이블 {y['L_min'] * 1000:.0f}~{(y['L_min'] + y['stroke_used']) * 1000:.0f} mm, 0~150° 에 줄어드는 길이 {y['stroke_used'] * 1000:.0f} mm, 최소 모멘트암 {y['r_min'] * 1000:.1f} mm")
    for S in sorted({S_MAIN, best}):
        for n in CABLES:
            o = res['opt'][f'{n}|{S * 1000:.0f}']; x = np.array(o['x']) * 1000
            print(f"SMA 반영 배치 S {S * 1000:.0f} mm {KOR[n]}: 위팔 점 ({x[0]:.0f}, {x[1]:.0f}, {x[2]:.0f}) 아래팔 점 ({x[3]:.0f}, {x[4]:.0f}, {x[5]:.0f}) mm, "
                  f"케이블 {o['L_min'] * 1000:.0f} mm~, 줄어드는 길이 {o['stroke_used'] * 1000:.0f} mm, 최소 모멘트암 {o['r_min'] * 1000:.1f} mm")
    print(f"\n두 근육 합이 가장 작은 스트로크: {best * 1000:.0f} mm.  걸린 시간 {time.time() - t0:.0f} 초. 결과: {json_path(req)}")
    figure(G, yubin, res)


def figure(G, yubin, res):
    plt = font(); fig = plt.figure(figsize=(16, 4.9)); gs = fig.add_gridspec(1, 3, width_ratios=[1.05, 1.0, 1.0])
    C = {'elbow_M0_r': '#C0392B', 'elbow_M1_r': '#2B6CB0'}; nt = G['nt']; SB = res['best_S']
    # (가) 스트로크 vs 필요한 힘
    ax = fig.add_subplot(gs[0]); Sx = np.array(STROKES) * 1000
    ax.add_patch(plt.Rectangle((25, 150), 35, 300, fc='#9AA5B5', alpha=.25, ec='#9AA5B5')); ax.text(27, 470, 'M3 규격\n25~60 mm\n150~450 N', fontsize=8, color='#4A5568', va='bottom')
    ax.add_patch(plt.Rectangle((40, 1000), 80, 4000, fc='#E8925E', alpha=.15, ec='#E8925E')); ax.text(100, 1080, 'M1 규격\n40~120 mm\n1,000~5,000 N', fontsize=8, color='#B8501E', va='bottom')
    for n in CABLES:
        o = [res['opt'][f'{n}|{s:.0f}']['F_need'] for s in Sx]; ax.plot(Sx, o, 'o-', c=C[n], lw=2.4, ms=7, label=f'{KOR[n]} SMA 반영 배치')
        y = [res['yubin'][n]['F'][f'{s:.0f}'][0] for s in Sx]
        for s, v in zip(Sx, y):
            if np.isfinite(v): ax.plot(s, min(v, 4800), 's', c=C[n], ms=9, mfc='white', mew=2); ax.annotate(f'유빈 {v:,.0f} N', (s, min(v, 4800)), textcoords='offset points', xytext=(8, -3), fontsize=8, color=C[n])
            else: ax.plot(s, 4600 if n == 'elbow_M0_r' else 4350, 'x', c=C[n], ms=9, mew=2.2)
    ax.text(45, 4780, '× = 유빈 배치로는 스트로크 부족', fontsize=8)
    ax.set(xlabel='근육 스트로크 S (mm)  ·  근육 길이 = S / 0.33', ylabel='필요한 최대 근육 힘 (N)', title='(가) 근육이 얼마나 커야 하나', xlim=(22, 122), ylim=(0, 5000))
    ax.plot([], [], 's', c='#555', mfc='white', mew=2, label='유빈 배치 (□, 자리 조건 빼고)')
    ax.legend(fontsize=7.5, loc='center right'); ax.grid(alpha=.3)
    # (나) 모멘트암 + 필요한 토크
    ax = fig.add_subplot(gs[1]); qg = np.degrees(G['q'][nt:])
    for n in CABLES:
        ry = SIGN[n] * arm(G, np.array(res['yubin'][n]['x'])[None])[0][nt:]; ro = SIGN[n] * arm(G, np.array(res['opt'][f'{n}|{SB * 1000:.0f}']['x'])[None])[0][nt:]
        ax.plot(qg, ry * 1000, '--', c=C[n], lw=1.6, label=f'{KOR[n]} 유빈'); ax.plot(qg, ro * 1000, '-', c=C[n], lw=2.2, label=f'{KOR[n]} SMA 반영')
    ax2 = ax.twinx(); ax2.plot(np.degrees(G['q'][:nt]), G['tau'], '.', c='#888', ms=2)
    if 'hold' in G: ax2.plot(qg, MARGIN * G['hold'][0], ':', c='#C0392B', lw=1.2); ax2.plot(qg, -MARGIN * G['hold'][1], ':', c='#2B6CB0', lw=1.2)
    ax2.set_ylabel('필요한 토크 N·m (점 = 들기 경로, 점선 = 버티기)', color='#666', fontsize=8); ax2.axhline(0, c='#bbb', lw=.8)
    ax.set(xlabel='팔꿈치 각도 (°)', ylabel='모멘트암 (mm, 방향 맞춤)', title=f'(나) 모멘트암 — SMA 반영은 스트로크 {SB * 1000:.0f} mm 배치'); ax.legend(fontsize=7.5, loc='upper left'); ax.grid(alpha=.3)
    # (다) 옆에서 본 배치 (팔꿈치 90°)
    ax = fig.add_subplot(gs[2]); J = G['joint']
    c, s = math.cos(math.radians(90)), math.sin(math.radians(90)); Rz = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    R = J['R_p'] @ Rz @ J['R_c'].T; T0 = J['t_p'] - R @ J['t_c']
    to_h = lambda p: T0 + R @ p
    el = J['t_p']; wr = to_h(WRIST)
    ax.plot([0, el[0]], [0, el[1]], '-', c='#CBD5E0', lw=14, solid_capstyle='round'); ax.plot([el[0], wr[0]], [el[1], wr[1]], '-', c='#CBD5E0', lw=11, solid_capstyle='round')
    ax.add_patch(plt.Rectangle((BOX['humerus_r'][0][0], BOX['humerus_r'][1][0]), np.ptp(BOX['humerus_r'][0]), np.ptp(BOX['humerus_r'][1]), fill=False, ec='#999', ls=':'))
    ax.text(0.012, 0.012, '어깨', fontsize=8); ax.text(el[0] - 0.005, el[1] + 0.03, '팔꿈치', fontsize=8, ha='right'); ax.text(wr[0], wr[1] + 0.012, '손목', fontsize=8, ha='center')
    for n in CABLES:
        for key, ls, lw, lab in [('yubin', '--', 1.6, '유빈'), ('opt', '-', 2.4, 'SMA 반영')]:
            x = np.array(res['yubin'][n]['x'] if key == 'yubin' else res['opt'][f'{n}|{SB * 1000:.0f}']['x'])
            a, b = x[:3], to_h(x[3:]); ax.plot([a[0], b[0]], [a[1], b[1]], ls, c=C[n], lw=lw, marker='o', ms=5, label=f'{KOR[n]} {lab}')
    ax.set_aspect('equal'); ax.set(xlabel='앞 ← x (m) →   (팔꿈치 90°, 옆에서 봄)', ylabel='y (m)', title='(다) 옆에서 본 케이블 (점선 = 위팔 후보 범위)')
    ax.legend(fontsize=7.5, loc='lower left'); ax.grid(alpha=.2)
    ttl = '들기 경로 + 0~150° 어디서든 물병 버티기(×1.2)' if res['req'] == 'range' else '들기 경로만 (1차 시도)'
    fig.suptitle(f'SMA 성질(스트로크 · 줄수록 힘 빠짐 · 근육 들어갈 자리 · 뼈와 거리)을 넣고 팔꿈치 부착점 다시 고르기 — 조건: {ttl}', fontsize=11.5, fontweight='bold')
    plt.tight_layout(); f = os.path.join(OUT, f"fig_place_elbow_{res['req']}.png"); plt.savefig(f, dpi=140); print('그림:', f)


# ───────────────────────── lift / view (MuJoCo) ─────────────────────────


def make_xml(res, S):
    """SMA 반영 배치(스트로크 S 용)로 부착점만 바꾼 MuJoCo xml"""
    xml = open(os.path.join(HERE, 'finalA_elbow_r.xml'), encoding='utf-8').read()
    for n in CABLES:
        x = res['opt'][f'{n}|{S * 1000:.0f}']['x']
        xml = re.sub(rf'(<site name="{n}_p0" pos=")[^"]*', lambda m: m.group(1) + ' '.join(f'{v:.10g}' for v in x[:3]), xml)
        xml = re.sub(rf'(<site name="{n}_p1" pos=")[^"]*', lambda m: m.group(1) + ' '.join(f'{v:.10g}' for v in x[3:]), xml)
    path = os.path.join(HERE, f'finalA_elbow_r_sma_place_{S * 1000:.0f}.xml'); open(path, 'w', encoding='utf-8').write(xml)
    return path


def lift_ns():
    src = open(os.path.join(HERE, 'lift_finalA.py'), encoding='utf-8').read()
    ns = {'__file__': os.path.join(HERE, 'lift_finalA.py'), '__name__': 'lift'}; exec(src[:src.index('mode = sys.argv')], ns)
    sm = ns['sm']

    def use(xml, F, S):                                   # 배치(xml) + 근육(힘 F, 스트로크 S — 숫자 하나 또는 케이블별 dict)
        def params(module, cable):
            Sc = S[cable] if isinstance(S, dict) else S
            p = copy.deepcopy(sm.M3); p['epsL'] = EPSL; p['force_model'] = 't31'
            p['FMAX'] = F[cable]; p['L0'] = Sc / EPSL; p['A_wire'] = p['FMAX'] / 400e6; p['K_mod'] = p['FMAX'] / EPSL
            return p
        ns['params'] = params; ns['XML'] = xml
    ns['use'] = use
    return ns


def check_lengths(xml_path, X):
    """정적 계산(파이썬) 과 MuJoCo 케이블 길이가 같은지 확인"""
    import mujoco
    G, _ = task('path'); m = mujoco.MjModel.from_xml_string(open(xml_path, encoding='utf-8').read()); d = mujoco.MjData(m)
    qa = m.jnt_qposadr[m.joint('elbow_flex_r').id]; err = 0.0; J = G['joint']
    for qd in range(0, 151, 30):
        d.qpos[qa] = math.radians(qd); mujoco.mj_forward(m, d)
        c, s = math.cos(math.radians(qd)), math.sin(math.radians(qd)); Rz = np.array([[[c, -s, 0], [s, c, 0], [0, 0, 1]]])
        R = J['R_p'] @ Rz @ J['R_c'].T
        G1 = dict(R=R, T0=J['t_p'] - R @ J['t_c'], dR=np.zeros_like(R), dT0=np.zeros((1, 3)))
        for n in CABLES: err = max(err, abs(lengths(G1, np.array(X[n])[None])[0][0, 0] - d.ten_length[m.tendon(n).id]))
    return err


def cases(res):
    S_BEST = res['best_S']
    small = {n: round(res['opt'][f'{n}|{S_BEST * 1000:.0f}']['F_need'] * 1.1, -1) for n in CABLES}   # 가장 좋은 배치가 필요한 힘 + 10 %
    yub = os.path.join(HERE, 'finalA_elbow_r.xml')
    return [('① 유빈 배치 · M1 (86 mm, 설계 힘)', yub, DESIGN_F, S_MAIN, '#9AA5B5'),
            ('② SMA 반영 배치 · 같은 M1', make_xml(res, S_MAIN), DESIGN_F, S_MAIN, '#2E7D4F'),
            (f"③ SMA 반영 배치 · 계산된 필요 힘 +10 % ({S_BEST * 1000:.0f} mm, {small['elbow_M0_r']:,.0f}/{small['elbow_M1_r']:,.0f} N)", make_xml(res, S_BEST), small, S_BEST, '#D69E2E')]


def lift(req):
    res = json.load(open(json_path(req), encoding='utf-8')); ns = lift_ns(); REF = ns['prepare'](); C = cases(res)
    for S in sorted({S_MAIN, res['best_S']}):
        err = check_lengths(make_xml(res, S), {n: res['opt'][f'{n}|{S * 1000:.0f}']['x'] for n in CABLES})
        print(f"스트로크 {S * 1000:.0f} mm 배치: 정적 계산 vs MuJoCo 케이블 길이 차이 {err * 1000:.3f} mm (0 이면 같은 배치)")
    print(f"\n{'경우':44s} {'평균 오차':>8s} {'최대 오차':>8s} {'10° 안':>7s} {'최고 온도':>8s} {'최대 힘 굽힘/폄':>15s}")
    rows, out = {}, {}
    for lab, xml, F, S, _ in C:
        ns['use'](xml, F, S); R = ns['simulate']('custom', REF); rows[lab] = R; e = R[:, 2] - R[:, 1]
        out[lab] = dict(rms=float(np.sqrt(np.mean(e ** 2))), maxerr=float(np.abs(e).max()), within10=float(np.mean(np.abs(e) < 10)),
                        Tmax=float(R[:, 3:5].max()), Fmax=[float(R[:, 5].max()), float(R[:, 6].max())], F=F, S=S)
        print(f"{lab:44s} {out[lab]['rms']:7.1f}° {out[lab]['maxerr']:7.1f}° {out[lab]['within10'] * 100:5.0f} % {out[lab]['Tmax']:7.1f}°C {out[lab]['Fmax'][0]:7.0f} / {out[lab]['Fmax'][1]:4.0f} N")
    json.dump(out, open(os.path.join(OUT, f'place_lift_{req}.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    plt = font(); fig, ax = plt.subplots(2, 1, figsize=(11, 7.5), sharex=True, gridspec_kw=dict(height_ratios=[2, 1]))
    R0 = rows[C[0][0]]; ax[0].plot(R0[:, 0], R0[:, 1], 'k', lw=3, label='목표 (유빈 물병 들기)')
    for lab, _, _, _, c in C:
        R = rows[lab]; ax[0].plot(R[:, 0], R[:, 2], c=c, lw=2.2, label=f"{lab}  — 평균 오차 {out[lab]['rms']:.1f}°")
        ax[1].plot(R[:, 0], R[:, 2] - R[:, 1], c=c, lw=1.8)
    ttl = '경로 + 전 범위 버티기 조건' if req == 'range' else '들기 경로만 조건 (1차 시도)'
    ax[0].set(ylabel='팔꿈치 각도 (°)', title=f'물병 3.79 kg 들기 — 부착점만 바꿨을 때 (근육 모델·제어는 같음) · 배치 최적화: {ttl}'); ax[0].legend(fontsize=8.5, loc='upper left')
    ax[1].axhline(0, c='k', lw=1); ax[1].set(ylabel='오차 (°)', xlabel='시간 (초)')
    for a in ax: a.grid(alpha=.3)
    plt.tight_layout(); f = os.path.join(OUT, f'fig_place_lift_{req}.png'); plt.savefig(f, dpi=140); print('\n그림:', f)


def playback(ns, REF, lab):
    """먼저 끝까지 계산한 뒤 실제 속도(6.4 초)로 다시 보여 줌 → 녹화용. 옆에서 본 화면, 끝나면 1 초 쉬고 반복"""
    import mujoco, mujoco.viewer
    print(f'{lab} — 계산 중 (몇 초) ...'); R = ns['simulate']('custom', REF); e = R[:, 2] - R[:, 1]
    S = ns['setup']('custom', REF); m = S['m']; d = mujoco.MjData(m); qa = S['qa']
    t, q = R[:, 0], np.radians(R[:, 2])
    print(f'{lab} · 평균 오차 {np.sqrt(np.mean(e ** 2)):.1f}° — 실제 속도로 반복 재생. 창을 닫으면 끝 (녹화: 창 클릭 후 Win+Alt+R)')
    try: v_ctx = mujoco.viewer.launch_passive(m, d, show_left_ui=False, show_right_ui=False)
    except TypeError: v_ctx = mujoco.viewer.launch_passive(m, d)
    with v_ctx as v:
        v.cam.lookat[:] = [0.05, 0.0, 0.15]; v.cam.distance = 0.9; v.cam.azimuth = 90; v.cam.elevation = -10
        t0 = time.time(); nxt = 0.0
        while v.is_running():
            tt = (time.time() - t0) % (t[-1] + 1.0)
            k = min(int(np.searchsorted(t, tt)), len(t) - 1)
            d.qpos[qa] = q[k]; mujoco.mj_forward(m, d); v.sync()
            if tt < nxt - 1: nxt = 0.0
            if tt >= nxt and tt <= t[-1]:
                print(f"t={t[k]:4.1f}s 목표 {R[k, 1]:6.1f}° 지금 {R[k, 2]:6.1f}°"); nxt += 0.5
            time.sleep(1 / 60)


def view(req, k):
    res = json.load(open(json_path(req), encoding='utf-8')); ns = lift_ns(); REF = ns['prepare']()
    lab, xml, F, S0, _ = cases(res)[k]; ns['use'](xml, F, S0); playback(ns, REF, lab)


def view_opt():
    """opt_dyn.py (CMA-ES) 로 찾고 validate_opt.py 로 30번 확인한 최적화 배치 — 케이블마다 들어가는 가장 긴 근육"""
    p = os.path.join(OUT, 'opt_validate.json')
    if not os.path.exists(p): sys.exit('out\\opt_validate.json 이 없어요 (zip 의 out 폴더째 복사하거나 opt_dyn.py → validate_opt.py 실행)')
    V = json.load(open(p, encoding='utf-8')); o = V['최적화 1']; ns = lift_ns(); REF = ns['prepare']()
    xml = xml_from({n: np.array(o['x'][n]) for n in CABLES}, os.path.join(OUT, 'finalA_elbow_r_opt_best.xml'))
    ns['use'](xml, DESIGN_F, o['S'])
    playback(ns, REF, '최적화 배치 (근육 ' + ' / '.join(f"{o['S'][n] * 1000:.0f}" for n in CABLES) + ' mm 스트로크)')


def view_dyn():
    """dyn 으로 고른 배치 (발표 18장) — 70 mm 근육 · M1 설계 힘"""
    p = os.path.join(OUT, 'place_dyn.json')
    if not os.path.exists(p): sys.exit('먼저  python place_elbow.py dyn')
    best = json.load(open(p, encoding='utf-8'))['best']; ns = lift_ns(); REF = ns['prepare']()
    xml = xml_from({n: np.array(best['x'][n]) for n in CABLES}, os.path.join(OUT, 'finalA_elbow_r_dyn_best.xml'))
    ns['use'](xml, DESIGN_F, S_DYN)
    playback(ns, REF, f"2단계로 고른 배치 (굽힘 후보 {best['i'] + 1} × 폄 후보 {best['j'] + 1}, {S_DYN * 1000:.0f} mm 근육)")


# ───────────────────────── dyn: 2단계 (정적으로 후보 추리기 → MuJoCo 로 골라내기) ─────────────────────────
def xml_from(X, path):
    xml = open(os.path.join(HERE, 'finalA_elbow_r.xml'), encoding='utf-8').read()
    for n in CABLES:
        x = X[n]
        xml = re.sub(rf'(<site name="{n}_p0" pos=")[^"]*', lambda m: m.group(1) + ' '.join(f'{v:.10g}' for v in x[:3]), xml)
        xml = re.sub(rf'(<site name="{n}_p1" pos=")[^"]*', lambda m: m.group(1) + ' '.join(f'{v:.10g}' for v in x[3:]), xml)
    open(path, 'w', encoding='utf-8').write(xml); return path


def candidates(G, cable, S, K, n=120000, seed=0):
    """1단계: 무작위 후보 n 개 중 기하 조건(스트로크 · 모멘트암 · 뼈 거리 · 근육 자리)을 지키는 것만 남기고,
    정적 1등(필요 힘 최소)에서 시작해 서로 가장 멀리 떨어진 K 개를 고름 (모양이 다른 후보들 — 정적 점수로 미리 거르지 않음)"""
    rng = np.random.default_rng(seed); keep = []
    centers = [np.array(json.load(open(json_path(r), encoding='utf-8'))['opt'][f'{cable}|{S * 1000:.0f}']['x'])
               for r in ('range', 'path') if os.path.exists(json_path(r))]          # search 결과 근처를 더 촘촘히 (자리 조건 때문에 좁음)
    for k in range(0, n, 4000):
        X = LO + rng.random((4000, 6)) * (HI - LO)
        if centers:
            c = np.array(centers)[rng.integers(len(centers), size=3000)]
            X[:3000] = np.clip(c + rng.normal(0, 0.02, (3000, 6)), LO, HI)
        _, F, _, _, _ = evaluate(G, X, cable, S, True)
        ok = np.isfinite(F); keep += list(zip(F[ok], X[ok]))
    keep.sort(key=lambda a: a[0]); top = np.array([x for _, x in keep]); Ft = np.array([f for f, _ in keep])   # 정적 1등 + 조건을 지키는 전체에서 고르게
    pick = [0]; d = np.linalg.norm(top - top[0], axis=1)
    while len(pick) < min(K, len(top)):
        j = int(np.argmax(d)); pick.append(j); d = np.minimum(d, np.linalg.norm(top - top[j], axis=1))
    return [(float(Ft[j]), top[j]) for j in pick], len(keep)


S_DYN = 0.070      # 86 mm(근육 261 mm)는 설계 상자 안에 일자로 넣을 자리가 거의 없음(모서리 한 점) → 70 mm(212 mm) 로 비교


def dyn(K=7):
    """M1 설계 힘(2,327 / 1,835 N) + 스트로크 70 mm(근육 212 mm, 자리에 들어감) 근육으로,
    정적 조건을 지키는 여러 배치를 MuJoCo 물병 들기로 직접 비교. 유빈 배치는 원래대로 86 mm 근육(자리 조건은 못 지킴)"""
    G, yubin = task('path'); S = S_DYN; t0 = time.time(); os.makedirs(OUT, exist_ok=True)
    print(f'1단계: 기하 조건을 지키는 후보 찾기 (스트로크 {S * 1000:.0f} mm, 근육 길이 {S / EPSL * 1000:.0f} mm) ...')
    C = {}
    for n in CABLES:
        C[n], nf = candidates(G, n, S, K); print(f"  {KOR[n]}: 조건을 지키는 후보 {nf:,} 개 → 모양이 다른 {len(C[n])} 개 고름 (필요 힘 {C[n][0][0]:,.0f} ~ {max(f for f, _ in C[n]):,.0f} N)")
    ns = lift_ns(); REF = ns['prepare'](); tmp = os.path.join(OUT, '_tmp_place.xml'); rows = []

    def run(X, S_mus=S):
        ns['use'](xml_from(X, tmp), DESIGN_F, S_mus); R = ns['simulate']('custom', REF); e = R[:, 2] - R[:, 1]
        return R, float(np.sqrt(np.mean(e ** 2))), float(np.abs(e).max())
    Ry, rms_y, mx_y = run({n: yubin[n] for n in CABLES}, S_MAIN)
    print(f'\n2단계: MuJoCo 물병 들기 {len(C[CABLES[0]]) * len(C[CABLES[1]])} 번 (근육 = M1 설계 힘, 스트로크 {S * 1000:.0f} mm) ...  유빈 배치(86 mm 근육): 평균 오차 {rms_y:.1f}°')
    for i, (Ff, xf) in enumerate(C['elbow_M0_r']):
        for j, (Fe, xe) in enumerate(C['elbow_M1_r']):
            R, rms, mx = run({'elbow_M0_r': xf, 'elbow_M1_r': xe})
            rows.append(dict(i=i, j=j, F=[Ff, Fe], rms=rms, maxerr=mx, x={'elbow_M0_r': xf.tolist(), 'elbow_M1_r': xe.tolist()}))
            print(f'  굽힘 후보 {i + 1} × 폄 후보 {j + 1}: 정적 필요 힘 {Ff:5,.0f} / {Fe:5,.0f} N → 평균 오차 {rms:5.1f}°, 최대 {mx:5.1f}°', flush=True)
    rows.sort(key=lambda r: r['rms'])
    # 3단계: 이 시뮬레이션은 켜짐/꺼짐 스위치(냉각 전환 · 히스테리시스 갈래 · 잠열 보정) 때문에 부착점을 0.01 mm 만 바꿔도
    #        평균 오차가 1~2° 달라짐 → 상위 5 개와 유빈 배치를 0.01 mm 씩 흔들어 8 번씩 다시 돌려 평균으로 비교
    print(f'\n3단계: 상위 5 개 + 유빈 배치를 부착점 0.01 mm 흔들어 8 번씩 다시 (시뮬레이션 흔들림 확인) ...')
    rng = np.random.default_rng(1)

    def jitter(X, S_mus):
        v = []
        for _ in range(8):
            v.append(run({n: np.array(X[n]) + rng.normal(0, 1e-5, 6) for n in CABLES}, S_mus)[1])
        return float(np.mean(v)), float(np.min(v)), float(np.max(v))
    ym = jitter({n: yubin[n] for n in CABLES}, S_MAIN); rms_y = ym[0]
    print(f'  유빈 배치 (86 mm 근육): 평균 {ym[0]:.1f}°  ({ym[1]:.1f} ~ {ym[2]:.1f}°)')
    for r in rows[:5]:
        r['jit'] = jitter(r['x'], S); print(f"  굽힘 후보 {r['i'] + 1} × 폄 후보 {r['j'] + 1}: 평균 {r['jit'][0]:.1f}°  ({r['jit'][1]:.1f} ~ {r['jit'][2]:.1f}°)")
    top = sorted(rows[:5], key=lambda r: r['jit'][0]); best = top[0]; best['rms'] = best['jit'][0]
    Rb, _, _ = run({n: np.array(best['x'][n]) for n in CABLES})
    gy = {n: float(evaluate(G, yubin[n][None], n, S, False)[1][0]) for n in CABLES}
    json.dump(dict(yubin=dict(rms=rms_y, jit=ym, maxerr=mx_y, F_nofit=gy), best=best, rows=rows), open(os.path.join(OUT, 'place_dyn.json'), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
    xb = {n: np.array(best['x'][n]) * 1000 for n in CABLES}
    print(f"\n가장 잘 따라간 배치: 굽힘 후보 {best['i'] + 1} × 폄 후보 {best['j'] + 1} → 평균 오차 {best['jit'][0]:.1f}° ({best['jit'][1]:.1f}~{best['jit'][2]:.1f}), "
          f"유빈 {ym[0]:.1f}° ({ym[1]:.1f}~{ym[2]:.1f}). 범위가 겹치면 '비슷하다' 로 읽기")
    for n in CABLES: print(f"  {KOR[n]}: 위팔 점 ({xb[n][0]:.0f}, {xb[n][1]:.0f}, {xb[n][2]:.0f}) 아래팔 점 ({xb[n][3]:.0f}, {xb[n][4]:.0f}, {xb[n][5]:.0f}) mm")
    rr = np.array([r['rms'] if 'jit' not in r else r['jit'][0] for r in rows]); ff = np.array([max(r['F'][0] / DESIGN_F['elbow_M0_r'], r['F'][1] / DESIGN_F['elbow_M1_r']) for r in rows])
    print(f"정적 필요 힘과 실제 오차의 상관: {np.corrcoef(ff, rr)[0, 1]:+.2f} (1 에 가까우면 정적 계산만으로 순위를 맞출 수 있음)")
    print(f'걸린 시간 {time.time() - t0:.0f} 초')
    plt = font(); fig, ax = plt.subplots(1, 2, figsize=(15, 4.8), gridspec_kw=dict(width_ratios=[1, 1.6]))
    ax[0].scatter(ff, rr, s=40, c='#2E7D4F', alpha=.75, label='SMA 조건을 지키는 배치 (굽힘 × 폄 조합)')
    ax[0].axhspan(ym[1], ym[2], color='#9AA5B5', alpha=.3, label=f'유빈 배치 {ym[1]:.1f}~{ym[2]:.1f}° (흔들림 범위, 86 mm 근육 · 자리 조건 못 지킴)')
    kb = rows.index(best)
    ax[0].scatter([ff[kb]], [rr[kb]], s=160, facecolors='none', edgecolors='#D69E2E', linewidths=2.5, label=f"고른 배치 (흔들어 평균 {best['jit'][0]:.1f}°)")
    ax[0].set(xlabel='정적 계산: 필요한 힘 ÷ M1 설계 힘 (큰 쪽)', ylabel='MuJoCo 물병 들기 평균 오차 (°)', title='(가) 정적 점수 vs 실제로 따라간 정도'); ax[0].grid(alpha=.3); ax[0].legend(fontsize=8)
    ax[1].plot(Ry[:, 0], Ry[:, 1], 'k', lw=3, label='목표 (유빈 들기)')
    ax[1].plot(Ry[:, 0], Ry[:, 2], c='#9AA5B5', lw=2.2, label=f'유빈 배치 (86 mm 근육) — 평균 오차 {ym[0]:.1f}° ({ym[1]:.1f}~{ym[2]:.1f})')
    ax[1].plot(Rb[:, 0], Rb[:, 2], c='#D69E2E', lw=2.2, label=f"고른 배치 ({S * 1000:.0f} mm 근육) — 평균 오차 {best['jit'][0]:.1f}° ({best['jit'][1]:.1f}~{best['jit'][2]:.1f})")
    ax[1].set(xlabel='시간 (초)', ylabel='팔꿈치 각도 (°)', title=f'(나) 물병 3.79 kg 들기 — 같은 힘(M1 설계), 유빈 86 mm 근육 vs 고른 배치 {S * 1000:.0f} mm 근육'); ax[1].grid(alpha=.3); ax[1].legend(fontsize=9)
    fig.suptitle('2단계 배치 고르기 — ① 기하 조건으로 후보 추리기 → ② MuJoCo(열 지연 · 히스테리시스 포함)로 직접 돌려 골라내기', fontsize=12, fontweight='bold')
    plt.tight_layout(); f = os.path.join(OUT, 'fig_place_dyn.png'); plt.savefig(f, dpi=140); print('그림:', f)
    try: os.remove(tmp)
    except OSError: pass


if __name__ == '__main__':
    mode = sys.argv[1] if len(sys.argv) > 1 else 'search'
    args = sys.argv[2:]; req = 'path' if 'path' in args else 'range'
    if mode == 'search': search(req)
    elif mode == 'dyn': dyn()
    elif mode == 'view' and 'dyn' in args: view_dyn()
    elif mode == 'view' and 'opt' in args: view_opt()
    elif mode in ('lift', 'view'):
        if not os.path.exists(json_path(req)): sys.exit(f"먼저  python place_elbow.py search{' path' if req == 'path' else ''}")
        if mode == 'lift': lift(req)
        else: view(req, int(next((a for a in args if a.isdigit()), '3')) - 1)
    else: print(__doc__)
