"""설계 방향 8개 → 피험자3 뼈대 위의 케이블 경로 (부착점 최적화).

케이블 i 의 모멘트암 벡터 r(q) = −∂L/∂q (tendon excursion, 중앙차분) 가 걷기·스쿼트·앉았다서기의
모든 자세에서 설계 방향 v_i 와 나란하도록 부착점(시작·경유·끝)을 고른다. 길이 계산은 뼈 변환을
미리 저장해 두고 numpy 로만 하므로(마이크로초) 전역 방법(differential evolution)을 쓸 수 있다.

    python report/path_design.py            → report/newmus_paths_k8.json, report/moco/newmus_paths_subject3_model.osim,
                                              report/fig_path_design.png
    python report/path_design.py --space arm_r --tasks report/upper/tasks.csv --design report/newmus_V_design_arm_r_k6.npz \
                                 --model data/moore2015/subject3/subject3_scaled.osim
                                            → report/newmus_paths_arm_r_k6.json, report/moco/newmus_paths_arm_r_k6_model.osim
"""
import argparse, csv, glob, json, os, sys, time
import numpy as np
from scipy.optimize import differential_evolution, minimize
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import addb_io as io                                   # noqa: E402
import design_space as ds                              # noqa: E402
import timedomain as td                                # noqa: E402
LAB = os.path.join(os.path.dirname(ROOT), "LabValidation_withoutVideos")
DELTA = 0.01                                          # rad, 중앙차분
STROKE = 0.40                                         # 모듈 스트로크 (변형률)
STROKE_MAX = None                                     # 모듈 스트로크 상한 m (--stroke-max, 과제 계획서 모듈표: M1 0.12, M2·M6 0.08, M3 0.06). None 이면 끔 (2026-11-04)
SMA_EPS = None                                        # SMA 묶음 제약 (L_max − L_min)/L_min ≤ SMA_EPS (--sma-eps, park2019 수축률 0.5 ÷ 여유). None 이면 끔 (2026-11-04)
MARGIN = 1.1                                          # 용량 여유 10 %
SP = COORDS = EXTRA = BODIES = BOX = F_MODULE = None
WRAP = None                                           # 감는 구 (design_space 의 wrap), 없으면 None
# 모듈 먼저 (--modules, 2026-11-05): 케이블마다 미리 배정한 등급의 힘·스트로크 상한. None 이면 부위 전체 상한(F_MODULE, STROKE_MAX) — 기존 동작
#   dict(fcap=(k,) 설계 힘 상한 N, smax=(k,) 길이 변화 상한 m, fr2=Frames(길이 검사 자세) 또는 None)
MODCAPS = None
W_CAP_F, W_CAP_S = 5.0, 50.0                          # 케이블별 상한 벌점 가중치 = 부위 전체 상한(F_MODULE 5, STROKE_MAX 50)과 같게. 두 방식의 차이가 "상한이 케이블마다 다르다"뿐이도록
#   (가중치 1e4 벽으로 해 봤더니 2단계 Powell 의 상자 전체 선탐색이 시작점보다 나쁜 점을 돌려줘 2단계가 아무것도 못 했다 — 2026-11-05 v2)
#   벌점이 부드러워 생기는 작은 초과는 2단계 뒤 repair_caps 가 그 케이블만 국소 탐색해 없앤다
#   (설계 상한을 조여 2단계를 다시 돌리는 것은 안 됐다: 몸통 M0 힘 +2 % 가 세 번 조여도 그대로 — 가중치 5 에서 2 % 초과 벌점이 0.01 이라 Powell 이 안 움직인다, v3)
CAP_MARGIN = 0.01                                     # 첫 회 설계 상한 = 등급 상한 × (1 − 1 %). 모듈 json 의 design_scale_force/_stroke 가 있으면 그 값
DEGEN_MF = 1e7                                        # 모듈 먼저에서 퇴화 해(구 안·모멘트암 0) 비용. 기본 1e3 이면 상한을 못 지키는 케이블의 벌점이 1e3 을 넘어 DE 가 퇴화 해로 도망간다 (어깨 M1)
# 시간 영역 (--force-mode time, 2026-11-19): 힘을 "용량 ÷ 최소 모멘트암" 대신 프레임마다 그 순간의 장력 T(t) 로, 모듈은 같은 순간의 (수축, 장력) 쌍이
# 힘–수축 곡선 안에 드는 묶음 수로 정한다 (report/timedomain.py). 목적 = 모듈 일 용량(묶음 × 최대 힘 × 최대 스트로크) 최소.
FORCE_MODE = "amin"
TDS = None                                            # dict(curves, grades, nb, w_work, w_h, U (N,k) 프레임별 토크 몫, Wref (k,) 정규화 J)
BRACKET = None                                        # --bracket H: dict(H, bones={body index: (p0, p1, R)}, clear={body index: 거리}) — 부착점은 뼈 관 표면 ~ +H 안
W_BRACKET_WALL = 1e4                                  # 브래킷 범위를 벗어난 거리(m)² 벌점


def configure(space="lower"):
    """설계 공간을 고른다 — 모듈 전역 COORDS/EXTRA/BODIES/BOX/F_MODULE 이 바뀐다. 다른 스크립트는 import 뒤 이걸 부른다."""
    global SP, COORDS, EXTRA, BODIES, BOX, F_MODULE
    SP = ds.get(space); COORDS = SP["coords"]; EXTRA = SP["extra"]; BODIES = SP["bodies"]; BOX = SP["box"]; F_MODULE = SP["f_module"]
    global WRAP; WRAP = SP.get("wrap")
    return SP


configure("lower")
_FONTS = font_manager.findSystemFonts()
for _key in ("applesdgothicneo", "malgun.ttf"):                # 맥 → 윈도(맑은 고딕) 순서 (2026-11-05)
    _f = next((f for f in _FONTS if _key in os.path.basename(f).lower()), None)
    if _f: plt.rcParams["font.family"] = font_manager.FontProperties(fname=_f).get_name(); break
plt.rcParams["axes.unicode_minus"] = False
SUBJECT3_BODY = (63.5, 1.69)                          # LabValidation 피험자3 질량 kg·키 m (design_directions 기본 REF 와 같음). poses() 의 --ref 환산에 쓴다


def poses(n_walk=10, n_sq=60, seed=0, ref=None):
    """자세 표본 Q (N × 5: 4 설계 좌표 + hip_rotation, rad) 와 같은 프레임의 ID 토크 T (N × 4, Nm).
    걷기 6 시행 + 스쿼트 (모두 AddB 변환 짝이라 IK·ID 가 일관). 스쿼트는 무릎이 가장 굽은 프레임을 반드시 넣는다.
    ref=(mass, height) 를 주면 토크를 그 몸 크기로 맞춘다 (design_directions 와 같은 질량×키 비율). 2026-11-05 전에는 ref 를 안 받아
    `--ref 40 1.60` 을 줘도 하체 2단계가 피험자3 크기(63.5 kg·1.69 m) 토크를 맞췄다."""
    rng = np.random.default_rng(seed); Q, T = [], []
    cols = list(COORDS) + list(EXTRA)
    trials = [os.path.basename(f).replace("_ik.mot", "") for f in sorted(glob.glob(os.path.join(ROOT, "data/subject3/walking*_ik.mot")))] + ["squats1_segment_0"]
    for tr in trials:
        m = io.read_mot(os.path.join(ROOT, "data/subject3", f"{tr}_ik.mot")); idm = io.read_mot(os.path.join(ROOT, "out/real/subject3", tr, f"{tr}_id.sto"))
        ok = (m.time >= idm.time[0]) & (m.time <= idm.time[-1]); tsel = m.time[ok]
        t = rng.choice(tsel, n_walk if tr.startswith("walking") else n_sq, replace=False)
        if not tr.startswith("walking"): t = np.append(t, tsel[np.argmax(m.column("knee_angle_r")[ok])])
        Q.append(np.column_stack([np.interp(t, m.time, m.column(c)) for c in cols]))
        T.append(np.column_stack([np.interp(t, idm.time, idm.column(c + "_moment")) for c in COORDS]))
    if ref: return np.vstack(Q), np.vstack(T) * (ref[0] * ref[1]) / (SUBJECT3_BODY[0] * SUBJECT3_BODY[1])
    return np.vstack(Q), np.vstack(T)


def poses_from_tasks(path, n_per=None, seed=0, ref=None):
    """ref=(mass, height) 를 주면 토크를 그 몸 크기로 맞춘다 (design_directions 와 같은 질량×키 비율)."""
    """tasks.csv 의 trial 마다 n개 프레임 + 설계 좌표별 극값 프레임(최대·최소). n_per 기본: 걷기 5, 그 외 10."""
    rng = np.random.default_rng(seed); Q, T = [], []; cols = list(COORDS) + list(EXTRA)
    for row in csv.DictReader(open(path, encoding="utf-8")):
        if not row.get("id") or row["id"].startswith("#"): continue
        ik = row["ik"] if os.path.isabs(row["ik"]) else os.path.join(ROOT, row["ik"]); idp = row["id"] if os.path.isabs(row["id"]) else os.path.join(ROOT, row["id"])
        m = io.read_mot(ik); idm = io.read_mot(idp)
        if m.in_degrees: m.data = np.radians(m.data)          # CARRT IK(OpenSim IK 툴 출력)는 도 단위, 설계 좌표는 전부 회전
        if not all(c in m.names for c in cols) or not all(c + "_moment" in idm.names for c in COORDS): print(f"  [건너뜀] {row['trial']}: 열 없음"); continue
        ok = (m.time >= idm.time[0]) & (m.time <= idm.time[-1]); tsel = m.time[ok]
        n = n_per or (5 if row["task"].startswith("walk") else 10)
        t = rng.choice(tsel, min(n, len(tsel)), replace=False)
        for c in COORDS:                                       # 극값 자세는 반드시 넣는다 (깊은 스쿼트 교훈)
            q = m.column(c)[ok]; t = np.append(t, [tsel[np.argmax(q)], tsel[np.argmin(q)]])
        tscale = (ref[0] * ref[1]) / (float(row["mass_kg"]) * float(row["height_m"])) if ref else 1.0
        Q.append(np.column_stack([np.interp(t, m.time, m.column(c)) for c in cols]))
        T.append(np.column_stack([np.interp(t, idm.time, idm.column(c + "_moment")) for c in COORDS]) * tscale)
    Q, T = np.vstack(Q), np.vstack(T)
    for cname, (lo, hi) in SP.get("pose_limit", {}).items():        # 특이점 근처 자세 제외 (design_space 참고)
        j = cols.index(cname); keep = (Q[:, j] >= lo) & (Q[:, j] <= hi)
        print(f"  자세 제외 {int((~keep).sum())}개: {cname} 가 [{lo:.2f}, {hi:.2f}] rad 밖"); Q, T = Q[keep], T[keep]
    return Q, T


def poses_time_from_tasks(path, dt=0.02, ref=None):
    """시간 영역 (2026-11-19): tasks.csv 의 trial 마다 dt 간격 연속 프레임 + 설계 좌표별 관절각·토크 최대/최소 프레임.
    poses_from_tasks 는 trial 당 무작위 5~10 개라 토크 피크를 놓칠 수 있고 시간 순서가 없다. 반환 Q, T, meta(trial (N,), time (N,), names)."""
    Q, T, tri, tim, names = [], [], [], [], []; cols = list(COORDS) + list(EXTRA)
    for row in csv.DictReader(open(path, encoding="utf-8")):
        if not row.get("id") or row["id"].startswith("#"): continue
        ik = row["ik"] if os.path.isabs(row["ik"]) else os.path.join(ROOT, row["ik"]); idp = row["id"] if os.path.isabs(row["id"]) else os.path.join(ROOT, row["id"])
        m = io.read_mot(ik); idm = io.read_mot(idp)
        if m.in_degrees: m.data = np.radians(m.data)
        if not all(c in m.names for c in cols) or not all(c + "_moment" in idm.names for c in COORDS): print(f"  [건너뜀] {row['trial']}: 열 없음"); continue
        ok = (m.time >= idm.time[0]) & (m.time <= idm.time[-1]); tsel = m.time[ok]
        t = list(np.arange(tsel[0], tsel[-1] + 1e-9, dt))
        for c in COORDS:
            q = m.column(c)[ok]; tau = np.interp(tsel, idm.time, idm.column(c + "_moment"))
            t += [tsel[np.argmax(q)], tsel[np.argmin(q)], tsel[np.argmax(tau)], tsel[np.argmin(tau)]]
        t = np.unique(np.round(t, 6))
        tscale = (ref[0] * ref[1]) / (float(row["mass_kg"]) * float(row["height_m"])) if ref else 1.0
        Q.append(np.column_stack([np.interp(t, m.time, m.column(c)) for c in cols]))
        T.append(np.column_stack([np.interp(t, idm.time, idm.column(c + "_moment")) for c in COORDS]) * tscale)
        tri.append(np.full(len(t), len(names))); tim.append(t); names.append(f"{row.get('dataset', '')}/{row.get('subject', '')}/{row['trial']}")
    Q, T, tri, tim = np.vstack(Q), np.vstack(T), np.concatenate(tri), np.concatenate(tim)
    for cname, (lo, hi) in SP.get("pose_limit", {}).items():
        j = cols.index(cname); keep = (Q[:, j] >= lo) & (Q[:, j] <= hi)
        print(f"  자세 제외 {int((~keep).sum())}개: {cname} 가 [{lo:.2f}, {hi:.2f}] rad 밖"); Q, T, tri, tim = Q[keep], T[keep], tri[keep], tim[keep]
    return Q, T, dict(trial=tri, time=tim, names=names)


def poses_time(dt=0.02, ref=None):
    """하체 기본 자세(poses() 와 같은 LabValidation 피험자3 걷기 6 + 스쿼트)를 시간 영역으로: dt 간격 + 좌표별 관절각·토크 극값."""
    Q, T, tri, tim, names = [], [], [], [], []; cols = list(COORDS) + list(EXTRA)
    trials = [os.path.basename(f).replace("_ik.mot", "") for f in sorted(glob.glob(os.path.join(ROOT, "data/subject3/walking*_ik.mot")))] + ["squats1_segment_0"]
    for tr in trials:
        m = io.read_mot(os.path.join(ROOT, "data/subject3", f"{tr}_ik.mot")); idm = io.read_mot(os.path.join(ROOT, "out/real/subject3", tr, f"{tr}_id.sto"))
        ok = (m.time >= idm.time[0]) & (m.time <= idm.time[-1]); tsel = m.time[ok]
        t = list(np.arange(tsel[0], tsel[-1] + 1e-9, dt))
        for c in COORDS:
            q = m.column(c)[ok]; tau = np.interp(tsel, idm.time, idm.column(c + "_moment"))
            t += [tsel[np.argmax(q)], tsel[np.argmin(q)], tsel[np.argmax(tau)], tsel[np.argmin(tau)]]
        t = np.unique(np.round(t, 6))
        Q.append(np.column_stack([np.interp(t, m.time, m.column(c)) for c in cols]))
        T.append(np.column_stack([np.interp(t, idm.time, idm.column(c + "_moment")) for c in COORDS]))
        tri.append(np.full(len(t), len(names))); tim.append(t); names.append(f"LabValidation/subject3/{tr}")
    Q, T = np.vstack(Q), np.vstack(T)
    if ref: T = T * (ref[0] * ref[1]) / (SUBJECT3_BODY[0] * SUBJECT3_BODY[1])
    return Q, T, dict(trial=np.concatenate(tri), time=np.concatenate(tim), names=names)


class Frames:
    """자세 × (기준, 좌표 j ±δ) × 바디 의 지면 변환을 미리 계산해 둔다. perturb=False 면 기준 자세만 (길이만 잴 때)."""
    def __init__(self, model_path, Q, perturb=True):
        import opensim as osim
        osim.Logger.setLevelString("Error")
        proc = osim.ModelProcessor(model_path); proc.append(osim.ModOpRemoveMuscles()); model = proc.process(); state = model.initSystem()
        cs = model.getCoordinateSet(); names = list(COORDS) + list(EXTRA)
        nP = 1 + 2 * len(COORDS) if perturb else 1
        self.R = np.zeros((len(Q), nP, len(BODIES), 3, 3)); self.p = np.zeros((len(Q), nP, len(BODIES), 3))
        e = [osim.Vec3(1, 0, 0), osim.Vec3(0, 1, 0), osim.Vec3(0, 0, 1)]; o = osim.Vec3(0, 0, 0)
        for t, q in enumerate(Q):
            for k in range(nP):
                qq = q.copy()
                if k > 0: j = (k - 1) // 2; qq[j] += DELTA if k % 2 == 1 else -DELTA
                for name, v in zip(names, qq): cs.get(name).setValue(state, float(v), False)
                model.realizePosition(state)
                for b, bn in enumerate(BODIES):
                    body = model.getBodySet().get(bn)
                    p0 = np.array([body.findStationLocationInGround(state, o).get(i) for i in range(3)])
                    self.p[t, k, b] = p0
                    for a in range(3):
                        self.R[t, k, b, :, a] = np.array([body.findStationLocationInGround(state, e[a]).get(i) for i in range(3)]) - p0
        self.model = model


def wrap_segment(chain):
    """감는 구가 있고 chain 에 그 뼈가 둘째 이후로 있으면, 구를 감는 구간 번호(점 n-1 → n). 없으면 None."""
    if not WRAP: return None
    wb = BODIES.index(WRAP["body"])
    return chain.index(wb) - 1 if wb in chain and chain.index(wb) >= 1 else None


def sphere_path(A, B, C, R):
    """A→B 직선이 중심 C 반지름 R 인 구를 뚫으면 구 표면을 따라 도는 최단 경로 길이(A, B, C 가 만드는 평면 안의 접선+호).
    직선이 구를 비켜 가면 |AB|. 점이 구 안에 있으면 R 로 당겨 계산하고 inside 로 알린다. 모양은 (..., 3)."""
    CA, CB = A - C, B - C; a = np.linalg.norm(CA, axis=-1); b = np.linalg.norm(CB, axis=-1)
    inside = (a < R) | (b < R); a_ = np.maximum(a, R * 1.0001); b_ = np.maximum(b, R * 1.0001)
    cos_t = np.clip(np.einsum("...i,...i->...", CA, CB) / (a * b + 1e-12), -1.0, 1.0); theta = np.arccos(cos_t)
    alpha = np.arccos(R / a_); beta = np.arccos(R / b_); arc = theta - alpha - beta
    L_wrap = np.sqrt(a_ ** 2 - R ** 2) + np.sqrt(b_ ** 2 - R ** 2) + R * np.maximum(arc, 0.0)
    L_line = np.linalg.norm(A - B, axis=-1)
    return np.where(arc > 0, L_wrap, L_line), inside


def cylinder_path(A, B, C, a, R, s, with_theta=False):
    """A→B 가 축 a (단위), 축 위 한 점 C, 반지름 R 인 원통을 정해진 방향 s 로 감는 최단 경로 길이 (2026-11-19, 노뼈 회내외 드럼).
    s = +1 이면 A 에서 B 로 갈 때 축 a 둘레 반시계(오른손), −1 이면 시계. 축에 수직인 평면에 투영하면 접선 + 호 + 접선 (L2D),
    원통을 펼치면 경로가 직선이라 L = sqrt(L2D² + Δz²) (Δz = 축 방향 거리). 감김 각 θ = (A→B 의 s 방향 각) − 두 접선 각.
    θ ≤ 0 이면 원통에 안 닿아 직선 |AB|. 점이 원통 안이면 R 로 당겨 계산하고 inside 로 알린다. 모양은 (..., 3), a 는 (..., 3) 또는 (3,)."""
    dA, dB = A - C, B - C; zA = np.einsum("...i,...i->...", dA, a); zB = np.einsum("...i,...i->...", dB, a)
    pA = dA - zA[..., None] * a; pB = dB - zB[..., None] * a; rA = np.linalg.norm(pA, axis=-1); rB = np.linalg.norm(pB, axis=-1)
    inside = (rA < R) | (rB < R); rA_ = np.maximum(rA, R * 1.0001); rB_ = np.maximum(rB, R * 1.0001)
    ang = np.arctan2(np.einsum("...i,...i->...", np.cross(pA, pB), a), np.einsum("...i,...i->...", pA, pB))     # A→B 반시계 각 (−π, π]
    theta = np.mod(s * ang, 2 * np.pi) - np.arccos(R / rA_) - np.arccos(R / rB_)
    L2D = np.sqrt(rA_ ** 2 - R ** 2) + np.sqrt(rB_ ** 2 - R ** 2) + R * np.maximum(theta, 0.0)
    L = np.where(theta > 0, np.sqrt(L2D ** 2 + (zB - zA) ** 2), np.linalg.norm(A - B, axis=-1))
    return (L, inside, theta) if with_theta else (L, inside)


def cyl_segment(chain, wc):
    """케이블별 원통 감김 wc (json 의 wrap_cyl) 이 걸리는 구간 번호: chain 에서 원통이 붙은 뼈로 들어가는 구간 (점 n-1 → n)."""
    wb = BODIES.index(wc["body"])
    if wb not in chain or chain.index(wb) < 1: raise SystemExit(f"wrap_cyl 뼈 {wc['body']} 가 경로 {[BODIES[c] for c in chain]} 의 둘째 이후에 없다")
    return chain.index(wb) - 1


def path_lengths(fr, chain, P, with_inside=False, wrap_cyl=None, with_theta=False):
    """chain: 바디 인덱스 목록, P: (len(chain), 3) 바디 좌표. 반환 L (poses, nP). 감는 구가 걸린 구간은 구 표면 경로로 잰다.
    wrap_cyl (케이블별, json 의 wrap_cyl dict: body, center, axis, radius, sign) 을 주면 그 구간은 원통 감김으로 잰다 (2026-11-19). 없으면 지금과 같다.
    with_theta 면 원통 감김 각 (poses, nP) 도 돌려준다 (원통 없으면 None)."""
    pts = np.stack([np.einsum("tkij,j->tki", fr.R[:, :, b], P[n]) + fr.p[:, :, b] for n, b in enumerate(chain)], axis=2)  # (T, nP, n, 3)
    seg = np.linalg.norm(np.diff(pts, axis=2), axis=3)                                                             # (T, nP, n-1)
    inside = np.zeros(seg.shape[:2], bool); ws = wrap_segment(chain); theta = None
    if ws is not None:
        wb = BODIES.index(WRAP["body"]); C = np.einsum("tkij,j->tki", fr.R[:, :, wb], np.asarray(WRAP["center"], float)) + fr.p[:, :, wb]
        seg[:, :, ws], inside = sphere_path(pts[:, :, ws], pts[:, :, ws + 1], C, float(WRAP["radius"]))
    if wrap_cyl:
        wc = wrap_cyl; cs = cyl_segment(chain, wc); wb = BODIES.index(wc["body"])
        C = np.einsum("tkij,j->tki", fr.R[:, :, wb], np.asarray(wc["center"], float)) + fr.p[:, :, wb]
        ax = np.asarray(wc["axis"], float); ax = np.einsum("tkij,j->tki", fr.R[:, :, wb], ax / np.linalg.norm(ax))
        seg[:, :, cs], ins_c, theta = cylinder_path(pts[:, :, cs], pts[:, :, cs + 1], C, ax, float(wc["radius"]), float(wc["sign"]), with_theta=True)
        inside = inside | ins_c
    L = seg.sum(axis=2)
    out = (L, inside) if with_inside else (L,)
    if with_theta: out = out + (theta,)
    return out if len(out) > 1 else out[0]


def moment_arms(fr, chain, P, wrap_cyl=None):
    L = path_lengths(fr, chain, P, wrap_cyl=wrap_cyl)
    r = np.stack([-(L[:, 1 + 2 * j] - L[:, 2 + 2 * j]) / (2 * DELTA) for j in range(len(COORDS))], axis=1)  # (T, 4)
    return r, L[:, 0]


def chain_for(v, thr=0.12):
    return ds.chain_for(SP, v, thr)


def cap_penalty(i, force, L, chain, P):
    """모듈 먼저: 케이블 i 에 배정한 등급의 힘·길이 변화 상한 벌점과 길이 변화 ΔL (m). ΔL 은 설계 자세 + stroke_check 자세(fr2)에서 잰다."""
    if MODCAPS.get("fr2") is not None:
        L2 = path_lengths(MODCAPS["fr2"], chain, P)[:, 0]; dL = float(max(L.max(), L2.max()) - min(L.min(), L2.min()))
    else: dL = float(L.max() - L.min())
    return W_CAP_F * max(0.0, force / MODCAPS["fcap"][i] - 1) ** 2 + W_CAP_S * max(0.0, dL / MODCAPS["smax"][i] - 1) ** 2, dL


def cost(x, fr, chain, v, cap, detail=False, i=None):
    P = x.reshape(len(chain), 3); r, L = moment_arms(fr, chain, P)
    if WRAP is not None and wrap_segment(chain) is not None:                       # 부착점이 감는 구 안에 들어가면 안 된다
        _, inside = path_lengths(fr, chain, P, with_inside=True)
        if inside.any(): return (1e3 if MODCAPS is None else DEGEN_MF) + float(inside.mean()) if not detail else None
    nr = np.linalg.norm(r, axis=1); ok = nr > 1e-4
    if ok.sum() < len(r) * 0.5: return (1e3 if MODCAPS is None else DEGEN_MF) if not detail else None
    e2 = np.linalg.norm(r[ok] / nr[ok, None] - v, axis=1) ** 2
    dir_err = float(0.5 * e2.mean() + 0.5 * e2.max())          # 평균 + 최악 자세 (깊은 스쿼트에서 방향이 틀어지는 것을 막는다, 9/15)
    arm_eff = (r @ v)                                          # 설계 방향으로의 유효 모멘트암 (m)
    arm_min = float(arm_eff.min()); arm_mean = float(arm_eff.mean())
    stroke = float((L.max() - L.min()) / L.mean()); sma = float((L.max() - L.min()) / L.min())
    force = MARGIN * cap / max(arm_min, 0.005)                 # 설계 용량을 가장 불리한 자세에서 내는 데 필요한 힘 (모멘트암 하한 5 mm — 0 에 가까우면 힘 벌점이 폭발해 DE 가 퇴화 해로 도망간다, 2026-11-03)
    pen = 20 * max(0.0, stroke - STROKE) ** 2 + 5 * max(0.0, force / F_MODULE - 1) ** 2 + 1e4 * max(0.0, 0.015 - arm_min) ** 2 + (50 * max(0.0, sma - SMA_EPS) ** 2 if SMA_EPS else 0.0) + (50 * max(0.0, (L.max() - L.min()) / STROKE_MAX - 1) ** 2 if STROKE_MAX else 0.0)
    if MODCAPS is not None and i is not None:                  # 모듈 먼저 (--modules): 부위 전체 상한 대신 케이블별 상한
        pc, dLc = cap_penalty(i, force, L, chain, P); pen += pc
        if detail: return dict(dir_err=dir_err, arm_min=arm_min, arm_mean=arm_mean, stroke=stroke, sma=sma, force=force, L_mean=float(L.mean()), pen=pen, dL_all=dLc)
    if detail: return dict(dir_err=dir_err, arm_min=arm_min, arm_mean=arm_mean, stroke=stroke, sma=sma, force=force, L_mean=float(L.mean()), pen=pen)
    return dir_err + pen


# ─────────────────────────────── 브래킷 제약 (--bracket) ───────────────────────────────
def setup_bracket(H, box_scale=1.0, body_scale=None):
    """design_space.BONES 를 이 공간의 뼈 인덱스로. 축은 box_scale (그리고 마디별 추가 배율) 로 같이 줄인다 (관 반지름은 하드웨어 값이라 그대로)."""
    bones, clear = {}, {}; body_scale = body_scale or {}
    for b, spec in ds.BONES.items():
        if b not in BODIES: continue
        j = BODIES.index(b); f = np.array(body_scale.get(b, (1.0, 1.0, 1.0))) * box_scale
        if "p0" in spec: bones[j] = (np.array(spec["p0"]) * f, np.array(spec["p1"]) * f, float(spec["radius"]))
        if "clear" in spec: clear[j] = float(spec["clear"])
    return dict(H=float(H), bones=bones, clear=clear)


def setup_bracket_cad(path):
    """--bracket-cad (2026-11-19, 연구실 팔 시제품): tools/cad_arm_model.py 의 맞춤 json + tools/cad_sdf.py 의 부품 거리장.
    CAD 부품이 있는 뼈는 관 대신 실제 부품 표면에서의 높이 h 를 쓰고, 와이어 구간이 부품 속을 지나는지도 본다 (표본 프레임)."""
    sys.path.insert(0, os.path.join(ROOT, "tools")); import cad_sdf as cs
    fj = json.load(open(path if os.path.isabs(path) else os.path.join(ROOT, path), encoding="utf-8"))
    names = {"humerus": "humerus_r", "ulna": "ulna_r", "radius": "radius_r"}; bodies = {}
    for bone, parts in fj["parts"].items():
        bn = names[bone]
        if bn not in BODIES: continue
        bodies[BODIES.index(bn)] = dict(A=np.array(fj["A"][bn]), t=np.array(fj["t"][bn]), grids=[cs.load(p) for p in parts])
    return dict(bodies=bodies, E_C=np.array(fj["E_C"]), sdf=cs.sdf, json=os.path.relpath(path if os.path.isabs(path) else os.path.join(ROOT, path), ROOT), frames={})


def cad_height(b, p):
    """뼈 b 좌표 점 p (m, (..., 3)) 의 CAD 부품 표면 높이 (m, 음수 = 부품 속). 그 뼈의 부품들 중 가장 가까운 것."""
    cb = BRACKET["cad"]["bodies"][b]; x = BRACKET["cad"]["E_C"] + np.einsum("ij,...j->...i", np.linalg.inv(cb["A"]), np.asarray(p) - cb["t"])
    return np.min([BRACKET["cad"]["sdf"](g, x) for g in cb["grids"]], axis=0) / 1000.0


CAD_SEG_FRAMES = 96; CAD_SEG_SAMPLES = 16            # 와이어 관통 검사: 관절각 분위 프레임 수, 구간당 표본 점 수
W_CAD_SEG = 1e7                                      # 관통 깊이(m)² 합 벌점 — W_BRACKET_WALL(1e4) 이면 2.7 mm 관통이 0.07 로 일 용량 항(≈0.3)보다 작아 그냥 뚫었다 (2026-11-19)


def cad_seg_frames(fr):
    """와이어 관통 검사용 프레임: CAD 뼈들 사이 상대 회전각의 분위 (끝 자세 포함). fr 마다 한 번 계산해 둔다."""
    key = (id(fr), fr.R.shape[0]); cache = BRACKET["cad"]["frames"]
    if key in cache: return cache[key]
    bs = sorted(BRACKET["cad"]["bodies"]); idx = set()
    for a_, b_ in zip(bs[:-1], bs[1:]):
        rel = np.einsum("tji,tjk->tik", fr.R[:, 0, a_], fr.R[:, 0, b_]); ang = np.arccos(np.clip((np.trace(rel, axis1=1, axis2=2) - 1) / 2, -1, 1))
        o = np.argsort(ang); idx |= set(o[np.linspace(0, len(o) - 1, CAD_SEG_FRAMES // max(len(bs) - 1, 1)).astype(int)].tolist())
    cache[key] = np.array(sorted(idx)); return cache[key]


def cad_seg_depth(fr, chain, P, frames=None):
    """케이블 구간의 표본 점이 CAD 부품 속으로 들어간 깊이 (m, ≥ 0) 배열. frames 없으면 cad_seg_frames."""
    f = cad_seg_frames(fr) if frames is None else frames
    pts = np.stack([np.einsum("tij,j->ti", fr.R[f, 0, b], P[n]) + fr.p[f, 0, b] for n, b in enumerate(chain)], axis=1)   # (F, n, 3)
    w = np.linspace(0.0, 1.0, CAD_SEG_SAMPLES + 2)[1:-1]
    S_ = (pts[:, :-1, None, :] * (1 - w)[None, None, :, None] + pts[:, 1:, None, :] * w[None, None, :, None]).reshape(len(f), -1, 3)   # (F, m, 3) 지면
    dep = np.zeros(S_.shape[:2])
    for b in BRACKET["cad"]["bodies"]:
        loc = np.einsum("tji,tmj->tmi", fr.R[f, 0, b], S_ - fr.p[f, 0, b][:, None, :])      # 지면 → 뼈 b 좌표
        dep = np.maximum(dep, -np.minimum(cad_height(b, loc), 0.0))
    return dep


def bracket_terms(fr, chain, P):
    """부착점마다 뼈 관 표면에서의 높이 h (축 거리 − 반지름)와 벌점: h < 0 (뼈 속) 또는 h > H (브래킷보다 높음),
    관절을 건너는 구간이 그 관절 중심에 clear 보다 가까이 지나감 (설계 자세 기준 상태 전부). 반환 (벌점, [h …] 축 있는 뼈의 점만).
    --bracket-cad 면 CAD 부품이 있는 뼈는 실제 부품 표면 높이를 쓰고, 와이어가 부품 속을 지나는 깊이²도 벌점."""
    pen, hs = 0.0, []
    cad = BRACKET.get("cad")
    for n, b in enumerate(chain):
        if cad and b in cad["bodies"]:
            h = float(cad_height(b, P[n])); hs.append(h)
            pen += W_BRACKET_WALL * (max(0.0, -h) ** 2 + max(0.0, h - BRACKET["H"]) ** 2)
            if BRACKET.get("xband"):                                    # --cad-xband: 점의 CAD x (굽힘축 방향, mm) 를 띠 안에 — 와이어를 팔의 앞뒤 평면에 둔다
                cb = cad["bodies"][b]; xc = float((cad["E_C"] + np.linalg.inv(cb["A"]) @ (np.asarray(P[n]) - cb["t"]))[0]); lo, hi = BRACKET["xband"]
                pen += W_BRACKET_WALL * (max(0.0, lo - xc, xc - hi) / 1000.0) ** 2
            if n == 0 and BRACKET.get("Hmin"): pen += W_BRACKET_WALL * max(0.0, BRACKET["Hmin"] - h) ** 2      # --bracket-min: 가이드(첫 점)는 부품 표면에서 Hmin 이상
            continue
        if b not in BRACKET["bones"]: continue
        p0, p1, R = BRACKET["bones"][b]; ax = p1 - p0; tpar = np.clip((P[n] - p0) @ ax / (ax @ ax), 0.0, 1.0)
        h = float(np.linalg.norm(P[n] - (p0 + tpar * ax)) - R); hs.append(h)
        pen += W_BRACKET_WALL * (max(0.0, -h) ** 2 + max(0.0, h - BRACKET["H"]) ** 2)
    for n in range(len(chain) - 1):
        bc = chain[n + 1]
        if chain[n] == bc or bc not in BRACKET["clear"]: continue
        A = np.einsum("tij,j->ti", fr.R[:, 0, chain[n]], P[n]) + fr.p[:, 0, chain[n]]; B = np.einsum("tij,j->ti", fr.R[:, 0, bc], P[n + 1]) + fr.p[:, 0, bc]
        C = fr.p[:, 0, bc]; AB = B - A; u = np.clip(np.einsum("ti,ti->t", C - A, AB) / (np.einsum("ti,ti->t", AB, AB) + 1e-12), 0.0, 1.0)
        dmin = float(np.linalg.norm(A + u[:, None] * AB - C, axis=1).min())
        pen += W_BRACKET_WALL * max(0.0, BRACKET["clear"][bc] - dmin) ** 2
    if cad: pen += W_CAD_SEG * float((cad_seg_depth(fr, chain, P) ** 2).sum())
    return pen, hs


# ─────────────────────────────── 시간 영역 1단계 비용 (--force-mode time) ───────────────────────────────
W_ARM_MIN = 20.0


def arm_min_pen(arm_act):
    """--arm-min (2026-11-20, 연구실 시험용 설계): 케이블이 당기는 프레임에서 유효 모멘트암이 arm_min 보다 작으면 벌점 W × 평균 (1 − arm/arm_min)².
    시간 영역 비용은 모듈 용량만 작게 하므로 와이어를 축에 바짝 붙인다 — 연구실이 요구한 "뼈 축에서 40 mm" 처럼 와이어를 멀리 두려면 이 하한을 준다. 기본 None = 지금과 같다."""
    am = TDS.get("arm_min") if TDS else None
    if not am or not len(arm_act): return 0.0
    return W_ARM_MIN * float((np.maximum(0.0, 1.0 - arm_act / am) ** 2).mean())


def len_min_pen(L):
    """--len-min: 가이드–부착점 와이어 길이가 어느 자세에서든 len_min 보다 짧으면 벌점 (두 점이 겹치는 답을 막는다; 모멘트암 하한만 주면 폄 와이어가 L≈0 으로 간다)."""
    lm = TDS.get("len_min") if TDS else None
    pen = W_ARM_MIN * float(max(0.0, 1.0 - L.min() / lm) ** 2) if lm else 0.0
    lx = TDS.get("len_max") if TDS else None                        # --len-max (2026-11-20): 와이어가 너무 길면 (가이드가 어깨 쪽으로 달아나 모듈 놓을 자리가 없어짐) 벌점
    return pen + (W_ARM_MIN * float(max(0.0, L.max() / lx - 1.0) ** 2) if lx else 0.0)


def cost_time(x, fr, chain, v, i, detail=False):
    """케이블 i 하나: 방향오차 + w_work × 모듈 일 용량 / 기준 + 스트로크 비 벌점 + 브래킷. 장력 T(t) = 토크 몫 u_i(t) ÷ 그 순간 유효 모멘트암."""
    P = x.reshape(len(chain), 3); r, L = moment_arms(fr, chain, P)
    if WRAP is not None and wrap_segment(chain) is not None:
        _, inside = path_lengths(fr, chain, P, with_inside=True)
        if inside.any(): return DEGEN_MF + float(inside.mean()) if not detail else None
    nr = np.linalg.norm(r, axis=1); ok = nr > 1e-4
    if ok.sum() < len(r) * 0.5: return DEGEN_MF if not detail else None
    U = TDS["U"][:, i]; act = ok & (U > 0.01 * max(U.max(), 1e-9))           # 방향은 케이블이 당기는 프레임에서만 본다 (쉬는 순간 모멘트암이 0 을 지나 부호가 바뀌는 것은 문제 아님)
    if act.sum() < 3: act = ok
    e2 = np.linalg.norm(r[act] / nr[act, None] - v, axis=1) ** 2; dir_err = float(0.5 * e2.mean() + 0.5 * e2.max())
    arm_eff = r @ v
    Tt = np.where(U > 1e-9, U / np.maximum(arm_eff, 0.002), 0.0); s = L.max() - L
    z = td.size_cable(Tt, s, TDS["curves"], TDS["grades"], TDS["nb"])
    stroke = float((L.max() - L.min()) / L.mean()); sma = float((L.max() - L.min()) / L.min())
    # 길이 비 벌점(ΔL/L ≤ STROKE)은 시간 영역에서 끈다: 모듈은 뼈를 따라 놓고 선(힘줄)으로 잇기 때문에 모듈 길이가 케이블 직선 길이에 묶이지 않는다.
    # 필요한 모듈 능동 길이(수축 최대 ÷ 0.4)는 결과에 보고만 한다 (--sma-eps 를 주면 그 제약은 건다)
    pen = (50 * max(0.0, sma - SMA_EPS) ** 2 if SMA_EPS else 0.0)
    bp, hs = bracket_terms(fr, chain, P) if BRACKET else (0.0, [])
    pen += arm_min_pen(arm_eff[act]) + len_min_pen(L)
    work = TDS["w_work"] * z["obj"] / TDS["Wref"][i]; total = dir_err + work + pen + bp + TDS["w_h"] * sum(max(0.0, h) for h in hs)
    if detail: return dict(dir_err=dir_err, arm_min=float(arm_eff.min()), arm_mean=float(arm_eff.mean()), stroke=stroke, sma=sma, force=z["T_peak"], L_mean=float(L.mean()),
                           pen=pen + bp, size=z, heights=hs, work_term=work, total=total)
    return total


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--space", default="lower"); ap.add_argument("--tasks", default=None, help="자세 표본을 뽑을 trial 목록 CSV (없으면 하체 기본)")
    ap.add_argument("--design", default=None, help="방향 설계 npz (기본: report/newmus_V_design_k8.npz 또는 _<space>_k<k>)")
    ap.add_argument("--k", type=int, default=None); ap.add_argument("--model", default=None, help="뼈대 모델 .osim (기본 data/subject3)")
    ap.add_argument("--out", default=None, help="결과 json (기본 report/newmus_paths_k8.json 또는 _<space>_k<k>)")
    ap.add_argument("--refine-from-json", action="store_true", help="1단계 결과를 다시 쓰고 2단계만")
    ap.add_argument("--box-scale", type=float, default=1.0, help="부착 상자·감는 구를 이 배율로 (로봇 키 / 피험자3 키, 예: 1.35/1.62)")
    ap.add_argument("--ref", type=float, nargs=2, metavar=("MASS", "HEIGHT"), default=None, help="자세 표본의 토크를 이 몸 크기로 맞춘다 (예: 40 1.35)")
    ap.add_argument("--out-tag", default="", help="출력 파일 이름 뒤에 붙일 말 (예: _40kg)")
    ap.add_argument("--wrap-radius", type=float, default=None, help="감는 구 반지름 m (기본 = 공간 값 × box-scale). 로봇 어깨 외피 크기")
    ap.add_argument("--f-module", type=float, default=None, help="케이블 힘 상한 N (기본 = 공간의 f_module). 낮추면 모멘트암 큰 경로를 찾는다")
    ap.add_argument("--maxfev", type=int, default=20000)
    ap.add_argument("--n-per", type=int, default=None, help="trial 당 자세 수 (기본: 걷기 5, 그 외 10). 들기 자세가 모자라면 늘린다 (2026-11-09)")
    ap.add_argument("--stroke-max", type=float, default=None, help="모듈 스트로크 상한 m: 케이블 길이 변화 (L_max−L_min) ≤ 이 값 (계획서 모듈표 M1 0.12, M2·M6 0.08, M3 0.06)")
    ap.add_argument("--sma-eps", type=float, default=None, help="SMA 묶음 제약: 케이블 (L_max−L_min)/L_min 상한 (park2019 수축률 0.5 ÷ 여유 1.2 = 0.417). 기본 없음")
    ap.add_argument("--modules", default=None, help="모듈 먼저 (2026-11-05): report/module_first.py 가 만든 케이블별 등급 json. 주면 --f-module/--stroke-max 대신 "
                    "케이블마다 그 등급의 힘(묶음 수 × 힘 상한 ÷ 힘 배율)·스트로크 상한을 1 %% 안쪽에 세게 건다. ΔL 은 stroke_check 자세까지 합쳐 잰다")
    ap.add_argument("--seed", type=int, default=0, help="1단계 차분진화 시드 묶음 (시드 = 근육 번호 + 10×시도 + 1000×이 값). 0 이면 기존과 같다 — 반복 실행으로 결과 흔들림을 잴 때 (2026-11-05)")
    ap.add_argument("--refine-trust", type=float, default=None, help="2단계 Powell 을 현재 점 ± 이 값(부착점 m, 힘 배율은 ×10) 안에서 돌리고 옮겨 가며 반복. "
                    "기본(없음)은 상자 전체 범위 — 1단계 해가 상한 경계에 붙어 있으면 상자 전체 선탐색이 시작점보다 나쁜 점을 돌려줘 2단계가 아무것도 못 한다 (모듈 먼저, 2026-11-05)")
    ap.add_argument("--force-mode", choices=("amin", "time"), default="amin", help="amin(기본, 지금까지): 힘 = 1.1 × 용량 ÷ 최소 모멘트암. "
                    "time (2026-11-19): 프레임마다 장력 T(t)·수축 s(t) 로 모듈 일 용량 최소화 (report/timedomain.py)")
    ap.add_argument("--dt", type=float, default=0.02, help="time: 프레임 간격 s (+ 좌표별 관절각·토크 극값)")
    ap.add_argument("--grades", default="M1,M2,M3,M4,M5,M6", help="time: 허용 모듈 등급 (예: 팔꿈치 M3,M4,M5)")
    ap.add_argument("--curves", default=None, help="time: 힘–수축 곡선 json (기본 report/module_curves.json, 임시)")
    ap.add_argument("--beta", type=float, default=None, help="time: 곡선의 수축 감소 β 를 모든 등급에 (기본 json 값)")
    ap.add_argument("--max-bundles", type=int, default=None, help=f"time: 케이블당 묶음 상한 (기본 {td.MAX_BUNDLES})")
    ap.add_argument("--w-work", type=float, default=0.1, help="time: 일 용량 항 가중치 (일 용량 ÷ 1.1·용량·관절 각도 범위)")
    ap.add_argument("--bracket-cad", default=None, help="--bracket 와 함께: CAD 맞춤 json (tools/cad_arm_model.py). CAD 부품이 있는 뼈는 관 대신 실제 부품 표면 높이 + 와이어 관통 벌점 (tools/cad_sdf.py 거리장, 연구실 팔 시제품, 2026-11-19). 없으면 지금과 같다")
    ap.add_argument("--bracket", type=float, default=None, help="브래킷 높이 상한 H (m): 축이 있는 뼈(design_space.BONES)의 부착점을 뼈 관 표면 ~ +H 안에, "
                    "관절을 건너는 구간은 관절 중심에서 clear 이상. 없으면 지금처럼 상자만")
    ap.add_argument("--cad-xband", default=None, help="--bracket-cad 와 함께: 모든 점의 CAD x 좌표 범위 'lo,hi' (mm). 1자유도 팔 시험처럼 와이어를 앞뒤 평면에 둘 때")
    ap.add_argument("--len-max", type=float, default=None, help="가이드–부착점 와이어 길이 상한 (m). 가이드가 어깨 쪽으로 멀리 가서 모듈 놓을 자리가 없어지는 답을 막는다")
    ap.add_argument("--bracket-min", type=float, default=None, help="--bracket-cad 와 함께 (2026-11-20): 가이드(케이블 첫 점)의 부품 표면 높이 하한 (m). 모듈이 뼈대에서 떨어져 놓여야 할 때 (예: 봉 중심에서 40 mm = 반지름 10 mm 봉 표면에서 0.030)")
    ap.add_argument("--w-bracket", type=float, default=2.0, help="브래킷 높이 합(m) 가중치 (낮은 브래킷 선호)")
    ap.add_argument("--body-scale", action="append", default=[], metavar="BODY=X,Y,Z",
                    help="그 뼈의 상자·브래킷 축만 --box-scale 위에 더 곱한다 (rescale_body --segment-scale 와 같은 값, 예: humerus_r=1,1.0974,1). 2026-11-19")
    ap.add_argument("--alloc-iters", type=int, default=300, help="time 2단계: 프레임별 장력 분배 반복 수")
    ap.add_argument("--stage1-stride", type=int, default=1, help="time 1단계(케이블별 DE)에서 프레임을 이 간격으로 솎는다 (극값·몫 상위 프레임은 항상). 2단계는 전부")
    ap.add_argument("--w-dir-time", type=float, default=0.0, help="time 2단계 방향오차 가중치 (기본 0: 실제 모멘트암으로 장력을 나누므로 배분비 V 에 맞출 필요가 없다. "
                    "3 으로 두었더니 팔꿈치 폄 케이블이 몇 각도에서 굽힘을 돕는 것까지 벌점 2×3 을 받아 2단계가 1단계보다 나빠졌다, 2026-11-19)")
    ap.add_argument("--refine-from", default=None, help="time: 이 경로 json 을 읽어 2단계만 다시 (결과는 --out-tag 이름으로 새로 저장)")
    ap.add_argument("--robust-frac", type=float, default=0.0, help="time: 모든 자세에서 좌표마다 양방향으로 그 좌표 최대 토크 × 이 비율을 낼 수 있게 가상 요구 프레임을 더한다 "
                    "(0 = 끔). 2026-11-19: 0 으로 설계한 팔꿈치가 Moco 걷기에서 미달 — 거의 편 자세에서 폄 케이블이 굽힘 쪽으로 넘어가 있었다")
    ap.add_argument("--robust-every", type=int, default=5, help="time: 가상 요구 프레임을 몇 프레임마다 둘지")
    ap.add_argument("--arm-min", type=float, default=None, help="time: 케이블이 당기는 프레임의 유효 모멘트암 하한 m (벌점). 연구실 시험용 '와이어를 축에서 멀리' 설계 (2026-11-20). 기본 없음")
    ap.add_argument("--len-min", type=float, default=None, help="time: 가이드–부착점 와이어 길이 하한 m (벌점; 두 점이 겹치는 답 방지). 기본 없음")
    a = ap.parse_args(); configure(a.space)
    if a.modules and a.force_mode == "time": raise SystemExit("--modules 와 --force-mode time 은 같이 쓰지 않는다 (time 은 등급을 경로와 함께 정한다)")
    if a.modules and (a.f_module or a.stroke_max): print("  --modules 를 주면 --f-module/--stroke-max 는 쓰지 않는다 (케이블별 상한으로 바뀐다)"); a.f_module = a.stroke_max = None
    if a.f_module: global F_MODULE; F_MODULE = a.f_module
    if a.sma_eps: global SMA_EPS; SMA_EPS = a.sma_eps; print(f"SMA 제약 ΔL/L_min ≤ {SMA_EPS:.3f}")
    if a.stroke_max: global STROKE_MAX; STROKE_MAX = a.stroke_max; print(f"모듈 스트로크 제약 ΔL ≤ {STROKE_MAX*1000:.0f} mm")
    if a.box_scale != 1.0:                                      # 로봇 크기에 맞춰 상자·감는 구를 줄인다
        global BOX, WRAP; BOX = {b: tuple(tuple(v * a.box_scale for v in ax) for ax in box) for b, box in BOX.items()}
        if WRAP: WRAP = dict(WRAP, center=tuple(v * a.box_scale for v in WRAP["center"]), radius=WRAP["radius"] * a.box_scale)
    if a.wrap_radius and WRAP: WRAP = dict(WRAP, radius=a.wrap_radius); print(f"감는 구 반지름 {a.wrap_radius} m")
    body_scale = {s.split("=")[0]: tuple(float(v) for v in s.split("=")[1].split(",")) for s in a.body_scale}
    for b, f in body_scale.items():                             # 마디별 추가 배율 (CAD 길이 맞춤) — 그 뼈 상자만
        if b in BOX: BOX = dict(BOX, **{b: tuple(tuple(v * fi for v in ax) for ax, fi in zip(BOX[b], f))})
    if body_scale: print(f"마디별 추가 배율 (상자·브래킷 축): {body_scale}")
    legacy = (a.space == "lower" and not a.tasks and not a.out_tag)
    k = a.k or SP["k"]; tag = ("_k8" if legacy else f"_{a.space}_k{k}") + a.out_tag
    design = a.design or os.path.join(ROOT, f"report/newmus_V_design{tag}.npz")
    d = np.load(design); V, c = d["V"], d["c"]; k = len(V); tag = ("_k8" if legacy else f"_{a.space}_k{k}") + a.out_tag   # k 는 설계 파일이 정한다
    model_path = a.model or os.path.join(ROOT, "data/subject3/subject3_scaled.osim")
    out_json = a.out or os.path.join(ROOT, f"report/newmus_paths{tag}.json")
    tasks_path = (os.path.join(ROOT, a.tasks) if not os.path.isabs(a.tasks) else a.tasks) if a.tasks else None
    meta = None
    if a.force_mode == "time":
        Q, T, meta = poses_time_from_tasks(tasks_path, dt=a.dt, ref=a.ref) if a.tasks else poses_time(dt=a.dt, ref=a.ref)
        if a.robust_frac > 0:                                   # 강건성 요구: 이 자세들에서 좌표 j 로 ±ρ·최대|τ_j| 를 따로 낼 수 있어야 한다 (trial −1)
            pk0 = np.abs(T).max(axis=0); iv = np.arange(0, len(Q), a.robust_every); Qv, Tv = [], []
            for j in range(T.shape[1]):
                for sg in (1.0, -1.0):
                    e = np.zeros(T.shape[1]); e[j] = sg * a.robust_frac * pk0[j]; Qv.append(Q[iv]); Tv.append(np.tile(e, (len(iv), 1)))
            nv = sum(len(x) for x in Qv); Q = np.vstack([Q] + Qv); T = np.vstack([T] + Tv)
            meta = dict(trial=np.r_[meta["trial"], np.full(nv, -1)], time=np.r_[meta["time"], np.full(nv, np.nan)], names=meta["names"])
            print(f"강건성 요구: {len(iv)} 자세 × 좌표 {T.shape[1]} × 양방향, 크기 = 좌표별 최대 토크 × {a.robust_frac} = {np.round(a.robust_frac * pk0, 1)} N·m → 가상 프레임 {nv}개", flush=True)
    else:
        Q, T = poses_from_tasks(tasks_path, n_per=a.n_per, ref=a.ref) if a.tasks else poses(ref=a.ref)
    print(f"[{a.space}] 근육 {k}개, 자세 표본 {len(Q)}개, 모델 {model_path}  (2단계 토크 좌표별 최대 {np.round(np.abs(T).max(axis=0), 1)} N·m)", flush=True)
    tic = time.time(); fr = Frames(model_path, Q); print(f"뼈 변환 저장 {time.time()-tic:.0f}s", flush=True)
    out = {"space": a.space, "coords": COORDS, "bodies": BODIES, "model": os.path.relpath(model_path, ROOT), "muscles": []}
    global FORCE_MODE, TDS, BRACKET
    if a.bracket is not None:
        BRACKET = setup_bracket(a.bracket, a.box_scale, body_scale)
        if a.bracket_cad:
            BRACKET["cad"] = setup_bracket_cad(a.bracket_cad); BRACKET["Hmin"] = a.bracket_min; BRACKET["xband"] = [float(v) for v in a.cad_xband.split(",")] if a.cad_xband else None
            print(f"브래킷 CAD: {BRACKET['cad']['json']} — 부품 표면 높이 0~{a.bracket*1000:.0f} mm, 와이어 관통 검사 {[BODIES[b] for b in BRACKET['cad']['bodies']]}", flush=True)
        print(f"브래킷: 높이 ≤ {a.bracket*1000:.0f} mm, 관 반지름·축 {[(BODIES[b], round(v[2]*1000)) for b, v in BRACKET['bones'].items()]}, 관절 간격 {[(BODIES[b], round(v*1000)) for b, v in BRACKET['clear'].items()]} mm", flush=True)
    if a.force_mode == "time":
        FORCE_MODE = "time"; F_MODULE = np.inf; STROKE_MAX = None             # 힘·스트로크 상한은 등급 판정이 대신한다
        grades = [g.strip() for g in a.grades.split(",") if g.strip()]
        U = td.shares(V, c, MARGIN * T)                                       # 프레임별 토크 몫 (N·m), 1단계용
        dth = np.array([max(np.ptp(Q[:, :len(COORDS)] @ V[i]), 1e-3) for i in range(k)]); Wref = np.maximum(MARGIN * c * dth, 1e-3)
        TDS = dict(curves=td.load_curves(a.curves, beta=a.beta), grades=grades, nb=a.max_bundles, w_work=a.w_work, w_h=a.w_bracket, U=U, Wref=Wref, arm_min=a.arm_min, len_min=a.len_min, len_max=a.len_max)
        fr_s1 = fr
        if a.stage1_stride > 1:                                               # 1단계(케이블별 DE)만 솎은 프레임 + 극값 + 케이블별 몫 상위 30 프레임. 2단계는 전부
            ext = set()
            for j in range(len(COORDS)): ext |= {int(np.argmax(T[:, j])), int(np.argmin(T[:, j])), int(np.argmax(Q[:, j])), int(np.argmin(Q[:, j]))}
            for i in range(k): ext |= set(int(t) for t in np.argsort(U[:, i])[-30:])
            idx = np.unique(np.r_[np.arange(0, len(Q), a.stage1_stride), sorted(ext)]).astype(int)
            fr_s1 = type("FramesView", (), {})(); fr_s1.R = fr.R[idx]; fr_s1.p = fr.p[idx]; TDS["U"] = U[idx]
            print(f"1단계 프레임: {len(idx)}개 (간격 ×{a.stage1_stride} + 극값·몫 상위), 2단계는 {len(Q)}개 전부", flush=True)
        print(f"시간 영역: 프레임 {len(Q)}개 (dt {a.dt} s, 시행 {len(meta['names'])}개), 허용 등급 {grades}, 묶음 ≤ {a.max_bundles or td.MAX_BUNDLES}, "
              f"β {[TDS['curves'][g]['beta'] for g in grades][0]:g} (임시 곡선), 일 용량 기준 {np.round(Wref, 1)} J", flush=True)
    mods = None
    if a.modules:                                               # 모듈 먼저: 케이블별 상한 (부위 전체 상한은 끈다)
        global MODCAPS                                          # F_MODULE·STROKE_MAX 는 위에서 이미 global
        mods = json.load(open(a.modules if os.path.isabs(a.modules) else os.path.join(ROOT, a.modules), encoding="utf-8"))
        assert [m["name"] for m in mods["muscles"]] == [f"M{i}" for i in range(k)], "모듈 json 의 케이블 수·이름이 설계와 다르다"
        assert np.allclose([m["capacity_Nm"] for m in mods["muscles"]], c), "모듈 json 이 다른 설계(용량 c)에서 나왔다"
        s = float(mods.get("force_scale", 1.0)); F_MODULE = np.inf; STROKE_MAX = None
        fcap = np.array([m["force_cap_N"] / s * m.get("design_scale_force", 1 - CAP_MARGIN) for m in mods["muscles"]])
        smax = np.array([m["stroke_cap_m"] * m.get("design_scale_stroke", 1 - CAP_MARGIN) for m in mods["muscles"]])
        Q2 = poses_from_tasks(tasks_path, n_per=20)[0] if a.tasks else poses(n_walk=20, n_sq=60)[0]      # stroke_check 와 같은 검사 자세
        MODCAPS = dict(fcap=fcap, smax=smax, fr2=Frames(model_path, Q2, perturb=False),
                       fcap_true=np.array([m["force_cap_N"] / s for m in mods["muscles"]]), smax_true=np.array([m["stroke_cap_m"] for m in mods["muscles"]]))
        print(f"모듈 먼저: {os.path.basename(a.modules)} (힘 배율 {s:g}), 검사 자세 {len(Q2)}개 추가")
        for i, m in enumerate(mods["muscles"]): print(f"  M{i}: {m['module']} ×{m['bundles']}  설계 힘 상한 {fcap[i]:.0f} N  ΔL 상한 {smax[i]*1000:.1f} mm")
    if a.refine_from_json:                                      # 1단계 결과를 다시 쓰고 2단계만 돌린다
        out = json.load(open(out_json))
    if a.refine_from:                                           # time: 다른 json 에서 2단계만 (2026-11-19)
        out = json.load(open(a.refine_from if os.path.isabs(a.refine_from) else os.path.join(ROOT, a.refine_from), encoding="utf-8")); a.refine_from_json = True
        out["refine_from"] = os.path.relpath(a.refine_from if os.path.isabs(a.refine_from) else os.path.join(ROOT, a.refine_from), ROOT)
    for i in range(k if not a.refine_from_json else 0):
        v = V[i]; chain = chain_for(v)
        bounds = [b for n in chain for b in BOX[BODIES[n]]]
        tic = time.time()
        if FORCE_MODE == "time":                                   # 시간 영역: 같은 DE·Powell, 비용만 cost_time
            fun, fargs = cost_time, (fr_s1, chain, v, i)
            fdet = lambda xx: cost_time(xx, fr_s1, chain, v, i, detail=True)   # noqa: E731
        else:
            fun, fargs = cost, (fr, chain, v, c[i], False, i)
            fdet = lambda xx: cost(xx, fr, chain, v, c[i], detail=True, i=i)   # noqa: E731
        for attempt in range(4):                                   # 퇴화(모멘트암 0 / 구 안) 해가 나오면 시드를 바꿔 다시
            res = differential_evolution(fun, bounds, args=fargs, seed=i + 10 * attempt + 1000 * a.seed, maxiter=400, popsize=25, tol=1e-7, polish=False)
            res2 = minimize(fun, res.x, args=fargs, method="Powell", bounds=bounds)
            x = res2.x if (res2.fun < res.fun and fdet(res2.x) is not None) else res.x
            det = fdet(x)
            if det is not None: break
            print(f"  M{i}: 시도 {attempt} 퇴화 해 (cost {res.fun:.1f}) — 시드 바꿔 다시", flush=True)
        if det is None: raise SystemExit(f"M{i}: 경로를 못 찾았다 — 상자/감는 구를 확인")
        r, _ = moment_arms(fr, chain, x.reshape(len(chain), 3)); rm = r.mean(axis=0)
        ang = np.degrees(np.arccos(np.clip(rm @ v / np.linalg.norm(rm), -1, 1)))
        print(f"M{i}: 바디 {[BODIES[n] for n in chain]}  방향오차 {det['dir_err']:.3f} (평균 모멘트암과 설계 방향 각도 {ang:.1f}°)  유효 모멘트암 {det['arm_min']*100:.1f}~{det['arm_mean']*100:.1f} cm  스트로크 {det['stroke']:.2f}  필요 힘 {det['force']:.0f} N  길이 {det['L_mean']*100:.0f} cm  벌점 {det['pen']:.3f}  {time.time()-tic:.0f}s", flush=True)
        print(f"     설계 v = {np.round(v, 2)}   평균 r/|r| = {np.round(rm / np.linalg.norm(rm), 2)}   r 범위(cm) = {np.round(r.min(axis=0)*100, 1)} ~ {np.round(r.max(axis=0)*100, 1)}", flush=True)
        if mods: print(f"     모듈 먼저: ΔL(검사 자세 포함) {det['dL_all']*1000:.1f} / 상한 {MODCAPS['smax'][i]*1000:.1f} mm   힘 {det['force']:.0f} / 상한 {MODCAPS['fcap'][i]:.0f} N", flush=True)
        if FORCE_MODE == "time": print(f"     시간 영역: {td.fmt_size(det['size'])}  브래킷 높이 {[round(h*1000, 1) for h in det['heights']]} mm", flush=True)
        out["muscles"].append(dict(name=f"M{i}", v=v.tolist(), capacity_Nm=float(c[i]), chain=[BODIES[n] for n in chain],
                                   points=x.reshape(len(chain), 3).tolist(), force_N=float(det["force"]), arm_min=det["arm_min"], arm_mean=det["arm_mean"],
                                   stroke=det["stroke"], dir_err=det["dir_err"], angle_deg=float(ang), r_mean=rm.tolist()))
        if FORCE_MODE == "time": out["muscles"][-1]["T_peak_N"] = float(det["size"]["T_peak"])
    if mods:                                                    # 1단계를 새로 돌렸든 json 에서 읽었든, 이번 회의 상한을 기록한다
        out["module_first"] = dict(modules=os.path.relpath(a.modules if os.path.isabs(a.modules) else os.path.join(ROOT, a.modules), ROOT), force_scale=float(mods.get("force_scale", 1.0)),
                                   w_cap_force=W_CAP_F, w_cap_stroke=W_CAP_S, n_check_poses=int(len(MODCAPS["fr2"].p)), refine_only=bool(a.refine_from_json))
        for i, (m, mm) in enumerate(zip(out["muscles"], mods["muscles"])):
            m["module"] = dict(name=mm["module"], bundles=mm["bundles"], force_cap_N=mm["force_cap_N"], stroke_cap_m=mm["stroke_cap_m"], design_force_cap_N=float(MODCAPS["fcap"][i]), design_stroke_cap_m=float(MODCAPS["smax"][i]))
    if FORCE_MODE == "time":
        ts = refine_time(out, fr, T, maxfev=a.maxfev, trust=a.refine_trust, alloc_iters=a.alloc_iters, w_dir=a.w_dir_time)
        for m in out["muscles"]: m["active_length_min_m"] = m["s_max_m"] / 0.40      # 모듈 능동 길이 하한 (최대 수축률 40 %) — 뼈를 따라 놓을 자리가 있는지 확인용
        out["refine"].update(dt=a.dt, n_trials=len(meta["names"]), grades=TDS["grades"], max_bundles=TDS["nb"] or td.MAX_BUNDLES, beta=TDS["curves"][TDS["grades"][0]]["beta"],
                             curves=os.path.relpath(a.curves or td.CURVES_DEFAULT, ROOT), w_work=a.w_work, bracket_H=a.bracket, w_bracket=a.w_bracket if a.bracket is not None else None, arm_min=a.arm_min, len_min=a.len_min, **({'bracket_cad': BRACKET['cad']['json']} if BRACKET and BRACKET.get('cad') else {}))
        np.savez(out_json.replace(".json", "_timeseries.npz"), T=ts["T"], s=ts["s"], L=ts["Ls"], tau=T, Q=Q, trial=meta["trial"], time=meta["time"], names=np.array(meta["names"]))
    else:
        refine(out, fr, T, maxfev=a.maxfev, trust=a.refine_trust)
    if a.refine_trust: out["refine"]["trust"] = a.refine_trust
    if a.seed: out["refine"]["seed"] = a.seed
    json.dump(out, open(out_json, "w"), indent=1, ensure_ascii=False); print("저장", out_json)
    model_out = os.path.join(ROOT, "report/moco/newmus_paths_subject3_model.osim") if legacy else os.path.join(ROOT, f"report/moco/newmus_paths{tag}_model.osim")
    build_model(out, model_path, model_out); plot(out, fr, V, os.path.join(ROOT, f"report/fig_path_design{'' if legacy else tag}.png"))


def batch_boxlsq(A, b, iters=400):
    """프레임마다 min ||A_t u − b_t||², 0 ≤ u ≤ 1 을 한꺼번에 (FISTA, 벡터화). A (T,4,8), b (T,4) → u (T,8)."""
    L = np.linalg.svd(A, compute_uv=False)[:, 0] ** 2 + 1e-9; step = 1.0 / L
    u = np.zeros((len(A), A.shape[2])); y = u.copy(); tk = 1.0
    AtA = np.einsum("tij,tik->tjk", A, A); Atb = np.einsum("tij,ti->tj", A, b)
    for _ in range(iters):
        g = np.einsum("tjk,tk->tj", AtA, y) - Atb
        un = np.clip(y - step[:, None] * g, 0.0, 1.0)
        tn = (1 + np.sqrt(1 + 4 * tk * tk)) / 2; y = un + ((tk - 1) / tn) * (un - u); u, tk = un, tn
    return u


TRUST_ROUNDS = 8                                     # --refine-trust 반복 상한


def refine(out, fr, T, maxfev=20000, w_dir=3.0, trust=None):
    """2단계: 8개 경로를 한꺼번에 조정해, 표본 프레임 전부에서 ID 토크×여유를 u∈[0,1] 로 재현하도록 한다.
    1단계는 근육 하나씩 '방향'만 맞추므로 깊은 스쿼트처럼 모멘트암이 크게 변하는 자세에서 부족분이 생겼다 (9/15).
    변수 = 부착점 + 힘 배율 f_i ≥ 1 (힘 = f_i × 여유 × 용량 / 최소 유효 모멘트암, 즉 설계 용량은 항상 보장).
    비용 = 최악 잔차 + 0.3 평균 잔차 + w_dir × 방향오차(설계 방향 유지) + 작은 힘 벌점 + 스트로크·힘 상한 벌점."""
    chains = [[BODIES.index(b) for b in m["chain"]] for m in out["muscles"]]; sizes = [len(c) * 3 for c in chains]; k = len(chains)
    V = np.array([m["v"] for m in out["muscles"]]); cap = np.array([m["capacity_Nm"] for m in out["muscles"]])
    x0 = np.concatenate([np.array(m["points"]).ravel() for m in out["muscles"]] + [np.ones(k)])
    bounds = [b for m, c in zip(out["muscles"], chains) for n in c for b in BOX[BODIES[n]]] + [(1.0, 3.0)] * k
    pk = np.abs(T).max(axis=0); Tm = MARGIN * T
    offs = np.concatenate([[0], np.cumsum(sizes)[:-1]]).astype(int)    # 케이블별 부착점 변수 시작 위치 (모듈 먼저에서만 쓴다)
    def arms(x):
        R, Ls = [], []; off = 0
        for c, sz in zip(chains, sizes):
            r, l = moment_arms(fr, c, x[off:off + sz].reshape(-1, 3)); R.append(r); Ls.append(l); off += sz
        return np.stack(R, axis=2), np.stack(Ls, axis=1)              # (T,4,8), (T,8)
    def cost_joint(x, detail=False):
        R, Ls = arms(x); f = x[-k:]
        nr = np.linalg.norm(R, axis=1) + 1e-9                            # (T,8)
        e2 = ((R / nr[:, None, :] - V.T[None]) ** 2).sum(axis=1)         # (T,8) 방향오차
        dir_err = 0.5 * e2.mean(axis=0) + 0.5 * e2.max(axis=0)
        arm_eff = np.einsum("tji,ij->ti", R, V); arm_min = np.maximum(arm_eff.min(axis=0), 0.005)
        F = f * MARGIN * cap / arm_min; A = R * F[None, None, :]
        u = batch_boxlsq(A, Tm); res = (np.einsum("tij,tj->ti", A, u) - Tm) / pk
        stroke = (Ls.max(axis=0) - Ls.min(axis=0)) / Ls.mean(axis=0); sma = (Ls.max(axis=0) - Ls.min(axis=0)) / Ls.min(axis=0)
        worst = np.abs(res).max(); mean = np.abs(res).mean()
        pen = 20 * (np.maximum(0.0, stroke - STROKE) ** 2).sum() + 5 * (np.maximum(0.0, F / F_MODULE - 1) ** 2).sum() + 1e4 * (np.maximum(0.0, 0.015 - arm_min) ** 2).sum() + 0.02 * (f - 1).mean() + (50 * (np.maximum(0.0, sma - SMA_EPS) ** 2).sum() if SMA_EPS else 0.0) + (50 * (np.maximum(0.0, (Ls.max(axis=0) - Ls.min(axis=0)) / STROKE_MAX - 1) ** 2).sum() if STROKE_MAX else 0.0)
        dL_all = None
        if MODCAPS is not None:                                          # 모듈 먼저: 케이블별 힘·ΔL 상한 (ΔL 은 stroke_check 자세까지)
            pcs = [cap_penalty(i, F[i], Ls[:, i], c, x[o:o + sz].reshape(-1, 3)) for i, (c, sz, o) in enumerate(zip(chains, sizes, offs))]
            pen = pen + sum(p for p, _ in pcs); dL_all = np.array([dl for _, dl in pcs])
        if detail: return dict(worst_pct=worst * 100, worst_coord=np.abs(res).max(axis=0) * 100, mean_pct=mean * 100, stroke=stroke, sma=sma, dL=(Ls.max(axis=0) - Ls.min(axis=0)), F=F, u_max=u.max(axis=0), R=R, Ls=Ls, arm_min=arm_min, arm_mean=arm_eff.mean(axis=0), dir_err=dir_err, f=f, dL_all=dL_all,
                               total=float(worst + 0.3 * mean + w_dir * dir_err.sum() + pen))
        return worst + 0.3 * mean + w_dir * dir_err.sum() + pen
    d0 = cost_joint(x0, detail=True); print(f"2단계 시작: 최악 잔차 {d0['worst_pct']:.1f} % (좌표별 {np.round(d0['worst_coord'],1)})  평균 {d0['mean_pct']:.2f} %  방향오차합 {d0['dir_err'].sum():.3f}", flush=True)
    tic = time.time(); n = [0]
    def cb(xk): n[0] += 1; d = cost_joint(xk, detail=True); print(f"  반복 {n[0]}: 최악 {d['worst_pct']:.1f} %  평균 {d['mean_pct']:.2f} %  방향오차합 {d['dir_err'].sum():.3f}  힘배율 {np.round(d['f'],2)}  {time.time()-tic:.0f}s", flush=True)
    if trust is None:                                            # 기존: 상자 전체 범위 Powell
        res = minimize(cost_joint, x0, method="Powell", bounds=bounds, callback=cb, options=dict(maxfev=maxfev, xtol=1e-4, ftol=1e-6))
        print(f"  Powell: fun={res.fun:.4f} (시작 {cost_joint(x0):.4f})  nfev={res.nfev}  {res.message}", flush=True)
        x = res.x if res.fun < cost_joint(x0) else x0
    else:                                                        # --refine-trust: 현재 점 ±trust(부착점 m, 힘 배율 ×10) 안에서 Powell, 나아지면 옮겨 반복 (2026-11-05)
        x = x0.copy(); fx = cost_joint(x); nP = len(x) - k
        for rr in range(TRUST_ROUNDS):
            dl = np.r_[np.full(nP, trust), np.full(k, 10 * trust)]
            tb = [(max(b[0], xi - di), min(b[1], xi + di)) for b, xi, di in zip(bounds, x, dl)]
            res = minimize(cost_joint, x, method="Powell", bounds=tb, callback=cb, options=dict(maxfev=maxfev // 4, xtol=1e-4, ftol=1e-6))
            print(f"  Powell 신뢰 구간 {rr + 1}회: fun={res.fun:.4f} (시작 {fx:.4f})  nfev={res.nfev}", flush=True)
            if res.fun < fx - 1e-6: x, fx = res.x, float(res.fun)
            else: break
    repair_log = None
    if MODCAPS is not None:                                  # 모듈 먼저: 벌점이 남긴 작은 초과분을 그 케이블만 국소 탐색으로 상한 안에 넣는다 (수선)
        x, repair_log = repair_caps(x, cost_joint, bounds, offs, sizes, k)
    d = cost_joint(x, detail=True); print(f"2단계 끝: 최악 잔차 {d['worst_pct']:.1f} % (좌표별 {np.round(d['worst_coord'],1)})  평균 {d['mean_pct']:.2f} %  {time.time()-tic:.0f}s", flush=True)
    if repair_log is not None: out["repair"] = repair_log
    off = 0
    for i, (m, c, sz) in enumerate(zip(out["muscles"], chains, sizes)):
        P = x[off:off + sz].reshape(-1, 3); off += sz; r = d["R"][:, :, i]; v = V[i]; rm = r.mean(axis=0)
        m.update(points=P.tolist(), force_N=float(d["F"][i]), force_mult=float(d["f"][i]), arm_min=float(d["arm_min"][i]), arm_mean=float(d["arm_mean"][i]), stroke=float(d["stroke"][i]), sma_ratio=float(d["sma"][i]), dL_m=float(d["dL"][i]), L_min_m=float(d["Ls"][:, i].min()), L_max_m=float(d["Ls"][:, i].max()),
                 dir_err=float(d["dir_err"][i]), angle_deg=float(np.degrees(np.arccos(np.clip(rm @ v / np.linalg.norm(rm), -1, 1)))), r_mean=rm.tolist(), u_max_design=float(d["u_max"][i]))
        print(f"  M{i}: 힘 {m['force_N']:.0f} N (배율 {m['force_mult']:.2f})  유효 모멘트암 {m['arm_min']*100:.1f}~{m['arm_mean']*100:.1f} cm  스트로크 {m['stroke']:.2f}  ΔL {m['dL_m']*1000:.0f} mm (ΔL/L_min {m['sma_ratio']:.2f})  각도차 {m['angle_deg']:.1f}°  u최대 {m['u_max_design']:.2f}", flush=True)
        if d["dL_all"] is not None:
            m["dL_all_m"] = float(d["dL_all"][i])
            print(f"      모듈 먼저: ΔL(검사 자세 포함) {m['dL_all_m']*1000:.1f} / 상한 {MODCAPS['smax'][i]*1000:.1f} mm   힘 {m['force_N']:.0f} / 상한 {MODCAPS['fcap'][i]:.0f} N", flush=True)
    out["refine"] = dict(worst_pct=float(d["worst_pct"]), worst_coord=d["worst_coord"].tolist(), mean_pct=float(d["mean_pct"]), n_frames=int(len(T)), margin=MARGIN, w_dir=w_dir)



def refine_time(out, fr, T, maxfev=20000, w_dir=3.0, trust=None, alloc_iters=300):
    """2단계 (시간 영역): 모든 케이블의 부착점을 한꺼번에 조정. 프레임마다 실제 모멘트암으로 장력을 나누고 (timedomain.allocate,
    min ‖R T − 1.1τ‖² + λ Σ(T/C)², T ≥ 0, C = 1단계 장력 최대), 같은 순간의 수축과 함께 모듈 일 용량을 잰다.
    비용 = 최악 잔차 + 0.3 평균 잔차 + w_dir × 방향오차 + w_work × Σ 일 용량/기준 + 스트로크 비 벌점 + 브래킷. 힘 배율 변수는 없다 (힘 = 장력 최대)."""
    chains = [[BODIES.index(b) for b in m["chain"]] for m in out["muscles"]]; sizes = [len(c) * 3 for c in chains]; k = len(chains)
    V = np.array([m["v"] for m in out["muscles"]])
    C = np.array([m.get("T_peak_N") or m["force_N"] for m in out["muscles"]])
    x0 = np.concatenate([np.array(m["points"]).ravel() for m in out["muscles"]])
    bounds = [b for m, c in zip(out["muscles"], chains) for n in c for b in BOX[BODIES[n]]]
    pk = np.abs(T).max(axis=0); Tm = MARGIN * T
    def arms(x):
        R, Ls = [], []; off = 0
        for c, sz in zip(chains, sizes):
            r, l = moment_arms(fr, c, x[off:off + sz].reshape(-1, 3)); R.append(r); Ls.append(l); off += sz
        return np.stack(R, axis=2), np.stack(Ls, axis=1)
    def cost_joint(x, detail=False, iters=None):
        R, Ls = arms(x)
        Ta, res, _ = td.allocate(R, Tm, C, iters=iters or alloc_iters)     # 매번 0 에서 (이어 풀면 같은 점의 비용이 호출마다 달라 Powell 이 헷갈린다)
        nr = np.linalg.norm(R, axis=1) + 1e-9; e2 = ((R / nr[:, None, :] - V.T[None]) ** 2).sum(axis=1)
        act = Ta > 0.01 * np.maximum(Ta.max(axis=0), 1e-9)[None, :]            # 방향오차는 그 케이블이 당기는 프레임에서만
        dir_err = np.array([0.5 * e2[act[:, i], i].mean() + 0.5 * e2[act[:, i], i].max() if act[:, i].sum() >= 3 else 0.5 * e2[:, i].mean() + 0.5 * e2[:, i].max() for i in range(k)])
        res = res / pk; s = Ls.max(axis=0) - Ls
        zs = [td.size_cable(Ta[:, i], s[:, i], TDS["curves"], TDS["grades"], TDS["nb"]) for i in range(k)]
        work = TDS["w_work"] * sum(z["obj"] / TDS["Wref"][i] for i, z in enumerate(zs))
        stroke = (Ls.max(axis=0) - Ls.min(axis=0)) / Ls.mean(axis=0); sma = (Ls.max(axis=0) - Ls.min(axis=0)) / Ls.min(axis=0)
        pen = (50 * (np.maximum(0.0, sma - SMA_EPS) ** 2).sum() if SMA_EPS else 0.0)           # 길이 비 벌점 없음 (cost_time 주석)
        if TDS.get("arm_min"): pen += sum(arm_min_pen((R[:, :, i] @ V[i])[act[:, i]]) for i in range(k))
        if TDS.get("len_min") or TDS.get("len_max"): pen += sum(len_min_pen(Ls[:, i]) for i in range(k))
        bp, hs = 0.0, []
        if BRACKET:
            off = 0
            for c, sz in zip(chains, sizes):
                p_, h_ = bracket_terms(fr, c, x[off:off + sz].reshape(-1, 3)); bp += p_; hs.append(h_); off += sz
        worst = np.abs(res).max(); mean = np.abs(res).mean()
        total = worst + 0.3 * mean + w_dir * dir_err.sum() + work + pen + bp + TDS["w_h"] * sum(max(0.0, h) for hh in hs for h in hh)
        if detail: return dict(worst_pct=worst * 100, worst_coord=np.abs(res).max(axis=0) * 100, mean_pct=mean * 100, R=R, Ls=Ls, T=Ta, s=s, sizes=zs, dir_err=dir_err,
                               stroke=stroke, sma=sma, heights=hs, work=work, total=float(total))
        return total
    d0 = cost_joint(x0, detail=True)
    print(f"2단계(시간 영역) 시작: 최악 잔차 {d0['worst_pct']:.1f} %  평균 {d0['mean_pct']:.2f} %  일 용량 항 {d0['work']:.3f}  방향오차합 {d0['dir_err'].sum():.3f}  "
          f"모듈 {[(z['grade'], z['n']) for z in d0['sizes']]}", flush=True)
    tic = time.time(); n = [0]
    def cb(xk):
        n[0] += 1; d = cost_joint(xk, detail=True)
        print(f"  반복 {n[0]}: 최악 {d['worst_pct']:.1f} %  평균 {d['mean_pct']:.2f} %  일 용량 항 {d['work']:.3f}  방향오차합 {d['dir_err'].sum():.3f}  "
              f"장력 최대 {[round(z['T_peak']) for z in d['sizes']]}  {time.time()-tic:.0f}s", flush=True)
    if trust is None:
        res = minimize(cost_joint, x0, method="Powell", bounds=bounds, callback=cb, options=dict(maxfev=maxfev, xtol=1e-4, ftol=1e-6))
        f0 = cost_joint(x0); print(f"  Powell: fun={res.fun:.4f} (시작 {f0:.4f})  nfev={res.nfev}  {res.message}", flush=True)
        x = res.x if res.fun < f0 else x0
    else:
        x = x0.copy(); fx = cost_joint(x)
        for rr in range(TRUST_ROUNDS):
            tb = [(max(b[0], xi - trust), min(b[1], xi + trust)) for b, xi in zip(bounds, x)]
            res = minimize(cost_joint, x, method="Powell", bounds=tb, callback=cb, options=dict(maxfev=maxfev // 4, xtol=1e-4, ftol=1e-6))
            print(f"  Powell 신뢰 구간 {rr + 1}회: fun={res.fun:.4f} (시작 {fx:.4f})  nfev={res.nfev}", flush=True)
            if res.fun < fx - 1e-6: x, fx = res.x, float(res.fun)
            else: break
    d = cost_joint(x, detail=True, iters=1500)
    print(f"2단계(시간 영역) 끝: 최악 잔차 {d['worst_pct']:.1f} % (좌표별 {np.round(d['worst_coord'],1)})  평균 {d['mean_pct']:.2f} %  {time.time()-tic:.0f}s", flush=True)
    off = 0
    for i, (m, c, sz) in enumerate(zip(out["muscles"], chains, sizes)):
        P = x[off:off + sz].reshape(-1, 3); off += sz; r = d["R"][:, :, i]; v = V[i]; rm = r.mean(axis=0); z = d["sizes"][i]; arm_eff = r @ v
        m.update(points=P.tolist(), force_N=float(z["T_peak"]), T_peak_N=float(z["T_peak"]), force_mode="time", arm_min=float(arm_eff.min()), arm_mean=float(arm_eff.mean()),
                 stroke=float(d["stroke"][i]), sma_ratio=float(d["sma"][i]), dL_m=float(d["Ls"][:, i].max() - d["Ls"][:, i].min()), s_max_m=float(z["s_max"]),
                 L_min_m=float(d["Ls"][:, i].min()), L_max_m=float(d["Ls"][:, i].max()), dir_err=float(d["dir_err"][i]),
                 angle_deg=float(np.degrees(np.arccos(np.clip(rm @ v / np.linalg.norm(rm), -1, 1)))), r_mean=rm.tolist(),
                 module=dict(name=z["grade"], bundles=z["n"], n_cont=z["n_cont"], work_J=z["work_J"], mass_park_g=z["mass_park_g"], feasible=z["feasible"]))
        if BRACKET: m["bracket_h_m"] = [float(h) for h in d["heights"][i]]
        print(f"  M{i}: {td.fmt_size(z)}  유효 모멘트암 {m['arm_min']*100:.1f}~{m['arm_mean']*100:.1f} cm  각도차 {m['angle_deg']:.1f}°"
              + (f"  브래킷 {[round(h*1000, 1) for h in m['bracket_h_m']]} mm" if BRACKET else ""), flush=True)
    out["refine"] = dict(worst_pct=float(d["worst_pct"]), worst_coord=d["worst_coord"].tolist(), mean_pct=float(d["mean_pct"]), n_frames=int(len(T)), margin=MARGIN, w_dir=w_dir, force_mode="time")
    return d


REPAIR_TOL = 0.15        # 모듈 먼저 수선: 이보다 많이 넘는 케이블은 그 등급으로는 기하학적으로 안 되는 것 → 손대지 않고 그대로 보고
REPAIR_STEP = (0.002, 0.02)   # 수선 국소 탐색(Nelder–Mead) 시작 크기: 부착점 m, 힘 배율
REPAIR_W = 1e4           # 수선 목적함수의 상한 벽 (국소 탐색이라 벽이 가팔라도 된다 — Powell 상자 전체 선탐색과 다름)


def repair_caps(x, cost_joint, bounds, offs, sizes, k):
    """2단계 해 x 에서 배정 등급 상한(힘: 등급 상한 ÷ 힘 배율, ΔL: 설계+검사 자세)을 조금(≤ REPAIR_TOL) 넘은 케이블마다
    그 케이블의 부착점과 힘 배율만 국소 탐색해 상한 안(× 0.999)으로 넣는다. 목적 = 2단계 비용 + 벽. 반환 (x, 기록)."""
    lo = np.array([b[0] for b in bounds]); hi = np.array([b[1] for b in bounds])
    fT, sT = MODCAPS["fcap_true"] * 0.999, MODCAPS["smax_true"] * 0.999
    def viol(xx):
        d = cost_joint(xx, detail=True); return d, d["F"] / fT - 1, d["dL_all"] / sT - 1
    d0, vf, vs = viol(x); log = dict(before_worst_pct=float(d0["worst_pct"]), before_mean_pct=float(d0["mean_pct"]), cables=[])
    for i in range(k):
        if max(vf[i], vs[i]) <= 0: continue
        e = dict(name=f"M{i}", force_before_N=float(d0["F"][i]), dL_before_mm=float(d0["dL_all"][i] * 1000), force_over_pct=float(vf[i] * 100), dL_over_pct=float(vs[i] * 100))
        if max(vf[i], vs[i]) > REPAIR_TOL:
            e["result"] = "수선 안 함 (15 % 넘게 초과)"; log["cables"].append(e); print(f"  수선 M{i}: 힘 {vf[i]*100:+.1f} %, ΔL {vs[i]*100:+.1f} % — 15 % 넘게 넘어 손대지 않음", flush=True); continue
        idx = np.r_[offs[i]:offs[i] + sizes[i], len(x) - k + i]
        def obj(p, i=i, idx=idx):
            xx = x.copy(); xx[idx] = np.clip(p, lo[idx], hi[idx]); d, a_, b_ = viol(xx)
            return d["total"] + REPAIR_W * (max(0.0, a_[i]) ** 2 + max(0.0, b_[i]) ** 2)
        p0 = x[idx].copy(); step = np.array([REPAIR_STEP[0]] * sizes[i] + [REPAIR_STEP[1]]); sim = [p0]
        for j in range(len(p0)):
            q = p0.copy(); q[j] = q[j] + step[j] if q[j] + step[j] <= hi[idx][j] else q[j] - step[j]; sim.append(q)
        r = minimize(obj, p0, method="Nelder-Mead", bounds=list(zip(lo[idx], hi[idx])), options=dict(initial_simplex=np.array(sim), maxfev=3000, xatol=1e-6, fatol=1e-7))
        x = x.copy(); x[idx] = np.clip(r.x, lo[idx], hi[idx]); d1, vf, vs = viol(x)
        ok = bool(d1["F"][i] <= MODCAPS["fcap_true"][i] and d1["dL_all"][i] <= MODCAPS["smax_true"][i])     # 판정은 실제 등급 상한 (목표 × 0.999 가 아니라)
        e.update(force_after_N=float(d1["F"][i]), dL_after_mm=float(d1["dL_all"][i] * 1000), ok=ok, nfev=int(r.nfev), worst_after_pct=float(d1["worst_pct"]), mean_after_pct=float(d1["mean_pct"]))
        e["result"] = "상한 안으로 넣음" if e["ok"] else "다 못 넣음"; log["cables"].append(e)
        print(f"  수선 M{i}: 힘 {e['force_before_N']:.0f}→{e['force_after_N']:.0f} N, ΔL {e['dL_before_mm']:.1f}→{e['dL_after_mm']:.1f} mm  {e['result']}  (잔차 최악 {d1['worst_pct']:.1f} %, 평균 {d1['mean_pct']:.2f} %, {r.nfev}회)", flush=True)
    d2 = cost_joint(x, detail=True); log.update(after_worst_pct=float(d2["worst_pct"]), after_mean_pct=float(d2["mean_pct"]))
    return x, log


def remove_design_actuators(model):
    """모델에 원래 있던 설계 좌표의 토크 모터(shoulder_flex_r, lumbar_ext …)를 뗀다 — 남겨 두면 새 근육 대신 그게 토크를 낸다.
    하체 모델에는 다리 토크 모터가 없어서 아무것도 안 뗀다."""
    import opensim as osim
    targets = {ds.full(b, s) for b in SP["base"] for s in SP["sides"]}
    fs = model.updForceSet(); removed = []
    for i in reversed(range(fs.getSize())):
        ca = osim.CoordinateActuator.safeDownCast(fs.get(i))
        if ca is not None and ca.get_coordinate() in targets:
            removed.append(ca.getName()); fs.remove(i)
    if removed: print(f"  설계 좌표의 토크 모터 제거: {removed}")
    model.finalizeConnections(); model.initSystem()
    return removed


def euler_xyz_body(Rm):
    """회전 행렬 R = Rx(a)·Ry(b)·Rz(c) → (a, b, c) (OpenSim xyz_body_rotation 과 같은 몸체 고정 X-Y-Z 순서)."""
    Rm = np.asarray(Rm, float); b = np.arcsin(np.clip(Rm[0, 2], -1.0, 1.0))
    return float(np.arctan2(-Rm[1, 2], Rm[2, 2])), float(b), float(np.arctan2(-Rm[0, 1], Rm[0, 0]))


def add_wrap_cylinder(model, wc, name, side, sgn):
    """json 의 wrap_cyl 을 OpenSim WrapCylinder 로 그 뼈(왼쪽은 거울상)에 붙이고 돌려준다. 원통 좌표계 = wc["frame"] (열 = x, y, z, z = 축).
    quadrant (기본 "+x") 로 감는 쪽을 고정한다 — 설계 계산(cylinder_path)의 감김 방향 s 와 같은 쪽이 원통 좌표계 +x 반쪽에 오도록 frame 의 x 를 정한다.
    왼쪽: 점처럼 z 부호를 뒤집는다 (R' = M R M, M = diag(1, 1, −1) — 오른손 좌표계가 유지되고 x 반쪽도 거울상으로 간다)."""
    import opensim as osim
    Mz = np.diag([1.0, 1.0, sgn]); Rm = Mz @ np.asarray(wc["frame"], float) @ Mz; c = Mz @ np.asarray(wc["center"], float)
    w = osim.WrapCylinder(); w.setName(name); w.set_radius(float(wc["radius"])); w.set_length(float(wc.get("length", 0.06)))
    w.set_translation(osim.Vec3(*c)); w.set_xyz_body_rotation(osim.Vec3(*euler_xyz_body(Rm))); w.set_quadrant(wc.get("quadrant", "+x"))
    model.getBodySet().get(ds.mirror_body(SP, wc["body"], side)).addWrapObject(w)
    return w


def build_model(out, model_path=None, out_path=None):
    """뼈대 모델(근육 제거)에 PathActuator k개 (× 좌우, mirror 공간이면 왼쪽은 z 부호 반전)."""
    import opensim as osim
    model_path = model_path or os.path.join(ROOT, "data/subject3/subject3_scaled.osim")
    out_path = out_path or os.path.join(ROOT, "report/moco/newmus_paths_subject3_model.osim")
    proc = osim.ModelProcessor(model_path); proc.append(osim.ModOpRemoveMuscles()); model = proc.process()
    remove_design_actuators(model)
    colors = [(0.9, 0.2, 0.2), (0.2, 0.6, 0.9), (0.2, 0.8, 0.3), (0.9, 0.6, 0.1), (0.6, 0.3, 0.8), (0.1, 0.7, 0.7), (0.8, 0.3, 0.6), (0.5, 0.5, 0.1), (0.3, 0.3, 0.3), (0.7, 0.5, 0.3)]
    sides = (("r", 1.0), ("l", -1.0)) if SP["mirror"] else (("", 1.0),)
    for side, sgn in sides:
        wrap_obj = None
        if WRAP:                                                                  # 감는 구를 그 뼈에 붙인다 (좌우 하나씩)
            wbody = model.getBodySet().get(ds.mirror_body(SP, WRAP["body"], side)); c = WRAP["center"]
            wrap_obj = osim.WrapSphere(); wrap_obj.setName(f"wrap_{WRAP['body'][:-2]}_{side}" if side else f"wrap_{WRAP['body']}")
            wrap_obj.set_radius(float(WRAP["radius"])); wrap_obj.set_translation(osim.Vec3(c[0], c[1], sgn * c[2])); wbody.addWrapObject(wrap_obj)
        for i, m in enumerate(out["muscles"]):
            nm = f"{m['name']}_{side}" if side else m["name"]
            pa = osim.PathActuator(); pa.setName(nm); pa.setOptimalForce(m["force_N"]); pa.setMinControl(0.0); pa.setMaxControl(1.0)
            for kk, (bn, pt) in enumerate(zip(m["chain"], m["points"])):
                pa.addNewPathPoint(f"{nm}_p{kk}", model.getBodySet().get(ds.mirror_body(SP, bn, side)), osim.Vec3(pt[0], pt[1], sgn * pt[2]))
            if wrap_obj is not None and wrap_segment([BODIES.index(b) for b in m["chain"]]) is not None:
                pa.updGeometryPath().addPathWrap(wrap_obj)
            if m.get("wrap_cyl"):                                                 # 케이블별 원통 감김 (2026-11-19, 노뼈 회내외 드럼): 그 케이블만 감는다
                pa.updGeometryPath().addPathWrap(add_wrap_cylinder(model, m["wrap_cyl"], f"wrapcyl_{nm}", side, sgn))
            pa.updGeometryPath().setDefaultColor(osim.Vec3(*colors[i % len(colors)]))
            model.addForce(pa)
    model.finalizeConnections(); model.initSystem()
    model.printToXML(out_path); print("모델 저장", out_path)


def plot(out, fr, V, png=None):
    k = len(out["muscles"]); ncol = 4 if k > 4 else k; nrow = int(np.ceil(k / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.5 * nrow), squeeze=False)
    labels = list(SP["labels"])
    for i, (ax, m) in enumerate(zip(axes.ravel(), out["muscles"])):
        chain = [BODIES.index(b) for b in m["chain"]]; r, _ = moment_arms(fr, chain, np.array(m["points"]), wrap_cyl=m.get("wrap_cyl"))
        x = np.arange(len(COORDS)); w = 0.38
        ax.bar(x - w / 2, np.array(m["v"]) * m["arm_mean"] * 100, w, color="0.6", label="설계 방향 × 평균 모멘트암")
        ax.bar(x + w / 2, r.mean(axis=0) * 100, w, color="#1e8449", label="경로의 평균 모멘트암")
        ax.errorbar(x + w / 2, r.mean(axis=0) * 100, yerr=[r.mean(axis=0) * 100 - r.min(axis=0) * 100, r.max(axis=0) * 100 - r.mean(axis=0) * 100], fmt="none", ecolor="k", lw=0.8)
        ax.axhline(0, color="k", lw=0.6); ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=8)
        ax.set_title(f"{m['name']}: {'→'.join(b.replace('_r','') for b in m['chain'])}  각도차 {m['angle_deg']:.0f}°  {m['force_N']:.0f} N  스트로크 {m['stroke']:.2f}", fontsize=9)
        if i == 0: ax.legend(fontsize=7); ax.set_ylabel("모멘트암 (cm)")
    fig.suptitle("케이블 경로 설계 — 설계 방향 vs 경로가 실제로 내는 모멘트암 (막대 = 자세 평균, 선 = 자세에 따른 범위)", fontsize=11)
    fig.tight_layout(); fig.savefig(png or os.path.join(ROOT, "report/fig_path_design.png"), dpi=150); print("그림 저장")


if __name__ == "__main__":
    main()
