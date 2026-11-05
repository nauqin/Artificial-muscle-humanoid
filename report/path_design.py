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


def path_lengths(fr, chain, P, with_inside=False):
    """chain: 바디 인덱스 목록, P: (len(chain), 3) 바디 좌표. 반환 L (poses, nP). 감는 구가 걸린 구간은 구 표면 경로로 잰다."""
    pts = np.stack([np.einsum("tkij,j->tki", fr.R[:, :, b], P[n]) + fr.p[:, :, b] for n, b in enumerate(chain)], axis=2)  # (T, nP, n, 3)
    seg = np.linalg.norm(np.diff(pts, axis=2), axis=3)                                                             # (T, nP, n-1)
    inside = np.zeros(seg.shape[:2], bool); ws = wrap_segment(chain)
    if ws is not None:
        wb = BODIES.index(WRAP["body"]); C = np.einsum("tkij,j->tki", fr.R[:, :, wb], np.asarray(WRAP["center"], float)) + fr.p[:, :, wb]
        seg[:, :, ws], inside = sphere_path(pts[:, :, ws], pts[:, :, ws + 1], C, float(WRAP["radius"]))
    L = seg.sum(axis=2)
    return (L, inside) if with_inside else L


def moment_arms(fr, chain, P):
    L = path_lengths(fr, chain, P)
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
    ap.add_argument("--stroke-max", type=float, default=None, help="모듈 스트로크 상한 m: 케이블 길이 변화 (L_max−L_min) ≤ 이 값 (계획서 모듈표 M1 0.12, M2·M6 0.08, M3 0.06)")
    ap.add_argument("--sma-eps", type=float, default=None, help="SMA 묶음 제약: 케이블 (L_max−L_min)/L_min 상한 (park2019 수축률 0.5 ÷ 여유 1.2 = 0.417). 기본 없음")
    ap.add_argument("--modules", default=None, help="모듈 먼저 (2026-11-05): report/module_first.py 가 만든 케이블별 등급 json. 주면 --f-module/--stroke-max 대신 "
                    "케이블마다 그 등급의 힘(묶음 수 × 힘 상한 ÷ 힘 배율)·스트로크 상한을 1 % 안쪽에 세게 건다. ΔL 은 stroke_check 자세까지 합쳐 잰다")
    ap.add_argument("--seed", type=int, default=0, help="1단계 차분진화 시드 묶음 (시드 = 근육 번호 + 10×시도 + 1000×이 값). 0 이면 기존과 같다 — 반복 실행으로 결과 흔들림을 잴 때 (2026-11-05)")
    ap.add_argument("--refine-trust", type=float, default=None, help="2단계 Powell 을 현재 점 ± 이 값(부착점 m, 힘 배율은 ×10) 안에서 돌리고 옮겨 가며 반복. "
                    "기본(없음)은 상자 전체 범위 — 1단계 해가 상한 경계에 붙어 있으면 상자 전체 선탐색이 시작점보다 나쁜 점을 돌려줘 2단계가 아무것도 못 한다 (모듈 먼저, 2026-11-05)")
    a = ap.parse_args(); configure(a.space)
    if a.modules and (a.f_module or a.stroke_max): print("  --modules 를 주면 --f-module/--stroke-max 는 쓰지 않는다 (케이블별 상한으로 바뀐다)"); a.f_module = a.stroke_max = None
    if a.f_module: global F_MODULE; F_MODULE = a.f_module
    if a.sma_eps: global SMA_EPS; SMA_EPS = a.sma_eps; print(f"SMA 제약 ΔL/L_min ≤ {SMA_EPS:.3f}")
    if a.stroke_max: global STROKE_MAX; STROKE_MAX = a.stroke_max; print(f"모듈 스트로크 제약 ΔL ≤ {STROKE_MAX*1000:.0f} mm")
    if a.box_scale != 1.0:                                      # 로봇 크기에 맞춰 상자·감는 구를 줄인다
        global BOX, WRAP; BOX = {b: tuple(tuple(v * a.box_scale for v in ax) for ax in box) for b, box in BOX.items()}
        if WRAP: WRAP = dict(WRAP, center=tuple(v * a.box_scale for v in WRAP["center"]), radius=WRAP["radius"] * a.box_scale)
    if a.wrap_radius and WRAP: WRAP = dict(WRAP, radius=a.wrap_radius); print(f"감는 구 반지름 {a.wrap_radius} m")
    legacy = (a.space == "lower" and not a.tasks and not a.out_tag)
    k = a.k or SP["k"]; tag = ("_k8" if legacy else f"_{a.space}_k{k}") + a.out_tag
    design = a.design or os.path.join(ROOT, f"report/newmus_V_design{tag}.npz")
    d = np.load(design); V, c = d["V"], d["c"]; k = len(V); tag = ("_k8" if legacy else f"_{a.space}_k{k}") + a.out_tag   # k 는 설계 파일이 정한다
    model_path = a.model or os.path.join(ROOT, "data/subject3/subject3_scaled.osim")
    out_json = a.out or os.path.join(ROOT, f"report/newmus_paths{tag}.json")
    tasks_path = (os.path.join(ROOT, a.tasks) if not os.path.isabs(a.tasks) else a.tasks) if a.tasks else None
    Q, T = poses_from_tasks(tasks_path, ref=a.ref) if a.tasks else poses(ref=a.ref)
    print(f"[{a.space}] 근육 {k}개, 자세 표본 {len(Q)}개, 모델 {model_path}  (2단계 토크 좌표별 최대 {np.round(np.abs(T).max(axis=0), 1)} N·m)", flush=True)
    tic = time.time(); fr = Frames(model_path, Q); print(f"뼈 변환 저장 {time.time()-tic:.0f}s", flush=True)
    out = {"space": a.space, "coords": COORDS, "bodies": BODIES, "model": os.path.relpath(model_path, ROOT), "muscles": []}
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
    for i in range(k if not a.refine_from_json else 0):
        v = V[i]; chain = chain_for(v)
        bounds = [b for n in chain for b in BOX[BODIES[n]]]
        tic = time.time()
        for attempt in range(4):                                   # 퇴화(모멘트암 0 / 구 안) 해가 나오면 시드를 바꿔 다시
            res = differential_evolution(cost, bounds, args=(fr, chain, v, c[i], False, i), seed=i + 10 * attempt + 1000 * a.seed, maxiter=400, popsize=25, tol=1e-7, polish=False)
            res2 = minimize(cost, res.x, args=(fr, chain, v, c[i], False, i), method="Powell", bounds=bounds)
            x = res2.x if (res2.fun < res.fun and cost(res2.x, fr, chain, v, c[i], detail=True, i=i) is not None) else res.x
            det = cost(x, fr, chain, v, c[i], detail=True, i=i)
            if det is not None: break
            print(f"  M{i}: 시도 {attempt} 퇴화 해 (cost {res.fun:.1f}) — 시드 바꿔 다시", flush=True)
        if det is None: raise SystemExit(f"M{i}: 경로를 못 찾았다 — 상자/감는 구를 확인")
        r, _ = moment_arms(fr, chain, x.reshape(len(chain), 3)); rm = r.mean(axis=0)
        ang = np.degrees(np.arccos(np.clip(rm @ v / np.linalg.norm(rm), -1, 1)))
        print(f"M{i}: 바디 {[BODIES[n] for n in chain]}  방향오차 {det['dir_err']:.3f} (평균 모멘트암과 설계 방향 각도 {ang:.1f}°)  유효 모멘트암 {det['arm_min']*100:.1f}~{det['arm_mean']*100:.1f} cm  스트로크 {det['stroke']:.2f}  필요 힘 {det['force']:.0f} N  길이 {det['L_mean']*100:.0f} cm  벌점 {det['pen']:.3f}  {time.time()-tic:.0f}s", flush=True)
        print(f"     설계 v = {np.round(v, 2)}   평균 r/|r| = {np.round(rm / np.linalg.norm(rm), 2)}   r 범위(cm) = {np.round(r.min(axis=0)*100, 1)} ~ {np.round(r.max(axis=0)*100, 1)}", flush=True)
        if mods: print(f"     모듈 먼저: ΔL(검사 자세 포함) {det['dL_all']*1000:.1f} / 상한 {MODCAPS['smax'][i]*1000:.1f} mm   힘 {det['force']:.0f} / 상한 {MODCAPS['fcap'][i]:.0f} N", flush=True)
        out["muscles"].append(dict(name=f"M{i}", v=v.tolist(), capacity_Nm=float(c[i]), chain=[BODIES[n] for n in chain],
                                   points=x.reshape(len(chain), 3).tolist(), force_N=float(det["force"]), arm_min=det["arm_min"], arm_mean=det["arm_mean"],
                                   stroke=det["stroke"], dir_err=det["dir_err"], angle_deg=float(ang), r_mean=rm.tolist()))
    if mods:                                                    # 1단계를 새로 돌렸든 json 에서 읽었든, 이번 회의 상한을 기록한다
        out["module_first"] = dict(modules=os.path.relpath(a.modules if os.path.isabs(a.modules) else os.path.join(ROOT, a.modules), ROOT), force_scale=float(mods.get("force_scale", 1.0)),
                                   w_cap_force=W_CAP_F, w_cap_stroke=W_CAP_S, n_check_poses=int(len(MODCAPS["fr2"].p)), refine_only=bool(a.refine_from_json))
        for i, (m, mm) in enumerate(zip(out["muscles"], mods["muscles"])):
            m["module"] = dict(name=mm["module"], bundles=mm["bundles"], force_cap_N=mm["force_cap_N"], stroke_cap_m=mm["stroke_cap_m"], design_force_cap_N=float(MODCAPS["fcap"][i]), design_stroke_cap_m=float(MODCAPS["smax"][i]))
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
            pa.updGeometryPath().setDefaultColor(osim.Vec3(*colors[i % len(colors)]))
            model.addForce(pa)
    model.finalizeConnections(); model.initSystem()
    model.printToXML(out_path); print("모델 저장", out_path)


def plot(out, fr, V, png=None):
    k = len(out["muscles"]); ncol = 4 if k > 4 else k; nrow = int(np.ceil(k / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4 * ncol, 3.5 * nrow), squeeze=False)
    labels = list(SP["labels"])
    for i, (ax, m) in enumerate(zip(axes.ravel(), out["muscles"])):
        chain = [BODIES.index(b) for b in m["chain"]]; r, _ = moment_arms(fr, chain, np.array(m["points"]))
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
