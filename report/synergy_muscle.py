"""상체 근육 시너지 분석 → 관절 배분비 → 직접 최적화 설계와 비교 (2026-11-06).

로봇 뼈대(robot40_h160)에는 상체 근육이 없어서, 근육 경로만 공개 모델에서 빌린다. 토크는 로봇 뼈대의 ID 그대로다.
  몸통: data/synergy_models/gait2392_simbody.osim 의 ercspn·intobl·extobl 좌우 6개 (OpenSim 공식 모델, 좌표 이름 같음)
  팔  : data/synergy_models/MOBL_ARMS_fixed_41.osim 의 어깨·팔꿈치 근육 24개 (MoBL-ARMS, Inria opensim_models 사본). 손목·손가락 근육은 뺀다.
        어깨 좌표가 달라(elv_angle·shoulder_elv·shoulder_rot) 위팔 자세를 회전행렬로 바꿔 넣고, 모멘트암은 우리 좌표로 근육 길이를 중앙차분해 잰다 (r = −∂L/∂q).
  참고: --space arm_10 은 full_upper_body_marks.osim 의 오른팔 10개. 근육을 다 써도 어깨 굴곡·안쪽 회전 토크를 못 내고 이두근 굴곡 부호가 반대라 쓰지 않는다.

순서 (spanning_check.py 의 하체 방법과 같다)
  1. tasks 의 자세·ID 토크를 뽑는다. 왼쪽 좌표는 오른팔 자세로 바꿔 넣는다(design_directions 가 좌우를 같이 쌓는 것과 같다).
  2. 자세마다 min Σa² s.t. R(q)·(a·Fmax) = τ, 0 ≤ a ≤ 1 (정적 최적화). 힘-길이·힘-속도는 넣지 않는다(F = a·Fmax).
  3. 근육 활성 a 를 근육 피크로 나눠 NMF, A ≈ C·W. VAF 로 개수를 본다.
  4. 근육 힘을 시너지에 귀속: share_ij = C_i W_ij / Σ_k C_k W_kj. 시너지 i 의 평균 토크 벡터를 단위벡터로 = 관절 배분비.
  5. design_directions 와 같은 LP·같은 평가 표본으로 총 토크 용량을 잰다.

한계: 근육 수가 적고(팔 10, 몸통 6) 모델 크기가 로봇과 다르다. 배분비는 단위벡터라 크기는 영향이 작다.

    python -I report/synergy_muscle.py --space trunk --check-npz report/newmus_V_design_trunk_k4_40kg_h160.npz
    python -I report/synergy_muscle.py --space arm_mobl --also shoulder_r --check-also report/newmus_V_design_shoulder_r_k6_40kg_h160_v635.npz
→ report/synergy_muscle_<space>_40kg_h160.csv, report/newmus_V_design_<space>_k<k>_40kg_h160_musyn.npz, out/synergy_upper/<space>_poses.npz (자세·모멘트암 저장)
"""
import argparse, csv, os, sys, time
import numpy as np
from scipy.optimize import lsq_linear

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import addb_io as io                                   # noqa: E402
import design_directions as dd                         # noqa: E402
import design_space as ds                              # noqa: E402
import synergy as sy                                   # noqa: E402
import synergy_torque as st                            # noqa: E402

ARM = ("ant_delt_r", "medial_delt_r", "post_delt_r", "TRIlong", "TRIlat", "BIClong", "BRD", "PT", "LAT2", "PECM2")
TRUNK = ("ercspn_r", "ercspn_l", "intobl_r", "intobl_l", "extobl_r", "extobl_l")
MOBL = ("DELT1", "DELT2", "DELT3", "SUPSP", "INFSP", "SUBSC", "TMIN", "TMAJ", "PECM1", "PECM2", "PECM3", "LAT1", "LAT2", "LAT3", "CORB",
        "TRIlong", "TRIlat", "TRImed", "ANC", "BIClong", "BICshort", "BRA", "BRD", "PT")
CONFIG = {"trunk": dict(space="trunk", model="data/synergy_models/gait2392_simbody.osim", muscles=TRUNK),
          "arm_10": dict(space="arm_r", model="data/synergy_models/full_upper_body_marks.osim", muscles=ARM),
          "arm_mobl": dict(space="arm_r", model="data/synergy_models/MOBL_ARMS_fixed_41.osim", muscles=MOBL, mobl=True)}
_a = np.array([0.0048, 0.9991, 0.0424]); _a /= np.linalg.norm(_a)      # MoBL elv_angle·shoulder_rot 축 (땅, OpenSim 에서 잼)
_b = np.array([-0.9983, 0.0023, 0.0589]); _x = -(_b - (_b @ _a) * _a); _x /= np.linalg.norm(_x)  # shoulder_elv 축, 수직으로 맞춘 −x
MOBL_B = np.column_stack([_x, _a, np.cross(_x, _a)])                     # 이 기저에서 elv·rot = y 축, elev = −x 축
W_TAU = 1e3                                            # 토크 맞춤 가중치 (활성 제곱 합 대비)


def load_tau(space):
    dd.SP = ds.get(space); dd.COORDS = dd.SP["base"]
    data = dd.dataset_from_tasks(TASKS); tau_all = np.vstack(list(data.values()))
    return tau_all, dd.subsample(tau_all, 3000, seed=7)


def sample_poses(space, n_frames):
    """자세 Q (N, 설계 좌표 + extra, 오른쪽 이름 기준) 와 ID 토크 T (N, 설계 좌표). trial 마다 고르게 + 좌표별 극값."""
    SP = ds.get(space); base = SP["base"]; rows = [r for r in csv.DictReader(open(TASKS, encoding="utf-8")) if r.get("id") and not r["id"].startswith("#")]
    total = 0; items = []
    for row in rows:
        m = io.read_mot(os.path.join(ROOT, row["ik"])); idm = io.read_mot(os.path.join(ROOT, row["id"]))
        if m.in_degrees: m.data = np.radians(m.data)
        ok = (m.time >= idm.time[0]) & (m.time <= idm.time[-1])
        items.append((row, m, idm, ok)); total += int(ok.sum()) * len(SP["sides"])
    stride = max(1, total // n_frames); Q, T = [], []
    for row, m, idm, ok in items:
        tsel = m.time[ok]
        for side in SP["sides"]:
            cols = [ds.full(b, side) for b in base]; extra = [e.replace("_r", f"_{side}") if side else e for e in SP["extra"]]
            if not all(c in m.names for c in cols + extra) or not all(c + "_moment" in idm.names for c in cols): print(f"  [건너뜀] {row['trial']} {side}"); continue
            idx = list(range(0, len(tsel), stride))
            for c in cols:
                q = m.column(c)[ok]; idx += [int(np.argmax(q)), int(np.argmin(q))]
                tt = np.interp(tsel, idm.time, idm.column(c + "_moment")); idx += [int(np.argmax(tt)), int(np.argmin(tt))]
            t = tsel[sorted(set(idx))]
            Q.append(np.column_stack([np.interp(t, m.time, m.column(c)) for c in cols + extra]))
            T.append(np.column_stack([np.interp(t, idm.time, idm.column(c + "_moment")) for c in cols]))
    Q, T = np.vstack(Q), np.vstack(T)
    names = [ds.full(b, SP["design_side"]) for b in base] + list(SP["extra"])
    for cname, (lo, hi) in SP.get("pose_limit", {}).items():
        j = names.index(cname); keep = (Q[:, j] >= lo) & (Q[:, j] <= hi); print(f"  자세 제외 {int((~keep).sum())}개: {cname} 특이점 근처"); Q, T = Q[keep], T[keep]
    return Q, T, names


def moment_arms(model_path, muscles, names, Q, n_design):
    """R (N, n_design, m) 미터, Fmax (m,). OpenSim computeMomentArm."""
    import opensim as osim
    osim.Logger.setLevelString("Error")
    model = osim.Model(model_path); state = model.initSystem(); cs = model.getCoordinateSet(); ms = model.getMuscles()
    coords = [cs.get(n) for n in names]; mus = [ms.get(n) for n in muscles]; Fmax = np.array([mu.getMaxIsometricForce() for mu in mus])
    R = np.zeros((len(Q), n_design, len(mus))); tic = time.time()
    for t, q in enumerate(Q):
        for c, v in zip(coords, q): c.setValue(state, float(v), False)
        model.realizePosition(state)
        for i in range(n_design):
            for j, mu in enumerate(mus): R[t, i, j] = mu.computeMomentArm(state, coords[i])
        if t % 1000 == 0: print(f"  모멘트암 {t}/{len(Q)}  {time.time()-tic:.0f}s", flush=True)
    return R, Fmax


def rot(axis, a):
    c, s = np.cos(a), np.sin(a)
    return {"x": np.array([[1, 0, 0], [0, c, -s], [0, s, c]]), "y": np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]]), "z": np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])}[axis]


def to_mobl(flex, add, rotation):
    """우리 어깨 좌표 Rz(flex)·Rx(add)·Ry(rot) → MoBL (elv_angle, shoulder_elv, shoulder_rot). 두 모델 모두 0 = 팔을 내린 자세, x 앞·y 위·z 오른쪽.
    MoBL 위팔 회전 = Rot(a, e)·Rot(b, el)·Rot(a, −e)·Rot(a, r) 이다: shoulder1_r2 가 elv_angle 을 되돌리는 연동 좌표라 shoulder_rot 은 평범한 세 번째 오일러각이 아니다.
    a·b 는 OpenSim 에서 잰 실제 축 (이름값 y·−x 에서 2~3° 기울어짐). a·b 가 수직이 되게 기저 B 를 잡아 분해한다."""
    Rt = MOBL_B.T @ (rot("z", flex) @ rot("x", add) @ rot("y", rotation)) @ MOBL_B
    el = np.arccos(np.clip(Rt[1, 1], -1, 1)); s = np.sin(el)
    if s < 1e-6: return 0.0, el, np.arctan2(Rt[0, 2], Rt[0, 0])
    e = np.arctan2(-Rt[0, 1], -Rt[2, 1]); phi = np.arctan2(-Rt[1, 0], Rt[1, 2])
    return e, el, (phi + e + np.pi) % (2 * np.pi) - np.pi


def moment_arms_mobl(model_path, muscles, Q, delta=0.01):
    """R (N, 4, m): 우리 좌표(arm_flex, arm_add, arm_rot, elbow_flex)에 대한 모멘트암. 근육 길이 중앙차분. 어깨뼈 연동은 assemble 이 풀어 준다."""
    import opensim as osim
    osim.Logger.setLevelString("Error")
    model = osim.Model(model_path); state = model.initSystem(); cs = model.getCoordinateSet(); ms = model.getMuscles()
    cE, cS, cR, cB = (cs.get(n) for n in ("elv_angle", "shoulder_elv", "shoulder_rot", "elbow_flexion")); mus = [ms.get(n) for n in muscles]; Fmax = np.array([mu.getMaxIsometricForce() for mu in mus])
    for co in (cE, cS, cR, cB):
        co.setDefaultClamped(False)
    model.getCoordinateSet().get("pro_sup").setValue(state, 0.0, False)
    def lengths(fl, ad, ro, eb):
        e, el, r = to_mobl(fl, ad, ro)
        cE.setValue(state, e, False); cS.setValue(state, el, False); cR.setValue(state, r, False); cB.setValue(state, eb, False)
        model.assemble(state); model.realizePosition(state)
        return np.array([mu.getGeometryPath().getLength(state) for mu in mus])
    R = np.zeros((len(Q), 4, len(mus))); tic = time.time()
    for t, q in enumerate(Q):
        for i in range(4):
            p, m = q[:4].copy(), q[:4].copy(); p[i] += delta; m[i] -= delta
            R[t, i] = -(lengths(*p) - lengths(*m)) / (2 * delta)
        if t % 500 == 0: print(f"  모멘트암 {t}/{len(Q)}  {time.time()-tic:.0f}s", flush=True)
    return R, Fmax


def static_opt(R, Fmax, T):
    """프레임별 min Σa² + W·|R F a − τ|², 0 ≤ a ≤ 1. 반환 a (N, m), 잔차 (N, n)."""
    m = len(Fmax); a = np.zeros((len(T), m)); res = np.zeros_like(T); s = np.sqrt(W_TAU)
    for t in range(len(T)):
        M = R[t] * Fmax
        sol = lsq_linear(np.vstack([s * M, np.eye(m)]), np.r_[s * T[t], np.zeros(m)], bounds=(0.0, 1.0))
        a[t] = sol.x; res[t] = M @ sol.x - T[t]
    return a, res


def allocation(C, W, F, R):
    """spanning_check.torque_signatures 와 같은 귀속. 반환 V (k, n) 단위벡터, 시너지별 평균 토크 크기."""
    k = C.shape[1]; share = C[:, :, None] * W[None, :, :]; share /= np.maximum(share.sum(axis=1, keepdims=True), 1e-12)      # (T,k,m)
    V = np.zeros((k, R.shape[1]))
    for i in range(k):
        Ti = np.einsum("tm,tnm->tn", share[:, i, :] * F, R)                                                                 # (T,n) 시너지 i 몫의 토크
        V[i] = (Ti * C[:, i:i + 1]).sum(axis=0) / max(C[:, i].sum(), 1e-9)
    return dd.unit(V), np.linalg.norm(V, axis=1)


def report_k(label, space, k, V, tau_eval, tau_all, writer, vaf_total=None, vaf_min=None, save=True):
    ev = st.evaluate(V, tau_eval, tau_all)
    print(f"  [{label}] k={k}: 토크 용량 합 {ev['cap']:.1f}  LP목적 {ev['total']:.1f}  slack {ev['slack']:.2f}  못 내는 토크 {ev['worst']:.2f} N·m", flush=True)
    writer.writerow([label, k, "" if vaf_total is None else f"{vaf_total:.4f}", "" if vaf_min is None else f"{vaf_min:.4f}", f"{ev['cap']:.1f}", f"{ev['total']:.1f}", f"{ev['slack']:.3f}", f"{ev['worst']:.3f}"])
    if save and ev["c"] is not None:
        np.savez(os.path.join(ROOT, f"report/newmus_V_design_{space}_k{k}_40kg_h160_musyn.npz"), V=V, c=ev["c"], V_nm=V * ev["c"][:, None], coords=np.array(ds.get(space)["base"]), space=space)
    return ev


def main():
    global TASKS
    ap = argparse.ArgumentParser()
    ap.add_argument("--space", required=True, choices=list(CONFIG), help="trunk, arm_mobl, arm_10 (쓰지 않음)"); ap.add_argument("--tasks", default="report/upper/tasks_40kg_h160.csv")
    ap.add_argument("--ref-mass", type=float, default=40.0); ap.add_argument("--ref-height", type=float, default=1.60)
    ap.add_argument("--n-frames", type=int, default=4000); ap.add_argument("--kmax", type=int, default=None)
    ap.add_argument("--check-npz", default=None, help="이 공간의 직접 최적화 npz (여러 개면 쉼표)"); ap.add_argument("--also", default=None, help="배분비를 이 좌표만 잘라 따로 평가 (arm_r → shoulder_r)")
    ap.add_argument("--check-also", default=None, help="--also 공간의 직접 최적화 npz"); ap.add_argument("--restarts", type=int, default=20)
    a = ap.parse_args(); TASKS = os.path.join(ROOT, a.tasks); dd.REF = a.ref_mass * a.ref_height
    cfg = CONFIG[a.space]; key = a.space; a.space = cfg["space"]; SP = ds.get(a.space); n = len(SP["base"]); model_path = os.path.join(ROOT, cfg["model"])
    cache = os.path.join(ROOT, f"out/synergy_upper/{key}_poses.npz"); os.makedirs(os.path.dirname(cache), exist_ok=True)
    if os.path.exists(cache):
        d = np.load(cache); Q, T, R, Fmax = d["Q"], d["T"], d["R"], d["Fmax"]; print(f"저장된 자세·모멘트암 사용 {len(Q)}개", flush=True)
    else:
        Q, T, names = sample_poses(a.space, a.n_frames); print(f"[{a.space}] 자세 {len(Q)}개, 근육 {len(cfg['muscles'])}개, 좌표 {names}", flush=True)
        R, Fmax = moment_arms_mobl(model_path, cfg["muscles"], Q) if cfg.get("mobl") else moment_arms(model_path, cfg["muscles"], names, Q, n); np.savez(cache, Q=Q, T=T, R=R, Fmax=Fmax, names=np.array(names))
    print("좌표별 토크 피크 N·m:", np.round(np.abs(T).max(axis=0), 1), flush=True)
    tic = time.time(); act, res = static_opt(R, Fmax, T); F = act * Fmax
    print(f"정적 최적화 {time.time()-tic:.0f}s. 좌표별 못 낸 토크 최대 N·m {np.round(np.abs(res).max(axis=0),2)}  (피크 대비 95퍼센타일 % {np.round(np.percentile(np.abs(res),95,axis=0)/np.abs(T).max(axis=0)*100,2)})", flush=True)
    peak = act.max(axis=0); print("근육별 활성 피크:", dict(zip(cfg["muscles"], np.round(peak, 3))), flush=True)
    live = peak > 1e-3; A = act[:, live] / peak[live]; muscles = [m for m, l in zip(cfg["muscles"], live) if l]; Rl, Fl = R[:, :, live], F[:, live]
    print(f"NMF 근육 {len(muscles)}개: {muscles}", flush=True)

    tau_all, tau_eval = load_tau(a.space)
    also = None
    if a.also:
        also_tau = load_tau(a.also); dd.SP = ds.get(a.space); dd.COORDS = dd.SP["base"]; also = (a.also, len(ds.get(a.also)["base"]), *also_tau)
    out = os.path.join(ROOT, f"report/synergy_muscle_{key}_40kg_h160.csv"); fh = open(out, "w", newline="", encoding="utf-8"); w = csv.writer(fh)
    w.writerow(["design", "k", "vaf_total", "vaf_muscle_min", "capacity_sum_Nm", "lp_objective", "slack_max_Nm", "shortfall_max_Nm"])
    print("\n비교 기준:", flush=True)
    report_k("축 정렬", a.space, 2 * n, dd.axis_aligned(n), tau_eval, tau_all, w, save=False)
    for p in (a.check_npz.split(",") if a.check_npz else []):
        Vd = np.load(os.path.join(ROOT, p))["V"]; report_k("직접 최적화", a.space, len(Vd), Vd, tau_eval, tau_all, w, save=False)
    if also:
        report_k(f"축 정렬 {also[0]}", also[0], 2 * also[1], dd.axis_aligned(also[1]), also[3], also[2], w, save=False)
        if a.check_also:
            Vd = np.load(os.path.join(ROOT, a.check_also))["V"]; report_k(f"직접 최적화 {also[0]}", also[0], len(Vd), Vd, also[3], also[2], w, save=False)
    print("\n근육 시너지:", flush=True); rows = []
    for k in range(1, (a.kmax or len(muscles)) + 1):
        C, W = sy.nmf(A, k, seed=k, n_restarts=a.restarts); total, per = sy.vaf(A, C @ W); rows.append((k, total, float(per.min())))
        print(f" k={k}: 전체 VAF {total*100:.1f} %  최저 근육 VAF {per.min()*100:.1f} %", flush=True)
        V, size = allocation(C, W, Fl, Rl); print("   배분비:", " | ".join(f"{np.round(v,2).tolist()}" for v in V), flush=True)
        report_k("근육 시너지", a.space, k, V, tau_eval, tau_all, w, total, float(per.min()))
        if also: report_k(f"근육 시너지 {also[0]}", also[0], k, dd.unit(V[:, :also[1]]), also[3], also[2], w, total, float(per.min()))
        fh.flush()
    fh.close()
    gk = next((r[0] for r in rows if r[1] >= sy.VAF_TOTAL), None); dk = next((r[0] for r in rows if r[1] >= sy.VAF_TOTAL and r[2] >= sy.VAF_MUSCLE), None)
    print(f"\nVAF 로 고른 개수: global {gk}, dual {dk}.  최소 개수 n+1 = {n+1}\n→ {out}", flush=True)


if __name__ == "__main__":
    main()
