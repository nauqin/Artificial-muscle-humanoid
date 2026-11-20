"""(한 번만 돌리는 보정) SMPL 몸 표면 점을 OpenSim 뼈에 붙인다 → `tools/smpl_markers.json`. **`osim` 환경.**

왜: 영상 한 개(`video_pose.py`, NLF)는 프레임마다 SMPL 표면 6,890점을 준다. OpenSim IK 는 "뼈에 붙은 마커"가 필요하다.
OpenCap 공개 세션(4d5c3eb1…, 키 1.86 m 남성)은 같은 순간의 **영상**과 OpenCap 이 계산한 **관절각·스케일 모델**을 같이 준다.
그래서 그 영상을 NLF 로 돌린 뒤, 표면 점마다 "어느 뼈 좌표계에서 가장 안 움직이나"를 재서 그 뼈의 마커로 삼는다.

    python tools/smpl_calibrate.py --session data/opencap/opencap_sample --trials squat walk_1_25ms DJ
        (먼저: tools/opencap_import.py <세션> --videos, 그리고 영상마다 tools/video_pose.py → <세션>/pose/<trial>_cam<i>.npz)

방법
  1. trial·카메라마다 NLF 관절 ↔ OpenCap 영상 키포인트(무릎·발목·어깨·팔꿈치·손목)로 카메라→세계 닮음 변환(회전·크기·이동)을 하나 구하고,
     프레임마다 이동만 따로 맞춘다(카메라 한 대는 깊이가 프레임마다 2~4 cm 흔들린다).
  2. OpenCap 관절각으로 스케일 모델을 세워 프레임마다 뼈 22개의 위치·방향을 구한다.
  3. 표면 점 × 뼈 마다, 그 뼈 좌표계에서 본 점 위치의 표준편차(강체 정도)를 잰다. 가장 작은 뼈가 그 점의 주인이다.
  4. 뼈마다 강체 정도가 좋은 점 중에서 서로 멀리 떨어진 것을 고른다(farthest point sampling).
  5. 마커 위치를 그 사람 모델의 뼈 스케일로 나눠 **일반(스케일 전) 모델 좌표**로 저장한다. 새 사람은 SMPL 뼈 길이 비로 스케일한다.

한계: 한 사람(남, 1.86 m, 83 kg)에서 만든 것이다. 체형이 많이 다르면(살집·옷) 표면 점이 뼈에서 떨어진 거리가 달라진다.
기준인 OpenCap 관절각도 영상 추정값이다(실험실 마커 대비 평균 4~5° 오차, Uhlrich 2023).
"""
import argparse, glob, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "tools"))
import addb_io as io                                   # noqa: E402
from trc_io import read_trc                             # noqa: E402

JOINT_PAIRS = {4: "LKnee", 5: "RKnee", 7: "LAnkle", 8: "RAnkle", 16: "LShoulder", 17: "RShoulder", 18: "LElbow", 19: "RElbow", 20: "LWrist", 21: "RWrist"}
#: 뼈마다 고를 마커 수. 발가락(toes)·손(hand)은 작고 흔들려 적게.
N_PER_BODY = {"pelvis": 8, "torso": 10, "femur_r": 6, "femur_l": 6, "tibia_r": 6, "tibia_l": 6, "calcn_r": 5, "calcn_l": 5, "toes_r": 2, "toes_l": 2,
              "humerus_r": 5, "humerus_l": 5, "ulna_r": 3, "ulna_l": 3, "radius_r": 4, "radius_l": 4, "hand_r": 2, "hand_l": 2}
#: 뼈마다 마커를 고를 수 있는 영역 (일반 모델 좌표, m). 관절 바로 옆 피부는 두 뼈 사이에서 늘어나므로 뺀다 — 처음엔 영역 없이 골랐더니
#: 골반 마커가 사타구니·엉덩이 점으로 뽑혀 앉았다 일어나기의 깊은 고관절 굴곡을 10~20° 작게 냈다 (2026-11-19).
#: 일반 LaiUhlrich2022: 고관절 중심 pelvis (−0.056, −0.078, ±0.077), 무릎 femur y −0.41, 어깨 torso (0.003, 0.372, ±0.170), 팔꿈치 humerus y −0.286
_ACR = np.array([0.003, 0.372, 0.170])
ZONES = {"pelvis": lambda p: p[1] > -0.03,
         "torso": lambda p: p[1] > 0.06 and min(np.linalg.norm(p - _ACR), np.linalg.norm(p - _ACR * [1, 1, -1])) > 0.09,
         "femur": lambda p: -0.36 < p[1] < -0.12, "tibia": lambda p: p[1] < -0.06,
         "humerus": lambda p: -0.24 < p[1] < -0.07, "ulna": lambda p: p[1] < -0.05, "radius": lambda p: p[1] < -0.05}
#: IK 가중치 (OpenCap Setup_IK 의 비율을 따름: 골반·다리·발 크게, 팔 작게)
WEIGHT = {"pelvis": 25, "torso": 8, "femur": 12, "tibia": 20, "calcn": 30, "toes": 10, "humerus": 5, "ulna": 5, "radius": 5, "hand": 2}
#: SMPL 관절로 재는 길이 (새 사람 스케일용). (관절 a, 관절 b)
LENGTHS = {"femur_r": (2, 5), "femur_l": (1, 4), "tibia_r": (5, 8), "tibia_l": (4, 7), "hip_width": (1, 2), "torso": (0, 12), "shoulder_width": (16, 17),
           "humerus_r": (17, 19), "humerus_l": (16, 18), "forearm_r": (19, 21), "forearm_l": (18, 20), "foot_r": (8, 11), "foot_l": (7, 10)}
#: 키 추정용 사슬: 발목→무릎→엉덩이, 골반→척추1→2→3→목→머리 (좌우 평균)
CHAIN = [((7, 4), (8, 5)), ((4, 1), (5, 2)), ((0, 3),), ((3, 6),), ((6, 9),), ((9, 12),), ((12, 15),)]


def similarity(A, B):
    ma, mb = A.mean(0), B.mean(0); A0, B0 = A - ma, B - mb
    U, S, Vt = np.linalg.svd(B0.T @ A0 / len(A)); d = np.sign(np.linalg.det(U @ Vt)); D = np.diag([1, 1, d])
    R = U @ D @ Vt; s = np.trace(np.diag(S) @ D) / (A0 ** 2).sum() * len(A); return s, R, mb - s * R @ ma


def smpl_lengths(J):
    """관절 (T, 24, 3) → 길이 dict (프레임 중앙값)."""
    return {k: float(np.nanmedian(np.linalg.norm(J[:, a] - J[:, b], axis=1))) for k, (a, b) in LENGTHS.items()}


def chain_length(J):
    return float(sum(np.mean([np.nanmedian(np.linalg.norm(J[:, a] - J[:, b], axis=1)) for a, b in seg]) for seg in CHAIN))


def body_frames(model_path, mot_path, times):
    """OpenCap 스케일 모델 + 관절각(도) → 프레임마다 뼈의 지면 변환. 반환 (뼈 이름, R (T, B, 3, 3), p (T, B, 3), 스케일 dict)."""
    import opensim as osim
    osim.Logger.setLevelString("Error")
    proc = osim.ModelProcessor(os.path.abspath(model_path)); proc.append(osim.ModOpRemoveMuscles()); model = proc.process(); state = model.initSystem()
    m = io.read_mot(mot_path); cs = model.getCoordinateSet(); bs = model.getBodySet(); names = [bs.get(i).getName() for i in range(bs.getSize())]
    R = np.zeros((len(times), len(names), 3, 3)); P = np.zeros((len(times), len(names), 3)); e = [osim.Vec3(1, 0, 0), osim.Vec3(0, 1, 0), osim.Vec3(0, 0, 1)]; o = osim.Vec3(0, 0, 0)
    for ti, t in enumerate(times):
        row = m.at(t)[0] if hasattr(m, "at") else np.array([np.interp(t, m.time, m.data[:, j]) for j in range(len(m.names))])
        for j, cn in enumerate(m.names):
            if not cs.contains(cn): continue
            c = cs.get(cn); v = row[j]
            if m.in_degrees and c.getMotionType() == 1: v = np.radians(v)
            c.setValue(state, float(v), False)
        model.realizePosition(state)
        for b, bn in enumerate(names):
            body = bs.get(bn); p0 = np.array([body.findStationLocationInGround(state, o).get(i) for i in range(3)]); P[ti, b] = p0
            for a in range(3): R[ti, b, :, a] = np.array([body.findStationLocationInGround(state, e[a]).get(i) for i in range(3)]) - p0
    scale = {}
    for bn in names:
        body = bs.get(bn)
        try: sf = body.get_attached_geometry(0).get_scale_factors(); scale[bn] = [sf.get(i) for i in range(3)]
        except Exception: scale[bn] = [1.0, 1.0, 1.0]
    return names, R, P, scale


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--session", default=os.path.join(ROOT, "data/opencap/opencap_sample")); ap.add_argument("--trials", nargs="+", default=["squat", "walk_1_25ms", "DJ"])
    ap.add_argument("--step", type=int, default=3, help="프레임 건너뛰기"); ap.add_argument("--out", default=os.path.join(ROOT, "tools", "smpl_markers.json"))
    a = ap.parse_args(); S = a.session
    model = glob.glob(os.path.join(S, "raw", "**", "Model", "*.osim"), recursive=True)[0]
    meta = json.load(open(os.path.join(S, "session.json"), encoding="utf-8"))
    Vw_all, R_all, P_all, J_all = [], [], [], []; bnames = scale = None
    for tr in a.trials:
        names, t, X = read_trc(os.path.join(S, "raw", "MarkerData", f"{tr}.trc")); idx = {n: i for i, n in enumerate(names)}
        mot = os.path.join(S, "raw", "OpenSimData", "Kinematics", f"{tr}.mot")
        for pz in sorted(glob.glob(os.path.join(S, "pose", f"{tr}_cam*.npz"))):
            d = np.load(pz); V, J = d["verts"].astype(np.float64), d["joints"].astype(np.float64); T = min(len(V), len(X)); sel = np.arange(0, T, a.step)
            sel = sel[~np.isnan(V[sel, 0, 0])]; V, J = V[sel], J[sel]; Xs = X[sel]
            A = np.concatenate([J[:, j] for j in JOINT_PAIRS]); B = np.concatenate([Xs[:, idx[k]] for k in JOINT_PAIRS.values()]); ok = ~np.isnan(B).any(1) & ~np.isnan(A).any(1)
            s, R, tr_ = similarity(A[ok], B[ok])
            shift = np.nanmean(np.stack([Xs[:, idx[k]] - (s * J[:, j] @ R.T + tr_) for j, k in JOINT_PAIRS.items()], axis=1), axis=1)
            W = (s * V @ R.T + tr_ + shift[:, None]).astype(np.float32)
            bnames, Rb, Pb, scale = body_frames(model, mot, t[sel])
            Vw_all.append(W); R_all.append(Rb); P_all.append(Pb); J_all.append(J * s)
            print(f"{tr} {os.path.basename(pz)}: {len(sel)} 프레임, 크기비 {s:.3f}, 깊이 흔들림 {np.nanstd(shift @ (R @ [0, 0, 1.0])) * 100:.1f} cm", flush=True)
    W = np.concatenate(Vw_all); Rb = np.concatenate(R_all); Pb = np.concatenate(P_all); J = np.concatenate(J_all); print(f"합계 {len(W)} 프레임")
    use = [b for b in bnames if b in N_PER_BODY]; std = np.full((len(use), W.shape[1]), np.inf); mean = np.zeros((len(use), W.shape[1], 3))
    for k, bn in enumerate(use):
        b = bnames.index(bn); loc = np.einsum("tvi,tij->tvj", W - Pb[:, None, b].astype(np.float32), Rb[:, b].astype(np.float32))      # Rᵀ (x − p)
        mean[k] = loc.mean(0); std[k] = np.sqrt(((loc - mean[k]) ** 2).sum(-1).mean(0))
    owner = std.argmin(0); best = std.min(0)
    markers = []
    for k, bn in enumerate(use):
        sf = np.array(scale[bn]); zone = ZONES.get(bn.rsplit("_", 1)[0] if bn[-2:] in ("_r", "_l") else bn)
        cand = np.where(owner == k)[0]
        if zone is not None:
            g = mean[k][cand] / sf; zc = np.array([zone(x * ([1, 1, -1] if bn.endswith("_l") else 1)) for x in g], bool)   # 왼쪽은 z 를 뒤집어 오른쪽 기준으로
            if zc.sum() >= N_PER_BODY[bn]: cand = cand[zc]
            else: print(f"[주의] {bn}: 허용 영역 안 점이 {int(zc.sum())}개뿐 — 영역 무시")
        if len(cand) < N_PER_BODY[bn]: print(f"[주의] {bn}: 주인인 점이 {len(cand)}개뿐");
        thr = max(np.percentile(best[cand], 40), 0.012) if len(cand) else 0; cand = cand[best[cand] <= thr]
        pick = [int(cand[np.argmin(best[cand])])]
        while len(pick) < min(N_PER_BODY[bn], len(cand)):            # 이미 고른 점들에서 가장 먼 점
            dmin = np.min(np.linalg.norm(mean[k][cand][:, None] - mean[k][pick][None], axis=2), axis=1); pick.append(int(cand[np.argmax(dmin)]))
        for v in pick:
            markers.append(dict(name=f"smpl_{v}", vert=v, body=bn, loc_generic=[round(float(x), 5) for x in mean[k][v] / sf], std_cm=round(float(best[v]) * 100, 2),
                                weight=WEIGHT[bn.rsplit("_", 1)[0] if bn[-2:] in ("_r", "_l") else bn]))
        print(f"{bn:10s} 주인 점 {int((owner == k).sum()):4d}개, 고른 {len(pick)}개 강체 정도 {np.mean([best[v] for v in pick]) * 100:.1f} cm (최대 {np.max([best[v] for v in pick]) * 100:.1f})")
    out = dict(_about="SMPL 표면 점 → LaiUhlrich2022 일반 모델 뼈 좌표. tools/smpl_calibrate.py (OpenCap 공개 세션 4d5c3eb1, NLF 0.3.2)",
               sample=dict(height_m=meta["height_m"], mass_kg=meta["mass_kg"], trials=a.trials, scale=scale, lengths=smpl_lengths(J), chain_m=chain_length(J)),
               stature_per_chain=meta["height_m"] / chain_length(J), markers=markers)
    json.dump(out, open(a.out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print(f"저장 {a.out}: 마커 {len(markers)}개, 키/사슬 길이 = {out['stature_per_chain']:.3f}")


if __name__ == "__main__":
    main()
