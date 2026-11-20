"""`video_pose.py` 의 결과(카메라 좌표의 몸 표면 점) → 세계 좌표 마커 → 스케일 → IK → 우리 형식 관절각. **`osim` 환경.**

    python tools/video_ik.py --pose data/custom/<묶음>/video/<이름>_pose.npz --height 1.70 --mass 60 --out data/custom/<묶음>/video/<이름>
→ <out>_markers.trc, <out>_scaled.osim (LaiUhlrich2022 를 그 사람 크기로), <out>_ik.mot (라디안, 우리 좌표 이름), <out>_video.json (진단)

보통은 직접 부르지 않는다 — `tools/add_motion.py --video ...` 가 부른다.

단계
  1. 크기: SMPL 관절 사슬 길이(발목→무릎→엉덩이→척추→목→머리)로 키를 추정하고, 입력한 키(--height)에 맞게 전체를 늘이거나 줄인다.
     (카메라 한 대는 절대 크기를 모른다. 키를 안 주면 NLF 가 추정한 크기를 그대로 쓴다.)
  2. 위쪽 방향: 카메라가 옆으로 안 기울었다고 보고(삼각대), 앞뒤 기울기는 "발 위에 무게중심이 평균적으로 있다"로 정한다.
     여러 걸음을 걸으면 발 디딘 자리들이 한 평면(바닥)에 있다는 것도 쓴다.
  3. 바닥 높이 = 발 표면 점의 가장 낮은 높이(5퍼센타일). 앞쪽(x) = 골반이 향한 쪽의 중앙값.
  4. 깊이 평활: 카메라 시선 방향의 몸 전체 이동만 1 Hz 로 거른다 (카메라 한 대는 거리가 프레임마다 2~4 cm 출렁인다 → 지면반력이 몇 배로 부푼다).
  5. 표면 점 중 `tools/smpl_markers.json` 의 것만 마커로 써서 .trc 저장.
  6. 스케일: 뼈마다 (보정 표본의 스케일) × (SMPL 뼈 길이 비). 질량은 --mass.
  7. IK (OpenSim InverseKinematicsTool) → 도 단위 .mot → 라디안·우리 좌표 이름으로.
  8. 발 기준 위치: 디딘 발이 미끄러지지 않게 골반 위치를 관절각에서 다시 만든다 (`root_from_contacts`). 트레드밀처럼 디딘 발이 움직이면 자동으로 꺼진다.

가정: 카메라 고정, 한 사람, 전신이 화면에 들어옴, 평평한 바닥.
"""
import argparse, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "tools"))
import addb_io as io                                   # noqa: E402
from trc_io import write_trc                            # noqa: E402
import smpl_calibrate as sc                             # noqa: E402

GENERIC = os.path.join(ROOT, "tools", "opencap_core", "LaiUhlrich2022.osim")
MARKERS = os.path.join(ROOT, "tools", "smpl_markers.json")
FOOT_JOINTS = {"r": (8, 11), "l": (7, 10)}             # SMPL 발목, 발(앞꿈치) 관절
#: 뼈 → 스케일에 쓸 SMPL 길이 (축마다). None 은 표본 스케일 그대로.
SCALE_BY = {"pelvis": ("hip_width",) * 3, "torso": ("shoulder_width", "torso", "shoulder_width"),
            **{f"femur_{s}": (f"femur_{s}",) * 3 for s in "rl"}, **{f"tibia_{s}": (f"tibia_{s}",) * 3 for s in "rl"}, **{f"patella_{s}": (f"tibia_{s}",) * 3 for s in "rl"},
            **{f"{b}_{s}": (f"foot_{s}",) * 3 for s in "rl" for b in ("talus", "calcn", "toes")},
            **{f"humerus_{s}": (f"humerus_{s}",) * 3 for s in "rl"}, **{f"{b}_{s}": (f"forearm_{s}",) * 3 for s in "rl" for b in ("ulna", "radius", "hand")}}


def fill_gaps(X, max_gap=15):
    """NaN 프레임을 선형보간으로 메운다 (X: (T, ...)). 너무 긴 공백은 경고."""
    bad = np.isnan(X.reshape(len(X), -1)).any(1)
    if not bad.any(): return X, 0
    good = np.where(~bad)[0]; flat = X.reshape(len(X), -1).copy()
    for j in range(flat.shape[1]): flat[bad, j] = np.interp(np.where(bad)[0], good, flat[good, j])
    return flat.reshape(X.shape), int(bad.sum())


def lowpass(X, fs, fc):
    from scipy.signal import butter, filtfilt
    if len(X) < 30: return X
    b, a = butter(4, fc / (fs / 2.0)); return filtfilt(b, a, X, axis=0, padlen=min(len(X) - 1, int(fs)))


def world_frame(V, J, fs, depth_lowpass=1.0, log=print):
    """카메라 좌표 (x 오른쪽, y 아래, z 앞) → 세계 (x 앞, y 위, z 오른쪽), 바닥 y=0. 반환 (Vw, Jw, 진단 dict)."""
    com = V.mean(axis=1)                                         # 표면 점 평균 — 무게중심 대용
    feet = 0.5 * (J[:, 7] + J[:, 8]); a = (com - feet); up = a.mean(0); up /= np.linalg.norm(up)
    up_bal = up.copy(); cam_x = np.array([1.0, 0.0, 0.0]); roll = np.degrees(np.arcsin(np.clip(up @ cam_x, -1, 1)))
    up = up - (up @ cam_x) * cam_x; up /= np.linalg.norm(up)       # 옆 기울기 0 가정
    # 여러 걸음: 발 디딘 자리의 주된 방향 d 에 수직이 되게
    ank = np.concatenate([J[:, 7], J[:, 8]]); span = ank - ank.mean(0); u, s, vt = np.linalg.svd(span - np.outer(span @ up, up), full_matrices=False)
    travel = float(np.ptp((J[:, 0] - J[0, 0]) @ vt[0])); method = "균형"
    if travel > 1.0:
        Jl = lowpass(J, fs, 4.0); sp = {k: np.linalg.norm(np.gradient(Jl[:, j], 1.0 / fs, axis=0), axis=1) for k, j in (("r", 8), ("l", 7))}
        pts = np.concatenate([Jl[sp["r"] < 0.4, 8], Jl[sp["l"] < 0.4, 7]])
        if len(pts) > 20:
            c = pts - pts.mean(0); _, s2, vt2 = np.linalg.svd(c, full_matrices=False); d = vt2[0]
            if s2[0] / np.sqrt(len(pts)) > 0.25: up = up - (up @ d) * d; up /= np.linalg.norm(up); method = "균형 + 걸은 방향"
    pitch = np.degrees(np.arccos(np.clip(up @ np.array([0.0, -1.0, 0.0]), -1, 1)))
    right = (J[:, 2] - J[:, 1]); right = right - np.outer(right @ up, up); fwd = np.cross(up, right); fwd = np.median(fwd / np.linalg.norm(fwd, axis=1, keepdims=True), axis=0)
    fwd = fwd - (fwd @ up) * up; fwd /= np.linalg.norm(fwd); z = np.cross(fwd, up)
    R = np.stack([fwd, up, z])                                     # 세계 = R · 카메라
    Vw = V @ R.T; Jw = J @ R.T
    # 발바닥 높이: 발목·앞꿈치 관절 아래의 표면 점 대신, 전체 표면 점 중 가장 낮은 것
    low = Vw[:, :, 1].min(axis=1); y0 = float(np.percentile(low, 5)); origin = np.array([Jw[0, 0, 0], y0, Jw[0, 0, 2]])
    Vw = Vw - origin; Jw = Jw - origin
    diag = dict(up_method=method, camera_pitch_deg=round(float(pitch), 1), balance_roll_deg=round(float(roll), 1), travel_m=round(travel, 2))
    # 깊이 방향 평활: 카메라 한 대는 몸 전체의 거리(깊이)가 프레임마다 2~4 cm 출렁인다. 그 출렁임은 무게중심 가속도 → 지면반력을 몇 배로 부풀린다
    # (2026-11-19 검증: 하체 토크가 OpenCap 의 3~5배). 자세는 두고, 시선 방향의 몸 전체 이동만 DEPTH_LOWPASS 로 거른다.
    dz = R @ np.array([0.0, 0.0, 1.0]); root = Jw[:, 0] @ dz
    if depth_lowpass and len(root) > 30:
        d_s = lowpass(root[:, None], fs, depth_lowpass)[:, 0]; Vw = Vw + ((d_s - root)[:, None, None] * dz); Jw = Jw + ((d_s - root)[:, None, None] * dz)
        diag["depth_jitter_cm"] = round(float(np.std(root - d_s)) * 100, 1)
    diag["depth_dir_world"] = [round(float(x), 4) for x in dz]
    log(f"  세계 좌표: 위쪽 = {method}, 카메라 앞뒤 기울기 {pitch:.1f}°, 이동 거리 {travel:.2f} m, 깊이 출렁임 {diag.get('depth_jitter_cm', 0):.1f} cm → {depth_lowpass:g} Hz 로 거름")
    return Vw, Jw, diag


def root_from_contacts(ik_path, model_path, dz, fs, mode="auto", log=print):
    """IK 뒤: 디딘 발이 미끄러지지 않게 골반 위치(pelvis_tx/ty/tz)를 관절각으로 다시 만든다. 관절각은 그대로.
    발마다 발바닥 세 점(뒤꿈치·중족골두·발끝) 중 가장 낮은 점을 보고, (높이가 바닥 + 4 cm 안) 이고 (깊이 성분을 뺀 속도 < 0.35 m/s) 이면 디딤 시작,
    (높이 > 바닥 + 8 cm) 또는 (속도 > 0.8 m/s) 면 디딤 끝 (히스테리시스). 디딘 동안 그 점을 처음 자리에 묶는다. 묶는 점이 바뀌면(뒤꿈치 → 발끝)
    새 점을 지금 자리에 묶어 이어지게 한다. 둘 다 뜬 동안은 보정량을 그대로 둔다. 디딘 프레임이 절반이 안 되면(트레드밀) 하지 않는다."""
    import grf_estimate as ge
    mot = io.read_mot(ik_path); names = list(mot.names); Q = mot.data.copy(); T = len(Q)
    body = ge.Body(model_path); zero = np.zeros(len(names)); P = {f: np.zeros((T, 3, 3)) for f in body.feet}
    for i in range(T):
        body.set(names, Q[i], zero); sp = body.sole_points()
        for f in body.feet: P[f][i] = sp[f]
    floor = float(np.percentile(np.min([P[f][:, :, 1].min(1) for f in body.feet], axis=0), 5))
    dz = np.asarray(dz, float); state = {}; contact = {f: np.zeros(T, bool) for f in body.feet}
    for f in body.feet:
        k = P[f][:, :, 1].argmin(1); low = P[f][np.arange(T), k]; Ls = lowpass(P[f].reshape(T, -1), fs, 6.0).reshape(T, 3, 3)
        v = np.gradient(Ls, 1.0 / fs, axis=0); v = v - np.einsum("tpi,i->tp", v, dz)[..., None] * dz; spd = np.linalg.norm(v, axis=2)[np.arange(T), k]
        h = low[:, 1] - floor; on = False
        for i in range(T):
            if not on and h[i] < 0.04 and spd[i] < 0.35: on = True
            elif on and (h[i] > 0.08 or spd[i] > 0.8): on = False
            contact[f][i] = on
    frac = float(np.mean(np.any([contact[f] for f in body.feet], axis=0)))
    if mode == "off" or (mode == "auto" and frac < 0.5):
        log(f"  발 기준 위치: 안 함 ({'끔' if mode == 'off' else f'디딘 프레임 {frac * 100:.0f} % — 트레드밀이거나 발이 안 보인다'})"); return dict(root_lock="off", contact_pct=round(frac * 100, 1))
    corr = np.zeros((T, 3)); cur = np.zeros(3); anchor = {}
    for i in range(T):
        est = []
        for f in body.feet:
            if not contact[f][i]: anchor.pop(f, None); continue
            k = int(P[f][i, :, 1].argmin()); pos = P[f][i, k] + cur
            if f not in anchor or anchor[f][0] != k: anchor[f] = (k, pos.copy())
            est.append(anchor[f][1] - P[f][i, k])
        if est: cur = np.mean(est, axis=0)
        corr[i] = cur
    corr = lowpass(corr, fs, 6.0)                                   # 묶는 점이 바뀔 때의 작은 꺾임을 편다
    for j, n in enumerate(("pelvis_tx", "pelvis_ty", "pelvis_tz")): Q[:, names.index(n)] += corr[:, j]
    io.write_mot(ik_path, names, mot.time, Q, in_degrees=False)
    log(f"  발 기준 위치: 디딘 프레임 {frac * 100:.0f} %, 보정 범위 {np.ptp(corr, axis=0).max() * 100:.1f} cm")
    return dict(root_lock="on", contact_pct=round(frac * 100, 1), corr_range_cm=round(float(np.ptp(corr, axis=0).max()) * 100, 1))


def scale_model(J, mass, out_model, cal, log=print):
    """일반 LaiUhlrich2022 + 보정 마커 → 그 사람 크기 모델."""
    import opensim as osim
    osim.Logger.setLevelString("Error")
    L = sc.smpl_lengths(J); L0 = cal["sample"]["lengths"]; m = osim.Model(GENERIC); st = m.initSystem(); ss = osim.ScaleSet(); fac = {}
    for i in range(m.getBodySet().getSize()):
        bn = m.getBodySet().get(i).getName(); f = np.array(cal["sample"]["scale"].get(bn, [1, 1, 1]), float)
        if bn in SCALE_BY: f = f * np.array([L[k] / L0[k] for k in SCALE_BY[bn]])
        fac[bn] = f; s = osim.Scale(); s.setSegmentName(bn); s.setScaleFactors(osim.Vec3(*[float(x) for x in f])); s.setApply(True); ss.cloneAndAppend(s)
    for mk in cal["markers"]:                                      # 스케일 전에 일반 좌표로 붙이면 scale() 이 같이 옮긴다
        mm = osim.Marker(); mm.setName(mk["name"]); mm.setParentFrameName(f"/bodyset/{mk['body']}"); mm.set_location(osim.Vec3(*mk["loc_generic"])); m.addMarker(mm)
    m.finalizeConnections(); st = m.initSystem()
    assert m.scale(st, ss, True, float(mass)), "scale 실패"
    st = m.initSystem(); m.setName(os.path.splitext(os.path.basename(out_model))[0]); m.printToXML(out_model)
    log(f"  스케일: 대퇴 ×{fac['femur_r'][0]:.3f}, 정강이 ×{fac['tibia_r'][0]:.3f}, 몸통 높이 ×{fac['torso'][1]:.3f}, 위팔 ×{fac['humerus_r'][0]:.3f}, 질량 {m.getTotalMass(st):.1f} kg")
    return out_model


def run_ik(model_path, trc, out_mot, cal, log=print):
    import opensim as osim
    ik = osim.InverseKinematicsTool(); ik.setName(os.path.basename(out_mot)[:-4]); model = osim.Model(model_path); model.initSystem(); ik.setModel(model)
    ik.setMarkerDataFileName(trc); ik.setOutputMotionFileName(out_mot); md = osim.MarkerData(trc); ik.setStartTime(md.getStartFrameTime()); ik.setEndTime(md.getLastFrameTime())
    tasks = osim.IKTaskSet()
    for mk in cal["markers"]:
        tk = osim.IKMarkerTask(); tk.setName(mk["name"]); tk.setWeight(float(mk["weight"])); tk.setApply(True); tasks.cloneAndAppend(tk)
    ik.set_IKTaskSet(tasks); ik.set_report_errors(True); ik.setResultsDir(os.path.dirname(out_mot)); ik.run()
    err = os.path.join(os.path.dirname(out_mot), f"{ik.getName()}_ik_marker_errors.sto"); rms = mx = float("nan")
    if os.path.exists(err): e = io.read_mot(err); rms = float(np.mean(e.column("marker_error_RMS"))); mx = float(np.max(e.column("marker_error_max")))
    log(f"  IK: 마커 오차 RMS {rms * 100:.1f} cm, 최대 {mx * 100:.1f} cm")
    return rms, mx


def run(pose_npz, out_prefix, height=None, mass=70.0, footlock="auto", marker_lowpass=8.0, depth_lowpass=1.0, log=print):
    """반환 dict(ik=우리 형식 mot 경로, model, height_m, diag)."""
    if not os.path.exists(MARKERS): raise SystemExit(f"{MARKERS} 가 없다 — tools/smpl_calibrate.py 를 먼저 돌린다")
    cal = json.load(open(MARKERS, encoding="utf-8")); d = np.load(pose_npz); fs = float(d["fps"]); t = d["time"].astype(float)
    V, nbad = fill_gaps(d["verts"].astype(float)); J, _ = fill_gaps(d["joints"].astype(float))
    est_h = sc.chain_length(J) * cal["stature_per_chain"]; s = (height / est_h) if height else 1.0
    log(f"  키: 영상 추정 {est_h:.3f} m" + (f" → 입력 {height:.2f} m 에 맞춰 ×{s:.3f}" if height else " (입력 없음 — 추정값 그대로)") + (f", 사람 못 잡은 {nbad} 프레임은 보간" if nbad else ""))
    V, J = V * s, J * s
    Vw, Jw, diag = world_frame(V, J, fs, depth_lowpass, log); diag.update(height_est_m=round(est_h, 3), scale_applied=round(s, 3), frames=len(t), fps=fs, missing_frames=nbad)
    names = [m["name"] for m in cal["markers"]]; X = Vw[:, [m["vert"] for m in cal["markers"]]]
    X = lowpass(X, fs, marker_lowpass)                             # 프레임별 잡음을 IK 전에 한 번 거른다 (최종 컷오프는 add_motion 의 --lowpass)
    os.makedirs(os.path.dirname(os.path.abspath(out_prefix)), exist_ok=True)
    trc = write_trc(out_prefix + "_markers.trc", names, t - t[0], X, rate=fs)
    model = scale_model(Jw, mass, out_prefix + "_scaled.osim", cal, log); raw = out_prefix + "_ik_deg.mot"
    rms, mx = run_ik(model, trc, raw, cal, log); diag.update(ik_rms_cm=round(rms * 100, 2), ik_max_cm=round(mx * 100, 2))
    from opencap_import import convert_mot
    convert_mot(raw, out_prefix + "_ik.mot")
    diag.update(root_from_contacts(out_prefix + "_ik.mot", model, diag["depth_dir_world"], fs, footlock, log))
    if rms > 0.04: log("  [주의] IK 마커 오차가 4 cm 를 넘는다 — 사람이 가려졌거나 화면 밖으로 나갔는지 영상을 확인")
    json.dump(diag, open(out_prefix + "_video.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    return dict(ik=out_prefix + "_ik.mot", model=model, height_m=height or est_h, diag=diag)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--pose", required=True); ap.add_argument("--out", required=True, help="출력 접두어")
    ap.add_argument("--height", type=float, default=None); ap.add_argument("--mass", type=float, default=70.0); ap.add_argument("--footlock", default="auto", choices=("auto", "on", "off"))
    a = ap.parse_args(); run(a.pose, a.out, a.height, a.mass, a.footlock)


if __name__ == "__main__":
    main()
