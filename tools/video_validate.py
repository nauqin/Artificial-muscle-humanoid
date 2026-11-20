"""영상 한 개 경로(`video_pose` → `video_ik`)가 OpenCap(아이폰 3대)과 얼마나 같은 관절각·로봇 토크를 내나. **`osim` 환경.**

    python tools/video_validate.py --trials STS jump run_2_5ms          # 보정(smpl_calibrate)에 안 쓴 trial 로
→ report/custom/video_validation.csv, 화면 표

trial × 카메라마다
  1. 그 카메라 영상 하나로 video_ik.run (키 1.86 m·질량 83 kg 은 세션 값을 준다 — 실제 쓸 때처럼)
  2. OpenCap 관절각(카메라 3대 + 마커 보강 + IK)과 관절각 RMSE (도). "모양" 은 평균 차이를 뺀 RMSE.
  3. 둘 다 add_motion 과 같은 과정(로봇 뼈대 → 지면반력 추정 → ID)으로 로봇 토크를 내서 부위별 피크 토크·RMSE 비교.
기준인 OpenCap 도 영상 추정이다 (실험실 마커 대비 관절각 평균 오차 4~5°, Uhlrich 2023). 즉 이 표는 "OpenCap 대비" 차이다.
"""
import argparse, csv, glob, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "tools")); sys.path.insert(0, os.path.join(ROOT, "report"))
import addb_io as io                                   # noqa: E402
import grf_estimate as ge                               # noqa: E402

ANGLES = {"고관절 굴곡": "hip_flexion", "고관절 내전": "hip_adduction", "고관절 회전": "hip_rotation", "무릎": "knee_angle", "발목": "ankle_angle",
          "요추 굴곡": "lumbar_extension", "요추 옆굽힘": "lumbar_bending", "요추 회전": "lumbar_rotation", "어깨 굴곡": "arm_flex", "어깨 내전": "arm_add", "팔꿈치": "elbow_flex"}
TORQUE_GROUPS = {"하체": ("hip_flexion", "hip_adduction", "knee_angle", "ankle_angle"), "몸통": ("lumbar_extension", "lumbar_bending", "lumbar_rotation"),
                 "어깨": ("arm_flex", "arm_add", "arm_rot"), "팔꿈치": ("elbow_flex",)}


def cols_for(mot, base):
    return [c for c in (base, f"{base}_r", f"{base}_l") if c in mot.names]


def angle_err(ref_path, test_path, lowpass=4.0):
    a = io.read_mot(ref_path); b = io.read_mot(test_path)
    na, ta, A, _ = ge.smooth_resample(a, 100.0, lowpass); nb, tb, B, _ = ge.smooth_resample(b, 100.0, lowpass)
    t0, t1 = max(ta[0], tb[0]) + 0.2, min(ta[-1], tb[-1]) - 0.2; t = np.arange(t0, t1, 0.01); out = {}
    for label, base in ANGLES.items():
        e, es = [], []
        for c in cols_for(a, base):
            if c not in nb: continue
            x = np.degrees(np.interp(t, ta, A[:, na.index(c)])); y = np.degrees(np.interp(t, tb, B[:, nb.index(c)])); d = y - x
            e.append(np.sqrt(np.mean(d ** 2))); es.append(np.sqrt(np.mean((d - d.mean()) ** 2)))
        if e: out[label] = (float(np.mean(e)), float(np.mean(es)))
    return out


def torque_err(id_ref, id_test):
    a = io.read_mot(id_ref); b = io.read_mot(id_test); t = a.time[(a.time >= b.t0) & (a.time <= b.t1)]; out = {}
    for g, bases in TORQUE_GROUPS.items():
        pr, pt, rm, qr, qt = [], [], [], [], []
        for base in bases:
            for c in (f"{base}_moment", f"{base}_r_moment", f"{base}_l_moment"):
                if c in a.names and c in b.names:
                    x = np.interp(t, a.time, a.column(c)); y = np.interp(t, b.time, b.column(c)); pr.append(np.abs(x).max()); pt.append(np.abs(y).max()); rm.append(np.sqrt(np.mean((y - x) ** 2)))
                    qr.append(np.percentile(np.abs(x), 95)); qt.append(np.percentile(np.abs(y), 95))
        if pr: out[g] = (float(max(pr)), float(max(pt)), float(max(rm)), float(max(qr)), float(max(qt)))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--session", default=os.path.join(ROOT, "data/opencap/opencap_sample"))
    ap.add_argument("--trials", nargs="+", default=["STS", "jump", "run_2_5ms"]); ap.add_argument("--no-torque", action="store_true")
    a = ap.parse_args()
    import opensim as osim; osim.Logger.setLevelString("Error")
    import video_ik, add_motion
    S = a.session; meta = json.load(open(os.path.join(S, "session.json"), encoding="utf-8")); H, M = meta["height_m"], meta["mass_kg"]
    cal = json.load(open(os.path.join(ROOT, "tools", "smpl_markers.json"), encoding="utf-8")); held = [t for t in a.trials if t in cal["sample"]["trials"]]
    if held: print(f"[주의] {held} 는 보정에 쓴 trial 이다 — 검증으로는 낙관적이다")
    rows = []; vdir = os.path.join(S, "video_check"); os.makedirs(vdir, exist_ok=True)
    for tr in a.trials:
        ref_ik = os.path.join(S, f"{tr}_ik.mot")
        if not a.no_torque:
            r = add_motion.process(f"oc_{tr}", ref_ik, H, "videocheck", task=tr, lowpass=add_motion.VIDEO_LOWPASS, log=lambda *x: None); ref_id = r["row"]["id"]
        for pz in sorted(glob.glob(os.path.join(S, "pose", f"{tr}_cam*.npz"))):
            cam = os.path.splitext(os.path.basename(pz))[0].split("_")[-1]; print(f"\n=== {tr} {cam} (영상 한 개) ===")
            v = video_ik.run(pz, os.path.join(vdir, f"{tr}_{cam}"), H, M)
            ae = angle_err(ref_ik, v["ik"]); row = dict(trial=tr, camera=cam, ik_rms_cm=v["diag"]["ik_rms_cm"], root_lock=v["diag"].get("root_lock", ""), contact_pct=v["diag"].get("contact_pct", ""), depth_jitter_cm=v["diag"].get("depth_jitter_cm", ""), height_est_m=v["diag"]["height_est_m"])
            print("  관절각 RMSE (도, 괄호 = 평균 차이 뺀 모양 RMSE): " + ", ".join(f"{k} {e:.1f} ({s:.1f})" for k, (e, s) in ae.items()))
            for k, (e, s) in ae.items(): row[f"ang_{k}"] = round(e, 1); row[f"shape_{k}"] = round(s, 1)
            if not a.no_torque:
                rv = add_motion.process(f"vid_{tr}_{cam}", v["ik"], H, "videocheck", task=tr, lowpass=add_motion.VIDEO_LOWPASS, log=lambda *x: None)
                te = torque_err(os.path.join(ROOT, ref_id), os.path.join(ROOT, rv["row"]["id"]))
                print("  로봇 토크 N·m (OpenCap → 영상: 95퍼센타일 / 최대, 시간 RMSE): " + ", ".join(f"{g} {p95:.0f}→{q95:.0f} / {p:.0f}→{q:.0f} ({e:.1f})" for g, (p, q, e, p95, q95) in te.items()))
                for g, (p, q, e, p95, q95) in te.items():
                    row[f"tq_{g}_p95_ref"] = round(p95, 1); row[f"tq_{g}_p95_video"] = round(q95, 1); row[f"tq_{g}_max_ref"] = round(p, 1); row[f"tq_{g}_max_video"] = round(q, 1); row[f"tq_{g}_rmse"] = round(e, 1)
            rows.append(row)
    out = os.path.join(ROOT, "report", "custom", "video_validation.csv"); keys = sorted({k for r in rows for k in r}, key=lambda k: list(rows[0]).index(k) if k in rows[0] else 999)
    with open(out, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    print(f"\n저장 {out}")
    for g in TORQUE_GROUPS:
        r95 = [r[f"tq_{g}_p95_video"] / r[f"tq_{g}_p95_ref"] for r in rows if r.get(f"tq_{g}_p95_ref")]
        if r95: print(f"  {g:4s} 로봇 토크 95퍼센타일 비 (영상/OpenCap) 중앙값 {np.median(r95):.2f}  범위 {min(r95):.2f}~{max(r95):.2f}  (trial {sorted({r['trial'] for r in rows})})")
    for k in ANGLES:
        v = [r[f"ang_{k}"] for r in rows if f"ang_{k}" in r]; s = [r[f"shape_{k}"] for r in rows if f"shape_{k}" in r]
        if v: print(f"  {k:8s} 관절각 RMSE 중앙값 {np.median(v):5.1f}° (모양 {np.median(s):4.1f}°)  범위 {min(v):.1f}~{max(v):.1f}")


if __name__ == "__main__":
    main()
