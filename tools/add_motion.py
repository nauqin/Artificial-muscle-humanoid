"""새 동작 하나(OpenCap 세션 · 영상 한 개 · 관절각 파일)를 로봇 뼈대의 관절 토크로 바꿔 설계 입력(tasks.csv)에 넣는다. **`osim` 환경.**

    # OpenCap 세션 (아이폰 2대 이상으로 찍은 것. 세션 번호·주소·zip·폴더 아무거나)
    python tools/add_motion.py --opencap 4d5c3eb1-1a59-4ea1-9178-d3634610561c --trials squat --task squat_deep --check
    # 영상 한 개 (폰 한 대, 삼각대 고정, 전신). 키·몸무게는 찍힌 사람 것
    python tools/add_motion.py --video ~/Videos/lift.mp4 --height 1.70 --mass 62 --task lift_box --hand-load 3 --bimanual --check
    # 이미 있는 관절각 (.mot, 우리 좌표 이름)
    python tools/add_motion.py --ik some_ik.mot --src-height 1.69 --name mytrial --task reach

한 trial 마다 하는 일
  1. 관절각을 100 Hz 로 맞추고 저역통과(영상 출처 기본 4 Hz, --ik 6 Hz)한 뒤 로봇 뼈대(40 kg · 1.60 m)로 옮긴다: 각도는 그대로, 골반 위치는 키 비로,
     로봇에서 잠긴 좌표(발목 아래·손목)는 0, 발바닥이 바닥(y=0)에 닿게 높이를 맞춘다.         → data/custom/<묶음>/<이름>_ik.mot
  2. 지면반력을 관절각에서 추정한다 (`grf_estimate.py`, 힘판 없음).                           → data/custom/<묶음>/<이름>_grf.mot
  3. 역동역학(ID) → 관절 토크. 손에 든 물건은 손 질량에 더한 모델로 (`--hand-load`).           → out/custom40_h160/<묶음>/<이름>/<이름>_id.sto, _id.png
  4. `report/custom/tasks_<묶음>.csv` 에 한 줄 (기존 `--tasks` 와 같은 열).
  5. `--check`: 지금 설계(최종 A안)가 이 동작을 낼 수 있는지 (`report/fit_check.py`).
  6. `--redesign lower,trunk`: 기존 과제 + 이 묶음으로 관절 배분비·토크 용량을 다시 설계 (`design_directions.py`), 경로 설계 명령을 알려 준다.

믿을 수 있는 정도 (자세한 표는 report/custom/README.md)
  - 상체(몸통·어깨·팔꿈치) 토크는 지면반력과 무관하다(관절 위쪽 분절과 손 하중만으로 정해진다). 관절각이 맞으면 맞는다.
  - 하체 토크는 추정 지면반력에 달렸다. 서서 하는 동작(스쿼트·들기)은 잘 맞고, 걷기는 두 발 지지 구간이 덜 맞는다.
  - 의자·난간·벽처럼 발 말고 다른 곳이 닿는 동작의 하체 토크는 틀리다 → `--upper-only`.
  - `--hand-load` 는 구간 내내 손에 들고 있다고 본다. 들기 전·내려놓은 뒤는 `--t0 --t1` 로 잘라낸다.
"""
import argparse, csv, json, os, re, subprocess, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "tools")); sys.path.insert(0, os.path.join(ROOT, "report"))
import _env  # noqa: F401,E402
import addb_io as io                                   # noqa: E402
import grf_estimate as ge                               # noqa: E402

ROBOT = dict(model="data/robot40_h160/robot40_h160.osim", mass=40.0, height=1.60, tag="40_h160")
TRANSL = ("pelvis_tx", "pelvis_ty", "pelvis_tz")
#: 영상 출처(OpenCap·카메라 한 대) 관절각의 기본 컷오프. 6 Hz 면 OpenCap 스쿼트에서 발별 수직력이 0.25~0.9 BW 로 출렁이고 골반 잔차 모멘트가
#: 28 N·m, 4 Hz 면 7 N·m (2026-11-19). 실험실 마커 데이터(--ik)는 지금까지처럼 6 Hz.
VIDEO_LOWPASS = 4.0
EDGE = 0.2                                              # 필터·미분이 망가지는 양 끝 (초) — ID 구간에서 뺀다
POSE_PY = os.environ.get("POSE_PYTHON") or os.path.join(os.path.expanduser("~"), "miniforge3", "envs", "pose", "python.exe")
BASE_TASKS = {"upper": "report/upper/tasks_40kg_h160.csv", "lower": "report/custom/tasks_lower_base_40kg_h160.csv"}
GROUPS = {"하체": ("hip_flexion", "hip_adduction", "knee_angle", "ankle_angle"), "몸통": ("lumbar_extension", "lumbar_bending", "lumbar_rotation"),
          "어깨": ("arm_flex", "arm_add", "arm_rot"), "팔꿈치": ("elbow_flex",)}


def safe(s):
    return re.sub(r"[^0-9A-Za-z_]+", "_", s).strip("_")


def robot_model_with_load(load_kg, bimanual, set_dir):
    """손 하중 모델. 이미 있는 것(0.52 양손, 0.9, 3.79)은 그대로 쓰고, 없으면 data/custom/<묶음>/models/ 에 만든다."""
    base = os.path.join(ROOT, ROBOT["model"])
    if not load_kg: return base
    name = f"{os.path.splitext(os.path.basename(base))[0]}_load{load_kg:g}{'_bi' if bimanual else ''}.osim"; p = os.path.join(os.path.dirname(base), name)
    if os.path.exists(p): return p
    from carrt_prepare import with_hand_load
    p = os.path.join(set_dir, "models", name); os.makedirs(os.path.dirname(p), exist_ok=True)
    if not os.path.exists(p): with_hand_load(base, p, load_kg, bimanual)
    return p


def retarget(src_ik, src_height, out_ik, *, lowpass=6.0, t0=None, t1=None, log=print):
    """사람 관절각 → 로봇 뼈대 관절각 (100 Hz, 필터, 라디안). 반환 (names, t, Q)."""
    mot = io.read_mot(src_ik)
    if t0 is not None or t1 is not None:
        keep = (mot.time >= (t0 if t0 is not None else -1e9)) & (mot.time <= (t1 if t1 is not None else 1e9)); mot.time = mot.time[keep]; mot.data = mot.data[keep]
    if len(mot.time) < 40: raise ValueError(f"구간이 너무 짧다 ({len(mot.time)} 프레임)")
    names, t, Q, Qd = ge.smooth_resample(mot, 100.0, lowpass)
    body = ge.Body(os.path.join(ROOT, ROBOT["model"])); cs = body.cs; ratio = ROBOT["height"] / src_height
    cn = [cs.get(i).getName() for i in range(cs.getSize()) if not cs.get(i).getName().endswith("_beta")]
    R = np.zeros((len(t), len(cn))); missing = []
    for j, n in enumerate(cn):
        if cs.get(n).get_locked(): continue                       # 로봇에서 잠긴 좌표는 0
        if n in names: R[:, j] = Q[:, names.index(n)] * (ratio if n in TRANSL else 1.0)
        else: missing.append(n)
    Rd = np.gradient(R, 0.01, axis=0); low = np.zeros(len(t))
    for i in range(len(t)):
        body.set(cn, R[i], Rd[i]); sp = body.sole_points(); low[i] = min(sp[f][:, 1].min() for f in sp)
    y0 = float(np.percentile(low, 5)); R[:, cn.index("pelvis_ty")] -= y0
    io.write_mot(out_ik, cn, t, R)
    log(f"  로봇 뼈대로: 키 비 ×{ratio:.3f}, {len(t)} 프레임 ({t[0]:.2f}~{t[-1]:.2f} s), 컷오프 {lowpass:g} Hz, 바닥 맞춤 {-y0 * 100:+.1f} cm" + (f", 없는 좌표 0 으로: {missing}" if missing else ""))
    return cn, t, R


def process(name, src_ik, src_height, set_name, *, task, hand_load=0.0, bimanual=False, lowpass=6.0, t0=None, t1=None, upper_only=False, dataset="Custom", note="", log=print):
    """trial 하나 끝까지. 반환 요약 dict (tasks.csv 한 줄 포함)."""
    import run_trial as rt, validate
    set_dir = os.path.join(ROOT, "data", "custom", set_name); out_dir = os.path.join(ROOT, "out", f"custom{ROBOT['tag']}", set_name, name)
    os.makedirs(set_dir, exist_ok=True); os.makedirs(out_dir, exist_ok=True)
    log(f"\n[{name}] 과제 {task}" + (f", 손 하중 {hand_load:g} kg{' (양손)' if bimanual else ''}" if hand_load else ""))
    ik = os.path.join(set_dir, f"{name}_ik.mot"); cn, t, R = retarget(src_ik, src_height, ik, lowpass=lowpass, t0=t0, t1=t1, log=log)
    model = robot_model_with_load(hand_load, bimanual, set_dir); a0, a1 = float(t[0] + EDGE), float(t[-1] - EDGE); summ = dict(name=name, task=task, t0=round(a0, 2), t1=round(a1, 2))
    if upper_only:
        s = rt.run_trial_no_grf(model, ik, out_dir, name, lowpass=-1, start=a0, end=a1, log=lambda *x: None); id_path = s["id_path"]; note = "[상체만] " + note
        log("  ID: 지면반력 없이 (상체 토크만 유효)")
    else:
        grf = os.path.join(set_dir, f"{name}_grf.mot"); r = ge.estimate(model, ik, grf, lowpass=0, floor=0.0, log=log); summ["grf"] = r["diag"]
        gm = io.read_mot(grf); plates = {k: p for k, p in io.find_plates(gm).items() if p.complete}
        ext = rt.write_external_loads(os.path.join(out_dir, f"{name}_external_loads.xml"), grf, plates, {k: {"body": k} for k in plates})
        id_path, _ = rt.run_inverse_dynamics(model, ik, ext, out_dir, name, a0, a1, lowpass=-1)
        idm = io.read_mot(id_path); W = r["diag"]["mass_kg"] * ge.G
        res = np.sqrt(sum(idm.column(f"pelvis_t{x}_force") ** 2 for x in "xyz")); resm = np.sqrt(sum(idm.column(f"pelvis_{x}_moment") ** 2 for x in ("tilt", "list", "rotation")))
        summ["residual_force_pctBW"] = round(float(np.percentile(res, 95) / W * 100), 1); summ["residual_moment_Nm"] = round(float(np.percentile(resm, 95)), 1)
        log(f"  ID: 골반 잔차 95퍼센타일 힘 {summ['residual_force_pctBW']} %BW, 모멘트 {summ['residual_moment_Nm']} N·m (추정 지면반력이라 0 에 가까워야 정상. 크면 공중·발 밖 구간)")
        try: validate.plot_id(os.path.join(out_dir, f"{name}_id.png"), idm, gm, plates, r["diag"]["mass_kg"], W, window=(a0, a1))
        except Exception as e: log(f"  [그림 생략] {type(e).__name__}: {e}")
    idm = io.read_mot(id_path); peaks = {}
    for gname, bases in GROUPS.items():
        if upper_only and gname == "하체": continue
        vals = []
        for b in bases:
            cols = [c for c in (f"{b}_moment", f"{b}_r_moment", f"{b}_l_moment") if c in idm.names]
            if cols: pk = max(float(np.abs(idm.column(c)).max()) for c in cols); peaks[b] = round(pk, 1); vals.append(f"{b.replace('_angle', '')} {pk:.1f}")
        log(f"  {gname} 피크 토크 (N·m): " + ", ".join(vals))
    summ["peaks_Nm"] = peaks
    row = dict(dataset=dataset, subject=set_name, trial=name, ik=os.path.relpath(ik, ROOT).replace("\\", "/"), id=os.path.relpath(id_path, ROOT).replace("\\", "/"),
               mass_kg=f"{ROBOT['mass']:g}", height_m=f"{ROBOT['height']:g}", task=task, hand_load_kg=f"{hand_load:g}",
               note=(note + f" 원본 키 {src_height:g} m → robot{ROBOT['tag']} 뼈대, 지면반력 {'없음' if upper_only else '추정'}, 컷오프 {lowpass:g} Hz" + ("_bi" if bimanual else "")).strip())
    json.dump(summ, open(os.path.join(out_dir, f"{name}_summary.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    summ["row"] = row; return summ


def update_tasks(set_name, rows):
    """report/custom/tasks_<묶음>.csv — 같은 trial 이름이면 바꿔 쓴다."""
    path = os.path.join(ROOT, "report", "custom", f"tasks_{set_name}.csv"); os.makedirs(os.path.dirname(path), exist_ok=True); old = []
    if os.path.exists(path): old = [r for r in csv.DictReader(open(path, encoding="utf-8")) if r["trial"] not in {x["trial"] for x in rows}]
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(old + rows)
    return path


def lower_base_tasks():
    """하체 재설계의 바탕 과제표 (이 컴퓨터에 있는 로봇 크기 ID 3건: 걷기·스쿼트·Moore 걷기)."""
    path = os.path.join(ROOT, BASE_TASKS["lower"])
    if os.path.exists(path): return path
    src = [("LabValidation", "subject3", "walking1_segment_0", "data/robot40_h160/lower/walking1_segment_0_ik.mot", "out/lower40_h160/subject3/walking1_segment_0_id.sto", "walk"),
           ("LabValidation", "subject3", "squats1_segment_0", "data/robot40_h160/lower/squats1_segment_0_ik.mot", "out/lower40_h160/subject3/squats1_segment_0_id.sto", "squat"),
           ("Moore2015", "subject3", "T047_segment_12", "data/robot40_h160/moore/T047_segment_12_ik.mot", "out/lower40_h160/moore/T047_segment_12_id.sto", "walk1.2")]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(["dataset", "subject", "trial", "ik", "id", "mass_kg", "height_m", "task", "hand_load_kg", "note"])
        for d, s, tr, ik, idp, task in src: w.writerow([d, s, tr, ik, idp, "40", "1.6", task, "0", "힘판 지면반력을 질량·키 비로 줄인 것 (rescale_body.py)"])
    return path


def redesign(set_name, spaces, tasks_path, maxfev=300, log=print):
    """기존 과제 + 새 묶음으로 관절 배분비·토크 용량을 다시 설계한다 (부위마다 수 분)."""
    import design_space as ds
    from fit_check import FINAL_A
    py = sys.executable; env = dict(os.environ, PYTHONUTF8="1", MKL_THREADING_LAYER="SEQUENTIAL"); new = list(csv.DictReader(open(tasks_path, encoding="utf-8")))
    for sp in spaces:
        base = lower_base_tasks() if sp == "lower" else os.path.join(ROOT, BASE_TASKS["upper"])
        rows = list(csv.DictReader(open(base, encoding="utf-8"))) + [r for r in new if not (sp == "lower" and r["note"].startswith("[상체만]"))]
        comb = os.path.join(ROOT, "report", "custom", f"tasks_{set_name}_plus_base_{'lower' if sp == 'lower' else 'upper'}.csv"); keys = list(rows[0])
        with open(comb, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore"); w.writeheader(); w.writerows(rows)
        k = int(np.load(os.path.join(ROOT, FINAL_A[sp]["design"]))["V"].shape[0]); tag = f"_40kg_h160_{set_name}"
        cmd = [py, os.path.join(ROOT, "report", "design_directions.py"), "--space", sp, "--tasks", os.path.relpath(comb, ROOT).replace(os.sep, "/"), "--k", str(k), "--ref-mass", "40", "--ref-height", "1.60", "--out-tag", tag, "--maxfev", str(maxfev)]
        log(f"\n=== 재설계 {sp} (근육 {k}개) ===\n  python report/design_directions.py " + " ".join(cmd[2:]))
        logp = os.path.join(ROOT, "report", "custom", f"design_{sp}_k{k}{tag}.log")
        with open(logp, "w", encoding="utf-8") as fh: rc = subprocess.run(cmd, cwd=ROOT, env=env, stdout=fh, stderr=subprocess.STDOUT).returncode
        tail = open(logp, encoding="utf-8", errors="replace").read().splitlines()
        for ln in [x for x in tail if ("총용량" in x or "전 프레임 검산" in x or "slack" in x)][-6:]: log("  " + ln.strip())
        if rc: log(f"  [실패] 로그 {logp}"); continue
        npz = f"report/newmus_V_design_{sp}_k{k}{tag}.npz"; old = np.load(os.path.join(ROOT, FINAL_A[sp]["design"]))["c"].sum(); newc = np.load(os.path.join(ROOT, npz))["c"].sum()
        log(f"  총 토크 용량: 최종 A안 {old:.1f} → 새 설계 {newc:.1f} N·m  ({(newc / old - 1) * 100:+.0f} %)   {npz}")
        stroke = FINAL_A[sp]["stroke"]; extra = "--ref 40 1.60 " if sp == "lower" else ""; comb_rel = os.path.relpath(comb, ROOT).replace(os.sep, "/")
        log(f"  다음(경로, 수 분~수십 분):\n    python report/path_design.py --space {sp} --tasks {comb_rel} --design {npz} --model {ROBOT['model']} --box-scale 0.9877 {extra}--out-tag {tag} --stroke-max {stroke} --refine-trust 0.01")


def main():
    ap = argparse.ArgumentParser(description="새 동작을 설계 입력으로")
    src = ap.add_argument_group("입력 (셋 중 하나)")
    src.add_argument("--opencap", help="OpenCap 세션 번호 / 주소 / zip / 폴더 / 이미 가져온 data/opencap/<이름>"); src.add_argument("--video", help="영상 파일 (카메라 한 대)")
    src.add_argument("--ik", help="관절각 .mot (우리 좌표 이름)")
    ap.add_argument("--trials", nargs="+", default=None, help="OpenCap: 이 trial 만"); ap.add_argument("--name", default=None, help="저장 이름 (기본: 세션·파일 이름)")
    ap.add_argument("--set", default="custom", help="묶음 이름 → report/custom/tasks_<묶음>.csv (기본 custom)")
    ap.add_argument("--task", default=None, help="과제 이름 (같은 이름끼리 한 과제로 묶인다. 기본: trial 이름)")
    ap.add_argument("--height", type=float, default=None, help="찍힌 사람 키 m (영상: 꼭 주는 게 좋다. OpenCap: 세션에 있으면 생략)"); ap.add_argument("--src-height", type=float, default=None, help="--ik 의 사람 키 m")
    ap.add_argument("--mass", type=float, default=None, help="찍힌 사람 몸무게 kg (영상의 사람 모델용. 로봇 토크에는 안 쓰인다)")
    ap.add_argument("--hand-load", type=float, default=0.0, help="손에 든 물건 kg (구간 내내)"); ap.add_argument("--bimanual", action="store_true", help="양손으로 든다 (반씩)")
    ap.add_argument("--t0", type=float, default=None); ap.add_argument("--t1", type=float, default=None, help="쓸 구간 (초)")
    ap.add_argument("--lowpass", type=float, default=None, help="관절각 저역통과 Hz (기본: OpenCap·영상 4, --ik 6. 점프·달리기처럼 빠른 동작은 8~10)")
    ap.add_argument("--upper-only", action="store_true", help="지면반력 없이 상체 토크만 (앉은 동작, 발이 안 보이는 영상)")
    ap.add_argument("--fov", type=float, default=55.0, help="영상 화각(도)"); ap.add_argument("--footlock", default="auto", choices=("auto", "on", "off"))
    ap.add_argument("--check", action="store_true", help="최종 A안이 이 동작을 낼 수 있나 검사"); ap.add_argument("--redesign", default=None, help="다시 설계할 부위 (쉼표: lower,trunk,shoulder_r,elbow_r)")
    ap.add_argument("--maxfev", type=int, default=300)
    ap.add_argument("--robot", nargs=3, metavar=("MODEL", "MASS", "HEIGHT"), default=None,
                    help="로봇 몸 (기본 data/robot40_h160/robot40_h160.osim 40 1.60). 예: data/robot40_h175sk/robot40_h175sk.osim 40 1.75 — 이때 --check/--redesign 은 160 cm 설계 기준이라 끈다")
    a = ap.parse_args()
    if a.robot:
        ROBOT.update(model=a.robot[0], mass=float(a.robot[1]), height=float(a.robot[2]), tag=os.path.splitext(os.path.basename(a.robot[0]))[0].replace("robot", ""))
        if a.check or a.redesign: print("[주의] --robot 을 바꾸면 --check/--redesign 은 하지 않는다 (최종 A안이 160 cm 기준)"); a.check = False; a.redesign = None
    if sum(bool(x) for x in (a.opencap, a.video, a.ik)) != 1: ap.error("--opencap / --video / --ik 중 하나만")
    import opensim as osim; osim.Logger.setLevelString("Error")
    jobs = []                                                    # (이름, 관절각, 키, dataset, 기본 과제 이름)
    if a.opencap:
        import opencap_import as oi
        sess = os.path.join(a.opencap, "session.json") if os.path.isdir(a.opencap) else None
        if sess and os.path.exists(sess): out_dir, meta = a.opencap, json.load(open(sess, encoding="utf-8"))
        else: out_dir, meta = oi.run(a.opencap, a.name, trials=a.trials)
        h = a.height or meta.get("height_m")
        if not h: ap.error("세션에 키가 없다 — --height 를 준다")
        for tr in meta["trials"]:
            if a.trials and tr not in [safe(x) for x in a.trials]: continue
            jobs.append((f"{meta['name']}_{tr}", os.path.join(out_dir, f"{tr}_ik.mot"), h, "OpenCap", tr))
    elif a.video:
        import video_ik
        name = safe(a.name or os.path.splitext(os.path.basename(a.video))[0]); vdir = os.path.join(ROOT, "data", "custom", a.set, "video"); os.makedirs(vdir, exist_ok=True)
        npz = os.path.join(vdir, f"{name}_pose.npz")
        if not os.path.exists(npz) or os.path.getmtime(npz) < os.path.getmtime(a.video):
            if not os.path.exists(POSE_PY): raise SystemExit(f"`pose` 환경 파이썬이 없다: {POSE_PY}\n  report/custom/README.md 의 '환경' 절대로 만들거나 POSE_PYTHON 환경변수로 알려 준다")
            print(f"[{name}] 영상에서 몸 표면 3D 점 추정 (NLF) …", flush=True)
            subprocess.run([POSE_PY, "-X", "utf8", os.path.join(ROOT, "tools", "video_pose.py"), "--video", a.video, "--out", npz, "--fov", str(a.fov)], check=True)
        if not a.height: print("[주의] --height 가 없다 — 영상이 추정한 키를 그대로 쓴다 (몇 % 틀릴 수 있다)")
        print(f"[{name}] 마커 → 스케일 → IK")
        v = video_ik.run(npz, os.path.join(vdir, name), a.height, a.mass or 70.0, a.footlock)
        jobs.append((name, v["ik"], v["height_m"], "Video", name))
    else:
        if not a.src_height: ap.error("--ik 에는 --src-height (그 사람 키) 가 필요하다")
        jobs.append((safe(a.name or os.path.splitext(os.path.basename(a.ik))[0].replace("_ik", "")), a.ik, a.src_height, "Custom", None))
    lp = a.lowpass or (6.0 if a.ik else VIDEO_LOWPASS)
    rows = []
    for name, ik, h, dset, default_task in jobs:
        try:
            s = process(name, ik, h, a.set, task=a.task or default_task or name, hand_load=a.hand_load, bimanual=a.bimanual, lowpass=lp, t0=a.t0, t1=a.t1, upper_only=a.upper_only, dataset=dset)
            rows.append(s["row"])
        except Exception as e:
            print(f"  [실패] {name}: {type(e).__name__}: {e}")
    if not rows: raise SystemExit("처리된 trial 이 없다")
    tasks = update_tasks(a.set, rows); rel = os.path.relpath(tasks, ROOT).replace("\\", "/")
    print(f"\n과제표 {rel} 에 {len(rows)}줄 (이 파일을 design_directions.py / path_design.py 의 --tasks 로 쓴다)")
    if a.check:
        import fit_check
        print("\n===== 최종 A안이 이 동작을 낼 수 있나 (정적 검사) =====")
        allr = []
        for r in rows:
            print(f"[{r['trial']}]"); allr += fit_check.check_trial(r, [s for s in fit_check.FINAL_A if not (s == "lower" and r["note"].startswith("[상체만]"))])
        outp = os.path.join(ROOT, "report", "custom", f"fit_check_{a.set}_last.csv")
        with open(outp, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=list(allr[0])); w.writeheader(); w.writerows(allr)
        print(f"저장 {os.path.relpath(outp, ROOT)}")
    if a.redesign: redesign(a.set, [s.strip() for s in a.redesign.split(",") if s.strip()], tasks, a.maxfev)
    else: print(f"다음: python report/fit_check.py --tasks {rel}     (검사)\n      python tools/add_motion.py ... --redesign lower,trunk,shoulder_r,elbow_r     (이 묶음을 넣어 다시 설계)")


if __name__ == "__main__":
    main()
