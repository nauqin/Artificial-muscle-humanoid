"""`.b3d` → `.osim` / `.mot` 변환. **`b3d` 환경에서 돌린다.**

    conda activate b3d
    python b3d_to_opensim.py raw_b3d --out data --walking-only

AddBiomechanics 공개 데이터셋은 subject 하나당 `.b3d` 바이너리 하나로 배포된다.
스케일된 모델과 모든 trial의 관절 각도·지면반력이 다 들어있지만 읽으려면
`nimblephysics`가 필요하고, 그 패키지는 OpenSim과 numpy 요구가 달라 환경을
따로 판다.

기본값은 `missingGRFReason`이 정상인 가장 긴 연속 구간만 남기는 것이다. GRF가
없는 프레임을 ID에 넣으면 체중 전체가 골반 잔차로 잡혀 결과가 통째로 무의미해진다.
"""

from __future__ import annotations

import argparse
import os
import re
import sys

import numpy as np

import addb_io as io   # numpy만 쓰므로 b3d 환경에서도 import 된다

#: 이 토큰이 이름에 있으면 보행/달리기로 본다. 판정은 trial 이름으로 한다 —
#: `DJ`(드롭점프), `STS`, `squats`, `static` 을 걸러내려는 것.
WALKING_TOKENS = ("walk", "run", "gait", "tread")
MIN_FRAMES = 50
AXES = ("x", "y", "z")

#: `--allow-manual-review` 로 눈감아 줄 사유. 프레임별 물리 판정이 아니라 trial
#: 전체에 균일하게 찍히는 "사람이 검토하지 않음" 표시라, 데이터가 나쁘다는 뜻이
#: 아니라 품질 보증이 안 됐다는 뜻이다. 실제로 subject10은 이 플래그가 붙었는데
#: AddBiomechanics 자체 잔차가 subject3보다도 낮았다. 물리로 판단하는 것이
#: 이 파이프라인의 방침이므로 허용하고 우리 검산기로 거른다.
SOFT_REASONS = ("manualReview",)


# ─────────────────────────────────────────────────────────────────────────────
# 열기
# ─────────────────────────────────────────────────────────────────────────────


def open_subject(path):
    import nimblephysics as nimble

    return nimble.biomechanics.SubjectOnDisk(path), nimble


def pick_pass(subject, nimble, requested=None) -> int:
    """쓸 processing pass. 기본은 마지막 `dynamics` pass."""
    n = subject.getNumProcessingPasses()
    if requested is not None:
        if not 0 <= requested < n:
            raise SystemExit(f"pass {requested}가 없다 (0~{n-1})")
        return requested
    dynamics = nimble.biomechanics.ProcessingPassType.DYNAMICS
    for i in reversed(range(n)):
        if subject.getProcessingPassType(i) == dynamics:
            return i
    return n - 1


def dof_names(subject, pass_idx) -> list:
    osim_file = subject.readOpenSimFile(pass_idx, "", True)
    skel = osim_file.skeleton
    return [skel.getDofByIndex(i).getName() for i in range(skel.getNumDofs())]


def count_muscles(osim_text: str) -> int:
    return len(re.findall(r"<\w*Muscle\s+name=", osim_text))


# ─────────────────────────────────────────────────────────────────────────────
# GRF 양호 구간
# ─────────────────────────────────────────────────────────────────────────────


def good_frames(subject, nimble, trial, allow_soft=False):
    """프레임별 GRF 양호 여부와 사유 집계."""
    reasons = subject.getMissingGRF(trial)
    ok_flag = nimble.biomechanics.MissingGRFReason.notMissingGRF
    tolerated = {getattr(nimble.biomechanics.MissingGRFReason, r)
                 for r in SOFT_REASONS} if allow_soft else set()
    mask = np.array([r == ok_flag or r in tolerated for r in reasons], dtype=bool)
    tally = {}
    for r, ok in zip(reasons, mask):
        if not ok:
            key = str(r).rsplit(".", 1)[-1]
            tally[key] = tally.get(key, 0) + 1
    return mask, tally


def longest_run(mask):
    """True 최장 연속 구간의 (시작, 길이). 없으면 (0, 0)."""
    best = (0, 0)
    start = None
    for i, v in enumerate(mask):
        if v and start is None:
            start = i
        elif not v and start is not None:
            if i - start > best[1]:
                best = (start, i - start)
            start = None
    if start is not None and len(mask) - start > best[1]:
        best = (start, len(mask) - start)
    return best


# ─────────────────────────────────────────────────────────────────────────────
# 변환
# ─────────────────────────────────────────────────────────────────────────────


def grf_columns(bodies) -> list:
    names = []
    for body in bodies:
        names += [f"ground_force_{body}_v{a}" for a in AXES]
        names += [f"ground_force_{body}_p{a}" for a in AXES]
        names += [f"ground_force_{body}_m{a}" for a in AXES]
    return names


def convert_trial(subject, nimble, trial, pass_idx, out_dir, subject_name, args, log=print):
    name = subject.getTrialName(trial) or f"trial{trial}"
    length = subject.getTrialLength(trial)
    dt = subject.getTrialTimestep(trial)

    if args.walking_only and not any(tok in name.lower() for tok in WALKING_TOKENS):
        return None

    mask, tally = good_frames(subject, nimble, trial, args.allow_manual_review)
    if args.keep_all_frames:
        start, count = 0, length
    else:
        start, count = longest_run(mask)

    if count < args.min_frames:
        detail = "  ".join(f"{k}×{v}" for k, v in sorted(tally.items())) or "-"
        log(f"  [건너뜀] {name}  GRF 양호 프레임 {int(mask.sum())} / {length} — {detail}")
        return None

    frames = subject.readFrames(trial, start, count, includeSensorData=False,
                                includeProcessingPasses=True)
    bodies = list(subject.getGroundForceBodies())
    dofs = dof_names(subject, pass_idx)

    times = np.empty(len(frames))
    pos = np.empty((len(frames), len(dofs)))
    grf = np.zeros((len(frames), len(bodies) * 9))

    for i, frame in enumerate(frames):
        p = frame.processingPasses[min(pass_idx, len(frame.processingPasses) - 1)]
        times[i] = (start + i) * dt
        pos[i] = np.asarray(p.pos).ravel()[:len(dofs)]
        force = np.asarray(p.groundContactForce).ravel()
        cop = np.asarray(p.groundContactCenterOfPressure).ravel()
        torque = np.asarray(p.groundContactTorque).ravel()
        for b in range(len(bodies)):
            grf[i, b * 9 + 0: b * 9 + 3] = force[b * 3: b * 3 + 3]
            grf[i, b * 9 + 3: b * 9 + 6] = cop[b * 3: b * 3 + 3]
            grf[i, b * 9 + 6: b * 9 + 9] = torque[b * 3: b * 3 + 3]

    os.makedirs(out_dir, exist_ok=True)
    ik_path = os.path.join(out_dir, f"{name}_ik.mot")
    grf_path = os.path.join(out_dir, f"{name}_grf.mot")
    # 라디안 그대로 쓴다 — inDegrees=no 로 표시하면 OpenSim이 알아서 읽는다.
    io.write_mot(ik_path, dofs, times, pos, name=f"{name}_ik", in_degrees=False)
    io.write_mot(grf_path, grf_columns(bodies), times, grf, name=f"{name}_grf")

    log(f"  trial {name}  {count} 프레임 ({times[0]:.2f}~{times[-1]:.2f}s, {count}/{length} 사용)")
    return {"name": name, "frames": count, "ik": ik_path, "grf": grf_path}


def convert_subject(path, out_root, args, log=print):
    subject, nimble = open_subject(path)
    subject_name = os.path.splitext(os.path.basename(path))[0]
    pass_idx = pick_pass(subject, nimble, args.pass_index)
    out_dir = os.path.join(out_root, subject_name)

    log(f"\n=== {subject_name} ===")
    log(f"  질량 {subject.getMassKg():.2f} kg, 키 {subject.getHeightM():.2f} m, "
        f"trial {subject.getNumTrials()}개, pass {pass_idx} "
        f"({str(subject.getProcessingPassType(pass_idx)).rsplit('.', 1)[-1]})")

    osim_text = subject.getOpensimFileText(pass_idx)
    n_muscles = count_muscles(osim_text)
    if n_muscles:
        log(f"  모델 근육 {n_muscles}개 — SO/시너지 가능")
    else:
        log(f"  [WARN] 모델에 근육이 없다 — ID까지만 된다")

    os.makedirs(out_dir, exist_ok=True)
    osim_path = os.path.join(out_dir, f"{subject_name}_scaled.osim")
    with open(osim_path, "w") as fh:
        fh.write(osim_text)

    made = []
    for trial in range(subject.getNumTrials()):
        try:
            result = convert_trial(subject, nimble, trial, pass_idx, out_dir,
                                   subject_name, args, log=log)
        except Exception as exc:
            log(f"  [ERROR] trial {trial}: {exc}")
            continue
        if result:
            made.append(result)
    log(f"  → {len(made)} trial 변환, {out_dir}")
    return made


def inspect_subject(path, args, log=print):
    subject, nimble = open_subject(path)
    name = os.path.splitext(os.path.basename(path))[0]
    pass_idx = pick_pass(subject, nimble, args.pass_index)
    log(f"\n=== {name} ===")
    log(f"  질량 {subject.getMassKg():.2f} kg  키 {subject.getHeightM():.2f} m  "
        f"DOF {subject.getNumDofs()}  trial {subject.getNumTrials()}개")
    log(f"  pass {pass_idx} / {subject.getNumProcessingPasses()}  "
        f"접지 바디 {list(subject.getGroundForceBodies())}")
    log(f"  근육 {count_muscles(subject.getOpensimFileText(pass_idx))}개")

    log(f"  {'trial':<26s} {'프레임':>7s} {'GRF양호':>8s} {'잔차 N (AddB)':>18s} {'마커mm':>7s}  사유")
    for trial in range(subject.getNumTrials()):
        tname = subject.getTrialName(trial) or f"trial{trial}"
        length = subject.getTrialLength(trial)
        mask, tally = good_frames(subject, nimble, trial)
        rtxt = _stat(subject.getTrialLinearResidualNorms, trial, pass_idx, "평균{:.1f} 최대{:.1f}")
        mtxt = _stat(subject.getTrialMarkerRMSs, trial, pass_idx, "{:.0f}", scale=1000.0, mean_only=True)
        detail = "  ".join(f"{k}×{v}" for k, v in sorted(tally.items())) or "-"
        log(f"  {tname:<26s} {length:7d} {int(mask.sum()):8d} {rtxt:>18s} {mtxt:>7s}  {detail}")


def _stat(getter, trial, pass_idx, fmt, *, scale=1.0, mean_only=False):
    """`getTrial*(trial, pass)` 결과를 요약. pass 인자를 빠뜨리면 조용히 '-'가 된다."""
    try:
        vals = np.asarray(getter(trial, pass_idx), dtype=float) * scale
        vals = vals[np.isfinite(vals)]
        if not vals.size:
            return "-"
        return fmt.format(vals.mean()) if mean_only else fmt.format(vals.mean(), vals.max())
    except Exception:
        return "-"


# ─────────────────────────────────────────────────────────────────────────────


def find_b3d(paths) -> list:
    found = []
    for path in paths:
        if os.path.isdir(path):
            for dirpath, _dirs, files in os.walk(path):
                found += [os.path.join(dirpath, f) for f in files if f.lower().endswith(".b3d")]
        elif path.lower().endswith(".b3d"):
            found.append(path)
    return sorted(set(found))


def main():
    ap = argparse.ArgumentParser(description=".b3d → .osim/.mot 변환 (b3d 환경)")
    ap.add_argument("inputs", nargs="+", help=".b3d 파일 또는 폴더")
    ap.add_argument("--out", default="data")
    ap.add_argument("--inspect", action="store_true", help="변환하지 않고 안을 훑기만")
    ap.add_argument("--walking-only", action="store_true", help="보행/달리기 trial만")
    ap.add_argument("--pass", dest="pass_index", type=int, default=None)
    ap.add_argument("--keep-all-frames", action="store_true", help="GRF 없는 구간도 그대로")
    ap.add_argument("--allow-manual-review", action="store_true",
                    help="manualReview 플래그를 무시한다 (물리 판정은 우리 검산기가 한다)")
    ap.add_argument("--min-frames", type=int, default=MIN_FRAMES)
    args = ap.parse_args()

    try:
        import nimblephysics  # noqa: F401
    except ImportError:
        raise SystemExit(
            "nimblephysics 가 없다. `conda activate b3d` 를 먼저 하거나,\n"
            "  conda create -n b3d python=3.11 \"numpy<2\" -y\n"
            "  conda activate b3d && python -m pip install \"nimblephysics==0.10.52.1\"\n"
            "※ 맥에서는 `pip`이 환경 밖 인터프리터를 가리킬 수 있다. `python -m pip`를 쓸 것.")

    files = find_b3d(args.inputs)
    if not files:
        raise SystemExit(f".b3d 를 찾지 못했다: {args.inputs}")
    print(f".b3d {len(files)}개")

    for path in files:
        try:
            if args.inspect:
                inspect_subject(path, args)
            else:
                convert_subject(path, args.out, args)
        except Exception as exc:
            print(f"[ERROR] {path}: {exc}", file=sys.stderr)


if __name__ == "__main__":
    main()
