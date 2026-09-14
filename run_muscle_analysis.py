"""모멘트암 추출 + 길항근 부호 검산.

부호 규약을 하드코딩하지 않는다. 모델마다 축 방향이 다를 수 있어서다. 대신
**길항근 쌍의 부호가 반대인지**와 크기가 0.02~0.08 m인지를 본다. 주동근끼리
부호가 갈리거나 길항근과 부호가 같으면 FAIL.
"""

from __future__ import annotations

import _env  # noqa: F401  환경 확인 — 다른 import 보다 먼저

import argparse
import glob
import os

import numpy as np

import addb_io as io
import validate

LOWPASS = 6.0

DEFAULT_COORDS = ("ankle_angle_r", "ankle_angle_l",
                  "knee_angle_r", "knee_angle_l",
                  "hip_flexion_r", "hip_flexion_l")

#: 좌표 -> (주동근 무리, 길항근 무리). 각 무리는 (전체, 대표)다.
#: 부호는 **전체**로 보고, 크기(0.02~0.08 m)는 **대표**로만 본다 — fdl·perbrev
#: 같은 작은 원위근은 모멘트암이 원래 1 cm 남짓이라 함께 재면 밴드를 벗어난다.
#: 이름은 Rajagopal 계열 기준이고, 모델에 없는 근육은 조용히 빠진다. 두 무리의
#: 부호가 서로 반대이기만 하면 되므로 어느 쪽을 '주동'이라 부르든 상관없다.
ANTAGONISTS = {
    "ankle_angle": (
        (("soleus", "gasmed", "gaslat", "tibpost", "fhl", "fdl", "perlong", "perbrev"),
         ("soleus", "gasmed", "gaslat")),
        (("tibant", "edl", "ehl"),
         ("tibant",)),
    ),
    "knee_angle": (
        (("vasmed", "vaslat", "vasint", "recfem"),
         ("vasmed", "vaslat", "vasint")),
        (("bflh", "semimem", "semiten", "bfsh"),
         ("bflh", "semimem", "semiten")),
    ),
    "hip_flexion": (
        (("psoas", "iliacus", "recfem", "sart", "tfl"),
         ("psoas", "iliacus")),
        (("glmax1", "glmax2", "glmax3", "semimem", "semiten", "bflh", "addmagIsch"),
         ("glmax1", "glmax2", "glmax3")),
    ),
}


def run_muscle_analysis(model_path, ik_path, out_dir, name, *, t0=None, t1=None,
                        lowpass=LOWPASS, keep_all=False, log=print):
    """MuscleAnalysis 실행. 결과 폴더 경로를 돌려준다."""
    import opensim as osim

    os.makedirs(out_dir, exist_ok=True)
    ik_mot = io.read_mot(ik_path)
    t0 = ik_mot.t0 if t0 is None else t0
    t1 = ik_mot.t1 if t1 is None else t1

    model = osim.Model(os.path.abspath(model_path))
    model.initSystem()

    tool = osim.AnalyzeTool()
    tool.setName(name)
    tool.setModel(model)
    tool.setModelFilename(os.path.abspath(model_path))
    tool.setCoordinatesFileName(os.path.abspath(ik_path))
    tool.setLowpassCutoffFrequency(lowpass)
    tool.setLoadModelAndInput(True)
    tool.setInitialTime(t0)
    tool.setFinalTime(t1)
    tool.setResultsDir(os.path.abspath(out_dir))

    ma = osim.MuscleAnalysis()
    ma.setName("MuscleAnalysis")
    ma.setStartTime(t0)
    ma.setEndTime(t1)
    tool.updAnalysisSet().cloneAndAppend(ma)

    setup_path = os.path.join(out_dir, f"{name}_ma_setup.xml")
    tool.printToXML(setup_path)
    # 프린트한 설정을 되읽어 돌린다 — 모델을 두 번 들지 않게 하려는 것.
    osim.AnalyzeTool(setup_path).run()
    if not keep_all:
        prune_outputs(out_dir, log=log)
    return out_dir


#: 남길 MuscleAnalysis 출력. 나머지는 안 쓰면서 trial마다 90개씩 쌓인다.
KEEP_PATTERNS = ("MomentArm_", "_Length.sto", "_TendonForce.sto")


def prune_outputs(out_dir, *, log=print):
    removed = 0
    for path in glob.glob(os.path.join(out_dir, "*MuscleAnalysis*.sto")):
        if not any(k in os.path.basename(path) for k in KEEP_PATTERNS):
            os.remove(path)
            removed += 1
    if removed:
        log(f"  MuscleAnalysis 출력 {removed}개 정리 (모멘트암·길이·건장력만 남김)")
    return removed


def load_moment_arms(ma_dir, name, coords=DEFAULT_COORDS) -> dict:
    """`*_MuscleAnalysis_MomentArm_<좌표>.sto` 들을 읽어 좌표별로."""
    out = {}
    for coord in coords:
        pattern = os.path.join(ma_dir, f"{name}_MuscleAnalysis_MomentArm_{coord}.sto")
        matches = glob.glob(pattern) or glob.glob(
            os.path.join(ma_dir, f"*MomentArm_{coord}.sto"))
        if matches:
            out[coord] = io.read_mot(matches[0])
    return out


def _peaks(ma, muscles, side) -> dict:
    """근육별 |모멘트암|이 최대인 시점의 부호 있는 값."""
    out = {}
    for base in muscles:
        col = f"{base}_{side}"
        if col not in ma.names:
            continue
        arm = ma.column(col)
        peak = float(arm[np.argmax(np.abs(arm))])
        if abs(peak) >= 1e-6:
            out[col] = peak
    return out


def _group_stats(ma, group, side):
    """(부호들, 대표 근육 크기 중앙값, 참여 근육 수). 없으면 None."""
    all_muscles, primary = group
    peaks = _peaks(ma, all_muscles, side)
    if not peaks:
        return None
    prim = _peaks(ma, primary, side) or peaks
    signs = [np.sign(v) for v in peaks.values()]
    return signs, float(np.median([abs(v) for v in prim.values()])), len(peaks)


def validate_moment_arms(ma_by_coord, coords=DEFAULT_COORDS) -> validate.Report:
    rep = validate.Report()
    for coord in coords:
        ma = ma_by_coord.get(coord)
        if ma is None:
            continue
        base, side = coord.rsplit("_", 1)
        pair = ANTAGONISTS.get(base)
        if pair is None:
            continue
        agon = _group_stats(ma, pair[0], side)
        antag = _group_stats(ma, pair[1], side)
        if agon is None or antag is None:
            continue

        a_signs, a_mag, a_n = agon
        b_signs, b_mag, b_n = antag
        problems = []
        if len(set(a_signs)) > 1:
            problems.append("주동근끼리 부호가 갈린다")
        if len(set(b_signs)) > 1:
            problems.append("길항근끼리 부호가 갈린다")
        if a_signs[0] == b_signs[0]:
            problems.append("길항근이 주동근과 같은 부호다")

        for label, mag, n in ((f"ma/{coord}+", a_mag, a_n), (f"ma/{coord}-", b_mag, b_n)):
            check = rep.add(label, mag, "moment_arm", note=f"n={n}")
            if problems:
                check.verdict = "FAIL"
                check.note = f"n={n}  " + " / ".join(problems)
    return rep


def main():
    ap = argparse.ArgumentParser(description="모멘트암 추출 + 부호 검산")
    ap.add_argument("--model", required=True)
    ap.add_argument("--ik", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--start", type=float, default=None)
    ap.add_argument("--end", type=float, default=None)
    ap.add_argument("--lowpass", type=float, default=LOWPASS)
    ap.add_argument("--coords", default=",".join(DEFAULT_COORDS))
    ap.add_argument("--keep-all", action="store_true", help="MuscleAnalysis 출력 전부 보관")
    args = ap.parse_args()

    import opensim as osim
    osim.Logger.setLevelString("Error")

    coords = tuple(c.strip() for c in args.coords.split(",") if c.strip())
    ma_dir = run_muscle_analysis(args.model, args.ik, args.out, args.name,
                                 t0=args.start, t1=args.end, lowpass=args.lowpass,
                                 keep_all=args.keep_all)
    ma = load_moment_arms(ma_dir, args.name, coords)
    print(f"  모멘트암 파일 {len(ma)}개 / 좌표 {len(coords)}개")
    rep = validate_moment_arms(ma, coords)
    print(rep.render())
    print(f"  판정 : {rep.verdict}")


if __name__ == "__main__":
    main()
