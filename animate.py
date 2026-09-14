"""Simbody 뷰어로 IK 모션 재생 — 눈으로 볼 때 쓴다.

파일로 남길 영상이 필요하면 `gait_video.py` 를 쓴다.

**OpenSim GUI의 체크무늬 지면은 신뢰하지 말 것.** 같은 모델·모션인데 GUI에서는
몸 전체가 지면 아래로 보이고 여기서는 정상으로 보인 적이 있다. 바닥 기준이
필요하면 이쪽을 쓴다 — 지면이 정확히 y=0에 그려진다.

재생 속도·일시정지는 뷰어 창 안에서 조절한다. 창을 닫으면 끝난다.
"""

from __future__ import annotations

import _env  # noqa: F401  환경 확인 — 다른 import 보다 먼저

import argparse
import os
import tempfile


def relax_muscles(model):
    """근육을 강체건·활성도 동역학 없음으로 바꾼다.

    운동학만 재생하는데도 뷰어는 Dynamics 단계까지 realize 한다. 그러면
    Millard 근육의 섬유 속도 Newton 반복이 수렴하지 않고 통째로 죽는다
    (`addlong_r Fiber velocity Newton method did not converge`). 재생에는
    근육 동역학이 필요 없으므로 꺼 버린다.
    """
    import opensim as osim

    muscles = model.updMuscles()
    for i in range(muscles.getSize()):
        muscle = muscles.get(i)
        muscle.set_ignore_tendon_compliance(True)
        muscle.set_ignore_activation_dynamics(True)
    return muscles.getSize()


def main():
    ap = argparse.ArgumentParser(description="IK 모션 재생 (Simbody 뷰어)")
    ap.add_argument("--model", required=True)
    ap.add_argument("--ik", required=True)
    ap.add_argument("--geometry", default=None, help="메시 폴더 (없으면 모델 옆을 찾는다)")
    args = ap.parse_args()

    import opensim as osim
    osim.Logger.setLevelString("Error")

    model_path = os.path.abspath(args.model)
    geometry = args.geometry or os.path.join(os.path.dirname(model_path), "Geometry")
    for cand in (geometry, os.path.join(os.path.dirname(os.path.dirname(model_path)), "Geometry")):
        if os.path.isdir(cand):
            osim.ModelVisualizer.addDirToGeometrySearchPaths(os.path.abspath(cand))

    model = osim.Model(model_path)
    relax_muscles(model)

    # 좌표를 라디안으로 맞춰 임시 파일로 넘긴다.
    # Storage 쪽 출력 API(`print`, `printResult`)는 이 빌드에서 세그폴트로 죽고
    # `exportToTable()` 은 SWIG가 감싸주지 않아 파이썬에서 못 쓴다. 그래서
    # 값을 numpy로 뽑아 직접 쓴다 — 파이프라인의 다른 곳과 같은 경로다.
    import addb_io as io

    pose = io.Pose(model, os.path.abspath(args.ik))

    with tempfile.TemporaryDirectory() as tmp:
        motion_path = os.path.join(tmp, "motion.mot")
        io.write_mot(motion_path, pose.names, pose.time, pose.data,
                     name="motion", in_degrees=False)
        table = osim.TimeSeriesTable(motion_path)
        print(f"  {table.getNumRows()} 프레임 — 뷰어 창에서 속도를 조절한다. 닫으면 끝난다.")
        osim.VisualizerUtilities.showMotion(model, table)


if __name__ == "__main__":
    main()
