"""환경 스모크 테스트 — 모델 로드, 체중, 좌표·바디 존재, 그림 저장.

`python check_env.py` 가 `RESULT: PASS` 를 찍으면 파이프라인을 돌릴 준비가 됐다.
"""

from __future__ import annotations

import os
import sys
import tempfile

DEFAULT_MODEL = os.path.join("dryrun", "Rajagopal", "Scale", "subject_scaled_walk.osim")
NEEDED_COORDS = ("ankle_angle_r", "ankle_angle_l", "knee_angle_r", "knee_angle_l",
                 "hip_flexion_r", "hip_flexion_l")
NEEDED_BODIES = ("calcn_r", "calcn_l", "pelvis")

failures = []


def check(label, fn):
    try:
        detail = fn()
        print(f"  [ OK ] {label}" + (f" — {detail}" if detail else ""))
    except Exception as exc:
        failures.append(label)
        print(f"  [FAIL] {label} — {exc}")


def main():
    model_path = sys.argv[1] if len(sys.argv) > 1 else DEFAULT_MODEL
    print(f"모델 : {model_path}")

    state = {}

    def load_opensim():
        import opensim as osim
        osim.Logger.setLevelString("Error")
        state["osim"] = osim
        return f"OpenSim {osim.__version__}"

    def load_numpy():
        import numpy
        return f"numpy {numpy.__version__}"

    def load_model():
        osim = state["osim"]
        if not os.path.exists(model_path):
            raise FileNotFoundError(model_path)
        model = osim.Model(os.path.abspath(model_path))
        model.initSystem()
        state["model"] = model
        import addb_io as io
        return f"체중 {io.total_mass(model):.2f} kg, 근육 {model.getMuscles().getSize()}개"

    def check_coords():
        coords = state["model"].getCoordinateSet()
        have = {coords.get(i).getName() for i in range(coords.getSize())}
        missing = [c for c in NEEDED_COORDS if c not in have]
        if missing:
            raise KeyError(f"없는 좌표: {', '.join(missing)}")
        return f"좌표 {len(have)}개"

    def check_bodies():
        bodies = state["model"].getBodySet()
        have = {bodies.get(i).getName() for i in range(bodies.getSize())}
        missing = [b for b in NEEDED_BODIES if b not in have]
        if missing:
            raise KeyError(f"없는 바디: {', '.join(missing)}")
        return f"바디 {len(have)}개"

    def check_plot():
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots()
        ax.plot([0, 1], [0, 1])
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "t.png")
            fig.savefig(path)
            size = os.path.getsize(path)
        plt.close(fig)
        return f"savefig {size} bytes"

    def check_modules():
        import addb_io, validate, run_trial, run_muscle_analysis  # noqa: F401
        import run_static_optimization, synergy, batch             # noqa: F401
        return "7개"

    check("opensim import", load_opensim)
    check("numpy import", load_numpy)
    check("파이프라인 모듈 import", check_modules)
    check("모델 로드", load_model)
    check("좌표 존재", check_coords)
    check("바디 존재", check_bodies)
    check("matplotlib 저장", check_plot)

    print(f"\nRESULT: {'FAIL — ' + ', '.join(failures) if failures else 'PASS'}")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
