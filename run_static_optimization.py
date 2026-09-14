"""근육별 힘 추정 + 보조 액추에이터 생성.

근육만으로 관절 회전력을 못 맞추면 최적화 문제가 아예 안 풀린다. 그래서 모든
자유 좌표에 `CoordinateActuator`를 하나씩 단다.

- `optimal_force = 1` — 활성도가 곧 실제 일반화 힘이 된다. 비용이 활성도
  제곱합이라 근육(활성도 ≤ 1)에 비해 수백~수만 배 비싸져서 최적화가 알아서 피한다.
- 한도 ±1000 — 그래도 정말 필요하면 쓸 수 있게 넉넉히.
- **종속 좌표는 제외** — knee_angle_r_beta 같은 구속된 좌표에 액추에이터를 달면
  SO가 깨진다.
"""

from __future__ import annotations

import _env  # noqa: F401  환경 확인 — 다른 import 보다 먼저

import argparse
import os

import addb_io as io
import run_muscle_analysis as rma
import validate

LOWPASS = 6.0
#: 활성도 1이 곧 1 N·m 이 되게 한다. 근육에 비해 압도적으로 비싸진다.
RESERVE_OPTIMAL_FORCE = 1.0
RESERVE_LIMIT = 1000.0
RESERVE_PREFIX = "reserve_"


def build_reserve_actuators(model_path, out_path, *, log=print):
    """모든 **비구속** 좌표에 CoordinateActuator를 하나씩. 경로와 개수를 반환."""
    import opensim as osim

    model = osim.Model(os.path.abspath(model_path))
    state = model.initSystem()

    forces = osim.ForceSet()
    coords = model.getCoordinateSet()
    skipped = []
    for i in range(coords.getSize()):
        coord = coords.get(i)
        name = coord.getName()
        if coord.isConstrained(state) or coord.get_locked():
            skipped.append(name)
            continue
        act = osim.CoordinateActuator(name)
        act.setName(f"{RESERVE_PREFIX}{name}")
        act.setOptimalForce(RESERVE_OPTIMAL_FORCE)
        act.setMinControl(-RESERVE_LIMIT)
        act.setMaxControl(RESERVE_LIMIT)
        forces.cloneAndAppend(act)

    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    forces.printToXML(out_path)
    log(f"  보조 액추에이터 {forces.getSize()}개"
        + (f" (종속·잠긴 좌표 {len(skipped)}개 제외: {', '.join(skipped[:4])}"
           + ("…" if len(skipped) > 4 else "") + ")" if skipped else ""))
    return out_path, forces.getSize()


def run_static_optimization(model_path, ik_path, ext_path, out_dir, name, t0, t1,
                            *, lowpass=LOWPASS, log=print):
    import opensim as osim

    os.makedirs(out_dir, exist_ok=True)
    reserve_path, _n = build_reserve_actuators(
        model_path, os.path.join(out_dir, f"{name}_reserve_actuators.xml"), log=log)

    tool = osim.AnalyzeTool()
    tool.setName(name)
    tool.setModelFilename(os.path.abspath(model_path))
    tool.setCoordinatesFileName(os.path.abspath(ik_path))
    tool.setExternalLoadsFileName(os.path.abspath(ext_path))
    tool.setLowpassCutoffFrequency(lowpass)
    tool.setInitialTime(t0)
    tool.setFinalTime(t1)
    tool.setResultsDir(os.path.abspath(out_dir))
    tool.setReplaceForceSet(False)

    force_files = osim.ArrayStr()
    force_files.append(os.path.abspath(reserve_path))
    tool.setForceSetFiles(force_files)

    so = osim.StaticOptimization()
    so.setName("StaticOptimization")
    so.setStartTime(t0)
    so.setEndTime(t1)
    so.setUseModelForceSet(True)
    tool.updAnalysisSet().cloneAndAppend(so)

    setup_path = os.path.join(out_dir, f"{name}_so_setup.xml")
    tool.printToXML(setup_path)
    osim.AnalyzeTool(setup_path).run()

    force_path = os.path.join(out_dir, f"{name}_StaticOptimization_force.sto")
    act_path = os.path.join(out_dir, f"{name}_StaticOptimization_activation.sto")
    if not os.path.exists(force_path):
        raise RuntimeError("StaticOptimization 출력이 없다")
    return force_path, act_path, reserve_path


def validate_run(id_path, force_path, ma_dir, ma_name, coords, log=print):
    """SO 결과를 ID 모멘트와 대조."""
    id_mot = io.read_mot(id_path)
    force_mot = io.read_mot(force_path)
    ma = rma.load_moment_arms(ma_dir, ma_name, coords)
    missing = [c for c in coords if c not in ma]
    if missing:
        log(f"  [WARN] 모멘트암이 없는 좌표: {', '.join(missing)}")
    rep = validate.validate_so(id_mot, force_mot, ma, coords,
                               reserve_prefix=RESERVE_PREFIX)
    return rep


def main():
    ap = argparse.ArgumentParser(description="Static Optimization + 검산")
    ap.add_argument("--model", required=True)
    ap.add_argument("--ik", required=True)
    ap.add_argument("--grf", required=True)
    ap.add_argument("--id", required=True, help="ID 결과 .sto")
    ap.add_argument("--out", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--ma-dir", default=None, help="생략하면 모멘트암을 직접 계산")
    ap.add_argument("--ext", default=None, help="ExternalLoads XML (생략하면 새로 만든다)")
    ap.add_argument("--start", type=float, default=None)
    ap.add_argument("--end", type=float, default=None)
    ap.add_argument("--lowpass", type=float, default=LOWPASS)
    ap.add_argument("--coords", default=",".join(rma.DEFAULT_COORDS))
    args = ap.parse_args()

    import opensim as osim
    osim.Logger.setLevelString("Error")

    import run_trial

    coords = tuple(c.strip() for c in args.coords.split(",") if c.strip())
    id_mot = io.read_mot(args.id)
    t0 = args.start if args.start is not None else id_mot.t0
    t1 = args.end if args.end is not None else id_mot.t1

    ext = args.ext
    if ext is None:
        grf_mot = io.read_mot(args.grf)
        plates = {k: p for k, p in io.find_plates(grf_mot).items() if p.complete}
        model = osim.Model(os.path.abspath(args.model))
        model.initSystem()
        mapping = io.assign_plates_to_feet(model, args.ik, grf_mot, plates)
        ext = run_trial.write_external_loads(
            os.path.join(args.out, f"{args.name}_external_loads.xml"),
            args.grf, plates, mapping)

    ma_dir = args.ma_dir
    if ma_dir is None:
        ma_dir = os.path.join(args.out, "ma")
        rma.run_muscle_analysis(args.model, args.ik, ma_dir, args.name,
                                t0=t0, t1=t1, lowpass=args.lowpass)

    force_path, act_path, _ = run_static_optimization(
        args.model, args.ik, ext, args.out, args.name, t0, t1, lowpass=args.lowpass)
    print(f"  근육 힘   : {force_path}")
    print(f"  활성도    : {act_path}")

    rep = validate_run(args.id, force_path, ma_dir, args.name, coords)
    print(rep.render())
    print(f"  판정 : {rep.verdict}")


if __name__ == "__main__":
    main()
