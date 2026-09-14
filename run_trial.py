"""trial 1개: 정합 확인 → 매핑 → 구간 결정 → ID → 검산 → 그래프."""

from __future__ import annotations

import _env  # noqa: F401  환경 확인 — 다른 import 보다 먼저

import argparse
import os

import numpy as np

import addb_io as io
import validate

#: 좌표 저역통과 컷오프. 6 Hz가 옳다 — 15 Hz면 잔차 22.9 %, 무필터면 98.1 %.
LOWPASS = 6.0
#: ID에서 제외할 힘. 근육이 관절 모멘트를 상쇄해버리면 안 된다.
EXCLUDE = ("Muscles", "Actuators")
#: 기본 접지 바디
FEET = ("calcn_r", "calcn_l")


def resolve_window(ik_mot, grf_mot, plates, weight_n, *, lowpass=LOWPASS,
                   start=None, end=None, log=print):
    """세 가지를 교집합해 분석 구간을 정한다.

      1. IK ∩ GRF 시간 교집합    — 한쪽에만 데이터가 있는 구간을 뺀다
      2. 지지 구간               — 발이 힘판 위에 없으면 체중 전체가 잔차로 잡힌다
      3. 양 끝 마진 = 1/컷오프    — 필터 후 미분한 가속도가 경계에서 망가진다

    셋을 다 적용하면 사람이 손으로 좋은 구간을 고른 것과 같은 결과가 나온다
    (드라이런 실측: 전체 93.6 % → 지지구간 24.4 % → 마진까지 7.8 %).
    """
    t0 = max(ik_mot.t0, grf_mot.t0)
    t1 = min(ik_mot.t1, grf_mot.t1)

    total = io.total_vertical(grf_mot, plates)
    s0, s1 = io.grf_support_window(grf_mot.time, total, weight_n)
    t0, t1 = max(t0, s0), min(t1, s1)

    margin = 1.0 / lowpass if lowpass and lowpass > 0 else 0.0
    t0, t1 = t0 + margin, t1 - margin

    if start is not None:
        t0 = start
    if end is not None:
        t1 = end
    if t1 <= t0:
        raise ValueError(f"분석 구간이 비었다 (t0={t0:.3f}, t1={t1:.3f})")
    log(f"  구간 : {t0:.3f} ~ {t1:.3f} s  (지지 {s0:.3f}~{s1:.3f}, 마진 {margin:.3f})")
    return t0, t1


def write_external_loads(path, grf_path, plates, mapping):
    """어느 힘판을 어느 발에 붙였는지를 XML로 남긴다."""
    import opensim as osim

    loads = osim.ExternalLoads()
    loads.setName("external_loads")
    loads.setDataFileName(os.path.abspath(grf_path))
    for key, info in mapping.items():
        plate = plates[key]
        force = osim.ExternalForce()
        force.setName(f"{key}_GRF")
        force.set_applied_to_body(info["body"])
        force.set_force_expressed_in_body("ground")
        force.set_point_expressed_in_body("ground")
        force.set_force_identifier(plate.identifier("force"))
        force.set_point_identifier(plate.identifier("point"))
        force.set_torque_identifier(plate.identifier("torque") or "")
        loads.cloneAndAppend(force)
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    loads.printToXML(path)
    return path


def run_inverse_dynamics(model_path, ik_path, ext_path, out_dir, name, t0, t1,
                         lowpass=LOWPASS, log=print):
    import opensim as osim

    os.makedirs(out_dir, exist_ok=True)
    tool = osim.InverseDynamicsTool()
    tool.setName(name)
    tool.setModelFileName(os.path.abspath(model_path))
    tool.setCoordinatesFileName(os.path.abspath(ik_path))
    tool.setExternalLoadsFileName(os.path.abspath(ext_path))
    tool.setLowpassCutoffFrequency(lowpass)
    tool.setStartTime(t0)
    tool.setEndTime(t1)
    tool.setResultsDir(os.path.abspath(out_dir))
    tool.setOutputGenForceFileName(f"{name}_id.sto")

    excluded = osim.ArrayStr()
    for e in EXCLUDE:
        excluded.append(e)
    tool.setExcludedForces(excluded)

    setup_path = os.path.join(out_dir, f"{name}_id_setup.xml")
    tool.printToXML(setup_path)
    if not tool.run():
        raise RuntimeError("InverseDynamicsTool 실패")
    return os.path.join(out_dir, f"{name}_id.sto"), setup_path


def summarize(model_path, ik_path, grf_path, id_path, out_dir, name, *, feet=FEET,
              lowpass=LOWPASS, start=None, end=None, log=print):
    """이미 있는 ID 결과로 요약·검산만 다시 만든다. ID 툴은 돌리지 않는다.

    배치가 "결과 있으면 건너뜀"으로 넘어갈 때 이걸 쓴다. 안 그러면 `summary.csv`
    의 그 줄이 통째로 비어 요약 표에 구멍이 난다.
    """
    import opensim as osim

    model = osim.Model(os.path.abspath(model_path))
    model.initSystem()
    mass = io.total_mass(model)
    weight_n = io.model_weight_n(model)

    grf_mot = io.read_mot(grf_path)
    plates = {k: p for k, p in io.find_plates(grf_mot).items() if p.complete}
    mapping = io.assign_plates_to_feet(model, ik_path, grf_mot, plates,
                                       feet=feet, weight_n=weight_n, log=lambda *a: None)
    id_mot = io.read_mot(id_path)
    t0, t1 = id_mot.t0, id_mot.t1

    mask = (grf_mot.time >= t0) & (grf_mot.time <= t1)
    total = io.total_vertical(grf_mot, plates)
    grf_peak_n = float(total[mask].max()) if mask.any() else float(total.max())
    cop_dist = max(info["dist"] for info in mapping.values()) if mapping else None

    rep = validate.validate_id(id_mot, mass_kg=mass, grf_peak_bw=grf_peak_n / weight_n,
                               grf_peak_n=grf_peak_n, cop_dist=cop_dist, window=(t0, t1))
    summary = {
        "model": os.path.basename(model_path),
        "mass_kg": round(mass, 3),
        "t0": round(t0, 3), "t1": round(t1, 3),
        "grf_peak": round(grf_peak_n / weight_n, 3),
        "cop_dist": round(cop_dist, 3) if cop_dist is not None else None,
        "id_verdict": rep.verdict,
        "id_path": id_path,
        "external_loads": os.path.join(out_dir, f"{name}_external_loads.xml"),
    }
    for check in rep.checks:
        if check.band in ("ankle", "knee", "hip"):
            summary[check.name] = round(check.value, 3)
        elif check.name == "residual":
            summary["residual"] = round(check.value, 3)
    return summary


def run_trial(model_path, ik_path, grf_path, out_dir, name, *, feet=FEET,
              lowpass=LOWPASS, start=None, end=None, plot=True, log=print):
    """trial 하나를 끝까지. 요약 dict를 돌려준다."""
    import opensim as osim

    model = osim.Model(os.path.abspath(model_path))
    model.initSystem()
    mass = io.total_mass(model)
    weight_n = io.model_weight_n(model)
    log(f"모델 : {model_path}  (체중 {mass:.2f} kg)")

    ik_mot = io.read_mot(ik_path)
    grf_mot = io.read_mot(grf_path)
    plates = {k: p for k, p in io.find_plates(grf_mot).items() if p.complete}
    if not plates:
        raise ValueError(f"GRF 힘판을 찾지 못했다: {grf_path}")

    mapping = io.assign_plates_to_feet(model, ik_path, grf_mot, plates,
                                       feet=feet, weight_n=weight_n, log=log)
    if not mapping:
        raise ValueError("힘판을 발에 배정하지 못했다")

    t0, t1 = resolve_window(ik_mot, grf_mot, plates, weight_n,
                            lowpass=lowpass, start=start, end=end, log=log)

    os.makedirs(out_dir, exist_ok=True)
    ext_path = write_external_loads(
        os.path.join(out_dir, f"{name}_external_loads.xml"), grf_path, plates, mapping)
    id_path, setup_path = run_inverse_dynamics(
        model_path, ik_path, ext_path, out_dir, name, t0, t1, lowpass=lowpass, log=log)

    id_mot = io.read_mot(id_path)

    window_mask = (grf_mot.time >= t0) & (grf_mot.time <= t1)
    total = io.total_vertical(grf_mot, plates)
    grf_peak_n = float(total[window_mask].max()) if window_mask.any() else float(total.max())
    cop_dist = max(info["dist"] for info in mapping.values())

    rep = validate.validate_id(id_mot, mass_kg=mass, grf_peak_bw=grf_peak_n / weight_n,
                               grf_peak_n=grf_peak_n, cop_dist=cop_dist, window=(t0, t1))
    log(rep.render())
    log(f"  판정 : {rep.verdict}")

    png = None
    if plot:
        png = os.path.join(out_dir, f"{name}_id.png")
        validate.plot_id(png, id_mot, grf_mot, plates, mass, weight_n, window=(t0, t1))

    summary = {
        "model": os.path.basename(model_path),
        "mass_kg": round(mass, 3),
        "t0": round(t0, 3),
        "t1": round(t1, 3),
        "grf_peak": round(grf_peak_n / weight_n, 3),
        "cop_dist": round(cop_dist, 3),
        "id_verdict": rep.verdict,
        "id_path": id_path,
        "external_loads": ext_path,
        "setup": setup_path,
        "png": png,
        "mapping": {k: v["body"] for k, v in mapping.items()},
    }
    for check in rep.checks:
        if check.band in ("ankle", "knee", "hip"):
            summary[check.name] = round(check.value, 3)
        elif check.name == "residual":
            summary["residual"] = round(check.value, 3)
    return summary


def main():
    ap = argparse.ArgumentParser(description="trial 1개 역동역학")
    ap.add_argument("--model", required=True)
    ap.add_argument("--ik", required=True)
    ap.add_argument("--grf", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--name", required=True)
    ap.add_argument("--lowpass", type=float, default=LOWPASS)
    ap.add_argument("--start", type=float, default=None)
    ap.add_argument("--end", type=float, default=None)
    ap.add_argument("--feet", default=",".join(FEET))
    ap.add_argument("--no-plot", action="store_true")
    args = ap.parse_args()

    import opensim as osim
    osim.Logger.setLevelString("Error")

    run_trial(args.model, args.ik, args.grf, args.out, args.name,
              feet=tuple(args.feet.split(",")), lowpass=args.lowpass,
              start=args.start, end=args.end, plot=not args.no_plot)


if __name__ == "__main__":
    main()
