"""Phase 0 — 일-용량 조건으로 모듈 개수의 하한을 낸다.

당기기만 하는 액추에이터 하나가 ROM 전체(각도폭 `ΔΘ`)를 덮으면서 어느 자세에서든
요구 토크 `τ_req` 를 내려면 모멘트암 `r` 이

    τ_req / (n·F_max)  ≤  r  ≤  stroke / ΔΘ

를 만족해야 한다. 좌변은 힘 조건, 우변은 변위 조건이다. 둘이 만나려면
`n·F_max·stroke ≥ τ_req·ΔΘ` 여야 하고, 그 비가 곧 **같은 방향으로 병렬로 필요한
모듈 개수의 하한**이다.

이건 자릿수를 보는 스크리닝이다. 실제 판정은 Phase 4 의 실현가능성 판정기가 한다.
여기서 쓰는 `τ_req` 는 관절별 **피크**이고 ROM 끝에서 나는 값이 아니므로,
"피크 토크를 ROM 양 끝에서도 내야 한다"는 보수적 가정이 들어가 있다.
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os

import _env  # noqa: F401
import numpy as np

import actuator_model as am
import addb_io as io
import models

ROBOT_MASSES_KG = (40.0, 55.0, 70.0)


# ─────────────────────────────────────────────────────────────────────────────
# 요구 토크 읽기
# ─────────────────────────────────────────────────────────────────────────────


def peaks_from_id_dir(id_dir, data_root, exclude=(), log=print) -> dict:
    """`*_id.sto` 를 훑어 좌표·방향별 요구 토크 피크 (Nm/kg).

    좌표 6개를 전부 얻으려면 ID 결과를 직접 봐야 한다 — `summary.csv` 에는
    발목·무릎·고관절 굴곡 3개만 들어 있다.
    """
    import opensim as osim
    osim.Logger.setLevelString("Error")

    masses = {}
    out = {}
    paths = sorted(glob.glob(os.path.join(id_dir, "*", "*", "*_id.sto")))
    if not paths:
        raise SystemExit(f"ID 결과가 없다: {id_dir}")

    for path in paths:
        trial_dir = os.path.dirname(path)
        subject = os.path.basename(os.path.dirname(trial_dir))
        trial = os.path.basename(trial_dir)
        if subject in exclude:
            continue
        if subject not in masses:
            model_path = io.pick_model(os.path.join(data_root, subject), log=lambda *a: None)
            model = osim.Model(model_path)
            model.initSystem()
            masses[subject] = io.total_mass(model)
        mass = masses[subject]

        mot = io.read_mot(path)
        for coord in models.LOWER_LIMB_COORDS + models.OPTIONAL_COORDS:
            for side in models.SIDES:
                col = f"{coord}_{side}_moment"
                if col not in mot.names:
                    continue
                v = mot.column(col) / mass
                for sign, value in (("+", float(v.max())), ("-", float(-v.min()))):
                    key = (coord, sign)
                    if value > out.get(key, (0.0, None, None))[0]:
                        out[key] = (value, subject, trial)

    log(f"  ID {len(paths)}개 / subject {len(masses)}명에서 요구 토크 봉투를 뽑았다")
    log(f"  (지금 입력은 **보행·보행수정만**이다. squats·STS·DJ 는 Phase 1 이후 추가된다)")
    return out


def peaks_from_summary(csv_path, log=print) -> dict:
    """`summary.csv` 에서. 발목·무릎·고관절 3개뿐이고 방향 구분이 없다."""
    out = {}
    mapping = {"ankle": "ankle_angle", "knee": "knee_angle", "hip": "hip_flexion"}
    with open(csv_path) as fh:
        for row in csv.DictReader(fh):
            for key, coord in mapping.items():
                for side in models.SIDES:
                    raw = row.get(f"{key}_{side}")
                    if not raw:
                        continue
                    value = float(raw)
                    for sign in ("+", "-"):
                        k = (coord, sign)
                        if value > out.get(k, (0.0, None, None))[0]:
                            out[k] = (value, row.get("subject"), row.get("trial"))
    log(f"  [주의] summary.csv 에는 좌표 3개뿐이고 **방향 구분이 없다** — "
        f"굴곡·신전 모두에 |피크|를 넣었다. 6 DOF 가 필요하면 --id-dir 을 쓴다")
    return out


# ─────────────────────────────────────────────────────────────────────────────
# 스크리닝
# ─────────────────────────────────────────────────────────────────────────────


def screen(tau_nm, span_rad, module, bound):
    """(r_min(n=1), r_max, n_min). 힘·변위 조건이 만나는 데 필요한 모듈 개수."""
    f_max, stroke = module.spec(bound)
    r_max = stroke / span_rad
    r_min = tau_nm / f_max
    n_min = max(1, math.ceil(tau_nm * span_rad / (f_max * stroke) - 1e-12))
    return r_min, r_max, n_min


def human_check(tau_nm, span_rad, arm_m, module, bound):
    """인체 부착점(모멘트암 고정)으로 갈 때 필요한 변위·힘."""
    f_max, stroke = module.spec(bound)
    travel = arm_m * span_rad                 # ROM 전체를 덮는 데 필요한 이동량
    force = tau_nm / arm_m                    # 그 모멘트암으로 τ 를 내는 데 필요한 힘
    return travel, force, travel <= stroke, max(1, math.ceil(force / f_max - 1e-12))


def smallest_module(tau_nm, span_rad, arm_m, bound):
    """인체 모멘트암 기준으로 변위·힘을 **하나로** 만족하는 가장 작은 등급."""
    for module in am.MODULES:
        travel, force, fits, n = human_check(tau_nm, span_rad, arm_m, module, bound)
        if fits and n == 1:
            return module.name
    return None


# ─────────────────────────────────────────────────────────────────────────────
# 보고
# ─────────────────────────────────────────────────────────────────────────────


def render(peaks, out_dir, bound="low", log=print):
    lines = []

    def emit(text=""):
        log(text)
        lines.append(text)

    f_note = "하한(보수적)" if bound == "low" else "상한"
    emit(f"\n## 요구 토크 봉투 (보행·보행수정 {len(peaks)}개 좌표·방향)\n")
    emit(f"  {'좌표':<16s} {'방향':<12s} {'Nm/kg':>8s} "
         + " ".join(f"{f'@{m:.0f}kg':>9s}" for m in ROBOT_MASSES_KG) + "   최대가 나온 곳")
    rows = []
    for (coord, sign), (value, subject, trial) in sorted(peaks.items()):
        if coord not in models.ROM_TARGET_DEG:
            continue
        label = models.DIRECTIONS[coord][sign]
        emit(f"  {coord:<16s} {label:<12s} {value:8.3f} "
             + " ".join(f"{value*m:9.1f}" for m in ROBOT_MASSES_KG)
             + f"   {subject}/{trial}")
        rows.append((coord, sign, label, value))

    emit(f"\n## 인체 부착점 그대로 가능한가  (모듈 {f_note} 규격, 로봇 70 kg)\n")
    emit("  ROM 전체를 덮는 데 필요한 이동량 Δl = r·ΔΘ, 그 모멘트암으로 τ 를 내는 데")
    emit("  필요한 힘 F = τ/r. 둘 다 만족하는 가장 작은 모듈 등급을 찾는다.\n")
    emit(f"  {'좌표':<16s} {'방향':<10s} {'대표근육':<22s} {'r(cm)':>6s} "
         f"{'Δl(mm)':>8s} {'F(N)':>8s} {'등급':>6s}  판정")
    verdicts = []
    for coord, sign, label, value in rows:
        rep = models.REPRESENTATIVE_ARM_M.get((coord, sign))
        span = math.radians(models.rom_span_deg(coord))
        tau = value * 70.0
        if rep is None:
            emit(f"  {coord:<16s} {label:<10s} {'모름':<22s} {'-':>6s} {'-':>8s} "
                 f"{'-':>8s} {'-':>6s}  대표 모멘트암 미측정 → Phase 2")
            continue
        name, arm = rep
        grade = smallest_module(tau, span, arm, bound)
        travel = arm * span
        force = tau / arm
        if grade:
            verdict, ok = f"{grade} 1개로 가능", True
        else:
            best = min(am.MODULES, key=lambda m: human_check(tau, span, arm, m, bound)[3])
            t, f, fits, n = human_check(tau, span, arm, best, bound)
            verdict = (f"불가 — 변위 {travel*1000:.0f} mm 가 최대 스트로크 "
                       f"{max(m.spec(bound)[1] for m in am.MODULES)*1000:.0f} mm 초과"
                       if not any(human_check(tau, span, arm, m, bound)[2] for m in am.MODULES)
                       else f"{best.name} {n}개 병렬 필요")
            ok = False
        emit(f"  {coord:<16s} {label:<10s} {name:<22s} {arm*100:6.1f} "
             f"{travel*1000:8.0f} {force:8.0f} {grade or '—':>6s}  {verdict}")
        verdicts.append((coord, sign, ok, verdict))

    emit(f"\n## 모듈 등급별 필요 개수 하한  (일-용량 조건, 로봇 70 kg, {f_note} 규격)\n")
    emit(f"  {'좌표':<16s} {'방향':<10s} {'ΔΘ(도)':>7s} {'r_max(cm)':>10s} "
         + " ".join(f"{m.name:>4s}" for m in am.MODULES) + "   ← 필요 모듈 수")
    for coord, sign, label, value in rows:
        span_deg = models.rom_span_deg(coord)
        span = math.radians(span_deg)
        tau = value * 70.0
        counts, r_max_show = [], None
        for module in am.MODULES:
            _r_min, r_max, n_min = screen(tau, span, module, bound)
            counts.append(n_min)
            if module.name == "M1":
                r_max_show = r_max
        emit(f"  {coord:<16s} {label:<10s} {span_deg:7.0f} {r_max_show*100:10.1f} "
             + " ".join(f"{n:4d}" for n in counts))
    emit("\n  r_max 는 M1 기준(스트로크가 가장 길다). 등급마다 다르다 — CSV 참고.")

    emit(f"\n## 이론적 하한\n")
    emit(f"  당김 전용 액추에이터로 6 DOF 토크 공간을 양의 결합으로 덮으려면")
    emit(f"  **다리당 최소 {am.positive_spanning_floor(6)}개, 양다리 {am.positive_spanning_floor(6)*2}개**가 필요하다.")
    emit(f"  위 표의 '필요 개수'는 방향별 힘·변위 조건일 뿐이고, 이 하한은 그것과 별개로 걸린다.")

    os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "screening.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["coord", "sign", "direction", "tau_nm_per_kg", "span_deg",
                    "module", "bound", "f_max_n", "stroke_mm", "r_min_cm", "r_max_cm",
                    "n_min", "arm_human_cm", "travel_mm", "force_n"])
        for coord, sign, label, value in rows:
            span = math.radians(models.rom_span_deg(coord))
            rep = models.REPRESENTATIVE_ARM_M.get((coord, sign))
            for module in am.MODULES:
                for b in ("low", "high"):
                    tau = value * 70.0
                    r_min, r_max, n_min = screen(tau, span, module, b)
                    f_max, stroke = module.spec(b)
                    row = [coord, sign, label, f"{value:.4f}", f"{models.rom_span_deg(coord):.0f}",
                           module.name, b, f_max, stroke * 1000,
                           f"{r_min*100:.2f}", f"{r_max*100:.2f}", n_min]
                    if rep:
                        travel, force, _fits, _n = human_check(tau, span, rep[1], module, b)
                        row += [f"{rep[1]*100:.1f}", f"{travel*1000:.0f}", f"{force:.0f}"]
                    else:
                        row += ["", "", ""]
                    w.writerow(row)
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser(description="Phase 0 — 모듈 개수 하한 스크리닝")
    ap.add_argument("--id-dir", default=None, help="배치 결과 폴더 (6 DOF 전부 얻는다)")
    ap.add_argument("--summary", default=None, help="summary.csv (3 좌표만)")
    ap.add_argument("--data", default=None, help="subject 모델이 있는 폴더 (--id-dir 과 함께)")
    ap.add_argument("--out", default="out/lowerlimb/phase0")
    ap.add_argument("--exclude-subject", action="append", default=[],
                    help="봉투에서 뺄 subject. 잔차가 높은 피험자의 영향을 보는 데 쓴다")
    ap.add_argument("--bound", default="both", choices=("low", "high", "both"))
    args = ap.parse_args()

    if args.id_dir:
        if not args.data:
            raise SystemExit("--id-dir 을 쓰려면 --data 도 필요하다 (체중을 알아야 Nm/kg 이 된다)")
        peaks = peaks_from_id_dir(args.id_dir, args.data, exclude=tuple(args.exclude_subject))
    elif args.summary:
        peaks = peaks_from_summary(args.summary)
    else:
        raise SystemExit("--id-dir 또는 --summary 중 하나가 필요하다")

    os.makedirs(args.out, exist_ok=True)
    bounds = ("low", "high") if args.bound == "both" else (args.bound,)
    chunks = []
    for bound in bounds:
        print(f"\n{'='*78}\n모듈 규격 {bound.upper()} 쪽\n{'='*78}")
        chunks.append(render(peaks, args.out, bound))

    header = (f"# Phase 0 — 모듈 개수 하한 스크리닝\n\n"
              f"입력: 보행·보행수정 ID 결과. **squats·STS·DJ 미포함** (Phase 1 이후 재실행).\n\n"
              f"파라미터: `ε_max={am.EPS_MAX}`, `ε_p={am.EPS_PASSIVE}`, "
              f"`ρ_F={am.RHO_F_TARGET} N/g`, 수동힘 상한 `{am.PASSIVE_CAP_RATIO}·F_max`\n\n"
              f"목표 ROM 은 AAOS 기준 임시값이며 **뼈대 팀과 확정 전**이다.\n")
    body = "\n".join(f"\n---\n\n# 모듈 규격 {b.upper()} 쪽\n\n```\n{c}\n```\n"
                     for b, c in zip(bounds, chunks))
    with open(os.path.join(args.out, "summary.md"), "w") as fh:
        fh.write(header + body)
    print(f"\n요약 : {os.path.join(args.out, 'summary.md')}")


if __name__ == "__main__":
    main()
