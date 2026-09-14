"""subject × trial 전수 처리, `summary.csv` 집계.

trial 하나가 죽어도 나머지는 계속 돈다. 실패 이유는 `error` 열에 남는다.
"""

from __future__ import annotations

import _env  # noqa: F401  환경 확인 — 다른 import 보다 먼저

import argparse
import csv
import os
import traceback

import addb_io as io
import run_muscle_analysis as rma
import run_static_optimization as rso
import run_trial as rt
import synergy as syn
import validate

#: 이보다 짧은 파일은 정적 자세로 보고 걸러낸다
MIN_DURATION = 0.5

COLUMNS = [
    "subject", "trial", "model", "mass_kg", "t0", "t1",
    "ankle_r", "ankle_l", "knee_r", "knee_l", "hip_r", "hip_l",
    "grf_peak", "residual", "cop_dist",
    "so_recon_max", "so_reserve_max",
    "syn_k_r", "syn_k_l", "syn_vaf_r", "syn_vaf_l",
    "syn_kglobal_r", "syn_kglobal_l", "syn_kdual_r", "syn_kdual_l",
    "id_verdict", "ma_verdict", "so_verdict", "syn_verdict", "error",
]


def condition_of(trial_name: str) -> str:
    """`walking3` → `walking`. 조건별로 이어붙일 때 쓴다."""
    return trial_name.rstrip("0123456789") or trial_name


def process_trial(subject, trial, model_path, out_root, args, log=print):
    row = {"subject": subject, "trial": trial.name}
    out_dir = os.path.join(out_root, subject, trial.name)
    id_path = os.path.join(out_dir, f"{trial.name}_id.sto")

    if os.path.exists(id_path) and not args.force:
        log(f"  건너뜀 (결과 있음) — 다시 하려면 --force")
        cached = True
    else:
        cached = False

    summary = None
    if not cached:
        summary = rt.run_trial(model_path, trial.ik, trial.grf, out_dir, trial.name,
                               lowpass=args.lowpass, start=args.start, end=args.end,
                               log=log)
    else:
        # ID 툴은 다시 안 돌리되 요약·검산은 다시 만든다 — 안 그러면 그 줄이 빈다.
        summary = rt.summarize(model_path, trial.ik, trial.grf, id_path, out_dir,
                               trial.name, lowpass=args.lowpass, log=log)
    row.update({k: v for k, v in summary.items() if k in COLUMNS})
    t0, t1 = summary["t0"], summary["t1"]
    ext_path = summary.get("external_loads") or os.path.join(
        out_dir, f"{trial.name}_external_loads.xml")

    coords = tuple(c.strip() for c in args.coords.split(",") if c.strip())

    ma_dir = None
    if not args.no_muscle:
        ma_dir = os.path.join(out_dir, "ma")
        marker = os.path.join(ma_dir, f"{trial.name}_MuscleAnalysis_MomentArm_{coords[0]}.sto")
        if not os.path.exists(marker) or args.force:
            rma.run_muscle_analysis(model_path, trial.ik, ma_dir, trial.name,
                                    t0=t0, t1=t1, lowpass=args.lowpass, log=log)
        ma = rma.load_moment_arms(ma_dir, trial.name, coords)
        ma_rep = rma.validate_moment_arms(ma, coords)
        log(ma_rep.render())
        row["ma_verdict"] = ma_rep.verdict

    if args.so or args.synergy:
        so_dir = os.path.join(out_dir, "so")
        force_path = os.path.join(so_dir, f"{trial.name}_StaticOptimization_force.sto")
        act_path = os.path.join(so_dir, f"{trial.name}_StaticOptimization_activation.sto")
        if not os.path.exists(force_path) or args.force:
            rso.run_static_optimization(model_path, trial.ik, ext_path, so_dir,
                                        trial.name, t0, t1, lowpass=args.lowpass, log=log)
        so_rep = rso.validate_run(id_path, force_path, ma_dir, trial.name, coords, log=log)
        log(so_rep.render())
        row["so_verdict"] = so_rep.verdict
        recon = validate.max_of(so_rep, "recon/")
        reserve = validate.max_of(so_rep, "reserve/")
        row["so_recon_max"] = round(recon, 3) if recon is not None else None
        row["so_reserve_max"] = round(reserve, 3) if reserve is not None else None

        if args.synergy:
            verdicts = []
            for side in _sides(args.side):
                res = syn.analyze([act_path], so_dir, trial.name, side,
                                  kmax=args.kmax, criterion=args.criterion,
                                  seed=args.seed, model_path=model_path, log=log)
                row[f"syn_k_{side}"] = res["k"]
                row[f"syn_vaf_{side}"] = res["vaf"]
                row[f"syn_kglobal_{side}"] = res["k_global"]
                row[f"syn_kdual_{side}"] = res["k_dual"]
                verdicts.append(res["verdict"])
            row["syn_verdict"] = validate.worst(verdicts)
        row["_activation"] = act_path
    return row


def _sides(side):
    return ("r", "l") if side == "both" else (side,)


def run_concat(subject, rows, out_root, model_path, args, log=print):
    """같은 조건의 trial 활성도를 이어붙여 시너지를 한 번 더 분해."""
    groups = {}
    for row in rows:
        act = row.get("_activation")
        if act and os.path.exists(act):
            groups.setdefault(condition_of(row["trial"]), []).append(act)

    for cond, acts in sorted(groups.items()):
        if len(acts) < 2:
            continue
        out_dir = os.path.join(out_root, subject, f"_concat_{cond}")
        log(f"\n[{subject}] 이어붙임 '{cond}' — trial {len(acts)}개")
        for side in _sides(args.side):
            try:
                syn.analyze(sorted(acts), out_dir, cond, side, kmax=args.kmax,
                            criterion=args.criterion, seed=args.seed,
                            model_path=model_path, log=log)
            except Exception as exc:
                log(f"  [ERROR] {cond} {side}: {exc}")


def main():
    ap = argparse.ArgumentParser(description="subject × trial 전수 처리")
    ap.add_argument("--data", default="data", help="subject 폴더들의 상위 폴더")
    ap.add_argument("--out", default="out", help="결과를 쌓을 곳")
    ap.add_argument("--subject", action="append", default=None, help="특정 subject만")
    ap.add_argument("--model", default=None, help="모델을 강제 지정")
    ap.add_argument("--force", action="store_true", help="기존 결과도 다시 계산")
    ap.add_argument("--no-muscle", action="store_true", help="MuscleAnalysis 생략, ID만")
    ap.add_argument("--so", action="store_true", help="Static Optimization까지")
    ap.add_argument("--synergy", action="store_true", help="NMF 시너지 (--so 포함)")
    ap.add_argument("--concat", action="store_true", help="같은 조건 trial을 이어붙여 재분해")
    ap.add_argument("--side", default="both", choices=("r", "l", "both"))
    ap.add_argument("--kmax", type=int, default=syn.K_MAX)
    ap.add_argument("--criterion", default="dual", choices=("dual", "global"))
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--coords", default=",".join(rma.DEFAULT_COORDS))
    ap.add_argument("--start", type=float, default=None)
    ap.add_argument("--end", type=float, default=None)
    ap.add_argument("--lowpass", type=float, default=rt.LOWPASS)
    ap.add_argument("--min-duration", type=float, default=MIN_DURATION)
    args = ap.parse_args()

    import opensim as osim
    osim.Logger.setLevelString("Error")

    if args.synergy:
        args.so = True
    if args.so and args.no_muscle:
        raise SystemExit("--so 에는 모멘트암이 필요하다. --no-muscle 과 함께 쓸 수 없다.")

    subjects = io.find_subjects(args.data)
    if args.subject:
        wanted = set(args.subject)
        subjects = [s for s in subjects if os.path.basename(s) in wanted]
    if not subjects:
        raise SystemExit(f"subject 폴더를 찾지 못했다: {args.data}")

    os.makedirs(args.out, exist_ok=True)
    all_rows = []

    for subject_dir in subjects:
        subject = os.path.basename(subject_dir)
        print(f"\n=== {subject} ===")
        model_path = args.model or io.pick_model(subject_dir)
        if not model_path:
            print(f"  [WARN] 모델이 없다 — 건너뜀")
            continue

        trials = io.find_trials(subject_dir, min_duration=args.min_duration)
        print(f"  trial: {len(trials)}개  {[t.name for t in trials]}")

        rows = []
        for trial in trials:
            print(f"\n[{subject}/{trial.name}]")
            try:
                rows.append(process_trial(subject, trial, model_path, args.out, args))
            except Exception as exc:
                traceback.print_exc()
                rows.append({"subject": subject, "trial": trial.name,
                             "model": os.path.basename(model_path), "error": str(exc)})

        if args.concat and args.synergy:
            run_concat(subject, rows, args.out, model_path, args)
        all_rows.extend(rows)

    csv_path = os.path.join(args.out, "summary.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)

    print(f"\n요약 : {csv_path}  ({len(all_rows)} trial)")
    for stage in ("id_verdict", "ma_verdict", "so_verdict", "syn_verdict"):
        counts = {}
        for row in all_rows:
            v = row.get(stage)
            if v:
                counts[v] = counts.get(v, 0) + 1
        if counts:
            print(f"  {stage:12s} " + "  ".join(f"{k} {v}" for k, v in sorted(counts.items())))
    failed = [r for r in all_rows if r.get("error")]
    if failed:
        print(f"  실패 {len(failed)}건: " + ", ".join(f"{r['subject']}/{r['trial']}" for r in failed))


if __name__ == "__main__":
    main()
