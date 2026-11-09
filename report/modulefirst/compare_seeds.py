"""모듈 먼저 실험 2차 비교표 (2026-11-05): 방식 4개 × 시드 3개 + 추가 검증 동작.

칸마다 "시드0 / 시드1 / 시드2 …"(어깨 현재 + 신뢰 구간은 시드 5 까지). 숫자는 결과 파일에서 읽는다:
  경로 json(2단계 잔차) · stroke_check csv(사후 등급) · 모듈 json(사전 등급) · results_paths.csv(Moco 보조 95퍼센타일, 좌표 중 최악) ·
  Moco 해 파일(어깨는 오른팔 = 설계 쪽만의 최악도 따로).
추정 질량 = Σ 묶음 수 × 등급 힘 상한 / 15 N/g (actuator_model.RHO_F_TARGET, 냉각·구조 제외). 배정 불가 케이블은 빼고 따로 센다.

    python report/modulefirst/compare_seeds.py   → 화면 마크다운, report/modulefirst/compare_seeds.csv
"""
import csv, json, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report")); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from actuator_model import MODULES, RHO_F_TARGET          # noqa: E402
import batch2 as b2                                         # noqa: E402
FCAP = {m.name: m.force_n[1] for m in MODULES}
METH = {"cur": "현재 방식", "curtr": "현재 + 신뢰 구간", "new": "모듈 먼저", "cal": "보정한 모듈 먼저"}
RNAME = {"lower": "하체", "trunk": "몸통", "elbow": "팔꿈치", "shoulder": "어깨"}
WINDOW = {"walk": (1, 4), "squat": (2.6, 4.4), "lift": (0.3, 5.5)}
IDS = {"walk": "out/upper40_h160/moore2015/subject3/T047_segment_12/T047_segment_12_id.sto", "lift": "out/upper40_h160/carrt/subject1/PFWJ2/PFWJ2_id.sto",
       "walk2": "out/upper40_h160/moore2015/subject3/T048_segment_5/T048_segment_5_id.sto", "lift2": "out/upper40_h160/carrt/subject1/PB1/PB1_id.sto"}
WINDOW.update(walk2=(1, 4), lift2=(0.3, 4.8))


def moco_rows():
    return {r[0]: r for r in csv.reader(open(os.path.join(ROOT, "report/moco/results_paths.csv"), encoding="utf-8")) if r}


def right_arm_r95(name, task):
    """어깨 Moco 해에서 오른팔(설계 쪽) 세 좌표의 보조 95퍼센타일 중 최악 (%)."""
    import opensim as osim, addb_io as io
    p = os.path.join(ROOT, "report/moco", f"{name}_solution.sto")
    if not os.path.exists(p): return None
    tr = osim.MocoTrajectory(p); tt = np.asarray(tr.getTimeMat()).ravel(); t0, t1 = WINDOW[task]; inner = (tt >= t0 + 0.05) & (tt <= t1 - 0.05)
    idm = io.read_mot(os.path.join(ROOT, IDS[task])); names = list(tr.getControlNames()); worst = 0.0
    for c in ("arm_flex_r", "arm_add_r", "arm_rot_r"):
        rn = [n for n in names if "reserve" in n and n.endswith("_" + c)][0]; r = np.abs(np.asarray(tr.getControlMat(rn)).ravel()[inner])
        worst = max(worst, float(np.percentile(r, 95)) / float(np.abs(idm.column(f"{c}_moment")).max()) * 100)
    return worst


def one(reg, meth, sd, mo):
    space, k, base, npz, common, curcap, mods, s = b2.REG[reg]; jf = os.path.join(ROOT, b2.json_of(reg, meth, sd)); j0 = os.path.join(ROOT, b2.json_of(reg, meth, sd, final=False))
    nm = os.path.basename(jf).replace("newmus_paths_", "").replace(".json", ""); sc = os.path.join(ROOT, f"report/modulefirst/stroke_check_{nm}.csv")
    if not (os.path.exists(j0) and os.path.exists(sc)): return None
    P = json.load(open(j0, encoding="utf-8")); S = list(csv.DictReader(open(sc, encoding="utf-8")))
    out = dict(region=reg, method=meth, seed=sd, worst=P["refine"]["worst_pct"], mean=P["refine"]["mean_pct"])
    for g in ("M1", "M2", "M3", "M4", "M6"): out[g] = sum(r["module"] == g for r in S)
    out["none"] = sum(not r["module"] for r in S)
    out["mass_g"] = sum(int(r["bundles"]) * FCAP[r["module"]] / RHO_F_TARGET for r in S if r["module"])
    out["viol"] = None
    if meth in ("new", "cal"):
        Mo = json.load(open(os.path.join(ROOT, f"report/modulefirst/{mods}{'_cal' if meth == 'cal' else ''}.json"), encoding="utf-8"))["muscles"]
        out["viol"] = sum(float(r["force_N"]) > m["force_cap_N"] + 1e-6 or float(r["dL_cm"]) / 100 > m["stroke_cap_m"] + 1e-9 for r, m in zip(S, Mo))
    grp = "lower" if reg == "lower" else "upper"
    for task in b2.TASKS[grp] + b2.EXTRA[grp]:
        name = f"mf_{reg}_{meth}{f'_sd{sd}' if sd else ''}_{task}"; row = mo.get(name)
        out[task] = float(row[3]) if row else None
        if reg == "shoulder" and row and task in ("walk", "lift", "walk2", "lift2"): out[task + "_right"] = right_arm_r95(name, task)
    return out


def fmt(vals, d=1):
    return " / ".join("—" if v is None else (f"{v:.{d}f}" if isinstance(v, float) else str(v)) for v in vals)


def main():
    mo = moco_rows(); rows = []
    for reg in b2.REG:
        for meth in METH:
            for sd in range(6):                                 # 어깨 현재 + 신뢰 구간만 시드 3~5 가 더 있다
                r = one(reg, meth, sd, mo)
                if r: rows.append(r)
    with open(os.path.join(ROOT, "report/modulefirst/compare_seeds.csv"), "w", newline="", encoding="utf-8") as fh:
        keys = sorted({k for r in rows for k in r}, key=lambda x: (x not in ("region", "method", "seed"), x))
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
    for reg in b2.REG:
        grp = "lower" if reg == "lower" else "upper"; t1, t2 = b2.TASKS[grp]
        print(f"\n### {RNAME[reg]} (칸 = 시드 0 / 1 / 2 …, 있는 것만)\n")
        print(f"| 방식 | M1 수 | 배정 불가 | 사전 등급 못 지킴 | 추정 질량 g | 2단계 잔차 최악 % | {t1} % | {t2} % |" + (" 걷기 오른팔 % | 들기 오른팔 % |" if reg == "shoulder" else ""))
        print("|---|---|---|---|---|---|---|---|" + ("---|---|" if reg == "shoulder" else ""))
        for meth, lab in METH.items():
            R = sorted([r for r in rows if r["region"] == reg and r["method"] == meth], key=lambda r: r["seed"])
            if not R: continue
            line = f"| {lab} | {fmt([r['M1'] for r in R])} | {fmt([r['none'] for r in R])} | {fmt([r['viol'] for r in R])} | {fmt([r['mass_g'] for r in R], 0)} | {fmt([r['worst'] for r in R])} | {fmt([r.get(t1) for r in R])} | {fmt([r.get(t2) for r in R])} |"
            if reg == "shoulder": line += f" {fmt([r.get('walk_right') for r in R])} | {fmt([r.get('lift_right') for r in R])} |"
            print(line)
    print("\n### 추가 검증 동작 (칸 = 시드 0 / 1 / 2 …, 있는 것만)\n")
    print("| 부위 | 방식 | 하체 LabValidation 걷기 · 상체 1.6 m/s 걷기 % | 상체 상자 들기 % | 어깨 오른팔만: 1.6 m/s 걷기 % | 어깨 오른팔만: 상자 % |\n|---|---|---|---|---|---|")
    for reg in b2.REG:
        for meth, lab in METH.items():
            R = sorted([r for r in rows if r["region"] == reg and r["method"] == meth], key=lambda r: r["seed"])
            if not R: continue
            a1 = [r.get("walklab") if reg == "lower" else r.get("walk2") for r in R]; a2 = [r.get("lift2") for r in R]
            sh = (fmt([r.get("walk2_right") for r in R]), fmt([r.get("lift2_right") for r in R])) if reg == "shoulder" else ("", "")
            print(f"| {RNAME[reg]} | {lab} | {fmt(a1)} | {fmt(a2) if reg != 'lower' else ''} | {sh[0]} | {sh[1]} |")
    print("\n→ report/modulefirst/compare_seeds.csv")


if __name__ == "__main__":
    main()
