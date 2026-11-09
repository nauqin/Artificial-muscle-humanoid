"""모듈 먼저 실험 2차 (2026-11-05): 시드 반복 + 보정한 모듈 먼저 + 검증 동작 추가를 한 번에.

방식 4개 × 부위 4개 × 시드:
  cur   = 현재 방식 (TUTORIAL 9-6 조건)                      시드 0 은 1차에서 함, 여기서 1·2
  curtr = 현재 방식 + --refine-trust 0.01                     시드 1·2
  new   = 모듈 먼저 (modules_<부위>.json, --refine-trust 0.01) 시드 1·2
  cal   = 보정한 모듈 먼저 (modules_<부위>_cal.json: W × 현재+신뢰 구간 시드 0 의 실제/예측 비)  시드 0·1·2
단계: 경로 설계 → 후처리(어깨 힘 ×1.5, 손 하중 모델 3.79·0.52_bi, stroke_check, 예측 비교) → Moco.
Moco: 모든 설계에 기본 두 동작(하체 걷기·스쿼트, 상체 걷기·찬 물병 들기). 시드 0 설계에는 추가 동작(하체 LabValidation 걷기, 상체 1.6 m/s 걷기·상자 들기).
이미 있는 결과(json, csv, results_paths.csv 의 이름)는 건너뛴다 → 다시 돌려도 이어서 한다.

    python report/modulefirst/batch2.py [--jobs 12] [--moco-jobs 8] [--only design|post|moco]
"""
import argparse, csv, os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); MF = os.path.join(ROOT, "report/modulefirst"); LOG = os.path.join(MF, "logs")
PY = sys.executable; ENV = dict(os.environ, PYTHONUTF8="1", MKL_THREADING_LAYER="SEQUENTIAL")
M = "data/robot40_h160/robot40_h160.osim"; T = "report/upper/tasks_40kg_h160.csv"; B = "0.9877"
REG = {   # 부위: (space, k, 기본 꼬리표, 설계 npz, path_design 공통 인자, 현재 방식 상한 인자, 모듈 json 이름, 힘 배율)
    "lower":    ("lower", 8, "_40kg_h160", "report/newmus_V_design_lower_k8_40kg_h160.npz", ["--model", M, "--box-scale", B, "--ref", "40", "1.60"], ["--stroke-max", "0.12"], "modules_lower_k8_40kg_h160", 1.0),
    "trunk":    ("trunk", 4, "_40kg_h160", "report/newmus_V_design_trunk_k4_40kg_h160.npz", ["--tasks", T, "--model", M, "--box-scale", B], ["--f-module", "1800", "--stroke-max", "0.08"], "modules_trunk_k4_40kg_h160", 1.0),
    "elbow":    ("elbow_r", 2, "_40kg_h160", "report/newmus_V_design_elbow_r_k2_40kg_h160.npz", ["--tasks", T, "--model", M, "--box-scale", B], ["--stroke-max", "0.08"], "modules_elbow_r_k2_40kg_h160", 1.0),
    "shoulder": ("shoulder_r", 6, "_40kg_h160_v635", "report/newmus_V_design_shoulder_r_k6_40kg_h160_v635.npz", ["--tasks", T, "--model", M, "--box-scale", B], ["--stroke-max", "0.08"], "modules_shoulder_r_k6_40kg_h160_v635_f15", 1.5),
}
TRUST = ["--refine-trust", "0.01"]
SEEDS = {"cur": (0, 1, 2), "curtr": (0, 1, 2), "new": (0, 1, 2), "cal": (0, 1, 2)}
TASKS = {"lower": ("walk", "squat"), "upper": ("walk", "lift")}; EXTRA = {"lower": ("walklab",), "upper": ("walk2", "lift2")}
PREFILTER = {"walk", "walk2"}                                  # 걷기만 --prefilter (맥과 같은 조건). LabValidation 걷기는 안 붙인 것이 맥 줄과 같다


def tag_of(reg, meth, seed):
    space, k, base, *_ = REG[reg]
    m = {"cur": "_s120_ref" if reg == "lower" else "_s80_win", "curtr": ("_s120_ref" if reg == "lower" else "_s80_win") + "_tr", "new": "_mf", "cal": "_cal"}[meth]
    return base + m + (f"_sd{seed}" if seed else "")


def json_of(reg, meth, seed, final=True):
    space, k = REG[reg][:2]; f15 = "_f15" if (final and REG[reg][7] != 1.0) else ""
    return f"report/newmus_paths_{space}_k{k}{tag_of(reg, meth, seed)}{f15}.json"


def run(cmd, log):
    with open(os.path.join(LOG, log), "w", encoding="utf-8") as fh:
        return subprocess.run(cmd, cwd=ROOT, env=ENV, stdout=fh, stderr=subprocess.STDOUT).returncode


def design_jobs():
    jobs = []
    for reg, (space, k, base, npz, common, curcap, mods, s) in REG.items():
        for meth, seeds in SEEDS.items():
            for sd in seeds:
                if os.path.exists(os.path.join(ROOT, json_of(reg, meth, sd, final=False))): continue
                args = ["--space", space, "--design", npz] + common + ["--out-tag", tag_of(reg, meth, sd)] + (["--seed", str(sd)] if sd else [])
                if meth == "cur": args += curcap
                elif meth == "curtr": args += curcap + TRUST
                elif meth == "new": args += ["--modules", f"report/modulefirst/{mods}.json"] + TRUST
                else: args += ["--modules", f"report/modulefirst/{mods}_cal.json"] + TRUST
                jobs.append((reg, meth, sd, [PY, "-u", "report/path_design.py"] + args))
    return jobs


def post(reg, meth, sd):
    space, k, base, npz, common, curcap, mods, s = REG[reg]; j0 = json_of(reg, meth, sd, final=False); jf = json_of(reg, meth, sd)
    nm = os.path.basename(jf).replace("newmus_paths_", "").replace(".json", ""); lg = f"post_{reg}_{meth}_sd{sd}.log"
    if space != "lower" or s != 1.0:
        cmd = [PY, "report/modulefirst/prepare_models.py", space, j0, "--loads", "3.79,0.52_bi"] + (["--force-scale", str(s)] if s != 1.0 else [])
        run(cmd, lg)
    sc = f"report/modulefirst/stroke_check_{nm}.csv"
    if not os.path.exists(os.path.join(ROOT, sc)):
        run([PY, "report/stroke_check.py", "--space", space, "--paths", jf, "--model", M, "--out", sc] + (["--tasks", T, "--box-scale", B] if space != "lower" else []), "stroke_" + lg)
    mj = f"report/modulefirst/{mods}{'_cal' if meth == 'cal' else ''}.json"
    run([PY, "report/module_first.py", "--space", space, "--design", "x", "--out", mj, "--check", jf, sc, "--check-out", f"report/modulefirst/check_{nm}.csv"], "check_" + lg)


def moco_jobs(done, extra_all=False):
    jobs = []
    for reg, (space, k, *_r) in REG.items():
        grp = "lower" if reg == "lower" else "upper"
        for meth, seeds in SEEDS.items():
            for sd in seeds:
                tasks = TASKS[grp] + (EXTRA[grp] if (sd == 0 or extra_all) else ())
                jf = json_of(reg, meth, sd); stem = os.path.basename(jf).replace(".json", "")
                for task in tasks:
                    name = f"mf_{reg}_{meth}{f'_sd{sd}' if sd else ''}_{task}"
                    if name in done: continue
                    model = f"report/moco/{stem}" + {"lift": "_load3.79", "lift2": "_load0.52_bi"}.get(task, "") + "_model.osim"
                    cmd = ["bash", "report/modulefirst/run_moco.sh", space, jf, model, task, name] + (["--prefilter"] if task in PREFILTER else [])
                    jobs.append((name, cmd))
    return jobs


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--jobs", type=int, default=12); ap.add_argument("--moco-jobs", type=int, default=8); ap.add_argument("--only", default=None); ap.add_argument("--extra-all", action="store_true", help="추가 검증 동작을 모든 시드에")
    a = ap.parse_args(); tic = time.time()
    if a.only in (None, "design"):
        jobs = design_jobs(); print(f"[설계] {len(jobs)}건", flush=True)
        def f(j):
            reg, meth, sd, cmd = j; rc = run(cmd, f"path2_{reg}_{meth}_sd{sd}.log")
            last = [l for l in open(os.path.join(LOG, f"path2_{reg}_{meth}_sd{sd}.log"), encoding="utf-8") if "2단계 끝" in l]
            print(f"[설계 끝] {reg} {meth} 시드{sd} rc={rc} {last[-1].strip() if last else '2단계 줄 없음'}  ({(time.time()-tic)/60:.0f}분)", flush=True)
        jobs.sort(key=lambda j: 0 if j[0] == "lower" else 1)    # 오래 걸리는 하체 먼저
        with ThreadPoolExecutor(a.jobs) as ex: list(ex.map(f, jobs))
    if a.only in (None, "post"):
        items = [(reg, meth, sd) for reg in REG for meth, seeds in SEEDS.items() for sd in seeds]
        def g(it):
            if os.path.exists(os.path.join(ROOT, json_of(*it, final=False))): post(*it)
        with ThreadPoolExecutor(a.jobs) as ex: list(ex.map(g, items))
        print(f"[후처리 끝] {len(items)}건  ({(time.time()-tic)/60:.0f}분)", flush=True)
    if a.only in (None, "moco"):
        done = {r[0] for r in csv.reader(open(os.path.join(ROOT, "report/moco/results_paths.csv"), encoding="utf-8")) if r}
        jobs = moco_jobs(done, a.extra_all); print(f"[Moco] {len(jobs)}건", flush=True)
        def h(j):
            name, cmd = j; out = subprocess.run(cmd, cwd=ROOT, env=ENV, capture_output=True, text=True, encoding="utf-8", errors="replace").stdout.strip()
            print(f"[Moco 끝] {out or name + ' 결과 줄 없음 (로그 확인)'}", flush=True)
        with ThreadPoolExecutor(a.moco_jobs) as ex: list(ex.map(h, jobs))
    print(f"[모두 끝] {(time.time()-tic)/60:.0f}분", flush=True)


if __name__ == "__main__":
    main()
