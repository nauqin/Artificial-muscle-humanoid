"""어깨 현재 방식 + 신뢰 구간을 시드 3·4·5 로 더 (2026-11-05 2차). 시드 1·2 가 어깨 네 동작을 모두 통과해서, 통과율을 가늠하려고.
batch2.py 와 같은 단계(설계 → 후처리 → Moco 기본 2 + 추가 2)를 어깨·curtr 에만.

    python report/modulefirst/batch2_shoulder_seeds.py
"""
import csv, os, subprocess, sys
from concurrent.futures import ThreadPoolExecutor
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import batch2 as b2                                           # noqa: E402

b2.REG = {"shoulder": b2.REG["shoulder"]}; b2.SEEDS = {"curtr": (3, 4, 5)}
jobs = b2.design_jobs(); print(f"[설계] {len(jobs)}건", flush=True)
with ThreadPoolExecutor(6) as ex:
    list(ex.map(lambda j: print(f"[설계 끝] 시드{j[2]} rc={b2.run(j[3], f'path2_shoulder_curtr_sd{j[2]}.log')}", flush=True), jobs))
for sd in b2.SEEDS["curtr"]: b2.post("shoulder", "curtr", sd)
done = {r[0] for r in csv.reader(open(os.path.join(b2.ROOT, "report/moco/results_paths.csv"), encoding="utf-8")) if r}
mj = b2.moco_jobs(done, extra_all=True); print(f"[Moco] {len(mj)}건", flush=True)
def h(j):
    out = subprocess.run(j[1], cwd=b2.ROOT, env=b2.ENV, capture_output=True, text=True, encoding="utf-8", errors="replace").stdout.strip()
    print(f"[Moco 끝] {out or j[0]}", flush=True)
with ThreadPoolExecutor(8) as ex: list(ex.map(h, mj))
print("[모두 끝]", flush=True)
