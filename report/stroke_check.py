"""케이블 스트로크·힘 → 과제 계획서 모듈 등급 배정 (+ SMA 수축률 참고 검사). (2026-11-04, 교수님 결정 8번 "근육 길이 단축")

모듈은 힘(N) 범위와 스트로크(mm) 범위만 본다 (대응 사람 근육 열은 무시, 11/4 사용자). 케이블마다 길이 변화 ΔL 과 필요 힘 F 를 재서,
스트로크 상한 ≥ ΔL 인 등급 중 묶음 수 n = ceil(F / 힘 상한) 이 가장 작은(같으면 작은 등급) 것을 고른다. ΔL 이 M1 의 120 mm 도 넘으면
어떤 모듈로도 못 만든다 → 경로를 바꿔 ΔL 을 줄여야 한다 (`path_design.py --stroke-max`).

설계 때 쓴 자세 표본(tasks.csv 의 trial 마다 n개 + 설계 좌표별 극값 자세)에서 케이블 길이 L 을 재고, 케이블마다
  L_min·L_max, ΔL = L_max − L_min, ΔL/L_mean(path_design 의 '스트로크' 지표, 상한 STROKE), ΔL/L_max,
  SMA 능동 길이 = ΔL / ε × 여유 (ε = SMA 수축률, 기본 0.5·여유 1.2) 와 나머지 수동 케이블 길이 = L_min − 능동 길이
를 표·CSV 로 낸다. 능동 길이가 L_min 보다 길면(수동 길이 음수) 그 케이블은 SMA 로 못 만든다 → 경로를 바꿔야 한다.

    python report/stroke_check.py --space lower --paths report/newmus_paths_k8.json --model data/subject3/subject3_scaled.osim
    python report/stroke_check.py --space shoulder_r --paths report/newmus_paths_shoulder_r_k5.json --tasks report/upper/tasks.csv \\
        --model data/moore2015/subject3/subject3_scaled.osim
    로봇 크기면 path_design 과 같은 --model/--box-scale/--ref 를 준다.
→ report/stroke_check_<json 이름>.csv
"""
import argparse, csv, json, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import path_design as pd                                     # noqa: E402
from actuator_model import MODULES                           # noqa: E402  과제 계획서 모듈표 (힘 N, 스트로크 mm)


def assign_module(force_N, dL_mm):
    """스트로크가 맞는 등급 중 묶음 수가 최소인 것. 반환 (등급 이름, 묶음 수) 또는 (None, None)."""
    cands = []
    for m in MODULES:
        if dL_mm > m.stroke_mm[1]: continue
        n = max(1, int(np.ceil(force_N / m.force_n[1])))
        cands.append((n, m.force_n[1], m.name))                 # 묶음 수 최소 → 같으면 힘 상한이 작은(가벼운) 등급
    if not cands: return None, None
    n, _, name = min(cands); return name, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--space", default="lower"); ap.add_argument("--paths", required=True, help="경로 json")
    ap.add_argument("--model", default=None, help="뼈대 .osim (기본 json 의 model)"); ap.add_argument("--tasks", default=None)
    ap.add_argument("--box-scale", type=float, default=1.0); ap.add_argument("--ref", type=float, nargs=2, default=None)
    ap.add_argument("--n-per", type=int, default=20, help="trial 당 자세 수 (극값 자세는 따로 항상 들어간다)")
    ap.add_argument("--eps", type=float, default=0.5, help="SMA 수축률 (park2019: 50 % 이상)"); ap.add_argument("--margin", type=float, default=1.2, help="능동 길이 여유 배율")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(); pd.configure(a.space)
    if a.box_scale != 1.0 and pd.WRAP: pd.WRAP = dict(pd.WRAP, center=tuple(v * a.box_scale for v in pd.WRAP["center"]), radius=pd.WRAP["radius"] * a.box_scale)
    pj = a.paths if os.path.isabs(a.paths) else os.path.join(ROOT, a.paths); paths = json.load(open(pj))
    model = a.model or paths.get("model", "data/subject3/subject3_scaled.osim"); model = model if os.path.isabs(model) else os.path.join(ROOT, model)
    if a.tasks: Q, T = pd.poses_from_tasks(a.tasks if os.path.isabs(a.tasks) else os.path.join(ROOT, a.tasks), n_per=a.n_per, ref=a.ref)
    else: Q, T = pd.poses(n_walk=a.n_per, n_sq=3 * a.n_per)
    fr = pd.Frames(model, Q)
    print(f"[{a.space}] {os.path.basename(pj)}  자세 {len(Q)}개  모델 {os.path.relpath(model, ROOT)}  SMA 수축률 {a.eps:.2f} × 여유 {a.margin}")
    print(f"{'근육':5s} {'바디':28s} {'L_min':>6s} {'L_max':>6s} {'ΔL':>5s}  {'힘 N':>6s}  {'모듈':8s} | {'ΔL/L평균':>8s} {'능동':>5s} {'수동':>5s} SMA   (cm)")
    rows = []
    for m in paths["muscles"]:
        chain = [pd.BODIES.index(b) for b in m["chain"]]; L = pd.path_lengths(fr, chain, np.array(m["points"]))[:, 0]
        lmin, lmax, lmean = float(L.min()), float(L.max()), float(L.mean()); dl = lmax - lmin
        act = dl / a.eps * a.margin; passive = lmin - act; ok = passive >= 0; F = float(m.get("force_N") or 0); mod, n = assign_module(F, dl * 1000)
        mod_s = f"{mod} ×{n}" if mod else "없음(ΔL 초과)"
        rows.append(dict(name=m["name"], chain="→".join(b.replace("_r", "") for b in m["chain"]), L_min_cm=lmin * 100, L_max_cm=lmax * 100, dL_cm=dl * 100, force_N=F,
                         module=mod, bundles=n, ratio_mean=dl / lmean, ratio_max=dl / lmax, active_cm=act * 100, passive_cm=passive * 100, sma_ok=ok))
        print(f"{m['name']:5s} {rows[-1]['chain']:28s} {lmin*100:6.1f} {lmax*100:6.1f} {dl*100:5.1f}  {F:6.0f}  {mod_s:8s} | {dl/lmean:8.2f} {act*100:5.1f} {passive*100:5.1f} {'가능' if ok else '불가'}")
    out = a.out or os.path.join(ROOT, "report", f"stroke_check_{os.path.basename(pj).replace('newmus_paths_', '').replace('.json', '')}.csv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    n_none = sum(r["module"] is None for r in rows); n_bad = sum(not r["sma_ok"] for r in rows)
    print(f"→ {out}   모듈 배정 불가 {n_none}개 / {len(rows)}개 (ΔL 최대 {max(r['dL_cm'] for r in rows)*10:.0f} mm)   SMA 참고: 불가 {n_bad}개")


if __name__ == "__main__":
    main()
