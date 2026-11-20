"""Moco 해가 실제로 쓴 (수축, 장력) 이 배정한 모듈 곡선 안에 드는지 (2026-11-19).

장력 T(t) = Moco 신호 u × force_N, 수축 s(t) = 설계 때의 최대 길이 L_max − 그 순간 길이 (같은 시각의 IK 로 잰다).
판정: T(t) ≤ 묶음 × 등급 g 한 묶음이 수축 s 에서 낼 수 있는 힘 (report/module_curves.json, 임시 β).

    python report/td/moco_curve_check.py --space elbow_r --paths report/newmus_paths_elbow_r_k2_40kg_h175sk_td_w0.json \
        --solution report/moco/td175_elbow_td_w0_lift_solution.sto --ik data/carrt/subject1/PFWJ2_ik.mot --model data/robot40_h175sk/robot40_h175sk.osim --box-scale 1.0802
여러 개: --batch report/td/moco_curve_batch.csv (열: space,paths,solution,ik,model,box_scale) → report/td/moco_curve_check.csv
"""
import argparse, csv, json, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import addb_io as io                                         # noqa: E402
import path_design as pd                                     # noqa: E402
import timedomain as td                                      # noqa: E402


def check(space, paths, solution, ik, model, box_scale, beta=None):
    pd.configure(space)
    if box_scale != 1.0 and pd.WRAP: pd.WRAP = dict(pd.WRAP, center=tuple(v * box_scale for v in pd.WRAP["center"]), radius=pd.WRAP["radius"] * box_scale)
    P = json.load(open(os.path.join(ROOT, paths), encoding="utf-8")); sol = io.read_mot(os.path.join(ROOT, solution)); m = io.read_mot(os.path.join(ROOT, ik))
    if m.in_degrees: m.data = np.radians(m.data)
    t = sol.time; cols = list(pd.COORDS) + list(pd.EXTRA)
    Q = np.column_stack([np.interp(t, m.time, m.column(c)) for c in cols]); fr = pd.Frames(os.path.join(ROOT, model), Q, perturb=False)
    curves = td.load_curves(beta=beta); out = []
    for mm in P["muscles"]:
        name = f"/forceset/{mm['name']}_r" if f"/forceset/{mm['name']}_r" in sol.names else f"/forceset/{mm['name']}"
        u = sol.column(name); T = u * mm["force_N"]
        L = pd.path_lengths(fr, [pd.BODIES.index(b) for b in mm["chain"]], np.array(mm["points"]))[:, 0]; s = mm["L_max_m"] - L
        g = (mm.get("module") or {}).get("name"); n = (mm.get("module") or {}).get("bundles")
        if not g: out.append(dict(name=mm["name"], grade=None)); continue
        cap = n * td.f_avail(curves[g], np.maximum(s, 0.0)); ratio = T / np.maximum(cap, 1e-6)
        out.append(dict(name=mm["name"], grade=g, bundles=n, T_max_N=float(T.max()), s_max_mm=float(s.max() * 1000), s_min_mm=float(s.min() * 1000),
                        ratio_max=float(ratio.max()), frac_over=float((ratio > 1.0 + 1e-6).mean()), stretched_beyond_design_mm=float(max(0.0, -s.min()) * 1000)))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--space"); ap.add_argument("--paths"); ap.add_argument("--solution"); ap.add_argument("--ik"); ap.add_argument("--model")
    ap.add_argument("--box-scale", type=float, default=1.0); ap.add_argument("--beta", type=float, default=None); ap.add_argument("--batch", default=None)
    a = ap.parse_args()
    jobs = list(csv.DictReader(open(os.path.join(ROOT, a.batch), encoding="utf-8"))) if a.batch else [dict(space=a.space, paths=a.paths, solution=a.solution, ik=a.ik, model=a.model, box_scale=a.box_scale)]
    rows = []
    for j in jobs:
        if not os.path.exists(os.path.join(ROOT, j["solution"])): print("  없음", j["solution"]); continue
        for r in check(j["space"], j["paths"], j["solution"], j["ik"], j["model"], float(j["box_scale"]), a.beta):
            r.update(solution=os.path.basename(j["solution"])); rows.append(r)
            if r.get("grade"): print(f"  {r['solution']:45s} {r['name']}: {r['grade']} ×{r['bundles']}  장력 최대 {r['T_max_N']:.0f} N  수축 {r['s_min_mm']:.0f}~{r['s_max_mm']:.0f} mm  "
                                   f"곡선 대비 최대 {r['ratio_max']*100:.0f} %  넘은 프레임 {r['frac_over']*100:.1f} %")
    if a.batch and rows:
        out = os.path.join(ROOT, "report/td/moco_curve_check.csv"); keys = sorted({k for r in rows for k in r})
        with open(out, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)
        print("→", os.path.relpath(out, ROOT))


if __name__ == "__main__":
    main()
