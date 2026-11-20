"""경로 json 을 시간 영역으로 평가한다 (2026-11-19): 프레임마다 실제 모멘트암으로 장력 T(t) 를 나누고, 같은 순간의 수축 s(t) 와 함께
모듈 등급의 힘–수축 곡선(report/module_curves.json, 임시)에 넣어 묶음 수를 정한다. 경로는 바꾸지 않는다.

지금 설계의 힘(force_N = 1.1 × 토크 용량 ÷ 최소 모멘트암, `amin`)과 나란히 보여 준다.

    python report/td_eval.py --space elbow_r --paths report/newmus_paths_elbow_r_k2_40kg_h160_s80.json \
        --tasks report/upper/tasks_40kg_h160.csv --model data/robot40_h160/robot40_h160.osim --box-scale 0.9877 --grades M3,M4,M5
→ report/td/td_<json 이름><tag>.csv (케이블 × β), _timeseries.npz, fig_td_<…>.png
"""
import argparse, csv, json, os, sys, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import path_design as pd                                     # noqa: E402
import timedomain as td                                      # noqa: E402
from stroke_check import assign_module                       # noqa: E402
from actuator_model import BY_NAME                           # noqa: E402
import matplotlib.pyplot as plt                              # noqa: E402  (path_design 이 Agg·글꼴을 이미 정했다)

ALL = ["M1", "M2", "M3", "M4", "M5", "M6"]


def evaluate(paths, fr, T, curves, grades, norm="force", max_bundles=None):
    """경로 json + 프레임 → dict(R, Ls, T_alloc, s, res, sizes)."""
    R, Ls = [], []
    for m in paths["muscles"]:
        r, L = pd.moment_arms(fr, [pd.BODIES.index(b) for b in m["chain"]], np.array(m["points"])); R.append(r); Ls.append(L)
    R = np.stack(R, axis=2); Ls = np.stack(Ls, axis=1)
    C = np.array([m["force_N"] for m in paths["muscles"]]) if norm == "force" else np.ones(R.shape[2]) * 1000.0
    Ta, res, _ = td.allocate(R, pd.MARGIN * T, C, iters=800); s = td.contraction(Ls)
    sizes = [td.size_cable(Ta[:, i], s[:, i], curves, grades, max_bundles) for i in range(R.shape[2])]
    return dict(R=R, Ls=Ls, T=Ta, s=s, res=res, sizes=sizes)


def plot_ts(paths, ev, curves, grades, meta, title, png):
    k = len(paths["muscles"]); ncol = min(4, k); nrow = int(np.ceil(k / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 3.6 * nrow), squeeze=False)
    for i, (ax, m) in enumerate(zip(axes.ravel(), paths["muscles"])):
        z = ev["sizes"][i]; s = ev["s"][:, i] * 1000; Tt = ev["T"][:, i]
        ax.scatter(s, Tt, s=3, c=meta["trial"], cmap="tab20", alpha=0.6, lw=0)
        if z["feasible"]:
            cv = curves[z["grade"]]; ss = np.linspace(0, cv["S"], 200); ax.plot(ss * 1000, z["n"] * td.f_avail(cv, ss), "k-", lw=1.2, label=f"{z['grade']} ×{z['n']} 곡선")
        ax.axhline(m["force_N"], color="r", ls="--", lw=1, label=f"json 설계 힘 (force_N) {m['force_N']:.0f} N")
        ax.set_title(f"{m['name']}: 장력 최대 {z['T_peak']:.0f} N, 수축 최대 {z['s_max']*1000:.0f} mm", fontsize=9)
        ax.set_xlabel("수축 s (mm)"); ax.set_ylabel("장력 T (N)"); ax.legend(fontsize=7)
    for ax in axes.ravel()[k:]: ax.axis("off")
    fig.suptitle(title, fontsize=10); fig.tight_layout(); fig.savefig(png, dpi=140); plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--space", required=True); ap.add_argument("--paths", required=True); ap.add_argument("--tasks", required=True)
    ap.add_argument("--model", default=None); ap.add_argument("--box-scale", type=float, default=1.0); ap.add_argument("--ref", type=float, nargs=2, default=None)
    ap.add_argument("--dt", type=float, default=0.02); ap.add_argument("--grades", default=",".join(ALL), help="허용 등급 (예: M3,M4,M5)")
    ap.add_argument("--curves", default=None); ap.add_argument("--betas", default=None, help="β 목록 (기본 module_curves.json 의 beta_sweep)")
    ap.add_argument("--force-scale", type=float, default=1.0, help="json force_N 에 곱해 둔 배율 (예: _f15 면 1.5) — 지금 설계 힘 비교용 표시에만")
    ap.add_argument("--max-bundles", type=int, default=None, help=f"케이블당 묶음 상한 (기본 {td.MAX_BUNDLES})"); ap.add_argument("--out-tag", default="")
    a = ap.parse_args(); pd.configure(a.space)
    if a.box_scale != 1.0 and pd.WRAP: pd.WRAP = dict(pd.WRAP, center=tuple(v * a.box_scale for v in pd.WRAP["center"]), radius=pd.WRAP["radius"] * a.box_scale)
    pj = a.paths if os.path.isabs(a.paths) else os.path.join(ROOT, a.paths); paths = json.load(open(pj, encoding="utf-8"))
    model = a.model or paths["model"]; model = model if os.path.isabs(model) else os.path.join(ROOT, model)
    tasks = a.tasks if os.path.isabs(a.tasks) else os.path.join(ROOT, a.tasks)
    Q, T, meta = pd.poses_time_from_tasks(tasks, dt=a.dt, ref=a.ref)
    tic = time.time(); fr = pd.Frames(model, Q); print(f"[{a.space}] {os.path.basename(pj)}  프레임 {len(Q)}개 (dt {a.dt} s + 극값, 시행 {len(meta['names'])}개)  뼈 변환 {time.time()-tic:.0f}s", flush=True)
    grades = [g.strip() for g in a.grades.split(",") if g.strip()]
    cjs = json.load(open(a.curves or td.CURVES_DEFAULT, encoding="utf-8")); betas = [float(b) for b in a.betas.split(",")] if a.betas else cjs.get("beta_sweep", [cjs.get("beta", 0.0)])
    pk = np.abs(T).max(axis=0); rows = []; stem = os.path.basename(pj).replace("newmus_paths_", "").replace(".json", "") + a.out_tag
    od = os.path.join(ROOT, "report", "td"); os.makedirs(od, exist_ok=True)
    for bi, beta in enumerate(betas):
        curves = td.load_curves(a.curves, beta=beta); ev = evaluate(paths, fr, T, curves, grades, max_bundles=a.max_bundles)
        res_pct = np.abs(ev["res"]).max(axis=0) / (pd.MARGIN * pk) * 100; res95 = np.percentile(np.abs(ev["res"]), 95, axis=0) / (pd.MARGIN * pk) * 100
        print(f"\n β = {beta:g} (허용 등급 {','.join(grades)}, 묶음 ≤ {a.max_bundles or td.MAX_BUNDLES})   분배 잔차 최대 {np.round(res_pct, 1)} % / 95퍼센타일 {np.round(res95, 1)} % (1.1×ID 피크 대비)")
        print(f"  {'근육':4s} {'지금 힘':>8s} {'ΔL':>6s} {'지금 등급':>9s} | {'장력 최대':>8s} {'비':>5s} {'수축 최대':>8s}  시간 영역 등급")
        for i, (m, z) in enumerate(zip(paths["muscles"], ev["sizes"])):
            F0 = m["force_N"]; dL = float(ev["Ls"][:, i].max() - ev["Ls"][:, i].min()); g0, n0 = assign_module(F0, dL * 1000)
            print(f"  {m['name']:4s} {F0:8.0f} {dL*1000:5.0f}mm {(f'{g0} ×{n0}' if g0 else '불가'):>9s} | {z['T_peak']:8.0f} {z['T_peak']/F0:5.2f} {z['s_max']*1000:7.0f}mm  {td.fmt_size(z)}")
            gall = td.size_cable(ev["T"][:, i], ev["s"][:, i], curves, ALL, a.max_bundles)
            rows.append(dict(beta=beta, name=m["name"], force_amin_N=F0, force_amin_noscale_N=F0 / a.force_scale, dL_mm=dL * 1000, grade_amin=g0, n_amin=n0,
                             T_peak_N=z["T_peak"], T_peak_over_amin=z["T_peak"] / F0, s_max_mm=z["s_max"] * 1000, grades_allowed=",".join(grades),
                             grade=z["grade"], bundles=z["n"], n_cont=z["n_cont"], work_J=z["work_J"], mass_park_g=z["mass_park_g"],
                             grade_any=gall["grade"], bundles_any=gall["n"], work_any_J=gall["work_J"], mass_any_park_g=gall["mass_park_g"],
                             arm_min_cm=float(np.einsum("tj,j->t", ev["R"][:, :, i], np.array(m["v"])).min() * 100),
                             res_max_pct=float(res_pct.max()), res95_pct=float(res95.max())))
        if bi == 0 or beta == cjs.get("beta"):
            np.savez(os.path.join(od, f"td_{stem}_b{beta:g}_timeseries.npz"), T=ev["T"], s=ev["s"], L=ev["Ls"], tau=T, Q=Q, trial=meta["trial"], time=meta["time"], names=np.array(meta["names"]), beta=beta, grades=np.array(grades))
            plot_ts(paths, ev, curves, grades, meta, f"{stem}: 시간 영역 (s, T) — β={beta:g} (임시 곡선), 허용 {','.join(grades)}", os.path.join(od, f"fig_td_{stem}_b{beta:g}.png"))
    out = os.path.join(od, f"td_{stem}.csv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("→", os.path.relpath(out, ROOT))


if __name__ == "__main__":
    main()
