"""모듈 먼저 정하기 (2026-11-05) — 경로를 찾기 전에 케이블별 SMA 모듈 등급·묶음 수와 필요한 모멘트암 범위를 정한다.

입력은 1단계 결과(V 단위 방향, c 용량 N·m)와 자세 표본뿐이다. 케이블 i 가 해야 하는 일
    W_i = s × 1.1 × c_i × Δθ_i   (J = N·m)
      1.1  : path_design 의 용량 여유 MARGIN
      s    : 힘 배율 (어깨는 경로 설계 뒤 힘 ×1.5 를 하므로 1.5, 나머지 1)
      Δθ_i : 자세 표본에서 v_i·q 가 움직이는 범위 (rad) = max_t(v_i·q_t) − min_t(v_i·q_t)
모멘트암이 자세에 따라 변하지 않고(크기 ρ) 방향이 v_i 와 같다고 보면 힘 F = s·1.1·c/ρ, 케이블 길이 변화 ΔL = ρ·Δθ 이므로 F × ΔL = W.
등급 g(힘 상한 F_g, 스트로크 상한 S_g)를 n 묶음 쓰면 F ≤ n·F_g, ΔL ≤ S_g 이어야 하고
    ρ ∈ [ s·1.1·c / (n·F_g) ,  S_g / Δθ ]
가 비어 있지 않으려면 n ≥ W / (F_g·S_g). 그래서 등급마다 n_g = max(1, ceil(W / (F_g·S_g))).
등급 고르기는 stroke_check.assign_module 과 같은 규칙: 묶음 수가 가장 적은 등급, 같으면 힘 상한이 작은 등급.
(ρ 상한 S_g/Δθ 가 path_design 의 모멘트암 하한 1.5 cm 보다 작은 등급은 못 쓴다 — M5 만 해당.)

이 예측은 낙관적이다. 실제 경로는 힘이 최소 모멘트암으로, 길이 변화가 평균 모멘트암으로 정해지고(F·ΔL = W × 평균/최소),
2단계가 힘 배율 f ≥ 1 을 올리기도 하며, 모멘트암이 v 와 꼭 나란하지도 않다. --check 로 실제 경로의 F·ΔL 과 비교한다.

자세 표본 = path_design 의 설계 자세 ∪ stroke_check 의 검사 자세(n_per 20) — 등급 판정(stroke_check)이 보는 자세를 빠뜨리지 않으려고.

    python report/module_first.py --space lower --design report/newmus_V_design_lower_k8_40kg_h160.npz --out report/modulefirst/modules_lower.json
    python report/module_first.py --space trunk --tasks report/upper/tasks_40kg_h160.csv --design report/newmus_V_design_trunk_k4_40kg_h160.npz --out ...
    python report/module_first.py --space shoulder_r --tasks ... --design ..._v635.npz --force-scale 1.5 --out ...
    비교: --check <경로 json (최종 힘)> <stroke_check csv> [--check-out csv]  (등급 json 은 --out 에 준 것을 읽는다)
"""
import argparse, csv, json, math, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import path_design as pd                                     # noqa: E402
from actuator_model import MODULES                           # noqa: E402  과제 계획서 모듈표 (힘 N, 스트로크 mm)

R_MIN = 0.015                                                # path_design 의 모멘트암 하한 벌점 기준 (m)


def grade_table(W, cs, dtheta):
    """등급마다 (묶음 수 n, 힘 상한 N, 이름, ρ 하한, ρ 상한, 쓸 수 있나). cs = s·1.1·c (N·m)."""
    rows = []
    for m in MODULES:
        Fg, Sg = float(m.force_n[1]), m.stroke_mm[1] / 1000.0
        n = max(1, int(math.ceil(W / (Fg * Sg) - 1e-12)))
        r_lo, r_hi = cs / (n * Fg), Sg / dtheta
        rows.append(dict(module=m.name, bundles=n, force_cap_N=n * Fg, stroke_cap_m=Sg, r_lo_m=r_lo, r_hi_m=r_hi, ok=bool(r_hi >= max(r_lo, R_MIN))))
    return rows


def assign(rows):
    """stroke_check.assign_module 과 같은 기준: 묶음 수 최소 → 같으면 힘 상한(1묶음)이 작은 등급."""
    ok = [r for r in rows if r["ok"]]
    return min(ok, key=lambda r: (r["bundles"], r["force_cap_N"] / r["bundles"], r["module"])) if ok else None


def pose_sample(space, tasks):
    pd.configure(space)
    if tasks:
        tp = tasks if os.path.isabs(tasks) else os.path.join(ROOT, tasks)
        Qd = pd.poses_from_tasks(tp)[0]; Qs = pd.poses_from_tasks(tp, n_per=20)[0]       # path_design 설계 자세, stroke_check 검사 자세
    else:
        Qd = pd.poses()[0]; Qs = pd.poses(n_walk=20, n_sq=60)[0]
    return np.vstack([Qd, Qs]), len(Qd), len(Qs)


def predict(a):
    Q, nd, ns = pose_sample(a.space, a.tasks)
    d = np.load(a.design if os.path.isabs(a.design) else os.path.join(ROOT, a.design)); V, c = d["V"], d["c"]
    th = Q[:, :len(pd.COORDS)] @ V.T                                                         # (자세, k) v·q
    dth = th.max(axis=0) - th.min(axis=0)
    cal = {}
    if a.calibrate:                                                                          # 보정: 케이블마다 W × (어떤 경로에서 잰 실제 F·ΔL / W) 로 등급을 고른다 (2026-11-05 2차)
        cal = {r["name"]: float(r["ratio"]) for r in csv.DictReader(open(a.calibrate if os.path.isabs(a.calibrate) else os.path.join(ROOT, a.calibrate), encoding="utf-8"))}
    out = dict(space=a.space, design=os.path.relpath(a.design if os.path.isabs(a.design) else os.path.join(ROOT, a.design), ROOT), tasks=a.tasks,
               force_scale=a.force_scale, margin=pd.MARGIN, n_poses=int(len(Q)), n_design_poses=nd, n_check_poses=ns, r_min_m=R_MIN, muscles=[])
    if a.calibrate: out["calibrated_from"] = a.calibrate
    print(f"[{a.space}] 자세 {len(Q)}개 (설계 {nd} + 검사 {ns}), 힘 배율 {a.force_scale:g}, 여유 {pd.MARGIN}" + (f", 보정 {a.calibrate}" if a.calibrate else ""))
    print(f"{'근육':4s} {'용량':>6s} {'Δθ':>5s} {'일 W':>6s}  {'등급':8s} {'힘 상한':>7s} {'ΔL 상한':>7s}  모멘트암 범위 (cm)   | 등급별 묶음 수")
    for i in range(len(V)):
        cs = a.force_scale * pd.MARGIN * c[i]; W = cs * dth[i]; k_cal = cal.get(f"M{i}", 1.0)
        rows = grade_table(W * k_cal, cs * k_cal, dth[i]); best = assign(rows)              # 보정하면 힘 쪽(ρ 하한)도 같은 비로 키운다: 실제 경로는 최소 모멘트암이 힘을 정하므로
        mm = dict(name=f"M{i}", v=V[i].tolist(), capacity_Nm=float(c[i]), dtheta_rad=float(dth[i]), work_J=float(W), calibration_ratio=k_cal, work_cal_J=float(W * k_cal),
                  module=best["module"] if best else None, bundles=best["bundles"] if best else None,
                  force_cap_N=best["force_cap_N"] if best else None, stroke_cap_m=best["stroke_cap_m"] if best else None,
                  r_lo_m=best["r_lo_m"] if best else None, r_hi_m=best["r_hi_m"] if best else None, grades=rows)
        out["muscles"].append(mm)
        per = "  ".join(f"{r['module']}×{r['bundles']}" + ("" if r["ok"] else "(불가)") for r in rows)
        if best: print(f"M{i:<3d} {c[i]:6.1f} {dth[i]:5.2f} {W:6.1f}" + (f" ×{k_cal:.2f}={W*k_cal:6.1f}" if a.calibrate else "") + f"  {best['module']} ×{best['bundles']:<3d} {best['force_cap_N']:7.0f} {best['stroke_cap_m']*1000:5.0f} mm  {best['r_lo_m']*100:5.2f} ~ {best['r_hi_m']*100:5.2f}      | {per}")
        else: print(f"M{i:<3d} {c[i]:6.1f} {dth[i]:5.2f} {W:6.1f}  배정 불가 | {per}")
    if a.out:
        op = a.out if os.path.isabs(a.out) else os.path.join(ROOT, a.out); os.makedirs(os.path.dirname(op), exist_ok=True)
        json.dump(out, open(op, "w", encoding="utf-8"), indent=1, ensure_ascii=False); print("→", op)
    return out


def check(mods, paths_json, stroke_csv, out_csv=None):
    """실제 경로와 예측 비교. F = 경로 json 의 최종 힘(어깨는 ×1.5 뒤), ΔL = stroke_check 가 잰 값.
    F·ΔL / W = f × (평균 모멘트암 / 최소 모멘트암) × (ΔL / (평균 모멘트암 × Δθ))  — 각각 2단계 힘 배율, 모멘트암 변화, 방향 어긋남·기타."""
    P = json.load(open(paths_json, encoding="utf-8")); S = {r["name"]: r for r in csv.DictReader(open(stroke_csv, encoding="utf-8"))}
    s = float(mods.get("force_scale", 1.0)); rows = []
    print(f"{'근육':4s} {'예측':10s} {'실제(사후)':10s} {'W J':>6s} {'F·ΔL J':>7s} {'비':>5s} = {'f':>5s} × {'평균/최소':>8s} × {'기타':>5s} | 모멘트암 예측 범위 → 실제 최소~평균 (cm)")
    for mm in mods["muscles"]:
        m = next(x for x in P["muscles"] if x["name"] == mm["name"]); st = S[mm["name"]]
        F = float(m["force_N"]); dL = float(st["dL_cm"]) / 100; arm_min = max(float(m["arm_min"]), 0.005); arm_mean = float(m["arm_mean"])
        ratio = F * dL / mm["work_J"]; f = float(m.get("force_mult", 1.0)); g1 = arm_mean / arm_min; g2 = dL / (arm_mean * mm["dtheta_rad"])
        post = f"{st['module']} ×{st['bundles']}" if st["module"] else "없음"
        pre = f"{mm['module']} ×{mm['bundles']}" if mm["module"] else "불가"
        rows.append(dict(name=mm["name"], predicted=pre, actual=post, work_J=mm["work_J"], F_N=F, dL_mm=dL * 1000, product_J=F * dL, ratio=ratio, force_mult=f,
                         arm_mean_over_min=g1, other=g2, r_lo_cm=(mm["r_lo_m"] or 0) * 100, r_hi_cm=(mm["r_hi_m"] or 0) * 100, arm_min_cm=arm_min * 100, arm_mean_cm=arm_mean * 100))
        print(f"{mm['name']:4s} {pre:10s} {post:10s} {mm['work_J']:6.1f} {F*dL:7.1f} {ratio:5.2f} = {f:5.2f} × {g1:8.2f} × {g2:5.2f} | {rows[-1]['r_lo_cm']:.2f}~{rows[-1]['r_hi_cm']:.2f} → {arm_min*100:.2f}~{arm_mean*100:.2f}")
    if out_csv:
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
        print("→", out_csv)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--space", default="lower"); ap.add_argument("--design", required=True, help="1단계 npz (V, c)")
    ap.add_argument("--tasks", default=None, help="자세 표본 trial 목록 (없으면 하체 기본 = path_design.poses)")
    ap.add_argument("--force-scale", type=float, default=1.0, help="경로 설계 뒤에 힘에 곱하는 배율 (어깨 1.5)")
    ap.add_argument("--calibrate", default=None, help="보정: --check 가 낸 csv(열 name, ratio). 케이블마다 W × ratio 로 등급을 고른다 (2026-11-05 2차)")
    ap.add_argument("--out", default=None, help="등급 json (path_design --modules 입력)")
    ap.add_argument("--check", nargs=2, default=None, metavar=("PATHS_JSON", "STROKE_CSV"), help="실제 경로와 예측 비교 (--out 의 등급 json 을 읽는다)")
    ap.add_argument("--check-out", default=None)
    a = ap.parse_args()
    if a.check:
        pd.configure(a.space); op = a.out if os.path.isabs(a.out) else os.path.join(ROOT, a.out)
        check(json.load(open(op, encoding="utf-8")), *a.check, out_csv=a.check_out)
    else:
        predict(a)


if __name__ == "__main__":
    main()
