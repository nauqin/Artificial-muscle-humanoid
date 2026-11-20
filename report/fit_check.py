"""새 동작을 **지금 설계(최종 A안)** 가 낼 수 있나 — Moco 없이 프레임마다 재는 빠른 검사 (수 초).

    python report/fit_check.py --tasks report/custom/tasks_custom.csv            # 표의 모든 trial
    python report/fit_check.py --tasks report/custom/tasks_custom.csv --trial opencap_sample_squat
→ 화면 표 + report/custom/fit_check_<tasks 이름>.csv

trial 마다, 부위(하체·몸통·어깨·팔꿈치)마다 세 가지를 잰다.
  1) 토크형: 관절 배분비 V·토크 용량 c 만으로 (0 ≤ u ≤ 1) ID 토크를 낼 수 있나 → 남는 토크. "근육 개수·용량이 모자란가"
  2) 케이블형: 실제 케이블 경로의 그 자세 모멘트암 r(q) × 케이블 힘으로 낼 수 있나 → 남는 토크. "그 자세에서 지렛대가 모자란가"
  3) 길이: 그 자세들에서 케이블 길이가 설계 때 범위를 얼마나 벗어나나 → 모듈 스트로크를 넘는가
`static_check.py` 를 설계 공간·tasks.csv 로 일반화한 것이다. 좌우가 있는 부위는 왼쪽 자세·토크를 오른쪽 케이블(거울 대칭)로 잰다.

판정: 케이블형 남는 토크의 95퍼센타일이 그 trial 피크 토크의 10 % 이하이고 길이 초과가 없으면 "가능".
이 검사는 정적이다(프레임 사이 연결·케이블 속도 없음). 통과해도 확정은 `moco_path_muscles.py` 로 한다. 미달이면 재설계 대상이다.
"""
import argparse, csv, json, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import addb_io as io                                   # noqa: E402
import design_space as ds                              # noqa: E402

#: 최종 A안 (2026-11-09 결정, report/upper/README.md 10-4): 하체 방식 3(_cal) · 몸통 방식 2(_mf) · 어깨 8개(_s80d_f15) · 팔꿈치(_s80)
FINAL_A = {
    "lower": dict(design="report/newmus_V_design_lower_k8_40kg_h160.npz", paths="report/newmus_paths_lower_k8_40kg_h160_cal.json", stroke=0.12),
    "trunk": dict(design="report/newmus_V_design_trunk_k4_40kg_h160.npz", paths="report/newmus_paths_trunk_k4_40kg_h160_mf.json", stroke=0.08),
    "shoulder_r": dict(design="report/newmus_V_design_shoulder_r_k8_40kg_h160.npz", paths="report/newmus_paths_shoulder_r_k8_40kg_h160_s80d_f15.json", stroke=0.08),
    "elbow_r": dict(design="report/newmus_V_design_elbow_r_k2_40kg_h160.npz", paths="report/newmus_paths_elbow_r_k2_40kg_h160_s80.json", stroke=0.08),
}
MODEL = "data/robot40_h160/robot40_h160.osim"; BOX_SCALE = 0.9877; REF = (40.0, 1.60)
PASS_PCT = 10.0


def side_name(name, side):
    """오른쪽 기준 좌표 이름을 side 쪽으로 ('hip_rotation_r' → 'hip_rotation_l')."""
    return name[:-2] + "_" + side if (side and name.endswith("_r")) else name


def box_lsq(A, b):
    from scipy.optimize import lsq_linear
    return lsq_linear(A, b, bounds=(0.0, 1.0), lsmr_tol="auto").x


def check_trial(row, spaces, max_frames=240, log=print):
    import path_design as pd
    ikp = row["ik"] if os.path.isabs(row["ik"]) else os.path.join(ROOT, row["ik"]); idp = row["id"] if os.path.isabs(row["id"]) else os.path.join(ROOT, row["id"])
    ik = io.read_mot(ikp); idm = io.read_mot(idp)
    if ik.in_degrees: ik.data = np.radians(ik.data)
    scale = (REF[0] * REF[1]) / (float(row["mass_kg"]) * float(row["height_m"]))
    sel = (ik.time >= idm.t0) & (ik.time <= idm.t1); tt = ik.time[sel]; tt = tt[::max(1, len(tt) // max_frames)]
    out = []
    for sp_name in spaces:
        cfg = FINAL_A[sp_name]; SP = pd.configure(sp_name)
        if pd.WRAP: pd.WRAP = dict(pd.WRAP, center=tuple(v * BOX_SCALE for v in pd.WRAP["center"]), radius=pd.WRAP["radius"] * BOX_SCALE)
        d = np.load(os.path.join(ROOT, cfg["design"])); Vnm = d["V_nm"]                       # (k, n)
        paths = json.load(open(os.path.join(ROOT, cfg["paths"]), encoding="utf-8")); mus = paths["muscles"]
        n = len(SP["base"]); res_t, res_c, pk_all, tau_all, Ls, n_skip = [], [], [], [], [], 0
        for side in SP["sides"]:
            qn = [side_name(c, side) for c in list(pd.COORDS) + list(pd.EXTRA)]; mn = [ds.full(b, side) + "_moment" for b in SP["base"]]
            if not all(c in ik.names for c in qn) or not all(c in idm.names for c in mn): continue
            Q = np.column_stack([np.interp(tt, ik.time, ik.column(c)) for c in qn]); tau = np.column_stack([np.interp(tt, idm.time, idm.column(c)) for c in mn]) * scale
            keep = np.ones(len(Q), bool)
            for cname, (lo, hi) in SP.get("pose_limit", {}).items():                           # 어깨 특이점 근처 자세는 설계에서도 뺐다
                j = list(pd.COORDS + pd.EXTRA).index(cname); keep &= (Q[:, j] >= lo) & (Q[:, j] <= hi)
            n_skip += int((~keep).sum()); Q, tau = Q[keep], tau[keep]
            if not len(Q): continue
            fr = pd.Frames(os.path.join(ROOT, MODEL), Q); A = np.zeros((len(Q), n, len(mus))); L = np.zeros((len(Q), len(mus)))
            for i, m in enumerate(mus):
                chain = [pd.BODIES.index(b) for b in m["chain"]]; r, L0 = pd.moment_arms(fr, chain, np.array(m["points"])); A[:, :, i] = float(m["force_N"]) * r; L[:, i] = L0
            rt_ = np.array([Vnm.T @ box_lsq(Vnm.T, tau[t]) - tau[t] for t in range(len(Q))]); rc_ = np.array([A[t] @ box_lsq(A[t], tau[t]) - tau[t] for t in range(len(Q))])
            res_t.append(rt_); res_c.append(rc_); tau_all.append(tau); Ls.append(L)
        if not tau_all: log(f"  {sp_name:11s} 좌표 없음 — 건너뜀"); continue
        res_t, res_c, tau_all, Ls = np.vstack(res_t), np.vstack(res_c), np.vstack(tau_all), np.vstack(Ls)
        pk = np.abs(tau_all).max(axis=0); ref = max(float(pk.max()), 1e-6)
        wt = np.abs(res_t).max(axis=1); wc = np.abs(res_c).max(axis=1); jworst = int(np.abs(res_c).max(axis=0).argmax())
        # 길이: 이 동작이 케이블 길이를 설계 때 범위(json 의 L_min_m/L_max_m) 밖으로 얼마나 넓히나, 그 결과 전체 ΔL 이 모듈 스트로크를 넘나.
        # 설계 범위 자체가 이미 상한을 몇 mm 넘는 케이블이 있어서(벌점 방식 — 팔꿈치 M0 83.7 mm vs 80), "1 mm 넘게 넓히고 그 결과 상한도 넘을 때" 만 초과로 본다.
        over, widen = [], []
        for i, m in enumerate(mus):
            lo, hi = float(Ls[:, i].min()), float(Ls[:, i].max()); cap = float((m.get("module") or {}).get("stroke_cap_m") or cfg["stroke"])
            if "L_min_m" in m:
                wd = max(0.0, hi - m["L_max_m"]) + max(0.0, m["L_min_m"] - lo); lo, hi = min(lo, m["L_min_m"]), max(hi, m["L_max_m"])
            else: wd = hi - lo
            widen.append(wd); over.append((hi - lo) - cap if wd > 0.001 else 0.0)
        over_mm = max(0.0, max(over)) * 1000; iw = int(np.argmax(over)) if over_mm > 0 else int(np.argmax(widen))
        r = dict(trial=row["trial"], task=row["task"], space=sp_name, frames=len(tau_all), skipped_pose=n_skip, peak_Nm=round(ref, 1),
                 peak_by_coord=" ".join(f"{p:.1f}" for p in pk), torque_type_p95_pct=round(float(np.percentile(wt, 95)) / ref * 100, 1), torque_type_max_Nm=round(float(wt.max()), 2),
                 cable_p95_pct=round(float(np.percentile(wc, 95)) / ref * 100, 1), cable_max_pct=round(float(wc.max()) / ref * 100, 1), cable_max_Nm=round(float(wc.max()), 2),
                 worst_coord=SP["base"][jworst], stroke_widen_mm=round(max(widen) * 1000, 1), stroke_over_mm=round(over_mm, 1), stroke_worst_cable=mus[iw]["name"], dL_has_design_range=("L_min_m" in mus[0]))
        r["verdict"] = "가능" if (r["cable_p95_pct"] <= PASS_PCT and over_mm <= 0) else ("토크 부족" if r["cable_p95_pct"] > PASS_PCT else "길이 초과")
        if r["verdict"] == "토크 부족" and r["torque_type_p95_pct"] <= PASS_PCT: r["verdict"] = "경로 부족 (용량은 됨)"
        out.append(r)
        log(f"  {sp_name:11s} 피크 {ref:6.1f} N·m | 토크형 남는 토크 95% {r['torque_type_p95_pct']:5.1f} % · 케이블형 95% {r['cable_p95_pct']:5.1f} % (최대 {r['cable_max_Nm']:.1f} N·m, {r['worst_coord']}) | "
            f"길이 범위 넓힘 {r['stroke_widen_mm']:.0f} mm·상한 초과 {r['stroke_over_mm']:.0f} mm{'' if r['dL_has_design_range'] else '*'} → {r['verdict']}" + (f"  [자세 제외 {n_skip}]" if n_skip else ""))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tasks", required=True); ap.add_argument("--trial", nargs="+", default=None, help="이 trial 만")
    ap.add_argument("--spaces", default="lower,trunk,shoulder_r,elbow_r"); ap.add_argument("--out", default=None)
    a = ap.parse_args(); tp = a.tasks if os.path.isabs(a.tasks) else os.path.join(ROOT, a.tasks); rows = []
    for row in csv.DictReader(open(tp, encoding="utf-8")):
        if not row.get("id") or row["id"].startswith("#") or (a.trial and row["trial"] not in a.trial): continue
        print(f"\n[{row['trial']}] 과제 {row['task']}, 손 하중 {row.get('hand_load_kg') or 0} kg")
        spaces = [s for s in a.spaces.split(",") if s]
        if (row.get("note") or "").startswith("[상체만]"): spaces = [s for s in spaces if s != "lower"]
        rows += check_trial(row, spaces)
    if not rows: raise SystemExit("검사할 trial 이 없다")
    out = a.out or os.path.join(os.path.dirname(tp), f"fit_check_{os.path.splitext(os.path.basename(tp))[0]}.csv")
    with open(out, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    print(f"\n저장 {out}\n(* = 경로 json 에 설계 때 길이 범위가 없어 이 trial 안의 길이 변화만 스트로크 상한과 비교)")


if __name__ == "__main__":
    main()
