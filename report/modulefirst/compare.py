"""모듈 먼저 실험(2026-11-05) 비교표: 부위별로 현재 방식(사후 배정)과 새 방식(모듈 먼저)을 나란히.

숫자는 전부 결과 파일에서 읽는다: 경로 json(힘·2단계 잔차), stroke_check csv(ΔL·사후 등급), 모듈 json(사전 배정·예측 일 W),
report/moco/results_paths.csv(Moco 보조 95퍼센타일 = 넷째 열).

    python report/modulefirst/compare.py            → 화면에 마크다운 표, report/modulefirst/compare_summary.csv · compare_cables.csv
"""
import csv, json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
MF = os.path.join(ROOT, "report/modulefirst"); GR = ("M1", "M2", "M3", "M4", "M5", "M6")

# (부위, 방식 라벨, 경로 json, stroke_check csv, 모듈 json 또는 None, Moco 이름 (걷기, 스쿼트/들기))
RUNS = [
    ("하체", "맥 채택 (수정 전, 사후)", "report/newmus_paths_lower_k8_40kg_h160_s120.json", "report/stroke_check_lower_k8_40kg_h160_s120.csv", None,
     ("winchk_lower_s120_walk_pf", "winchk_lower_s120_squat")),
    ("하체", "현재 방식 (토크 40 kg 수정, 사후)", "report/newmus_paths_lower_k8_40kg_h160_s120_ref.json", f"{MF}/stroke_check_lower_k8_40kg_h160_s120_ref.csv", None,
     ("mf_lower_cur_walk", "mf_lower_cur_squat")),
    ("하체", "현재 방식 + 신뢰 구간 (민감도)", "report/newmus_paths_lower_k8_40kg_h160_s120_ref_tr.json", f"{MF}/stroke_check_lower_k8_40kg_h160_s120_ref_tr.csv", None,
     ("mf_lower_curtr_walk", "mf_lower_curtr_squat")),
    ("하체", "새 방식 (모듈 먼저, 신뢰 구간)", "report/newmus_paths_lower_k8_40kg_h160_mf.json", f"{MF}/stroke_check_lower_k8_40kg_h160_mf.csv", f"{MF}/modules_lower_k8_40kg_h160.json",
     ("mf_lower_new_walk", "mf_lower_new_squat")),
    ("몸통", "맥 채택 (사후)", "report/newmus_paths_trunk_k4_40kg_h160_s80.json", "report/stroke_check_trunk_k4_40kg_h160_s80.csv", None,
     ("winchk_trunk_s80_walk", "winchk_trunk_s80_lift")),
    ("몸통", "현재 방식 (윈도 재실행, 사후)", "report/newmus_paths_trunk_k4_40kg_h160_s80_win.json", f"{MF}/stroke_check_trunk_k4_40kg_h160_s80_win.csv", None,
     ("mf_trunk_cur_walk", "mf_trunk_cur_lift")),
    ("몸통", "현재 방식 + 신뢰 구간 (민감도)", "report/newmus_paths_trunk_k4_40kg_h160_s80_win_tr.json", f"{MF}/stroke_check_trunk_k4_40kg_h160_s80_win_tr.csv", None,
     ("mf_trunk_curtr_walk", "mf_trunk_curtr_lift")),
    ("몸통", "새 방식 (모듈 먼저, 신뢰 구간)", "report/newmus_paths_trunk_k4_40kg_h160_mf.json", f"{MF}/stroke_check_trunk_k4_40kg_h160_mf.csv", f"{MF}/modules_trunk_k4_40kg_h160.json",
     ("mf_trunk_new_walk", "mf_trunk_new_lift")),
    ("팔꿈치", "맥 채택 (사후)", "report/newmus_paths_elbow_r_k2_40kg_h160_s80.json", "report/stroke_check_elbow_r_k2_40kg_h160_s80.csv", None,
     ("winchk_elbow_s80_walk", "winchk_elbow_s80_lift")),
    ("팔꿈치", "현재 방식 (윈도 재실행, 사후)", "report/newmus_paths_elbow_r_k2_40kg_h160_s80_win.json", f"{MF}/stroke_check_elbow_r_k2_40kg_h160_s80_win.csv", None,
     ("mf_elbow_cur_walk", "mf_elbow_cur_lift")),
    ("팔꿈치", "현재 방식 + 신뢰 구간 (민감도)", "report/newmus_paths_elbow_r_k2_40kg_h160_s80_win_tr.json", f"{MF}/stroke_check_elbow_r_k2_40kg_h160_s80_win_tr.csv", None,
     ("mf_elbow_curtr_walk", "mf_elbow_curtr_lift")),
    ("팔꿈치", "새 방식 (모듈 먼저, 신뢰 구간)", "report/newmus_paths_elbow_r_k2_40kg_h160_mf.json", f"{MF}/stroke_check_elbow_r_k2_40kg_h160_mf.csv", f"{MF}/modules_elbow_r_k2_40kg_h160.json",
     ("mf_elbow_new_walk", "mf_elbow_new_lift")),
    ("어깨", "맥 채택 (사후, 힘 ×1.5)", "report/newmus_paths_shoulder_r_k6_40kg_h160_v635_s80_f15.json", "report/stroke_check_shoulder_r_k6_40kg_h160_v635_s80_f15.csv", None,
     ("winchk_shoulder_s80f15_walk", "winchk_shoulder_s80f15_lift")),
    ("어깨", "현재 방식 (윈도 재실행, 사후, 힘 ×1.5)", "report/newmus_paths_shoulder_r_k6_40kg_h160_v635_s80_win_f15.json", f"{MF}/stroke_check_shoulder_r_k6_40kg_h160_v635_s80_win_f15.csv", None,
     ("mf_shoulder_cur_walk", "mf_shoulder_cur_lift")),
    ("어깨", "현재 방식 + 신뢰 구간 (민감도, 힘 ×1.5)", "report/newmus_paths_shoulder_r_k6_40kg_h160_v635_s80_win_tr_f15.json", f"{MF}/stroke_check_shoulder_r_k6_40kg_h160_v635_s80_win_tr_f15.csv", None,
     ("mf_shoulder_curtr_walk", "mf_shoulder_curtr_lift")),
    ("어깨", "새 방식 (모듈 먼저, 신뢰 구간, 힘 ×1.5)", "report/newmus_paths_shoulder_r_k6_40kg_h160_v635_mf_f15.json", f"{MF}/stroke_check_shoulder_r_k6_40kg_h160_v635_mf_f15.csv", f"{MF}/modules_shoulder_r_k6_40kg_h160_v635_f15.json",
     ("mf_shoulder_new_walk", "mf_shoulder_new_lift")),
]


def p(x): return x if os.path.isabs(x) else os.path.join(ROOT, x)


def moco_r95():
    out = {}
    for row in csv.reader(open(os.path.join(ROOT, "report/moco/results_paths.csv"), encoding="utf-8")):
        if len(row) > 3: out[row[0]] = float(row[3])                  # 같은 이름이 여러 번이면 마지막 줄
    return out


def main():
    mo = moco_r95(); summ, cab = [], []
    for region, label, pj, sc, mj, (mw, ml) in RUNS:
        if not (os.path.exists(p(pj)) and os.path.exists(p(sc))): print(f"  [없음] {region} {label}: {pj if not os.path.exists(p(pj)) else sc}"); continue
        P = json.load(open(p(pj), encoding="utf-8")); S = list(csv.DictReader(open(p(sc), encoding="utf-8"))); Mo = json.load(open(p(mj), encoding="utf-8")) if mj else None
        cnt = {g: sum(r["module"] == g for r in S) for g in GR}; none = sum(not r["module"] for r in S)
        viol = None
        if Mo:                                                          # 사전 배정 등급의 상한을 넘은 케이블 수 (힘 또는 ΔL)
            viol = sum(float(r["force_N"]) > m["force_cap_N"] + 1e-6 or float(r["dL_cm"]) / 100 > m["stroke_cap_m"] + 1e-9 for r, m in zip(S, Mo["muscles"]))
        rf = P.get("refine", {})
        summ.append(dict(region=region, method=label, **{g: cnt[g] for g in GR}, none=none, pre_violations=viol, worst_pct=rf.get("worst_pct"), mean_pct=rf.get("mean_pct"),
                         walk_r95=mo.get(mw), other_r95=mo.get(ml), moco_walk=mw, moco_other=ml))
        for i, r in enumerate(S):
            m = Mo["muscles"][i] if Mo else None
            cab.append(dict(region=region, method=label, cable=r["name"], dL_mm=float(r["dL_cm"]) * 10, force_N=float(r["force_N"]),
                            post=f"{r['module']} ×{r['bundles']}" if r["module"] else "없음", pre=(f"{m['module']} ×{m['bundles']}" if m else ""),
                            pre_cap=(f"{m['force_cap_N']:.0f} N / {m['stroke_cap_m']*1000:.0f} mm" if m else ""),
                            ok=("" if not m else ("지킴" if float(r["force_N"]) <= m["force_cap_N"] + 1e-6 and float(r["dL_cm"]) / 100 <= m["stroke_cap_m"] + 1e-9 else "넘음"))))
    fmt = lambda v, d=1: "—" if v is None else f"{v:.{d}f}"
    print("\n| 부위 | 방식 | M1 | M2 | M3 | M4 | M6 | 배정 불가 | 사전 등급 넘은 케이블 | 2단계 잔차 최악 % | 평균 % | Moco 걷기 95퍼센타일 % | 스쿼트/들기 % |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for s in summ:
        print(f"| {s['region']} | {s['method']} | {s['M1']} | {s['M2']} | {s['M3']} | {s['M4']} | {s['M6']} | {s['none']} | {'—' if s['pre_violations'] is None else s['pre_violations']} | "
              f"{fmt(s['worst_pct'])} | {fmt(s['mean_pct'], 2)} | {fmt(s['walk_r95'])} | {fmt(s['other_r95'])} |")
    for region in dict.fromkeys(s["region"] for s in summ):
        meths = [s["method"] for s in summ if s["region"] == region]
        print(f"\n**{region}** — 케이블별 ΔL mm / 힘 N / 등급 (새 방식은 사후 등급, 괄호 안이 사전 배정)\n")
        print("| 케이블 | " + " | ".join(meths) + " |"); print("|---|" + "---|" * len(meths))
        names = [c["cable"] for c in cab if c["region"] == region and c["method"] == meths[0]]
        for nm in names:
            cells = []
            for me in meths:
                c = next(c for c in cab if c["region"] == region and c["method"] == me and c["cable"] == nm)
                cells.append(f"{c['dL_mm']:.0f} / {c['force_N']:.0f} / {c['post']}" + (f" ({c['pre']}, {c['ok']})" if c["pre"] else ""))
            print(f"| {nm} | " + " | ".join(cells) + " |")
    for name, rows in (("compare_summary.csv", summ), ("compare_cables.csv", cab)):
        with open(os.path.join(MF, name), "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print("\n→", os.path.join(MF, "compare_summary.csv"), os.path.join(MF, "compare_cables.csv"))


if __name__ == "__main__":
    main()
