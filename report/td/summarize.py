"""시간 영역 설계 결과 모아 보기 (2026-11-19): 경로 json 의 케이블별 모듈(등급 × 묶음, 장력 최대, 수축 최대) + Moco 보조 95퍼센타일 (results_paths.csv).

    python report/td/summarize.py 160      → report/td/summary_160.md
    python report/td/summarize.py 175      → report/td/summary_175.md
지금 방식(amin) json 은 stroke_check 규칙(assign_module: 힘·ΔL 상한)으로 등급을 매긴다.
"""
import csv, json, os, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
from actuator_model import BY_NAME                           # noqa: E402
RHO = 130.0 / 12.0
SFX = os.environ.get("SFX", "_rb")

SETS = {
    "160": [("팔꿈치", "지금 방식 (_s80)", "elbow_r_k2_40kg_h160_s80", "td160_elbow_s80"),
            ("팔꿈치", "시간 영역, 상자", "elbow_r_k2_40kg_h160_td", "td160_elbow_td")]
           + [("팔꿈치", f"시간 영역, 브래킷 ≤ {h} mm", f"elbow_r_k2_40kg_h160_td_br{h}", f"td160_elbow_td_br{h}") for h in (0, 20, 30, 40, 50)]
           + [("어깨", "지금 방식 ×1.5 (_s80d_f15)", "shoulder_r_k8_40kg_h160_s80d_f15", "td160_shoulder_s80d_f15"),
              ("어깨", "지금 방식 ×1.5 없음 (_s80d)", "shoulder_r_k8_40kg_h160_s80d", "td160_shoulder_s80d"),
              ("어깨", "시간 영역, 상자", "shoulder_r_k8_40kg_h160_td", "td160_shoulder_td")]
           + [("어깨", f"시간 영역, 브래킷 ≤ {h} mm", f"shoulder_r_k8_40kg_h160_td_br{h}", f"td160_shoulder_td_br{h}") for h in (20, 40)],
    "175": [("하체", "시간 영역, 상자", "lower_k8_40kg_h175sk_td", "td175_lower_td")]
           + [("하체", f"시간 영역, 브래킷 ≤ {h} mm", f"lower_k8_40kg_h175sk_td_br{h}", f"td175_lower_td_br{h}") for h in (20, 40)]
           + [("몸통", "시간 영역", "trunk_k4_40kg_h175sk_td", "td175_trunk_td")]
           + [("어깨", "시간 영역, 상자", "shoulder_r_k8_40kg_h175sk_td", "td175_shoulder_td")]
           + [("어깨", f"시간 영역, 브래킷 ≤ {h} mm", f"shoulder_r_k8_40kg_h175sk_td_br{h}", f"td175_shoulder_td_br{h}") for h in (20, 40)]
           + [("팔꿈치", "시간 영역, 상자", "elbow_r_k2_40kg_h175sk_td", "td175_elbow_td")]
           + [("팔꿈치", f"시간 영역, 브래킷 ≤ {h} mm", f"elbow_r_k2_40kg_h175sk_td_br{h}", f"td175_elbow_td_br{h}") for h in (20, 30, 40, 50)]
           + [("팔꿈치", "CAD 기존 앵커 쌍", "elbow_r_k2_40kg_h175sk_cadanchor", "td175_elbow_cadanchor")],
}
TASKS = {"하체": [("walk", "Moore 걷기"), ("walklab", "LabValidation 걷기"), ("squat", "스쿼트")]}
UP_TASKS = [("walk", "걷기 s3"), ("walks9", "걷기 s9"), ("lift", "들기 s1"), ("lifts6", "들기 s6")]


def assign_amin(F, dL_mm):
    from stroke_check import assign_module
    return assign_module(F, dL_mm)


def moco_rows():
    rows = {}
    p = os.path.join(ROOT, "report/moco/results_paths.csv")
    for r in csv.reader(open(p, encoding="utf-8")):
        if len(r) > 5: rows[r[0]] = r                           # 같은 이름이면 마지막 행
    return rows


def main(which):
    mr = moco_rows(); lines = [f"# 시간 영역 설계 비교 — {which} cm (자동 생성: report/td/summarize.py)", "",
                               "모듈 = 등급 × 묶음. 일 용량 = 묶음 × 등급 최대 힘 × 최대 스트로크 (J). park 질량 = 묶음 × 최대 힘 ÷ 10.8 N/g. "
                               "Moco 칸 = 보조 액추에이터 95퍼센타일 ÷ ID 피크 (기준 10 %), 괄호 = 케이블 토크 vs ID 오차 95퍼센타일 (100 % 넘으면 순간 튐 의심). 힘–수축 곡선은 임시 (β 0.25)."]
    last = None; items = []
    for region, label, stem, moco in SETS[which]:                 # 시간 영역: 강건성 없음(팔꿈치 _w0, 나머지 원래 설계)과 강건성 0.25 (_rb) 두 줄
        if "_td" in stem:
            for sfx, tag in ((("_w0" if region == "팔꿈치" else ""), "강건성 없음"), ("_rb", "강건성 0.25")): items.append((region, f"{label} ({tag})", stem + sfx, moco + sfx))
        else: items.append((region, label, stem, moco))
    for region, label, stem, moco in items:
        pj = os.path.join(ROOT, f"report/newmus_paths_{stem}.json")
        if not os.path.exists(pj): lines.append(f"\n- {region} {label}: (결과 없음 `{stem}`)"); continue
        P = json.load(open(pj, encoding="utf-8")); tasks = TASKS.get(region, UP_TASKS)
        if region != last:
            lines += ["", f"## {region}", "", "| 설계 | 케이블별 모듈 (장력 최대 N / 수축 최대 mm) | M1 수 | 일 용량 합 J | park 질량 g | " + " | ".join(t[1] for t in tasks) + " |",
                      "| --- | --- | --- | --- | --- | " + " | ".join("---" for _ in tasks) + " |"]; last = region
        cells, work, mass, n_m1 = [], 0.0, 0.0, 0
        for m in P["muscles"]:
            if m.get("force_mode") == "time" or isinstance(m.get("module"), dict) and m.get("T_peak_N"):
                g, n = m["module"]["name"], m["module"]["bundles"]; Tpk = m.get("T_peak_N", m["force_N"]); s = m.get("s_max_m", m.get("dL_m", 0)) * 1000
            else:
                dL = m.get("dL_m") or (m["L_max_m"] - m["L_min_m"]); g, n = assign_amin(m["force_N"], dL * 1000); Tpk = m["force_N"]; s = dL * 1000
            if g:
                F, S = BY_NAME[g].force_n[1], BY_NAME[g].stroke_mm[1] / 1000; work += n * F * S; mass += n * F / RHO; n_m1 += (g == "M1")
            cells.append(f"{m['name']} {g or '불가'}{'×' + str(n) if g and n and n > 1 else ''} ({Tpk:.0f}/{s:.0f})")
        res = []
        for t, _ in tasks:
            r = mr.get(f"{moco}_{t}"); res.append(f"{float(r[3]):.1f} ({float(r[5]):.0f}) %" if r else "—")
        lines.append(f"| {label} | {', '.join(cells)} | {n_m1} | {work:.0f} | {mass:.0f} | " + " | ".join(res) + " |")
    out = os.path.join(ROOT, f"report/td/summary_{which}.md"); open(out, "w", encoding="utf-8").write("\n".join(lines) + "\n"); print("\n".join(lines)); print("→", os.path.relpath(out, ROOT))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "160")
