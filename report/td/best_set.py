"""부위마다 Moco 를 모두 통과(보조 95퍼센타일 < 10 %)한 설계 중 모듈 일 용량이 가장 작은 것을 고르고, 160 cm 최종 A안과 같은 기준으로 합계를 낸다 (2026-11-19).

    python report/td/best_set.py → report/td/best_set.md
모듈 수·일 용량·park 질량은 한쪽 × 2 (하체·어깨·팔꿈치는 좌우, 몸통은 하나).
"""
import csv, json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
from actuator_model import BY_NAME                           # noqa: E402
from stroke_check import assign_module                       # noqa: E402
RHO = 130.0 / 12.0
UP = ["walk", "walks9", "lift", "lifts6"]; LO = ["walk", "walklab", "squat"]
REG = {  # 부위: (좌우 배수, Moco 시행, [(라벨, json 줄기, moco 이름)])
    "하체": (2, LO, [(f"{v}{s}", f"lower_k8_40kg_h175sk_{v}{s}", f"td175_lower_{v}{s}") for v in ("td", "td_br20", "td_br40") for s in ("", "_rb")]),
    "몸통": (1, UP, [(f"td{s}", f"trunk_k4_40kg_h175sk_td{s}", f"td175_trunk_td{s}") for s in ("", "_rb")]),
    "어깨": (2, UP, [(f"{v}{s}", f"shoulder_r_k8_40kg_h175sk_{v}{s}", f"td175_shoulder_{v}{s}") for v in ("td", "td_br20", "td_br40") for s in ("", "_rb")]),
    "팔꿈치": (2, UP, [(f"{v}{s}", f"elbow_r_k2_40kg_h175sk_{v}{s}", f"td175_elbow_{v}{s}") for v in ("td", "td_br20", "td_br30", "td_br40", "td_br50") for s in ("_w0", "_rb")]),
}
A160 = {"하체": (2, "lower_k8_40kg_h160_cal"), "몸통": (1, "trunk_k4_40kg_h160_mf"), "어깨": (2, "shoulder_r_k8_40kg_h160_s80d_f15"), "팔꿈치": (2, "elbow_r_k2_40kg_h160_s80")}


def modules(stem):
    P = json.load(open(os.path.join(ROOT, f"report/newmus_paths_{stem}.json"), encoding="utf-8")); out = []
    for m in P["muscles"]:
        if m.get("force_mode") == "time": g, n = m["module"]["name"], m["module"]["bundles"]
        else: g, n = assign_module(m["force_N"], (m.get("dL_m") or (m["L_max_m"] - m["L_min_m"])) * 1000)
        out.append((g, n))
    return out


def totals(mods, mult):
    w = sum(n * BY_NAME[g].force_n[1] * BY_NAME[g].stroke_mm[1] / 1000 for g, n in mods if g) * mult
    kg = sum(n * BY_NAME[g].force_n[1] / RHO for g, n in mods if g) * mult / 1000
    m1 = sum(1 for g, n in mods if g == "M1") * mult; bad = sum(1 for g, n in mods if not g) * mult
    return w, kg, m1, bad, len(mods) * mult


def main():
    mr = {}
    for r in csv.reader(open(os.path.join(ROOT, "report/moco/results_paths.csv"), encoding="utf-8")):
        if len(r) > 5: mr[r[0]] = (float(r[3]), float(r[5]))      # (보조 95퍼센타일, 케이블 토크 vs ID 95퍼센타일)
    lines = ["# 175 cm 시간 영역 최종 조합 (자동 생성: report/td/best_set.py)", "", "부위마다 Moco 시행을 모두 통과(보조 95퍼센타일 < 10 %, 케이블 토크 vs ID 95퍼센타일 < 100 %)한 변형 중 일 용량이 가장 작은 것.", "",
             "| 부위 | 고른 변형 | 케이블별 모듈 | Moco 최악 보조 / 토크 오차 (95퍼센타일) | 후보 (통과/전체) |", "| --- | --- | --- | --- | --- |"]
    pick, T = {}, [0, 0, 0, 0, 0]
    for reg, (mult, tasks, cands) in REG.items():
        ok = []
        for lab, stem, mo in cands:
            if not os.path.exists(os.path.join(ROOT, f"report/newmus_paths_{stem}.json")): continue
            res = [mr.get(f"{mo}_{t}") for t in tasks]
            if None in res: continue
            mods = modules(stem); w = totals(mods, 1)[0]; r95 = max(x[0] for x in res); e95 = max(x[1] for x in res)
            # 통과 = 보조 < 10 % 이고 케이블 토크 vs ID < 100 % (윈도 어깨 6개를 3,930 % 순간 튐으로 안 쓴 것과 같은 기준, 맥 최종 어깨 63 %)
            ok.append((r95 < 10.0 and e95 < 100.0 and all(g for g, _ in mods), w, lab, stem, (r95, e95), mods))
        good = sorted([o for o in ok if o[0]], key=lambda o: o[1])
        if not good: lines.append(f"| {reg} | (통과한 것 없음 — 결과 {len(ok)}개) | | | |"); continue
        _, w, lab, stem, worst, mods = good[0]; pick[reg] = stem
        t = totals(mods, mult); T = [a + b for a, b in zip(T, t)]
        lines.append(f"| {reg} | `{lab}` | {', '.join(f'{g}×{n}' if n > 1 else g for g, n in mods)} | {worst[0]:.1f} / {worst[1]:.0f} % | {len(good)}/{len(ok)} |")
    A = [0, 0, 0, 0, 0]
    for reg, (mult, stem) in A160.items(): A = [a + b for a, b in zip(A, totals(modules(stem), mult))]
    lines += ["", "| 양쪽 합 | 케이블 | M1 | 배정 불가 | 일 용량 합 | park 질량 (10.8 N/g) |", "| --- | --- | --- | --- | --- | --- |",
              f"| 160 cm 최종 A안 (지금 방식, 같은 기준으로 다시 셈) | {A[4]} | {A[2]} | {A[3]} | {A[0]:,.0f} J | {A[1]:.1f} kg |",
              f"| 175 cm 시간 영역 조합 | {T[4]} | {T[2]} | {T[3]} | {T[0]:,.0f} J | {T[1]:.1f} kg |", "",
              "park 질량 = 묶음 × 등급 최대 힘 ÷ 10.8 N/g (발표의 '케이블 힘 합 ÷ 10.8' 과 다르다 — 등급 최대 힘으로 셈). 힘–수축 곡선 임시 (β 0.25)."]
    out = os.path.join(ROOT, "report/td/best_set.md"); open(out, "w", encoding="utf-8").write("\n".join(lines) + "\n"); print("\n".join(lines)); print("→", os.path.relpath(out, ROOT))
    json.dump(pick, open(os.path.join(ROOT, "report/td/best_set.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)


if __name__ == "__main__":
    main()
