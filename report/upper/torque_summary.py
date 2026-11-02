"""tasks.csv 의 trial 들에서 상체 좌표별 피크 토크(+/−)를 과제별로 정리한다 — 설계 좌표·길항쌍 판단 근거.

    python report/upper/torque_summary.py [report/upper/tasks.csv]
→ 표(stdout) + report/upper/torque_summary.csv + fig_torque_summary.png
토크는 REF(63.5 kg × 1.69 m)/(mass×height) 로 정규화해 피험자3 몸 크기로 맞춘 값(N·m)이다.
"""
import csv, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, ROOT)
import addb_io as io                                   # noqa: E402
REF = 63.5 * 1.69
COORDS = ["lumbar_extension", "lumbar_bending", "lumbar_rotation", "arm_flex_r", "arm_add_r", "arm_rot_r", "elbow_flex_r", "pro_sup_r", "wrist_flex_r", "wrist_dev_r"]
for f in font_manager.findSystemFonts():
    if "AppleSDGothicNeo" in f or "Malgun" in f: plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
plt.rcParams["axes.unicode_minus"] = False


def main():
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "report/upper/tasks.csv")
    per = {}                                           # (task, coord) -> list of (max, min)
    for row in csv.DictReader(open(path, encoding="utf-8")):
        if not row.get("id") or row["id"].startswith("#"): continue
        f = row["id"] if os.path.isabs(row["id"]) else os.path.join(ROOT, row["id"])
        if not os.path.exists(f): print("  없음:", f); continue
        mot = io.read_mot(f); scale = REF / (float(row["mass_kg"]) * float(row["height_m"]))
        for c in COORDS:
            for side_col in (c + "_moment", c.replace("_r", "_l") + "_moment"):
                if side_col in mot.names:
                    x = mot.column(side_col) * scale; per.setdefault((row["task"], c), []).append((x.max(), x.min()))
    tasks = sorted({k[0] for k in per}); out = []
    print(f"{'좌표':18s}" + "".join(f"{t:>22s}" for t in tasks) + "   |  전체 +피크 / −피크")
    for c in COORDS:
        line = f"{c:18s}"; allmax, allmin = [], []
        for t in tasks:
            v = per.get((t, c))
            if v: mx = max(a for a, _ in v); mn = min(b for _, b in v); line += f"{mx:9.1f} /{mn:8.1f}   "; allmax.append(mx); allmin.append(mn)
            else: line += f"{'-':>22s}"
            out.append(dict(task=t, coord=c, peak_pos=round(float(mx), 2) if v else "", peak_neg=round(float(mn), 2) if v else ""))
        if allmax: line += f"   |  {max(allmax):7.1f} / {min(allmin):7.1f}"
        print(line)
    with open(os.path.join(ROOT, "report/upper/torque_summary.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["task", "coord", "peak_pos", "peak_neg"]); w.writeheader(); w.writerows(out)
    fig, ax = plt.subplots(figsize=(11, 4.5)); x = np.arange(len(COORDS)); wdt = 0.8 / max(1, len(tasks))
    for i, t in enumerate(tasks):
        mx = [max((a for a, _ in per.get((t, c), [(0, 0)])), default=0) for c in COORDS]; mn = [min((b for _, b in per.get((t, c), [(0, 0)])), default=0) for c in COORDS]
        ax.bar(x + i * wdt - 0.4 + wdt / 2, mx, wdt, label=f"{t} (+)"); ax.bar(x + i * wdt - 0.4 + wdt / 2, mn, wdt, alpha=0.5)
    ax.axhline(0, color="k", lw=0.6); ax.set_xticks(x); ax.set_xticklabels(COORDS, rotation=30, ha="right", fontsize=8); ax.set_ylabel("N·m (피험자3 크기)"); ax.legend(fontsize=7, ncol=3)
    ax.set_title("상체 관절 토크 피크 — 과제별 (+ 위, − 아래 반투명)"); fig.tight_layout(); fig.savefig(os.path.join(ROOT, "report/upper/fig_torque_summary.png"), dpi=150); print("그림 저장")


if __name__ == "__main__":
    main()
