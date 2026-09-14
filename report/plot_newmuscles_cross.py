"""새 근육 8개 — 방향 V 를 어디서 만들었나에 따른 교차 검증 그림 (results_newmus.csv → fig_newmuscles_cross.png)."""
import csv, os
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
for f in font_manager.findSystemFonts():
    if "AppleSDGothicNeo" in f: plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
plt.rcParams["axes.unicode_minus"] = False

rows = {r[0]: float(r[5]) for r in csv.reader(open(os.path.join(ROOT, "report/moco/results_newmus.csv"))) if len(r) > 5}
TRIALS = [("subject3", "walking1"), ("subject3", "walking2"), ("subject3", "walking3"), ("subject3", "walkingTS2"),
          ("subject3", "walkingTS3"), ("subject3", "walkingTS4"), ("subject4", "walking1")]
SERIES = [  # (라벨, 태그 생성함수, 색)
    ("V = 자기 걸음 (보조 가중치 1)",        lambda s, t: f"newmus_{s}_{t}_k8_self" if not (s == "subject3" and t == "walking1") else "newmus_subject3_walking1_k8_w0_u10", "#bdc3c7"),
    ("V = subject3 walking1 고정 (가중치 1)", lambda s, t: f"newmus_{s}_{t}_k8_vW1", "#e59866"),
    ("V = subject3 walking1 고정 (가중치 100)", lambda s, t: f"newmus_{s}_{t}_k8_vW1_rw100" if not (s == "subject3" and t == "walking1") else "newmus_subject3_walking1_k8_self_rw100", "#d35400"),
    ("V = subject3 6걸음 합침 (가중치 100)",   lambda s, t: f"newmus_{s}_{t}_k8_{'vALL' if s == 'subject3' else 'vS3ALL'}_rw100", "#1e8449"),
]
fig, ax = plt.subplots(figsize=(12, 4.6))
x = np.arange(len(TRIALS)); w = 0.2
for i, (lab, fn, col) in enumerate(SERIES):
    vals = [rows.get(fn(s, t), np.nan) for s, t in TRIALS]
    b = ax.bar(x + (i - 1.5) * w, vals, w, label=lab, color=col)
    for xi, v in zip(x + (i - 1.5) * w, vals):
        if not np.isnan(v): ax.text(xi, v + 0.5, f"{v:.0f}" if v >= 1 else f"{v:.1f}", ha="center", fontsize=7)
ax.axhline(10, color="k", ls="--", lw=0.8); ax.text(len(TRIALS) - 0.5, 10.5, "10 % (재현 기준)", ha="right", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels([f"{s}\n{t}" for s, t in TRIALS], fontsize=9)
ax.set_ylabel("보조 액추에이터 최악 (% of ID 피크)"); ax.set_ylim(0, 40)
ax.set_title("새 근육 8개/다리 — 근육 방향 V 를 고정하고 다른 걸음을 걷게 했을 때 (MocoInverse, u ≤ 10)", fontsize=11)
ax.legend(fontsize=8, loc="upper left")
ax.text(0.99, 0.6, "k=5 (합친 V, 가중치 100): 97 %", transform=ax.transAxes, ha="right", fontsize=9, color="#922b21")
fig.tight_layout(); p = os.path.join(ROOT, "report/fig_newmuscles_cross.png"); fig.savefig(p, dpi=150); print("저장", p)
