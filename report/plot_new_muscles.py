"""새 근육 k개 MocoInverse 결과 그림 — 신호 u_i(t), 새 근육 토크(= ID − 보조), 보조 토크."""
import os, sys
import numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import addb_io as io
for f in font_manager.findSystemFonts():
    if "AppleSDGothicNeo" in f: plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
plt.rcParams["axes.unicode_minus"] = False
COORDS = ["hip_flexion_r", "hip_adduction_r", "knee_angle_r", "ankle_angle_r"]

def main(tag="_w0", ks=(8, 5), umax=3.0, out="fig_newmuscles_walk.png"):
    idm = io.read_mot(os.path.join(ROOT, "out/real/subject3/walking1_segment_0/walking1_segment_0_id.sto"))
    fig, axes = plt.subplots(len(ks), 3, figsize=(15, 3.8 * len(ks)), squeeze=False)
    for row, k in enumerate(ks):
        sol = io.read_mot(os.path.join(ROOT, f"report/moco/newmus_subject3_walking1_k{k}{tag}_solution.sto")); t = sol.time
        ax = axes[row][0]
        for n in [n for n in sol.names if "newmuscles_r" in n]:
            ax.plot(t, sol.column(n), label="S" + n.split("_")[-1])
        ax.axhline(umax, color="0.5", ls="--", lw=0.8); ax.text(t[0], umax * 0.97, f"한도 {umax:g}", fontsize=8, va="top", color="0.4")
        ax.set_title(f"새 근육 {k}개 (오른쪽) — 신호 u_i(t)"); ax.set_ylabel("u (1 = 그 시너지의 피크)"); ax.legend(fontsize=8, ncol=2); ax.set_ylim(0, umax * 1.05)
        ax = axes[row][1]; worst = 0
        for c in COORDS:
            res = sol.column([n for n in sol.names if "reserve" in n and n.endswith("_" + c)][0])
            idc = np.interp(t, idm.time, idm.column(c + "_moment"))
            ln, = ax.plot(t, idc - res, label=c.replace("_r", "")); ax.plot(t, idc, ":", color=ln.get_color(), lw=1)
            worst = max(worst, np.abs(res).max() / np.abs(idc).max() * 100)
        ax.set_title(f"새 근육 {k}개가 낸 토크(실선) vs ID 모멘트(점선)"); ax.set_ylabel("Nm"); ax.legend(fontsize=8)
        ax = axes[row][2]
        for c in COORDS:
            ax.plot(t, sol.column([n for n in sol.names if "reserve" in n and n.endswith("_" + c)][0]), label=c.replace("_r", ""))
        ax.set_title(f"보조 액추에이터가 떠맡은 토크 — 최악 {worst:.0f} % of ID 피크"); ax.set_ylabel("Nm"); ax.set_ylim(-30, 30); ax.legend(fontsize=8)
    for ax in axes[-1]: ax.set_xlabel("time [s]")
    fig.suptitle("subject3 walking1 — 사람 근육 80개를 떼고 새 근육 k개(토크 비율 고정)로 같은 걸음 재현 (MocoInverse, 운동학·GRF 실측 고정)", fontsize=12)
    fig.tight_layout(); p = os.path.join(ROOT, "report", out); fig.savefig(p, dpi=150); print("저장", p)

if __name__ == "__main__":
    tag = sys.argv[1] if len(sys.argv) > 1 else "_w0"; umax = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
    main(tag=tag, umax=umax, out=f"fig_newmuscles_walk{tag}.png")
