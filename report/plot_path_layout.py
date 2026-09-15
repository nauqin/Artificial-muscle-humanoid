"""케이블 경로를 자세 3개(서기·다리 흔들기·깊은 스쿼트)에서 옆에서 본 그림 — 최적화 결과는 눈으로 확인한다.
사용: python report/plot_path_layout.py  → report/fig_path_layout.png"""
import os, sys, json
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import path_design as pd

def main():
    plt.rcParams["axes.unicode_minus"] = False                # 한글 글꼴은 path_design 이 import 될 때 이미 잡혀 있다
    out = json.load(open(os.path.join(ROOT, "report/newmus_paths_k8.json")))
    Q, T = pd.poses()
    hip, knee = Q[:, 0], Q[:, 2]
    idx = [int(np.argmin(np.abs(hip) + np.abs(knee))), int(np.argmax(knee[:60] if len(knee) > 60 else knee)), int(np.argmax(knee))]
    names = ["서기 (걷기 중 가장 곧은 자세)", "걷기 — 다리 흔들기 (무릎 최대)", "깊은 스쿼트 (무릎 최대)"]
    fr = pd.Frames(os.path.join(ROOT, "data/subject3/subject3_scaled.osim"), Q[idx])
    colors = ["#e53935", "#3b82f6", "#22a04a", "#f39c12", "#8e44ad", "#17a2b8", "#d63384", "#8a8a1a"]
    fig, axes = plt.subplots(1, 3, figsize=(15, 6.5))
    for k, ax in enumerate(axes):
        P0 = fr.p[k, 0]                                              # 바디 원점 (골반, 대퇴=고관절, 정강이=무릎, 발꿈치)
        toe = fr.R[k, 0, 3] @ np.array([0.2, 0.0, 0.0]) + P0[3]
        bone = np.vstack([P0[0], P0[1], P0[2], P0[3], toe])
        ax.plot(bone[:, 0], bone[:, 1], "-", color="0.55", lw=6, solid_capstyle="round", zorder=1)
        ax.plot(bone[1:4, 0], bone[1:4, 1], "o", color="0.3", ms=7, zorder=2)
        for i, m in enumerate(out["muscles"]):
            chain = [pd.BODIES.index(b) for b in m["chain"]]
            pts = np.array([fr.R[k, 0, b] @ np.array(p) + fr.p[k, 0, b] for b, p in zip(chain, m["points"])])
            ax.plot(pts[:, 0], pts[:, 1], "-", color=colors[i], lw=2, label=f"{m['name']} {m['force_N']:.0f} N" if k == 0 else None, zorder=3)
            ax.plot(pts[:, 0], pts[:, 1], ".", color=colors[i], ms=6, zorder=4)
        ax.set_aspect("equal"); ax.set_title(f"{names[k]}\n고관절 {np.degrees(Q[idx[k],0]):.0f}°  무릎 {np.degrees(Q[idx[k],2]):.0f}°  발목 {np.degrees(Q[idx[k],3]):.0f}°", fontsize=10)
        ax.set_xlabel("앞 (m)"); ax.grid(alpha=0.3)
        if k == 0: ax.set_ylabel("위 (m)"); ax.legend(fontsize=7, loc="lower left")
    fig.suptitle("케이블 경로 8개 — 옆에서 본 모습 (회색 = 뼈: 골반–고관절–무릎–발꿈치–발끝, 점 = 부착·경유점)", fontsize=11)
    fig.tight_layout(); p = os.path.join(ROOT, "report/fig_path_layout.png"); fig.savefig(p, dpi=150); print("저장", p)

if __name__ == "__main__":
    main()
