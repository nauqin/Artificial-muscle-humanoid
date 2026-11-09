"""어깨 케이블이 팔을 뒤로 젖힌 자세에서 '앞으로 당기는 힘(굴곡 모멘트암)'을 잃는 것을 옆모습으로 보여 주는 그림·영상.

케이블 모델(.osim)의 실제 경로점(감는 구 포함)과 OpenSim 모멘트암을 그대로 쓴다. 옆모습 = 몸통 좌표의 x(앞)·y(위) 평면.
굴곡 축이 몸통 z 축이라, 이 평면에서 '케이블 직선과 어깨 중심 사이 거리'가 곧 굴곡 모멘트암이다.

    python report/plot_shoulder_flip.py --model report/moco/newmus_paths_shoulder_r_k6_40kg_h160_v635_s80_f15_model.osim \\
        --ik data/moore2015/subject3/T047_segment_12_ik.mot \\
        --id out/upper40_h160/moore2015/subject3/T047_segment_12/T047_segment_12_id.sto --name 40kg_h160_s80_f15
→ report/upper/fig_shoulder_flex_flip_<name>.png, report/moco/video_shoulder_flex_flip_<name>.mp4
"""
import argparse, os, re, shutil, sys
import numpy as np, opensim as osim
from scipy.optimize import lsq_linear
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib.patches import Circle, FancyBboxPatch
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import addb_io as io                                   # noqa: E402
plt.rcParams["font.family"] = ["Apple SD Gothic Neo", "AppleGothic"]; plt.rcParams["axes.unicode_minus"] = False
COORDS = ("arm_flex_r", "arm_add_r", "arm_rot_r"); CNAME = ("굴곡(앞뒤)", "내전(옆)", "회전(비틀기)")
POS, NEG, ZERO = "#1668b8", "#d1332e", "#8a8a8a"; EPS = 0.10      # cm: 이보다 작으면 '거의 0'


def sign_color(r_cm):
    return ZERO if abs(r_cm) < EPS else (POS if r_cm > 0 else NEG)


class Shoulder:
    def __init__(self, model, ik, idf):
        self.model = osim.Model(model); self.st = self.model.initSystem()
        self.cs = self.model.getCoordinateSet(); self.bs = self.model.getBodySet(); self.torso = self.bs.get("torso")
        self.ik = io.read_mot(ik); self.id = io.read_mot(idf); fs = self.model.getForceSet()
        self.names = sorted((fs.get(i).getName() for i in range(fs.getSize()) if re.fullmatch(r"M\d+_r", fs.get(i).getName())), key=lambda s: int(s[1:-2]))
        self.acts = [osim.PathActuator.safeDownCast(fs.get(n)) for n in self.names]
        self.F = np.array([a.getOptimalForce() for a in self.acts]); self.labels = [n[:-2] for n in self.names]
        wo = self.bs.get("humerus_r").getWrapObjectSet(); self.radius = osim.WrapSphere.safeDownCast(wo.get(0)).get_radius() * 100 if wo.getSize() else None

    def _in_torso(self, frame, v=(0.0, 0.0, 0.0)):
        p = frame.findStationLocationInAnotherFrame(self.st, osim.Vec3(*v), self.torso); return np.array([p.get(i) for i in range(3)])

    def at(self, t):
        st, cs = self.st, self.cs
        for n in self.ik.names:
            if cs.contains(n) and not cs.get(n).getDefaultLocked():
                v = float(np.interp(t, self.ik.time, self.ik.column(n)))
                if self.ik.in_degrees and cs.get(n).getMotionType() == osim.Coordinate.Rotational: v = np.radians(v)
                cs.get(n).setValue(st, v, False)
        self.model.assemble(st); self.model.realizePosition(st)
        c = self._in_torso(self.bs.get("humerus_r")); paths = []
        for a in self.acts:
            pts = a.getGeometryPath().getCurrentPath(st)
            paths.append(np.array([self._in_torso(pts.get(k).getParentFrame(), [pts.get(k).getLocation(st).get(j) for j in range(3)]) - c for k in range(pts.getSize())]) * 100)
        R = np.array([[a.getGeometryPath().computeMomentArm(st, cs.get(co)) for a in self.acts] for co in COORDS])
        tau = np.array([float(np.interp(t, self.id.time, self.id.column(co + "_moment"))) for co in COORDS])
        u = lsq_linear(R * self.F, tau, bounds=(0, 1)).x
        return dict(t=t, q=np.degrees([cs.get(co).getValue(st) for co in COORDS]), elbow=(self._in_torso(self.bs.get("ulna_r")) - c) * 100,
                    hand=(self._in_torso(self.bs.get("hand_r")) - c) * 100, paths=paths, R=R, tau=tau, u=u, ach=R * self.F @ u,
                    cap=(np.clip(R[0], 0, None) * self.F).sum())


def draw_side(ax, sh, d, title, big=True, mark=None):
    """옆모습(몸통 x-y). 원점 = 어깨 중심, 오른쪽 = 몸 앞. mark = 모멘트암을 검은 선으로 표시할 케이블 번호."""
    ax.clear(); ax.set_aspect("equal"); ax.set_xlim(-10.5, 13.5); ax.set_ylim(-19, 6.8)
    ax.add_patch(FancyBboxPatch((-9, -45), 20, 51, boxstyle="round,pad=0,rounding_size=5", fc="#f1f1f1", ec="#d5d5d5", lw=1, zorder=1))
    ax.text(-8.6, -18.6, "등 쪽", fontsize=10.5, color="#777", va="bottom"); ax.text(10.6, -18.6, "가슴 쪽", fontsize=10.5, color="#777", va="bottom", ha="right")
    arm = np.array([[0, 0], d["elbow"][:2], d["hand"][:2]])
    if not hasattr(ax, "_pt_per_cm"):                                                  # 위팔 굵기(7 cm)를 실제 축척으로
        ax.figure.canvas.draw(); ax._pt_per_cm = ax.get_window_extent().height / 25.8 * 72 / ax.figure.dpi
    ax.plot(arm[:, 0], arm[:, 1], color="#cdcdcd", lw=7.0 * ax._pt_per_cm, solid_capstyle="round", solid_joinstyle="round", zorder=2)
    if sh.radius: ax.add_patch(Circle((0, 0), sh.radius, fc="none", ec="#444", lw=1.1, ls=(0, (4, 3)), zorder=3))
    for i, (p, lab) in enumerate(zip(d["paths"], sh.labels)):
        r = d["R"][0, i] * 100; col = sign_color(r)
        ax.plot(p[:, 0], p[:, 1], color=col, lw=3.0 if abs(r) >= EPS else 2.0, solid_capstyle="round", zorder=5, alpha=0.95)
        ax.plot(*p[-1, :2], "o", ms=6, color=col, mec="white", mew=0.8, zorder=6, clip_on=True); ax.plot(*p[0, :2], "s", ms=5, color="#222", zorder=6)
    ax.plot(0, 0, "+", ms=17, mew=2.6, color="k", zorder=8)
    if mark is not None:                                                               # 어깨 중심에서 케이블 직선까지의 수직 거리 = 굴곡 모멘트암
        p = d["paths"][mark]; a, b = p[-2, :2], p[-1, :2]; v = (b - a) / max(np.linalg.norm(b - a), 1e-9); foot = a + v * float(-a @ v); r = d["R"][0, mark] * 100
        s0 = float(-a @ v)
        if s0 < 0: ax.plot([foot[0], a[0]], [foot[1], a[1]], color=sign_color(r), lw=1.1, ls=(0, (2, 2)), zorder=4)
        ax.plot([0, foot[0]], [0, foot[1]], color="k", lw=3.2, solid_capstyle="butt", zorder=7)
        ax.annotate(f"{sh.labels[mark]} 모멘트암 {abs(r):.1f} cm\n(어깨 중심 ~ 케이블 직선 거리)", xy=foot * 0.55, xytext=(3.2, 4.3), fontsize=10.5, fontweight="bold", va="center",
                    arrowprops=dict(arrowstyle="-", color="k", lw=0.9, shrinkA=2, shrinkB=1), zorder=9)
    x0, y0 = (-10.0, 6.3)
    ax.text(x0, y0, "굴곡 모멘트암", fontsize=10.5, va="top", color="#222", zorder=8)
    for i, lab in enumerate(sh.labels):
        r = d["R"][0, i] * 100
        ax.text(x0, y0 - 1.3 * (i + 1), f"{lab}  {r:+.1f} cm", fontsize=10.5, va="top", color=sign_color(r), fontweight="bold" if abs(r) >= EPS else "normal", zorder=8)
    ax.set_title(title, fontsize=14 if big else 12, fontweight="bold", pad=8); ax.set_xlabel("앞뒤 위치 (cm)   뒤 ←   → 앞"); ax.set_ylabel("높이 (cm)")
    ax.tick_params(labelsize=9); [s.set_visible(False) for s in (ax.spines["top"], ax.spines["right"])]


def legend_text(fig, y, x=0.05, fs=11.0, short=False):
    items = [("━ 앞으로 당기는 케이블", POS, True), ("━ 뒤로 당기는 케이블", NEG, True), ("━ 거의 0 (0.1 cm 미만)", ZERO, True)]
    if not short: items += [("✚ 어깨 중심", "k", True), ("■ 몸통 쪽 부착점   ● 위팔 쪽 부착점   점선 원 = 감는 구   굵은 회색 = 위팔", "#333", False)]
    for txt, col, bold in items:
        t = fig.text(x, y, txt, color=col, fontsize=fs, fontweight="bold" if bold else "normal")
        fig.canvas.draw(); x = t.get_window_extent().x1 / (fig.get_figwidth() * fig.dpi) + 0.018


def figure(sh, a, out):
    df, db = sh.at(a.t_front), sh.at(a.t_back); T = np.arange(a.plot0, a.plot1 + 1e-9, 0.01); S = [sh.at(t) for t in T]
    fig = plt.figure(figsize=(16, 11.6)); gs = fig.add_gridspec(4, 3, height_ratios=[3.5, 0.9, 1.25, 1.25], width_ratios=[1, 1, 0.95], hspace=0.42, wspace=0.22,
                                                                 left=0.05, right=0.985, top=0.875, bottom=0.05)
    fig.suptitle("팔을 뒤로 젖히면 어깨 케이블이 '앞으로 당기는 힘'을 잃는다", fontsize=19, fontweight="bold", y=0.992)
    fig.text(0.5, 0.943, "오른쪽 어깨를 옆에서 본 모습 · 로봇 160 cm·40 kg · 어깨 케이블 6개(스트로크 상한 설계) · 걷기 Moore T047_12", ha="center", fontsize=11.5, color="#444")
    km = int(np.argmax(df["R"][0])); legend_text(fig, 0.915)
    draw_side(fig.add_subplot(gs[0, 0]), sh, df, f"① 팔이 앞에 있을 때  ({df['q'][0]:+.0f}°, {a.t_front:.2f} s)", mark=km)
    draw_side(fig.add_subplot(gs[0, 1]), sh, db, f"② 팔을 뒤로 젖혔을 때  ({db['q'][0]:+.0f}°, {a.t_back:.2f} s)", mark=km)
    sub = gs[0, 2].subgridspec(2, 1, height_ratios=[1.5, 1], hspace=0.95)
    ax = fig.add_subplot(sub[0]); n = len(sh.labels); y = np.arange(n)[::-1]
    rf, rb = df["R"][0] * 100, db["R"][0] * 100
    ax.barh(y + 0.19, rf, 0.36, color=[sign_color(r) for r in rf], alpha=0.45, label="① 팔 앞")
    ax.barh(y - 0.19, rb, 0.36, color=[sign_color(r) for r in rb], label="② 팔 뒤")
    for yi, r1, r2 in zip(y, rf, rb):
        ax.text(r1 + (0.08 if r1 >= 0 else -0.08), yi + 0.19, f"{r1:+.1f}", va="center", ha="left" if r1 >= 0 else "right", fontsize=9, color="#555")
        ax.text(r2 + (0.08 if r2 >= 0 else -0.08), yi - 0.19, f"{r2:+.1f}", va="center", ha="left" if r2 >= 0 else "right", fontsize=9.5, fontweight="bold")
    ax.axvline(0, color="k", lw=1); ax.set_yticks(y); ax.set_yticklabels(sh.labels, fontsize=11); ax.set_xlim(-2.4, 3.4)
    ax.set_title("케이블별 굴곡 모멘트암 (cm)\n연한 막대 = ① 팔 앞, 진한 막대 = ② 팔 뒤", fontsize=11.5, fontweight="bold"); ax.set_xlabel("← 뒤로 당김      앞으로 당김 →", fontsize=10)
    [s.set_visible(False) for s in (ax.spines["top"], ax.spines["right"])]
    ax = fig.add_subplot(sub[1]); k = int(np.argmax(db["R"][0])); rr = db["R"][:, k] * 100
    ax.barh([2, 1, 0], rr, 0.55, color=[POS, "#6b6b6b", "#6b6b6b"])
    for yi, r in zip([2, 1, 0], rr): ax.text(r + (0.1 if r >= 0 else -0.1), yi, f"{r:+.1f}", va="center", ha="left" if r >= 0 else "right", fontsize=10.5, fontweight="bold")
    ax.axvline(0, color="k", lw=1); ax.set_yticks([2, 1, 0]); ax.set_yticklabels(CNAME, fontsize=10.5); ax.set_xlim(-4.6, 2.2)
    ax.set_title(f"② 에서 혼자 남은 {sh.labels[k]} 의 세 방향 모멘트암 (cm)\n앞으로보다 옆·비틀기로 {abs(rr[1] / rr[0]):.1f}배·{abs(rr[2] / rr[0]):.1f}배 세게 끈다", fontsize=11, fontweight="bold")
    [s.set_visible(False) for s in (ax.spines["top"], ax.spines["right"])]
    q = np.array([s["q"][0] for s in S]); Rf = np.array([s["R"][0] for s in S]) * 100; need = np.array([s["tau"][0] for s in S]); ach = np.array([s["ach"][0] for s in S])
    back = q < -17
    ax1 = fig.add_subplot(gs[1, :]); ax1.plot(T, q, color="k", lw=1.8); ax1.axhline(0, color="#999", lw=0.8)
    ax1.fill_between(T, -30, 35, where=back, color=NEG, alpha=0.13, lw=0); ax1.set_ylim(-30, 35); ax1.set_ylabel("팔 굴곡 각도 (°)\n+ 앞 / - 뒤", fontsize=10)
    ax1.set_title("걷는 동안의 변화 (붉은 띠 = 팔이 뒤로 17° 넘게 젖혀진 구간, 세로 점선 = 위 그림 ①·② 의 시각)", fontsize=12, fontweight="bold", loc="left")
    ax2 = fig.add_subplot(gs[2, :], sharex=ax1); pal = ["#b0b0b0", "#0b3d91", "#1a9e8f", "#c9c9c9", "#8e44ad", "#9a9a9a"]
    flex = [i for i in range(n) if Rf[:, i].max() > 0.5]
    for i in range(n):
        ax2.plot(T, Rf[:, i], color=pal[i % len(pal)], lw=2.4 if i in flex else 1.0, label=sh.labels[i] + (" (굴곡 담당)" if i in flex else ""), zorder=3 if i in flex else 2)
    ax2.axhline(0, color="k", lw=1.2); ax2.fill_between(T, -2.2, 3.2, where=back, color=NEG, alpha=0.13, lw=0); ax2.set_ylim(-2.2, 3.2)
    ax2.set_ylabel("굴곡 모멘트암 (cm)\n0 아래 = 뒤로 당김", fontsize=10); ax2.legend(ncol=n, fontsize=9.5, loc="upper right", frameon=False, bbox_to_anchor=(1.0, 1.22))
    ax3 = fig.add_subplot(gs[3, :], sharex=ax1); ax3.plot(T, need, color="k", lw=1.8, label="필요한 굴곡 토크 (역동역학)"); ax3.plot(T, ach, color=POS, lw=1.8, label="케이블 6개로 낼 수 있는 값 (세 방향을 같이 맞출 때)")
    ax3.fill_between(T, ach, need, where=(need - ach) > 0.05, color=NEG, alpha=0.55, lw=0, label="못 내는 양"); ax3.axhline(0, color="#999", lw=0.8)
    ax3.fill_between(T, -5, 5.5, where=back, color=NEG, alpha=0.13, lw=0); ax3.set_ylim(-5, 5.5); ax3.set_ylabel("굴곡 토크 (N·m)\n+ 앞으로", fontsize=10); ax3.set_xlabel("시간 (s)")
    ax3.legend(ncol=3, fontsize=9.5, loc="upper right", frameon=False, bbox_to_anchor=(1.0, 1.22)); ax1.set_xlim(T[0], T[-1])
    for ax in (ax1, ax2, ax3):
        [s.set_visible(False) for s in (ax.spines["top"], ax.spines["right"])]; ax.tick_params(labelsize=9)
        for tt, lab in ((a.t_front, "①"), (a.t_back, "②")):
            ax.axvline(tt, color="#333", lw=0.9, ls="--")
            if ax is ax1: ax.text(tt + 0.012, 31, lab, ha="left", va="top", fontsize=13, fontweight="bold")
    fig.savefig(out, dpi=150); plt.close(fig); print("그림", out)
    return df, db, T, S


def video(sh, a, out, T, S):
    fps, slow = 30, a.slow; tt = np.arange(a.t0, a.t1 + 1e-9, 1.0 / (fps * slow)); q_all = np.array([s["q"][0] for s in S])
    holds = {}                                                                         # 팔이 가장 뒤로 간 순간마다 멈춤
    for i in range(1, len(T) - 1):
        if a.t0 < T[i] < a.t1 and q_all[i] < -17 and q_all[i] <= q_all[i - 1] and q_all[i] < q_all[i + 1]: holds[int(np.argmin(abs(tt - T[i])))] = int(2.6 * fps)
    seq = [0] * fps
    for i in range(len(tt)): seq += [i] * (1 + holds.get(i, 0))
    D = [sh.at(t) for t in tt]; need = np.array([s["tau"][0] for s in S]); ach = np.array([s["ach"][0] for s in S]); n = len(sh.labels)
    fig = plt.figure(figsize=(12.8, 7.2), dpi=100); gs = fig.add_gridspec(2, 2, width_ratios=[1.12, 1], height_ratios=[1.15, 1], left=0.06, right=0.975, top=0.86, bottom=0.085, hspace=0.5, wspace=0.2)
    axs, axb, axt = fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])
    fig.suptitle("팔을 뒤로 젖히면 어깨 케이블이 '앞으로 당기는 힘'을 잃는다", fontsize=17, fontweight="bold", y=0.975)
    legend_text(fig, 0.905, x=0.06, fs=10.5, short=True); cap = fig.text(0.975, 0.905, "", ha="right", fontsize=11.5, fontweight="bold", color="#222")
    km = int(np.argmax(sh.at(a.t_front)["R"][0]))
    win = (T >= a.t0) & (T <= a.t1)
    axt.plot(T[win], need[win], color="k", lw=1.8, label="필요한 굴곡 토크"); axt.plot(T[win], ach[win], color=POS, lw=1.8, label="케이블로 낼 수 있는 값")
    axt.fill_between(T[win], ach[win], need[win], where=(need - ach)[win] > 0.05, color=NEG, alpha=0.55, lw=0, label="못 내는 양"); axt.axhline(0, color="#999", lw=0.8)
    axt.set_xlim(a.t0, a.t1); axt.set_ylim(-5, 6.6); axt.set_xlabel("시간 (s)"); axt.set_ylabel("굴곡 토크 (N·m)"); axt.legend(ncol=3, fontsize=9, frameon=False, loc="upper center", columnspacing=1.0, handlelength=1.4)
    axt.set_title("필요한 토크와 낼 수 있는 토크", fontsize=12, fontweight="bold", loc="left"); cur = axt.axvline(a.t0, color="#333", lw=1.2, ls="--")
    [s.set_visible(False) for s in (axt.spines["top"], axt.spines["right"])]

    def frame(j):
        d = D[seq[j]]; held = j > 0 and seq[j] == seq[j - 1] and seq[j] in holds
        draw_side(axs, sh, d, f"옆에서 본 오른쪽 어깨   팔 각도 {d['q'][0]:+.0f}°   ({d['t']:.2f} s, {slow}배 느리게)", big=False, mark=km)
        axb.clear(); r = d["R"][0] * 100; y = np.arange(n)[::-1]; axb.barh(y, r, 0.6, color=[sign_color(x) for x in r])
        for yi, x in zip(y, r): axb.text(x + (0.08 if x >= 0 else -0.08), yi, f"{x:+.1f}", va="center", ha="left" if x >= 0 else "right", fontsize=10.5, fontweight="bold")
        axb.axvline(0, color="k", lw=1); axb.set_yticks(y); axb.set_yticklabels(sh.labels, fontsize=11); axb.set_xlim(-2.4, 3.4); axb.set_xlabel("← 뒤로 당김      앞으로 당김 →", fontsize=10)
        axb.set_title("케이블별 굴곡 모멘트암 (cm)", fontsize=12, fontweight="bold", loc="left"); [s.set_visible(False) for s in (axb.spines["top"], axb.spines["right"])]
        cur.set_xdata([d["t"], d["t"]]); pos = [sh.labels[i] for i in range(n) if r[i] >= EPS]
        if held: cap.set_text(f"멈춤: 앞으로 당기는 케이블은 {'·'.join(pos) or '없음'} 뿐 — 필요 {d['tau'][0]:+.1f}, 낼 수 있는 값 {d['ach'][0]:+.1f} N·m"); cap.set_color(NEG)
        else: cap.set_text(f"앞으로 당기는 케이블 {len(pos)}개 ({'·'.join(pos) or '없음'})"); cap.set_color("#222")
        return []

    ff = shutil.which("ffmpeg")
    if not ff:
        import imageio_ffmpeg; ff = imageio_ffmpeg.get_ffmpeg_exe()
    plt.rcParams["animation.ffmpeg_path"] = ff
    animation.FuncAnimation(fig, frame, frames=len(seq), blit=False).save(out, writer=animation.FFMpegWriter(fps=fps, bitrate=2600, extra_args=["-pix_fmt", "yuv420p"]))
    plt.close(fig); print("영상", out, f"({len(seq) / fps:.1f} s, 멈춤 {len(holds)}회)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True); ap.add_argument("--ik", required=True); ap.add_argument("--id", required=True); ap.add_argument("--name", required=True)
    ap.add_argument("--t-front", type=float, default=1.00); ap.add_argument("--t-back", type=float, default=1.55)
    ap.add_argument("--plot0", type=float, default=1.0); ap.add_argument("--plot1", type=float, default=4.0)
    ap.add_argument("--t0", type=float, default=1.0, help="영상 시작 (s)"); ap.add_argument("--t1", type=float, default=2.6); ap.add_argument("--slow", type=int, default=4)
    ap.add_argument("--no-video", action="store_true"); a = ap.parse_args()
    p = lambda s: s if os.path.isabs(s) else os.path.join(ROOT, s)
    sh = Shoulder(p(a.model), p(a.ik), p(a.id))
    df, db, T, S = figure(sh, a, os.path.join(ROOT, f"report/upper/fig_shoulder_flex_flip_{a.name}.png"))
    for d in (df, db):
        print(f"t={d['t']:.2f} 굴곡 {d['q'][0]:+.1f}°  모멘트암(cm) " + " ".join(f"{x * 100:+.2f}" for x in d["R"][0]) + f"  앞으로 최대 {d['cap']:.1f} N·m  필요 {d['tau'][0]:+.2f}  낼 수 있는 값 {d['ach'][0]:+.2f}")
    if not a.no_video: video(sh, a, os.path.join(ROOT, f"report/moco/video_shoulder_flex_flip_{a.name}.mp4"), T, S)
