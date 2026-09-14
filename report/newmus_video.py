"""새 근육 8개/다리로 걷는 모습을 영상으로 — OpenSim GUI 없이.

왼쪽: 하지 스틱피겨(Moco 에 준 운동학) + 지면반력. 오른쪽 위: 이 순간 새 근육 8개의 신호
(막대, 오른다리 빨강·왼다리 파랑). 오른쪽 아래: 새 근육이 낸 토크(실선) vs ID(점선), 커서.

    python report/newmus_video.py --k 8 --tag _w0_u10 --slow 3 --loops 3
"""
import argparse, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
from matplotlib.collections import LineCollection

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import addb_io as io                                  # noqa: E402
import gait_video as gv                               # noqa: E402

for f in font_manager.findSystemFonts():
    if "AppleSDGothicNeo" in f:
        plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
plt.rcParams["axes.unicode_minus"] = False
try:
    import imageio_ffmpeg
    plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
except ImportError:
    pass

COORDS = ("hip_flexion", "hip_adduction", "knee_angle", "ankle_angle")
KLABEL = {"hip_flexion": "고관절 굴곡", "hip_adduction": "고관절 내전", "knee_angle": "무릎", "ankle_angle": "발목"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--subject", default="subject3"); ap.add_argument("--trial", default="walking1")
    ap.add_argument("--k", type=int, default=8); ap.add_argument("--tag", default="_w0_u10")
    ap.add_argument("--umax", type=float, default=10.0)
    ap.add_argument("--fps", type=int, default=30); ap.add_argument("--slow", type=float, default=3.0, help="배속 분모 (3 = 1/3 배속)")
    ap.add_argument("--loops", type=int, default=3); ap.add_argument("--out", default=None)
    ap.add_argument("--gif", action="store_true", help="mp4 대신 gif")
    a = ap.parse_args()
    import opensim as osim; osim.Logger.setLevelString("Error")

    name = f"newmus_{a.subject}_{a.trial}_k{a.k}{a.tag}"
    mdir = os.path.join(ROOT, "report", "moco")
    model = osim.Model(os.path.join(mdir, f"{name}_model.osim")); model.initSystem()
    weight_n = io.model_weight_n(model)
    ik_path = os.path.join(mdir, f"{name}_ik_full.mot")
    sol = io.read_mot(os.path.join(mdir, f"{name}_solution.sto"))
    idm = io.read_mot(os.path.join(ROOT, "out/real", a.subject, f"{a.trial}_segment_0", f"{a.trial}_segment_0_id.sto"))
    grf_mot = io.read_mot(os.path.join(ROOT, "data", a.subject, f"{a.trial}_segment_0_grf.mot"))

    t0, t1 = float(sol.time[0]), float(sol.time[-1])
    times = np.arange(t0, t1, 1.0 / (a.fps * a.slow))       # 느리게: 실제 1초를 slow초에
    n = len(times)

    # 신호 u_i(t), 토크
    U = {s: np.column_stack([np.interp(times, sol.time, sol.column(c)) for c in sol.names if f"newmuscles_{s}" in c]) for s in "rl"}
    tau_id, tau_new = {}, {}
    worst = 0.0
    for c in COORDS:
        for s in "rl":
            coord = f"{c}_{s}"
            res = np.interp(times, sol.time, sol.column([m for m in sol.names if "reserve" in m and m.endswith("_" + coord)][0]))
            idc = np.interp(times, idm.time, idm.column(coord + "_moment"))
            tau_id[coord], tau_new[coord] = idc, idc - res
            res_full = sol.column([m for m in sol.names if "reserve" in m and m.endswith("_" + coord)][0])   # 해의 원 표본으로 최악값
            worst = max(worst, np.abs(res_full).max() / np.abs(idm.column(coord + "_moment")).max() * 100)

    # 스틱피겨
    bodies = gv.collect_bodies(model)
    have = {model.getBodySet().get(i).getName() for i in range(model.getBodySet().getSize())}
    stations = [gv.TRUNK_TOP] if gv.TRUNK_TOP[0] in have else []
    pos = gv.sample_positions(model, ik_path, times, bodies, stations)
    labels = bodies + (["trunk_top"] if stations else [])
    segs = gv.segments(labels)
    plates = {k: p for k, p in io.find_plates(grf_mot).items() if p.complete}
    grf = []
    for key, plate in plates.items():
        cop = np.column_stack([np.interp(times, grf_mot.time, grf_mot.column(plate.point[ax])) for ax in "xyz"])
        force = np.column_stack([np.interp(times, grf_mot.time, grf_mot.column(plate.force[ax])) for ax in "xyz"])
        grf.append((key, cop, force))

    fig = plt.figure(figsize=(14, 7.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.15, 1], height_ratios=[1, 1.15], left=0.04, right=0.98, top=0.9, bottom=0.08, wspace=0.18, hspace=0.42)
    ax_fig = fig.add_subplot(gs[:, 0]); ax_bar = fig.add_subplot(gs[0, 1]); ax_tau = fig.add_subplot(gs[1, 1])
    fig.suptitle(f"새 근육 {a.k}개/다리로 걷기 — {a.subject} {a.trial}   (사람 근육 80개 제거, MocoInverse, 보조 최악 {worst:.0f} %)   1/{a.slow:g} 배속", fontsize=13)

    ax_fig.set_aspect("equal"); ax_fig.axhline(0, color="#7f8c8d", lw=1.5, zorder=1)
    ax_fig.set_ylim(-0.12, float(pos[:, :, 1].max()) + 0.3); ax_fig.set_yticks([0, 0.5, 1.0]); ax_fig.set_xlabel("[m]")
    for sp in ("top", "right"): ax_fig.spines[sp].set_visible(False)
    span = 1.7
    lines = LineCollection([], linewidths=4, zorder=3); ax_fig.add_collection(lines)
    joints = ax_fig.plot([], [], "o", color="#2c3e50", ms=5, zorder=4)[0]
    arrows = [ax_fig.annotate("", xy=(0, 0), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=gv.COLOR_GRF, lw=2.2)) for _ in grf]
    tlabel = ax_fig.text(0.02, 0.97, "", transform=ax_fig.transAxes, fontsize=11, va="top")
    ax_fig.text(0.02, 0.05, "빨강 = 오른다리, 파랑 = 왼다리, 초록 = 지면반력", transform=ax_fig.transAxes, fontsize=9, color="0.35")

    # 막대: 오른다리 k개, 왼다리 k개
    x = np.arange(a.k); w = 0.4
    br = ax_bar.bar(x - w / 2, np.zeros(a.k), w, color=gv.COLOR_R, label="오른다리")
    bl = ax_bar.bar(x + w / 2, np.zeros(a.k), w, color=gv.COLOR_L, label="왼다리")
    ax_bar.set_xticks(x); ax_bar.set_xticklabels([f"S{i}" for i in x]); ax_bar.set_ylim(0, max(U["r"].max(), U["l"].max()) * 1.15)
    ax_bar.axhline(1.0, color="0.5", ls="--", lw=0.8); ax_bar.text(a.k - 0.5, 1.02, "시너지 피크(=1)", fontsize=8, ha="right", va="bottom", color="0.4")
    ax_bar.set_title(f"이 순간 새 근육 {a.k}개의 신호 u_i", fontsize=11); ax_bar.set_ylabel("u"); ax_bar.legend(fontsize=8, loc="upper right")

    # 토크: 오른다리 4좌표
    for c in COORDS:
        coord = f"{c}_r"
        ln, = ax_tau.plot(times, tau_new[coord], lw=1.8, label=KLABEL[c]); ax_tau.plot(times, tau_id[coord], ":", color=ln.get_color(), lw=1.2)
    cursor = ax_tau.axvline(t0, color="k", lw=1)
    ax_tau.set_title("오른다리 관절 토크 — 새 근육이 낸 것(실선) vs 측정 ID(점선)", fontsize=11); ax_tau.set_ylabel("Nm"); ax_tau.set_xlabel("time [s]")
    ax_tau.legend(fontsize=8, ncol=4, loc="lower left"); ax_tau.set_xlim(t0, t1)
    ax_tau.text(0.99, 0.97, "점선이 안 보이면 실선에 완전히 덮인 것", transform=ax_tau.transAxes, fontsize=8, ha="right", va="top", color="0.4")

    def update(fi):
        i = fi % n
        p = pos[i]
        lines.set_segments([[(p[s0][0], p[s0][1]), (p[s1][0], p[s1][1])] for s0, s1, _ in segs]); lines.set_color([c for _a, _b, c in segs])
        joints.set_data(p[:, 0], p[:, 1])
        centre = float(p[0][0]); ax_fig.set_xlim(centre - span / 2, centre + span / 2)
        text = f"t = {times[i]:5.2f} s   (반복 {fi // n + 1}/{a.loops})"
        total = 0.0
        for arrow, (key, cop, force) in zip(arrows, grf):
            mag = float(np.linalg.norm(force[i])); total += mag
            if mag > gv.GRF_MIN_BW * weight_n:
                arrow.set_position((cop[i][0], cop[i][1])); arrow.xy = (cop[i][0] + force[i][0] / weight_n * gv.GRF_SCALE, cop[i][1] + force[i][1] / weight_n * gv.GRF_SCALE); arrow.set_visible(True)
            else: arrow.set_visible(False)
        tlabel.set_text(text + f"   vGRF = {total / weight_n:4.2f} BW")
        for b, v in zip(br, U["r"][i]): b.set_height(v)
        for b, v in zip(bl, U["l"][i]): b.set_height(v)
        cursor.set_xdata([times[i], times[i]])
        return [lines, joints, tlabel, cursor] + list(br) + list(bl) + arrows

    anim = FuncAnimation(fig, update, frames=n * a.loops, interval=1000 / a.fps, blit=False)
    out = a.out or os.path.join(ROOT, "report", f"video_{name}.{'gif' if a.gif else 'mp4'}")
    print(f"{n} 프레임 × {a.loops} 반복, {a.fps} fps → {out}", flush=True)
    if out.endswith(".mp4"):
        anim.save(out, writer=FFMpegWriter(fps=a.fps, bitrate=3000, extra_args=["-pix_fmt", "yuv420p"]), dpi=110)
    else:
        anim.save(out, writer=PillowWriter(fps=a.fps), dpi=80)
    print(f"저장 {out}  ({os.path.getsize(out) / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
