"""케이블 근육 8개로 움직이는 영상 — 뼈대 + 케이블(활성도로 색) + 지면반력 + 신호 막대 + 토크 추적.

    python report/paths_video.py --name paths_subject3_squats1 --grf <forces.mot> --id <id.sto> [--slow 2]
"""
import argparse, json, os, re, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.collections import LineCollection
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import addb_io as io, gait_video as gv, path_design as pd                       # noqa: E402
for f in font_manager.findSystemFonts():
    if "AppleSDGothicNeo" in f: plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
plt.rcParams["axes.unicode_minus"] = False
try:
    import imageio_ffmpeg; plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
except ImportError: pass
COORDS = ("hip_flexion", "hip_adduction", "knee_angle", "ankle_angle"); KL = {"hip_flexion": "고관절 굴곡", "hip_adduction": "고관절 내전", "knee_angle": "무릎", "ankle_angle": "발목"}
COLORS = ["#c0392b", "#2980b9", "#27ae60", "#e67e22", "#8e44ad", "#16a085", "#d35400", "#7f8c8d"]


def read_grf(path):
    m = io.read_mot(path); out = []
    for pre in sorted({re.match(r"(.*)_v[xyz]$", n).group(1) for n in m.names if re.match(r".*_v[xyz]$", n)}):
        ppre = pre.replace("force", "force") ; pp = [n for n in m.names if n.startswith(pre.rsplit("_", 1)[0]) and n.endswith("_px")]
        if not pp: continue
        pbase = pp[0][:-1]
        out.append((pre, m.time, np.column_stack([m.column(f"{pbase}{a}") for a in "xyz"]), np.column_stack([m.column(f"{pre}_v{a}") for a in "xyz"])))
    return out


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--name", required=True); ap.add_argument("--grf", required=True); ap.add_argument("--id", required=True)
    ap.add_argument("--title", default=""); ap.add_argument("--fps", type=int, default=30); ap.add_argument("--slow", type=float, default=2.0); ap.add_argument("--loops", type=int, default=2)
    a = ap.parse_args(); import opensim as osim; osim.Logger.setLevelString("Error")
    mdir = os.path.join(ROOT, "report/moco"); model_path = os.path.join(ROOT, "report/moco/newmus_paths_subject3_model.osim")
    sol = io.read_mot(os.path.join(mdir, f"{a.name}_solution.sto")); ik_path = os.path.join(mdir, f"{a.name}_ik_full.mot"); idm = io.read_mot(a.id)
    paths = json.load(open(os.path.join(ROOT, "report/newmus_paths_k8.json")))
    model = osim.Model(model_path); model.initSystem(); weight_n = io.model_weight_n(model)
    t0, t1 = float(sol.time[0]), float(sol.time[-1]); times = np.arange(t0, t1, 1.0 / (a.fps * a.slow)); n = len(times)
    U = np.column_stack([np.interp(times, sol.time, sol.column(f"/forceset/{m['name']}_r")) for m in paths["muscles"]])
    # 뼈대
    bodies = gv.collect_bodies(model); have = {model.getBodySet().get(i).getName() for i in range(model.getBodySet().getSize())}
    stations = [gv.TRUNK_TOP] if gv.TRUNK_TOP[0] in have else []
    pos = gv.sample_positions(model, ik_path, times, bodies, stations); segs = gv.segments(bodies + (["trunk_top"] if stations else []))
    # 케이블 점 (오른다리) — 자세별 뼈 변환
    ik = io.read_mot(ik_path); Q = np.column_stack([np.interp(times, ik.time, ik.column(c)) for c in list(pd.COORDS) + list(pd.EXTRA)])
    # 골반 자세도 필요: Frames 는 골반을 기준 자세로 두므로, 대신 모델 상태로 직접 station 위치를 잰다
    pose = io.Pose(model, ik_path); cab = []
    for m in paths["muscles"]:
        pts = np.zeros((n, len(m["chain"]), 3))
        for i, t in enumerate(times):
            pose.set(float(t))
            for k, (b, p) in enumerate(zip(m["chain"], m["points"])): pts[i, k] = pose.station_position(b, tuple(p))
        cab.append(pts)
    # 토크 (오른다리) = Σ F u r
    fr = pd.Frames(model_path, Q); tau = np.zeros((n, 4))
    for j, m in enumerate(paths["muscles"]):
        r, _ = pd.moment_arms(fr, [pd.BODIES.index(b) for b in m["chain"]], np.array(m["points"])); tau += (m["force_N"] * U[:, j])[:, None] * r
    grf = read_grf(a.grf)
    # 보조 액추에이터(새 근육이 못 낸 몫) — 오른다리 4좌표
    RES = np.column_stack([np.interp(times, sol.time, sol.column([k for k in sol.names if "reserve" in k and k.endswith(f"_{c}_r")][0])) for c in COORDS])
    fig = plt.figure(figsize=(14, 8.4)); gs = fig.add_gridspec(3, 2, width_ratios=[1.1, 1], height_ratios=[0.9, 1.1, 0.8], left=0.04, right=0.98, top=0.92, bottom=0.07, wspace=0.18, hspace=0.5)
    axf = fig.add_subplot(gs[:, 0]); axb = fig.add_subplot(gs[0, 1]); axt = fig.add_subplot(gs[1, 1]); axr = fig.add_subplot(gs[2, 1])
    fig.suptitle(f"케이블 근육 8개/다리 — {a.title or a.name}   (색 진하기 = 신호 u, MocoInverse, 운동학·지면반력 실측 고정)   1/{a.slow:g} 배속", fontsize=12)
    axf.set_aspect("equal"); axf.axhline(0, color="#7f8c8d", lw=1.5); axf.set_ylim(-0.12, float(pos[:, :, 1].max()) + 0.3); axf.set_yticks([0, 0.5, 1.0]); axf.set_xlabel("[m]")
    for sp in ("top", "right"): axf.spines[sp].set_visible(False)
    bones = LineCollection([], linewidths=5, colors="#bbbbbb", zorder=1); axf.add_collection(bones)
    cables = [axf.plot([], [], "-", lw=2.5, color=c, zorder=3)[0] for c in COLORS]
    arrows = [axf.annotate("", xy=(0, 0), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=gv.COLOR_GRF, lw=2)) for _ in grf]
    tl = axf.text(0.02, 0.97, "", transform=axf.transAxes, fontsize=11, va="top"); span = 1.6
    x = np.arange(8); bars = axb.bar(x, np.zeros(8), color=COLORS); axb.set_xticks(x); axb.set_xticklabels([f"M{i}" for i in x]); axb.set_ylim(0, 1.05)
    axb.set_title("이 순간 케이블 근육 8개의 신호 u (1 = 용량)", fontsize=11); axb.set_ylabel("u")
    for j, c in enumerate(COORDS):
        ln, = axt.plot(times, tau[:, j], lw=1.8, label=KL[c]); axt.plot(times, np.interp(times, idm.time, idm.column(f"{c}_r_moment")), ":", color=ln.get_color(), lw=1.2)
    cur = axt.axvline(t0, color="k", lw=1); axt.set_title("오른다리 관절 토크 — 케이블 근육이 낸 것(실선) vs 측정 ID(점선)", fontsize=11); axt.set_ylabel("Nm"); axt.legend(fontsize=8, ncol=4, loc="lower left"); axt.set_xlim(t0, t1)
    idpk = {c: float(np.abs(idm.column(f"{c}_r_moment")).max()) for c in COORDS}; worst = max(np.abs(RES[:, j]).max() / idpk[c] * 100 for j, c in enumerate(COORDS))
    for j, c in enumerate(COORDS): axr.plot(times, RES[:, j], lw=1.5, label=KL[c])
    cur2 = axr.axvline(t0, color="k", lw=1); lim = max(10.0, float(np.abs(RES).max()) * 1.2); axr.set_ylim(-lim, lim); axr.set_xlim(t0, t1)
    axr.set_title(f"보조 액추에이터가 떠맡은 토크 (케이블 근육이 못 낸 몫) — 최악 {worst:.0f} % of ID 피크", fontsize=11); axr.set_ylabel("Nm"); axr.set_xlabel("time [s]"); axr.legend(fontsize=8, ncol=4, loc="lower left")
    rtxt = axr.text(0.99, 0.95, "", transform=axr.transAxes, ha="right", va="top", fontsize=9)

    def update(fi):
        i = fi % n; p = pos[i]
        bones.set_segments([[(p[s0][0], p[s0][1]), (p[s1][0], p[s1][1])] for s0, s1, _ in segs])
        for j, ln in enumerate(cables): ln.set_data(cab[j][i, :, 0], cab[j][i, :, 1]); ln.set_alpha(0.15 + 0.85 * min(U[i, j], 1.0)); ln.set_linewidth(1.5 + 3 * min(U[i, j], 1.0))
        centre = float(p[0][0]); axf.set_xlim(centre - span / 2, centre + span / 2); tot = 0.0
        for arrow, (key, gt, cop, force) in zip(arrows, grf):
            f = np.array([np.interp(times[i], gt, force[:, k]) for k in range(3)]); c = np.array([np.interp(times[i], gt, cop[:, k]) for k in range(3)]); mag = np.linalg.norm(f); tot += mag
            if mag > gv.GRF_MIN_BW * weight_n: arrow.set_position((c[0], c[1])); arrow.xy = (c[0] + f[0] / weight_n * gv.GRF_SCALE, c[1] + f[1] / weight_n * gv.GRF_SCALE); arrow.set_visible(True)
            else: arrow.set_visible(False)
        tl.set_text(f"t = {times[i]:5.2f} s   vGRF = {tot / weight_n:4.2f} BW")
        for b, v in zip(bars, U[i]): b.set_height(v)
        cur.set_xdata([times[i], times[i]]); cur2.set_xdata([times[i], times[i]])
        rtxt.set_text("지금: " + "  ".join(f"{KL[c]} {RES[i, j]:+.1f}" for j, c in enumerate(COORDS)) + " Nm")
        return [bones, tl, cur, cur2, rtxt] + cables + list(bars) + arrows
    anim = FuncAnimation(fig, update, frames=n * a.loops, interval=1000 / a.fps, blit=False)
    out = os.path.join(ROOT, "report", f"video_{a.name}.mp4"); anim.save(out, writer=FFMpegWriter(fps=a.fps, bitrate=3000, extra_args=["-pix_fmt", "yuv420p"]), dpi=110)
    print(f"저장 {out}  ({os.path.getsize(out)/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
