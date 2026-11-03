"""상체 케이블 근육(몸통 + 어깨 + 팔꿈치)으로 움직이는 영상 — 뼈대 + 케이블(신호로 색) + 신호 막대 + 토크 추적 + 보조.

`paths_video.py`(하체)의 상체판. 케이블 위치는 `gui_export_upper.py` 가 만든 근육 모델에서 OpenSim 이 실제로 계산한 경로
(감는 구를 돌아가는 점 포함)를 읽는다. 신호·토크·보조는 설계 공간별 Moco 해에서.

    python report/paths_video_upper.py --name lift_s1_PFWJ2 --base data/carrt/subject1/subject1_load3.79.osim \\
        --paths trunk=report/newmus_paths_trunk_k4.json shoulder_r=report/newmus_paths_shoulder_r_k5_wrap.json elbow_r=report/newmus_paths_elbow_r_k2.json \\
        --sol trunk=paths_trunk_k4_lift_s1_PFWJ2 shoulder_r=paths_shoulder_k5wrap_lift_s1_PFWJ2 elbow_r=paths_elbow_k2_lift_s1_PFWJ2 \\
        --id out/upper/carrt/subject1/PFWJ2/PFWJ2_id.sto --title "찬 물병 3.79 kg 들기" --azimuth 90 --slow 1.5
→ report/video_upper_<name>.mp4   (`--png 2.0` 이면 그 시각 한 장만 report/frame_upper_<name>.png 으로)
--azimuth: 보는 방향. 0 = 월드 x 축이 화면 가로(옆), 90 = z 축이 가로(앞). 피험자가 어디를 보느냐에 따라 맞춘다.
"""
import argparse, json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.animation import FuncAnimation, FFMpegWriter
from matplotlib.collections import LineCollection
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import addb_io as io, gait_video as gv, path_design as pd, design_space as ds      # noqa: E402
for f in font_manager.findSystemFonts():
    if "AppleSDGothicNeo" in f: plt.rcParams["font.family"] = font_manager.FontProperties(fname=f).get_name(); break
plt.rcParams["axes.unicode_minus"] = False
try:
    import imageio_ffmpeg; plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
except ImportError: pass
GROUP_COLORS = {"trunk": "#c0392b", "shoulder_r": "#2980b9", "elbow_r": "#27ae60", "lower": "#8e44ad"}
KL = {"lumbar_extension": "요추 신전", "lumbar_bending": "요추 옆굽힘", "lumbar_rotation": "요추 회전",
      "arm_flex_r": "어깨 굴곡", "arm_add_r": "어깨 내전", "arm_rot_r": "어깨 회전", "elbow_flex_r": "팔꿈치",
      "hip_flexion_r": "고관절 굴곡", "hip_adduction_r": "고관절 내전", "knee_angle_r": "무릎", "ankle_angle_r": "발목"}
ARM_CHAIN = (("humerus", "ulna"), ("ulna", "hand"))


def kv(items): return dict(s.split("=", 1) for s in items)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True); ap.add_argument("--base", required=True); ap.add_argument("--paths", nargs="+", required=True); ap.add_argument("--sol", nargs="+", required=True)
    ap.add_argument("--id", required=True); ap.add_argument("--id-lower", default=None, help="하체 ID .sto (하체 공간이 있을 때, 상체 ID 와 다른 파일이면)"); ap.add_argument("--grf", default=None, help="지면반력 .mot (있으면 화살표)")
    ap.add_argument("--title", default=""); ap.add_argument("--azimuth", default="auto", help="보는 방향(도). auto = 골반이 향한 방향을 화면 가로로 (옆모습), 거기에 90 을 더하면 앞모습")
    ap.add_argument("--fps", type=int, default=30); ap.add_argument("--slow", type=float, default=1.5); ap.add_argument("--loops", type=int, default=1); ap.add_argument("--png", type=float, default=None); ap.add_argument("--tag", default="", help="출력 파일 이름 뒤에 붙일 말 (예: _side)"); ap.add_argument("--box-scale", type=float, default=1.0, help="감는 구 배율 (로봇 크기)")
    a = ap.parse_args(); import opensim as osim; osim.Logger.setLevelString("Error")
    mdir = os.path.join(ROOT, "report/moco"); paths_j = {k: (v if os.path.isabs(v) else os.path.join(ROOT, v)) for k, v in kv(a.paths).items()}; sols = kv(a.sol)
    gui_model = os.path.join(mdir, f"gui_upper_{a.name}_model.osim")
    if not os.path.exists(gui_model): raise SystemExit(f"{gui_model} 없음 — 먼저 report/gui_export_upper.py 를 돌린다")
    first = next(iter(sols.values())); ik_path = os.path.join(mdir, f"{first}_ik_full.mot"); idm = io.read_mot(a.id); idm_lower = io.read_mot(a.id_lower) if a.id_lower else idm
    sol = {sp: io.read_mot(os.path.join(mdir, f"{n}_solution.sto")) for sp, n in sols.items()}
    t0 = max(float(s.time[0]) for s in sol.values()); t1 = min(float(s.time[-1]) for s in sol.values())
    times = np.arange(t0, t1, 1.0 / (a.fps * a.slow)); n = len(times)
    if str(a.azimuth).startswith("auto"):                                      # 골반 x 축(앞)을 지면에 내린 방향 = 옆모습의 가로축
        mdl = osim.Model(a.base if os.path.isabs(a.base) else os.path.join(ROOT, a.base)); ps = io.Pose(mdl, ik_path); ps.set(float(times[n // 2]))
        fwd = ps.station_position("pelvis", (1.0, 0.0, 0.0)) - ps.body_position("pelvis"); az = np.arctan2(fwd[2], fwd[0]) + np.radians(float(a.azimuth[4:] or 0))
        print(f"azimuth auto: 골반 앞 방향 {np.degrees(np.arctan2(fwd[2], fwd[0])):.0f}° → 화면 가로 {np.degrees(az):.0f}°")
    else: az = np.radians(float(a.azimuth))
    hdir = np.array([np.cos(az), 0.0, np.sin(az)])
    # 근육 목록 (설계 공간 순서), 신호 U, 토크 tau, 보조 RES
    muscles = []; U = []; tau = {}; RES = {}; coords_all = []
    base_model_path = a.base if os.path.isabs(a.base) else os.path.join(ROOT, a.base)
    ik = io.read_mot(ik_path)
    for sp in sols:
        spd = pd.configure(sp); pj = json.load(open(paths_j[sp])); side = spd["design_side"]
        if pd.WRAP and a.box_scale != 1.0: pd.WRAP = dict(pd.WRAP, center=tuple(v * a.box_scale for v in pd.WRAP["center"]), radius=pd.WRAP["radius"] * a.box_scale)
        Q = np.column_stack([np.interp(times, ik.time, ik.column(c)) for c in list(pd.COORDS) + list(pd.EXTRA)]); fr = pd.Frames(base_model_path, Q)
        tt = np.zeros((n, len(pd.COORDS)))
        for m in pj["muscles"]:
            col = f"/forceset/{m['name']}" + (f"_{side}" if side else ""); u = np.interp(times, sol[sp].time, sol[sp].column(col)); U.append(u)
            muscles.append((sp, f"{sp.replace('_r', '')}_{m['name']}" + (f"_{side}" if side else "")))
            r, _ = pd.moment_arms(fr, [pd.BODIES.index(b) for b in m["chain"]], np.array(m["points"])); tt += (m["force_N"] * u)[:, None] * r
        for j, c in enumerate(pd.COORDS):
            tau[c] = tt[:, j]; cands = [k for k in sol[sp].names if "reserve" in k and k.endswith("_" + c)]
            RES[c] = np.interp(times, sol[sp].time, sol[sp].column(cands[0])) if cands else np.zeros(n); coords_all.append(c)
    U = np.column_stack(U); nm = len(muscles)
    # 뼈대 + 케이블 경로 (OpenSim 이 계산한 경로 점)
    model = osim.Model(gui_model); pose = io.Pose(model, ik_path); have = {model.getBodySet().get(i).getName() for i in range(model.getBodySet().getSize())}
    bodies = gv.collect_bodies(model) + [f"{b}_{s}" for s in "rl" for b in ("humerus", "ulna", "hand") if f"{b}_{s}" in have] + ["torso"]
    stations = [gv.TRUNK_TOP]; segs = gv.segments(bodies + ["trunk_top"]); idx = {b: i for i, b in enumerate(bodies + ["trunk_top"])}
    for s, col in (("r", gv.COLOR_R), ("l", gv.COLOR_L)):
        for up, low in ARM_CHAIN:
            if f"{up}_{s}" in idx and f"{low}_{s}" in idx: segs.append((idx[f"{up}_{s}"], idx[f"{low}_{s}"], col))
    if "humerus_r" in idx and "humerus_l" in idx: segs.append((idx["humerus_r"], idx["humerus_l"], gv.COLOR_TRUNK))
    if "torso" in idx: segs.append((idx["torso"], idx["trunk_top"], gv.COLOR_TRUNK))
    pos = np.empty((n, len(bodies) + 1, 3)); cab = [[] for _ in muscles]; fs = model.getForceSet()
    for i, t in enumerate(times):
        pose.set(float(t))
        for j, b in enumerate(bodies): pos[i, j] = pose.body_position(b)
        pos[i, -1] = pose.station_position(*gv.TRUNK_TOP)
        for j, (sp, mname) in enumerate(muscles):
            gp = osim.Muscle.safeDownCast(fs.get(mname)).getGeometryPath(); pts = gp.getCurrentPath(pose.state)
            P = np.array([[pts.get(k).getLocationInGround(pose.state).get(d) for d in range(3)] for k in range(pts.getSize())]); cab[j].append(P)
    H = lambda P: P @ hdir                                                       # 화면 가로 좌표
    has_lower = "lower" in sols; nrow = 5 if has_lower else 4
    fig = plt.figure(figsize=(15, 10.5 if has_lower else 9)); gs = fig.add_gridspec(nrow, 2, width_ratios=[1.0, 1.15], height_ratios=[0.8, 1, 1] + ([1] if has_lower else []) + [0.8], left=0.04, right=0.98, top=0.93, bottom=0.05, wspace=0.15, hspace=0.6)
    axf = fig.add_subplot(gs[:, 0]); axb = fig.add_subplot(gs[0, 1]); axt1 = fig.add_subplot(gs[1, 1]); axt2 = fig.add_subplot(gs[2, 1]); axt3 = fig.add_subplot(gs[3, 1]) if has_lower else None; axr = fig.add_subplot(gs[nrow - 1, 1])
    counts = {sp: sum(1 for s_, _ in muscles if s_ == sp) for sp in sols}
    fig.suptitle(f"케이블 근육 " + " + ".join(f"{ {'trunk': '몸통', 'shoulder_r': '어깨', 'elbow_r': '팔꿈치', 'lower': '다리'}[sp]} {n_}" for sp, n_ in counts.items()) + f" — {a.title or a.name}   (색 진하기 = 신호 u, MocoInverse, 운동학 실측 고정)   1/{a.slow:g} 배속", fontsize=12)
    axf.set_aspect("equal"); axf.axhline(0, color="#7f8c8d", lw=1.5); axf.set_ylim(-0.1, float(pos[:, :, 1].max()) + 0.25); axf.set_xlabel("[m]")
    for s_ in ("top", "right"): axf.spines[s_].set_visible(False)
    bones = LineCollection([], linewidths=4, colors=[c for _, _, c in segs], zorder=1); axf.add_collection(bones)
    cables = [axf.plot([], [], "-", lw=2.5, color=GROUP_COLORS[sp], zorder=3)[0] for sp, _ in muscles]
    tl = axf.text(0.02, 0.97, "", transform=axf.transAxes, fontsize=11, va="top"); span = 1.4
    x = np.arange(nm); bars = axb.bar(x, np.zeros(nm), color=[GROUP_COLORS[sp] for sp, _ in muscles]); axb.set_xticks(x); axb.set_xticklabels([m.replace("shoulder_", "어깨").replace("elbow_", "팔꿈치").replace("trunk_", "몸통").replace("lower_", "다리").replace("_r", "") for _, m in muscles], fontsize=7, rotation=30)
    axb.set_ylim(0, 1.05); axb.set_title("이 순간 케이블 근육의 신호 u (1 = 용량)", fontsize=11); axb.set_ylabel("u")
    LOW = ("hip_flexion_r", "hip_adduction_r", "knee_angle_r", "ankle_angle_r"); idsrc = lambda c: idm_lower if c in LOW else idm
    panels = [(axt1, [c for c in coords_all if c.startswith("lumbar")], "몸통 토크 — 케이블(실선) vs 측정 ID(점선)"), (axt2, [c for c in coords_all if c.startswith(("arm", "elbow"))], "오른팔 토크 — 케이블(실선) vs 측정 ID(점선)")]
    if has_lower: panels.append((axt3, [c for c in coords_all if c in LOW], "오른다리 토크 — 케이블(실선) vs 측정 ID(점선)"))
    for ax, cs, ttl in panels:
        for c in cs:
            ln, = ax.plot(times, tau[c], lw=1.8, label=KL.get(c, c)); ax.plot(times, np.interp(times, idsrc(c).time, idsrc(c).column(f"{c}_moment")), ":", color=ln.get_color(), lw=1.2)
        ax.set_title(ttl, fontsize=11); ax.set_ylabel("Nm"); ax.legend(fontsize=8, ncol=4, loc="lower left"); ax.set_xlim(t0, t1)
    curs = [ax.axvline(t0, color="k", lw=1) for ax, _, _ in panels]
    idpk = {c: float(np.abs(idsrc(c).column(f"{c}_moment")).max()) for c in coords_all}
    p95 = max(float(np.percentile(np.abs(RES[c]), 95)) / idpk[c] * 100 for c in coords_all)
    for c in coords_all: axr.plot(times, RES[c], lw=1.4, label=KL.get(c, c))
    curs.append(axr.axvline(t0, color="k", lw=1)); lim = max(5.0, float(max(np.abs(RES[c]).max() for c in coords_all)) * 1.2); axr.set_ylim(-lim, lim); axr.set_xlim(t0, t1)
    axr.set_title(f"보조 액추에이터가 떠맡은 토크 (케이블이 못 낸 몫) — 95퍼센타일 최악 {p95:.1f} % of ID 피크", fontsize=11); axr.set_ylabel("Nm"); axr.set_xlabel("time [s]"); axr.legend(fontsize=7, ncol=6, loc="lower left")
    grf = []; arrows = []; weight_n = 0.0
    if a.grf:
        import paths_video as pv; grf = pv.read_grf(a.grf if os.path.isabs(a.grf) else os.path.join(ROOT, a.grf)); weight_n = io.model_weight_n(model)
        arrows = [axf.annotate("", xy=(0, 0), xytext=(0, 0), arrowprops=dict(arrowstyle="-|>", color=gv.COLOR_GRF, lw=2)) for _ in grf]

    def update(fi):
        i = fi % n; p = pos[i]; h = H(p)
        bones.set_segments([[(h[s0], p[s0][1]), (h[s1], p[s1][1])] for s0, s1, _ in segs])
        for j, ln in enumerate(cables):
            P = cab[j][i]; ln.set_data(H(P), P[:, 1]); ln.set_alpha(0.35 + 0.65 * min(U[i, j], 1.0)); ln.set_linewidth(2.0 + 3 * min(U[i, j], 1.0))
        centre = float(h[0]); axf.set_xlim(centre - span / 2, centre + span / 2)
        tot = 0.0
        for arrow, (key, gt, cop, force) in zip(arrows, grf):
            f = np.array([np.interp(times[i], gt, force[:, k]) for k in range(3)]); cc = np.array([np.interp(times[i], gt, cop[:, k]) for k in range(3)]); mag = np.linalg.norm(f); tot += mag
            if mag > gv.GRF_MIN_BW * weight_n: arrow.set_position((float(cc @ hdir), cc[1])); arrow.xy = (float((cc + f / weight_n * gv.GRF_SCALE) @ hdir), cc[1] + f[1] / weight_n * gv.GRF_SCALE); arrow.set_visible(True)
            else: arrow.set_visible(False)
        tl.set_text(f"t = {times[i]:5.2f} s" + (f"   vGRF = {tot / weight_n:4.2f} BW" if grf else ""))
        for b, v in zip(bars, U[i]): b.set_height(v)
        for cur in curs: cur.set_xdata([times[i], times[i]])
        return [bones, tl] + curs + cables + list(bars) + arrows
    if a.png is not None:
        update(int(np.argmin(np.abs(times - a.png)))); out = os.path.join(ROOT, "report", f"frame_upper_{a.name}{a.tag}.png"); fig.savefig(out, dpi=110); print("저장", out); return
    anim = FuncAnimation(fig, update, frames=n * a.loops, interval=1000 / a.fps, blit=False)
    out = os.path.join(ROOT, "report", f"video_upper_{a.name}{a.tag}.mp4"); anim.save(out, writer=FFMpegWriter(fps=a.fps, bitrate=3000, extra_args=["-pix_fmt", "yuv420p"]), dpi=110)
    print(f"저장 {out}  ({os.path.getsize(out)/1e6:.1f} MB)")


if __name__ == "__main__":
    main()
