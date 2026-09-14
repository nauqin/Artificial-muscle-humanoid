"""보행을 영상으로 — 하지 스틱피겨 + 지면반력 벡터.

화면 없이 돌아간다. `animate.py`(Simbody 뷰어)는 창을 띄워 눈으로 볼 때 쓰고,
이쪽은 발표자료에 붙일 파일을 만들 때 쓴다.

    python gait_video.py --model <osim> --ik <ik.mot> --grf <grf.mot> --out walk.mp4

지면은 정확히 y=0 이다. 발이 살짝 잠겨 보이면 `check_foot_height.py`로 확인한다 —
이 피험자는 발 Y 스케일이 1.66배라 calcn 원점이 발 바닥보다 위에 있다.
"""

from __future__ import annotations

import _env  # noqa: F401  환경 확인 — 다른 import 보다 먼저

import argparse
import os

import numpy as np

import addb_io as io

#: 하지 골격. (윗점, 아랫점) 쌍을 이어 그린다. 바디 원점이 곧 관절 중심이다.
CHAIN = (
    ("pelvis", "femur"), ("femur", "tibia"), ("tibia", "talus"),
    ("talus", "calcn"), ("calcn", "toes"),
)
#: 몸통은 방향만 보이면 되므로 한 줄로. 끝점은 torso 좌표계에서 위로 잡는다 —
#: torso 바디 원점은 골반 바로 위라, 원점만 이으면 토막처럼 보인다.
TRUNK = ("pelvis", "torso")
TRUNK_TOP = ("torso", (0.0, 0.42, 0.0))

COLOR_R = "#c0392b"
COLOR_L = "#2471a3"
COLOR_TRUNK = "#555555"
COLOR_GRF = "#27ae60"

#: GRF 화살표 길이 = 힘 / 체중 × 이 값 (m)
GRF_SCALE = 0.35
#: 이보다 작은 힘은 안 그린다. 낮게 두면 화살촉만 남아 얼룩처럼 보인다.
GRF_MIN_BW = 0.05


def sample_positions(model, ik_path, times, bodies, stations=()):
    """각 시각의 위치. (T, n_bodies + n_stations, 3)"""
    pose = io.Pose(model, ik_path)
    out = np.empty((len(times), len(bodies) + len(stations), 3))
    for i, t in enumerate(times):
        pose.set(float(t))
        for j, body in enumerate(bodies):
            out[i, j] = pose.body_position(body)
        for j, (body, offset) in enumerate(stations):
            out[i, len(bodies) + j] = pose.station_position(body, offset)
    return out


def collect_bodies(model):
    have = {model.getBodySet().get(i).getName()
            for i in range(model.getBodySet().getSize())}
    bodies = ["pelvis"]
    for side in ("r", "l"):
        for _up, low in CHAIN:
            name = low if low == "pelvis" else f"{low}_{side}"
            if name in have and name not in bodies:
                bodies.append(name)
    return bodies


def segments(bodies):
    """(인덱스1, 인덱스2, 색) 목록."""
    idx = {b: i for i, b in enumerate(bodies)}
    out = []
    for side, color in (("r", COLOR_R), ("l", COLOR_L)):
        for up, low in CHAIN:
            a = up if up == "pelvis" else f"{up}_{side}"
            b = low if low == "pelvis" else f"{low}_{side}"
            if a in idx and b in idx:
                out.append((idx[a], idx[b], color))
    # 골반 폭 — 좌우 고관절을 잇는다
    if "femur_r" in idx and "femur_l" in idx:
        out.append((idx["femur_r"], idx["femur_l"], COLOR_TRUNK))
    if "pelvis" in idx and "trunk_top" in idx:
        out.append((idx["pelvis"], idx["trunk_top"], COLOR_TRUNK))
    return out


def render(model_path, ik_path, out_path, *, grf_path=None, fps=30, view="sagittal",
           follow=True, start=None, end=None, dpi=110, log=print):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, FFMpegWriter, PillowWriter
    from matplotlib.collections import LineCollection

    import opensim as osim
    osim.Logger.setLevelString("Error")

    model = osim.Model(os.path.abspath(model_path))
    model.initSystem()
    weight_n = io.model_weight_n(model)

    ik_mot = io.read_mot(ik_path)
    t0 = ik_mot.t0 if start is None else start
    t1 = ik_mot.t1 if end is None else end
    times = np.arange(t0, t1, 1.0 / fps)
    if times.size < 2:
        raise ValueError(f"프레임이 부족하다 ({times.size}개)")

    bodies = collect_bodies(model)
    have = {model.getBodySet().get(i).getName() for i in range(model.getBodySet().getSize())}
    stations = [TRUNK_TOP] if TRUNK_TOP[0] in have else []
    pos = sample_positions(model, ik_path, times, bodies, stations)
    labels = bodies + (["trunk_top"] if stations else [])
    segs = segments(labels)
    log(f"  {len(times)} 프레임 ({t0:.2f}~{t1:.2f}s, {fps} fps), 마디 {len(segs)}개")

    # 시상면은 x-y, 관상면은 z-y
    hi, vi = (0, 1) if view == "sagittal" else (2, 1)

    grf = None
    if grf_path:
        grf_mot = io.read_mot(grf_path)
        plates = {k: p for k, p in io.find_plates(grf_mot).items() if p.complete}
        grf = []
        for key, plate in plates.items():
            cop = np.column_stack([np.interp(times, grf_mot.time, grf_mot.column(plate.point[a]))
                                   for a in "xyz"])
            force = np.column_stack([np.interp(times, grf_mot.time, grf_mot.column(plate.force[a]))
                                     for a in "xyz"])
            grf.append((key, cop, force))
        log(f"  GRF 힘판 {len(grf)}개: {[g[0] for g in grf]}")

    fig, ax = plt.subplots(figsize=(9, 5))
    ax.set_aspect("equal")
    ax.axhline(0.0, color="#7f8c8d", lw=1.5, zorder=1)   # 지면 y=0
    ax.set_ylim(-0.15, float(pos[:, :, 1].max()) + 0.35)
    ax.set_xlabel("[m]")
    ax.set_yticks([0, 0.5, 1.0])
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    span = 1.6
    if not follow:
        lo, hi_x = float(pos[:, :, hi].min()), float(pos[:, :, hi].max())
        ax.set_xlim(lo - 0.4, hi_x + 0.4)

    lines = LineCollection([], linewidths=3.5, zorder=3)
    ax.add_collection(lines)
    joints = ax.plot([], [], "o", color="#2c3e50", ms=4.5, zorder=4)[0]
    arrows = [ax.annotate("", xy=(0, 0), xytext=(0, 0),
                          arrowprops=dict(arrowstyle="-|>", color=COLOR_GRF, lw=2.2))
              for _ in (grf or [])]
    label = ax.text(0.02, 0.95, "", transform=ax.transAxes, fontsize=10,
                    family="monospace", va="top")

    def update(i):
        p = pos[i]
        lines.set_segments([[(p[a][hi], p[a][vi]), (p[b][hi], p[b][vi])] for a, b, _ in segs])
        lines.set_color([c for _a, _b, c in segs])
        joints.set_data(p[:, hi], p[:, vi])
        if follow:
            centre = float(p[0][hi])
            ax.set_xlim(centre - span / 2, centre + span / 2)
        text = f"t = {times[i]:5.2f} s"
        for arrow, (key, cop, force) in zip(arrows, grf or []):
            mag = float(np.linalg.norm(force[i]))
            if mag > GRF_MIN_BW * weight_n:
                tail = (cop[i][hi], cop[i][vi])
                tip = (cop[i][hi] + force[i][hi] / weight_n * GRF_SCALE,
                       cop[i][vi] + force[i][vi] / weight_n * GRF_SCALE)
                arrow.set_position(tail)
                arrow.xy = tip
                arrow.set_visible(True)
            else:
                arrow.set_visible(False)
        if grf:
            total = sum(float(np.linalg.norm(f[i])) for _k, _c, f in grf)
            text += f"   vGRF = {total / weight_n:4.2f} BW"
        label.set_text(text)
        return [lines, joints, label] + arrows

    anim = FuncAnimation(fig, update, frames=len(times), interval=1000 / fps, blit=False)

    os.makedirs(os.path.dirname(os.path.abspath(out_path)) or ".", exist_ok=True)
    if out_path.lower().endswith(".mp4"):
        try:
            writer = FFMpegWriter(fps=fps, bitrate=2400)
            anim.save(out_path, writer=writer, dpi=dpi)
        except Exception as exc:
            out_path = os.path.splitext(out_path)[0] + ".gif"
            log(f"  [WARN] mp4 저장 실패({exc}) — GIF로 바꾼다")
            anim.save(out_path, writer=PillowWriter(fps=fps), dpi=dpi)
    else:
        anim.save(out_path, writer=PillowWriter(fps=fps), dpi=dpi)
    plt.close(fig)
    log(f"  → {out_path}  ({os.path.getsize(out_path) / 1e6:.1f} MB)")
    return out_path


def main():
    ap = argparse.ArgumentParser(description="보행 영상 (하지 스틱피겨 + GRF)")
    ap.add_argument("--model", required=True)
    ap.add_argument("--ik", required=True)
    ap.add_argument("--grf", default=None, help="있으면 지면반력 벡터를 함께 그린다")
    ap.add_argument("--out", required=True, help=".mp4 또는 .gif")
    ap.add_argument("--fps", type=int, default=30)
    ap.add_argument("--view", default="sagittal", choices=("sagittal", "frontal"))
    ap.add_argument("--no-follow", action="store_true", help="카메라 고정 (기본은 골반 추적)")
    ap.add_argument("--start", type=float, default=None)
    ap.add_argument("--end", type=float, default=None)
    args = ap.parse_args()

    render(args.model, args.ik, args.out, grf_path=args.grf, fps=args.fps,
           view=args.view, follow=not args.no_follow, start=args.start, end=args.end)


if __name__ == "__main__":
    main()
