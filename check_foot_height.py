"""발·골반 높이로 바닥 관통·단위를 확인한다.

이 피험자는 발 Y 스케일이 1.66배라 발 메시가 calcn 원점 아래로 내려간다.
GUI에서 발이 살짝 잠겨 보이는 건 그 때문이고 IK 결과는 정상이다.
힘판 표면은 y=0 이다.
"""

from __future__ import annotations

import _env  # noqa: F401  환경 확인 — 다른 import 보다 먼저

import argparse

import numpy as np

import addb_io as io


def main():
    ap = argparse.ArgumentParser(description="발·골반 높이 확인")
    ap.add_argument("--model", required=True)
    ap.add_argument("--ik", required=True)
    ap.add_argument("--bodies", default="calcn_r,calcn_l")
    ap.add_argument("--samples", type=int, default=60)
    args = ap.parse_args()

    import opensim as osim
    osim.Logger.setLevelString("Error")

    model = osim.Model(args.model)
    model.initSystem()
    pose = io.Pose(model, args.ik)
    times = np.linspace(pose.time[0], pose.time[-1], args.samples)
    bodies = [b.strip() for b in args.bodies.split(",") if b.strip()]

    heights = {b: [] for b in bodies}
    heights["pelvis"] = []
    markers = []
    marker_set = model.getMarkerSet()

    for t in times:
        state = pose.set(float(t))
        for b in list(heights):
            if b == "pelvis" and "pelvis" not in bodies:
                heights[b].append(pose.body_position("pelvis")[1])
            elif b in bodies:
                heights[b].append(pose.body_position(b)[1])
        lows = []
        for i in range(marker_set.getSize()):
            marker = marker_set.get(i)
            name = marker.getName().lower()
            if any(k in name for k in ("toe", "heel", "ankle", "mt5", "cal")):
                lows.append(marker.getLocationInGround(state).get(1))
        if lows:
            markers.append(min(lows))

    for name, vals in heights.items():
        vals = np.asarray(vals)
        print(f"  {name:10s} 최저 {vals.min()*1000:+7.1f} mm   최고 {vals.max()*1000:+7.1f} mm")
    if markers:
        m = np.asarray(markers)
        print(f"  {'발 마커':10s} 최저 {m.min()*1000:+7.1f} mm")

    pelvis = np.asarray(heights["pelvis"])
    print(f"\n  골반 높이 {pelvis.min():.3f} ~ {pelvis.max():.3f} m")
    if pelvis.min() < 0.4:
        print("  [WARN] 골반이 너무 낮다 — 단위(mm vs m)나 좌표계를 의심한다")
    else:
        print("  자릿수 정상 (힘판 표면 y=0 기준)")


if __name__ == "__main__":
    main()
