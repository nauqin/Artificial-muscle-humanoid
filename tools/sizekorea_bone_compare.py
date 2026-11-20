"""사이즈코리아 측정기준점과 '같은 정의·같은 바닥'으로 모델 치수를 비교 (2026-11-19).

sizekorea_compare.py 의 두 문제를 고친다 (그 파일과 결과 csv 는 그대로 둔다):
  ① 바닥: 모델 기본 자세에서 발이 y = 0 보다 아래에 있다 (175 cm 로 키우면서 다리가 길어졌는데 pelvis_ty 는 0.94 그대로).
     y = 0 을 바닥으로 잡아 높이가 모두 약 40 mm 낮게 나왔다.
     → 바닥 = 발뼈 메쉬 가장 아래 − 체중 실린 발뒤꿈치 살 11 mm. 확인: 머리 마커 − 바닥 − 마커 반지름 7 mm ≈ 키 1,750 mm.
  ② 재는 점: 어깨너비(겉 어깨점 vs 관절 중심), 무릎 높이(무릎뼈 가운데 vs 가쪽 위관절융기 마커) 처럼 정의가 달랐다.
     → 사이즈코리아 점과 같은 해부학 점을 쓴다: 피부 위 마커(마커 반지름 7 mm 보정) 또는 뼈 메쉬(무릎뼈).
     모델에 같은 점이 없는 것(겉 어깨점)은 비교하지 않고 그렇게 적는다.
    python tools/sizekorea_bone_compare.py → data/robot40_h175sk/sizekorea_compare_same.csv
"""
import argparse, csv, os
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEOM = os.path.join(ROOT, "report/moco/Geometry")
HEEL_PAD = 11.0          # mm, 체중이 실린 발뒤꿈치 살 두께 (가정, 대략 10~12 mm)
MARKER_R = 7.0           # mm, 마커 반지름 (지름 14 mm 가정)
TABLE = {"30대": dict(엉덩이너비=350.21, 발직선길이=256.29)}          # 8차_남성_모델링 치수.xlsx 30대 열


def mesh_ground(osim, m, s, body, name):
    b = m.getBodySet().get(body)
    for i in range(b.getPropertyByName("attached_geometry").size()):
        g = osim.Mesh.safeDownCast(b.get_attached_geometry(i))
        if g is None or g.get_mesh_file() != name: continue
        pm = osim.PolygonalMesh(); pm.loadFile(os.path.join(GEOM, name))
        V = np.array([[pm.getVertexPosition(k).get(j) for j in range(3)] for k in range(pm.getNumVertices())])
        sf = np.array([g.get_scale_factors().get(j) for j in range(3)]); T = g.getFrame().getTransformInGround(s)
        Rm = T.R().asMat33(); Rn = np.array([[Rm.get(r, c) for c in range(3)] for r in range(3)]); p = np.array([T.p().get(j) for j in range(3)])
        return (V * sf) @ Rn.T * 1000 + p * 1000                                             # mm, 지면 좌표 (x 앞, y 위, z 오른쪽)
    raise SystemExit(f"{body} 에 {name} 없음")


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--model", default="data/robot40_h175sk/robot40_h175sk.osim"); ap.add_argument("--age", default="30대")
    a = ap.parse_args()
    import opensim as osim
    osim.Logger.setLevelString("Error"); m = osim.Model(os.path.join(ROOT, a.model)); s = m.initSystem()
    M = {mk.getName(): mk.getLocationInGround(s).to_numpy() * 1000 for mk in m.getMarkerSet()}
    J = {j.getName(): j.getChildFrame().getPositionInGround(s).to_numpy() * 1000 for j in m.getJointSet()}
    L = np.loadtxt(os.path.join(ROOT, "data/사이즈코리아/표준체형_8차_남자_인체형상데이터", f"{a.age} 남성", f"{a.age}_남성.asc")); sk = lambda i: L[i - 1]   # noqa: E731
    T = TABLE[a.age]
    foot = mesh_ground(osim, m, s, "calcn_r", "r_foot.vtp"); toes = mesh_ground(osim, m, s, "toes_r", "r_bofoot.vtp")
    floor = foot[:, 1].min() - HEEL_PAD; H = lambda y: y - floor                              # noqa: E731
    pat = mesh_ground(osim, m, s, "patella_r", "r_patella.vtp")
    foot_bone = max(foot[:, 0].max(), toes[:, 0].max()) - foot[:, 0].min()
    rows = [("머리 위 (키)", sk(1)[1], H(M["THEAD"][1]) - MARKER_R, "머리 마커 − 반지름"),
            ("목뒤 높이 (C7)", sk(6)[1], H(M["C7"][1]), "C7 가시돌기 마커"),
            ("어깨점 높이", (sk(41)[1] + sk(42)[1]) / 2, H((M["RSHO"][1] + M["LSHO"][1]) / 2) - MARKER_R, "어깨뼈 봉우리 위 마커 − 반지름 (SK 는 봉우리 가장 가쪽)"),
            ("무릎뼈 가운데 높이", sk(26)[1], H((pat[:, 1].max() + pat[:, 1].min()) / 2), "무릎뼈 뼈 메쉬 위가·아래가 가운데"),
            ("가쪽 복사 높이", sk(28)[1], H(M["RLM"][1]), "가쪽 복사 마커"),
            ("엉덩이 너비", T["엉덩이너비"], float(np.linalg.norm(M["RGTRO"] - M["LGTRO"])) - 2 * MARKER_R, "큰돌기 마커 사이 − 반지름 2개"),
            ("발 길이", T["발직선길이"], foot_bone + 12.0, "발뼈 길이 + 뒤꿈치·발끝 살 약 12 mm"),
            ("어깨너비", float(np.linalg.norm(sk(41) - sk(42))), float("nan"), f"모델에 겉 어깨점 없음 (어깨 관절 중심 사이 {np.linalg.norm(J['acromial_r'] - J['acromial_l']):.0f} mm)")]
    out = os.path.join(ROOT, "data/robot40_h175sk/sizekorea_compare_same.csv")
    with open(out, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.writer(f); w.writerow(["치수", f"사이즈코리아_{a.age}_mm", "모델_mm", "차이_mm", "모델 점"])
        for n, v1, v2, note in rows:
            d = "" if np.isnan(v2) else f"{v2 - v1:+.1f}"
            w.writerow([n, f"{v1:.1f}", "" if np.isnan(v2) else f"{v2:.1f}", d, note]); print(f"{n:14s} SK {v1:7.1f}  모델 {v2:7.1f}  차이 {d:>7s}  ({note})")
    print(f"바닥 y = {floor:.1f} mm (발뼈 최저 {foot[:, 1].min():.1f} − {HEEL_PAD:.0f}) · 확인: 머리 마커 높이 {H(M['THEAD'][1]):.1f} − 반지름 {MARKER_R:.0f} = {H(M['THEAD'][1]) - MARKER_R:.1f} (키 1,750)")
    print(f"발뼈 길이 {foot_bone:.1f} mm · 무릎뼈 높이 범위 {pat[:, 1].max() - pat[:, 1].min():.1f} mm (SK 가운데–아래가 {sk(26)[1] - sk(39)[1]:.1f} mm) → {out}")


if __name__ == "__main__":
    main()
