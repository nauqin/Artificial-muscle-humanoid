"""사이즈코리아 8차 남성 표준체형(30대, 키 1,751 mm) 측정기준점과 로봇 모델 치수 비교 (2026-11-19).

data/사이즈코리아/표준체형_8차_남자_인체형상데이터/30대 남성/30대_남성.asc = 측정기준점 47개 (mm, 순서는 같은 폴더의 '8차 표준인체형상 측정기준점.xlsx',
y 위·z 앞·x− 오른쪽), 치수표 = data/사이즈코리아/8차_남성_모델링 치수.xlsx.
    python tools/sizekorea_compare.py [--model data/robot40_h175sk/robot40_h175sk.osim] [--age 30대]
→ data/robot40_h175sk/sizekorea_compare.csv

주의 (2026-11-19 밤): 높이를 y = 0 기준으로 쟀는데, 모델 기본 자세의 바닥은 약 −43 mm 다 (높이가 약 40 mm 낮게 나옴).
어깨너비·무릎 높이도 사이즈코리아와 재는 점이 다르다. 같은 바닥·같은 점 비교는 tools/sizekorea_bone_compare.py.
"""
import argparse, csv, os
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLE = {"30대": dict(엉덩이너비=350.21, 발직선길이=256.29)}          # 8차_남성_모델링 치수.xlsx 의 30대 열 (이 비교에 쓰는 것만)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--model", default="data/robot40_h175sk/robot40_h175sk.osim"); ap.add_argument("--age", default="30대")
    a = ap.parse_args()
    import opensim as osim
    L = np.loadtxt(os.path.join(ROOT, "data/사이즈코리아/표준체형_8차_남자_인체형상데이터", f"{a.age} 남성", f"{a.age}_남성.asc")); sk = lambda i: L[i - 1]   # noqa: E731
    osim.Logger.setLevelString("Error"); m = osim.Model(os.path.join(ROOT, a.model)); s = m.initSystem()
    M = {mk.getName(): mk.getLocationInGround(s).to_numpy() * 1000 for mk in m.getMarkerSet()}
    J = {j.getName(): j.getChildFrame().getPositionInGround(s).to_numpy() * 1000 for j in m.getJointSet()}
    d = lambda p, q: float(np.linalg.norm(p - q))   # noqa: E731
    T = TABLE[a.age]
    rows = [("키 (머리 위)", sk(1)[1], M["THEAD"][1], "머리마루점 / THEAD 마커"),
            ("목뒤 높이 (C7)", sk(6)[1], M["C7"][1], "목뒤점 / C7 마커"),
            ("어깨 높이", (sk(41)[1] + sk(42)[1]) / 2, (M["RSHO"][1] + M["LSHO"][1]) / 2, "어깨점 / RSHO·LSHO"),
            ("무릎 높이", sk(26)[1], M["RLEK"][1], "무릎뼈가운데점 / 가쪽 위관절융기(관절선) — 기준점이 20~25 mm 다름"),
            ("가쪽 복사 높이", sk(28)[1], M["RLM"][1], "가쪽복사점 / RLM"),
            ("엉덩이 너비", T["엉덩이너비"], d(M["RGTRO"], M["LGTRO"]), "치수표 / 큰돌기 마커 사이"),
            ("발 길이", T["발직선길이"], d(M["RHEE"], M["RTOE"]), "치수표 / 발뒤꿈치–발끝 마커"),
            ("위팔 (어깨점–팔꿈치)", d(sk(41), sk(33)), d(M["RSHO"], M["RLEE"]), "모델은 마디별 추가 배율 적용"),
            ("아래팔 (팔꿈치–손목)", d(sk(33), (sk(30) + sk(31)) / 2), d(M["RLEE"], (M["RMW"] + M["RLW"]) / 2), "모델은 마디별 추가 배율 적용"),
            ("어깨너비", d(sk(41), sk(42)), d(J["acromial_r"], J["acromial_l"]), "SK 겉 어깨점 / 모델 어깨 관절 중심 (정의 다름, 관절 중심은 겉에서 25~30 mm 안쪽)")]
    print(f"{'치수':22s} {'SK':>7s} {'모델':>7s} {'차이':>6s}  기준")
    for n, x, y, note in rows: print(f"{n:22s} {x:7.0f} {y:7.0f} {y - x:+6.0f}  {note}")
    out = os.path.join(ROOT, os.path.dirname(a.model), "sizekorea_compare.csv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh); w.writerow(["치수", f"사이즈코리아_{a.age}_mm", "모델_mm", "차이_mm", "기준"]); w.writerows([(n, round(x, 1), round(y, 1), round(y - x, 1), note) for n, x, y, note in rows])
    print("→", os.path.relpath(out, ROOT))


if __name__ == "__main__":
    main()
