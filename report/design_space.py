"""설계 공간 — "어느 좌표에, 어느 뼈 사이에, 어디까지 붙일 수 있나"를 한 곳에 모은다.

하체(`lower`)는 2026-09 설계 그대로이고, 상체는 `trunk`(요추 3자유도)와 `arm_r`(어깨 3 + 팔꿈치 1)이다.
`design_directions.py`, `path_design.py`, `moco_new_muscles.py`, `moco_path_muscles.py` 가 `--space` 로 고른다.

필드
  base       토크 열의 기본 이름 (side 없이). ID 열은 f"{base}_{side}_moment" 또는 side 가 "" 이면 f"{base}_moment"
  sides      토크를 모을 때 쌓는 side. 요추는 ("",) 하나.
  design_side 경로 설계·Moco 지표를 재는 side ("r" 또는 "")
  extra      자세에는 넣되 설계 방향에는 없는 좌표 (full name)
  bodies     케이블이 지날 수 있는 뼈, 몸통 쪽 → 말단 쪽 순서
  joint_of   설계 좌표 i 가 bodies[j]↔bodies[j+1] 사이 관절이면 j
  box        뼈별 부착 허용 상자 (뼈 좌표계, m). 오른쪽 기준 (z+ 바깥)
  mirror     왼쪽을 z 부호 반전으로 만들지 (요추는 False)
  k          기본 근육 수
  f_module   케이블 힘 상한 (N). 상체는 기계연 사양 확인 전 가정값
  labels     그림용 한글 라벨
"""

SPACES = {
    "lower": dict(
        base=("hip_flexion", "hip_adduction", "knee_angle", "ankle_angle"),
        sides=("r", "l"), design_side="r",
        extra=("hip_rotation_r",),
        bodies=("pelvis", "femur_r", "tibia_r", "calcn_r"),
        joint_of=(0, 0, 1, 2),
        box={"pelvis": ((-0.12, 0.08), (-0.12, 0.05), (0.02, 0.16)),
             "femur_r": ((-0.07, 0.07), (-0.42, -0.02), (-0.06, 0.07)),
             "tibia_r": ((-0.07, 0.06), (-0.40, -0.02), (-0.05, 0.05)),
             "calcn_r": ((-0.03, 0.18), (-0.02, 0.06), (-0.05, 0.05))},
        mirror=True, k=8, f_module=5000.0,
        labels=("고관절 굴곡", "고관절 내전", "무릎", "발목"),
    ),
    # 요추: back 관절은 pelvis (-0.114, 0.065, 0) = torso 원점. 어깨는 torso (0.004, 0.33, ±0.161).
    "trunk": dict(
        base=("lumbar_extension", "lumbar_bending", "lumbar_rotation"),
        sides=("",), design_side="",
        extra=(),
        bodies=("pelvis", "torso"),
        joint_of=(0, 0, 0),
        box={"pelvis": ((-0.18, 0.02), (-0.06, 0.09), (-0.14, 0.14)),
             "torso": ((-0.14, 0.10), (0.03, 0.32), (-0.17, 0.17))},
        mirror=False, k=6, f_module=5000.0,
        labels=("요추 굴곡(+)/신전", "요추 옆굽힘", "요추 회전"),
    ),
    # 오른팔: acromial_r 은 torso (0.004, 0.33, 0.161) = humerus_r 원점, elbow_r 은 humerus_r (0.013, -0.278, -0.01) = ulna_r 원점.
    "arm_r": dict(
        base=("arm_flex", "arm_add", "arm_rot", "elbow_flex"),
        sides=("r", "l"), design_side="r",
        extra=("pro_sup_r",),
        bodies=("torso", "humerus_r", "ulna_r"),
        joint_of=(0, 0, 0, 1),
        box={"torso": ((-0.13, 0.10), (0.08, 0.36), (0.00, 0.22)),
             "humerus_r": ((-0.07, 0.07), (-0.27, -0.02), (-0.07, 0.07)),
             "ulna_r": ((-0.05, 0.05), (-0.20, -0.02), (-0.04, 0.04))},
        mirror=True, k=6, f_module=1800.0,          # M3 450 N × 4 묶음 가정 (기계연 확인 전)
        # 경로 설계 자세에서 뺄 범위: 외전 75° 를 넘으면(팔을 머리 위로) 어깨 오일러각이 특이점이라 IK 가 뒤집힌다 (2026-11-03, CARRT subject6 찬 물병)
        pose_limit={"arm_add_r": (-1.31, 1.31)},
        labels=("어깨 굴곡", "어깨 내전", "어깨 회전", "팔꿈치 굴곡"),
    ),
    # 팔을 단관절로 나눈 것 (2026-11-03): 어깨·팔꿈치를 같이 당기는 직선 케이블은 어깨 굴곡 모멘트암이 안 나와 걷기(팔 흔들기)를 못 냈다.
    "shoulder_r": dict(
        base=("arm_flex", "arm_add", "arm_rot"),
        sides=("r", "l"), design_side="r",
        extra=("elbow_flex_r", "pro_sup_r"),
        bodies=("torso", "humerus_r"),
        joint_of=(0, 0, 0),
        box={"torso": ((-0.13, 0.10), (0.08, 0.36), (0.00, 0.22)),
             "humerus_r": ((-0.07, 0.07), (-0.27, -0.05), (-0.07, 0.07))},   # y ≤ −5 cm: 감는 구(반지름 4 cm) 밖
        mirror=True, k=4, f_module=1800.0,
        labels=("어깨 굴곡", "어깨 내전", "어깨 회전"),
        pose_limit={"arm_add_r": (-1.31, 1.31)},
        # 감는 면: 상완골두 자리(humerus_r 원점 = 어깨 중심)에 반지름 4 cm 구. 몸통→상완 케이블이 이 구를 감아 돈다 (2026-11-03).
        # 직선 케이블은 팔을 올리면 굴곡→신전으로 뒤집혀 걷기 팔 흔들기를 못 냈다. 사람 삼각근이 상완골두를 감는 것과 같은 해법.
        wrap=dict(body="humerus_r", center=(0.0, 0.0, 0.0), radius=0.04),
    ),
    "elbow_r": dict(
        base=("elbow_flex",),
        sides=("r", "l"), design_side="r",
        extra=("pro_sup_r",),
        bodies=("humerus_r", "ulna_r"),
        joint_of=(0,),
        # 전완 상자는 팔꿈치 뒤·위(주두 자리, x −8 / y +4 cm)까지 — 직선 케이블로 신전을 내려면 팔꿈치 뒤에 지렛대가 있어야 한다 (2026-11-03)
        box={"humerus_r": ((-0.07, 0.05), (-0.27, -0.02), (-0.05, 0.05)),
             "ulna_r": ((-0.08, 0.05), (-0.20, 0.04), (-0.04, 0.04))},
        mirror=True, k=2, f_module=1800.0,
        labels=("팔꿈치 굴곡",),
    ),
}


def get(name: str) -> dict:
    if name not in SPACES:
        raise SystemExit(f"설계 공간 {name!r} 없음 — {list(SPACES)}")
    sp = dict(SPACES[name]); sp["name"] = name
    sp["coords"] = tuple(full(b, sp["design_side"]) for b in sp["base"])
    return sp


def full(base: str, side: str) -> str:
    return f"{base}_{side}" if side else base


def chain_for(sp: dict, v, thr=0.12):
    """설계 방향 v 가 쓰는 관절들을 모두 덮는 연속 뼈 구간 (bodies 인덱스 목록)."""
    joints = [sp["joint_of"][i] for i in range(len(v)) if abs(v[i]) > thr]
    if not joints:
        joints = [sp["joint_of"][int(max(range(len(v)), key=lambda i: abs(v[i])))]]
    return list(range(min(joints), max(joints) + 2))


def mirror_body(sp: dict, body: str, side: str) -> str:
    """오른쪽 기준 뼈 이름을 side 쪽으로. 좌우가 없는 뼈(pelvis, torso)는 그대로."""
    if not sp["mirror"] or not side:
        return body
    return body[:-2] + "_" + side if body.endswith("_r") else body
