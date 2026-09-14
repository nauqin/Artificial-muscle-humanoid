"""모델별 설정 — 근육 이름·좌표 이름·ROM·대표 모멘트암을 여기 한 곳에만 둔다.

스크립트 본문에 `soleus_r` 같은 이름을 하드코딩하지 않는다. 나중에 상지로 확장할 때
(MoBL-ARMS 등) **설정 추가만으로** 끝나게 하려는 것이다. 근육 목록은 항상
`model.getMuscles()` 에서 읽고, 이 파일은 그 위에 얹는 주석(그룹·역할)이다.
"""

from __future__ import annotations

# ─────────────────────────────────────────────────────────────────────────────
# 하지 관절 세트 (다리당 6 DOF)
# ─────────────────────────────────────────────────────────────────────────────

#: mtp_angle 은 기본 제외(잠금). 실측 모멘트가 0.3 Nm 로 다른 좌표의 1/300 이라
#: 구동원을 배정할 이유가 없다 — Phase 0 에서 확인했다.
LOWER_LIMB_COORDS = (
    "hip_flexion", "hip_adduction", "hip_rotation",
    "knee_angle", "ankle_angle", "subtalar_angle",
)
OPTIONAL_COORDS = ("mtp_angle",)

SIDES = ("r", "l")

#: 좌표의 부호 규약. 일반화 모멘트가 **양수면 좌표를 키우는 방향**이다.
#: 실제 ID 결과에서 확인했다 (Phase 0): 발목은 −쪽이 189 Nm 로 압도적 = 저굴,
#: 무릎은 −쪽 = 신전, 고관절 외전이 −쪽.
DIRECTIONS = {
    "hip_flexion":     {"+": "굴곡", "-": "신전"},
    "hip_adduction":   {"+": "내전", "-": "외전"},
    "hip_rotation":    {"+": "내회전", "-": "외회전"},
    "knee_angle":      {"+": "굴곡", "-": "신전"},
    "ankle_angle":     {"+": "배굴", "-": "저굴"},
    "subtalar_angle":  {"+": "회외(추정)", "-": "회내(추정)"},
    "mtp_angle":       {"+": "신전", "-": "굴곡"},
}

# ─────────────────────────────────────────────────────────────────────────────
# 목표 ROM (도). AAOS 정상 가동범위 기준 — **뼈대 팀과 확정 전 임시값**
# ─────────────────────────────────────────────────────────────────────────────

ROM_TARGET_DEG = {
    "hip_flexion":    (-20.0, 120.0),
    "hip_adduction":  (-45.0,  30.0),
    "hip_rotation":   (-45.0,  45.0),
    "knee_angle":     (  0.0, 135.0),
    "ankle_angle":    (-50.0,  20.0),
    "subtalar_angle": (-15.0,  35.0),
}

# ─────────────────────────────────────────────────────────────────────────────
# 대표 근육과 실측 모멘트암 (m)
# ─────────────────────────────────────────────────────────────────────────────

#: NOTES.md 2절, 드라이런 Rajagopal(85.07 kg)에서 실측한 무리별 중앙값.
#: 부호는 뺀 크기다. Phase 2 가 ROM 격자에서 근육마다 다시 계산한다.
#: 값이 없는 방향은 **모름** — 지어내지 않는다.
REPRESENTATIVE_ARM_M = {
    ("ankle_angle", "-"): ("soleus/gasmed/gaslat", 0.051),
    ("ankle_angle", "+"): ("tibant", 0.050),
    ("knee_angle", "-"):  ("vasti", 0.045),
    ("knee_angle", "+"):  ("hamstrings", 0.053),
    ("hip_flexion", "+"): ("psoas/iliacus", 0.037),
    ("hip_flexion", "-"): ("glmax", 0.072),
}

#: 길항근 무리 — run_muscle_analysis.ANTAGONISTS 와 같은 정의를 쓴다.
#: 중복 정의를 만들지 않으려고 여기서는 import 로 참조만 한다.


def rom_span_deg(coord: str) -> float:
    lo, hi = ROM_TARGET_DEG[coord]
    return hi - lo


def coord_names(coord: str, sides=SIDES):
    return tuple(f"{coord}_{s}" for s in sides)
