"""인공근육 모듈의 개념 모델 — 모든 Phase가 여기서 파라미터를 가져간다.

**파라미터 정의가 두 곳에 있으면 안 된다.** 값을 바꾸려면 이 파일 상단만 고친다.

출처
----
- 모듈 규격표(`MODULES`) : 한국기계연구원 발표자료 p.35~36. 힘·변위는 **범위**만
  주어져 있어 하한(보수적)과 상한을 둘 다 돌린다.
- 개념 모델 파라미터    : 기계연 fabric muscle 논문(Park et al., IEEE TNSRE 33, 2025)
  에서 **형태와 자릿수만** 가져왔다. 실제 과제 근육은 더 발전된 소재(KAIST TC-μSMA)
  이므로 **수치는 논문이 아니라 발표자료의 목표값**을 쓴다.
- 힘밀도 목표 15 N/g   : 발표자료 p.14, p.20. 민감도로 논문값 5.4 N/g 도 본다.

이번 범위에서 다루지 않는 것: 속도·대역폭·온도. SMA는 가열/냉각이 비대칭이지만
그건 소재·냉각 팀 몫이고, 여기서는 **준정적(quasi-static)** 실현가능성만 본다.
"""

from __future__ import annotations

from dataclasses import dataclass

# ─────────────────────────────────────────────────────────────────────────────
# 파라미터
# ─────────────────────────────────────────────────────────────────────────────

#: 최대 수축 변형률. 스트로크는 모듈 활성 길이의 비율이다 (논문 39~43.5 %).
EPS_MAX = 0.40
#: 수동(이완) 강성을 정하는 변형률 — 자유길이에서 이만큼 늘면 F_max 만큼 저항한다.
#: 논문의 "98 N에 10 % 늘어남" 에서 왔다.
EPS_PASSIVE = 0.10
#: 힘밀도 (N/g). 질량 = F_max / RHO_F. 냉각·구조 질량은 제외한다.
RHO_F_TARGET = 15.0
RHO_F_PAPER = 5.4

#: 수동 힘 상한 (F_max 배수).
#: **모델 그대로 두면 발산한다** — EPS_PASSIVE=0.10 인데 스트로크가 EPS_MAX=0.40
#: 이므로, 완전 신장 자세에서 f_passive = 4·F_max 가 된다. 길항근이 자기 능동
#: 최대의 4배로 버티는 셈이라 어떤 배치도 불가 판정이 난다. 물리적으로도 SMA
#: 코일실은 그 전에 항복한다. 그래서 상한을 둔다. 이 값을 바꿔 민감도를 본다.
PASSIVE_CAP_RATIO = 1.0


@dataclass(frozen=True)
class Module:
    """표준 인공근육 모듈 한 등급. 힘·변위는 (하한, 상한)."""

    name: str
    label: str
    force_n: tuple          # (하한, 상한) N
    stroke_mm: tuple        # (하한, 상한) mm
    form: str
    example: str

    def spec(self, bound: str = "low"):
        """`low`(보수적) 또는 `high` 쪽 (F_max N, stroke m)."""
        i = 0 if bound == "low" else 1
        return self.force_n[i], self.stroke_mm[i] / 1000.0

    def mass_g(self, bound="low", rho_f=RHO_F_TARGET) -> float:
        return self.spec(bound)[0] / rho_f


#: 표준 모듈 규격 (발표자료 p.35~36)
MODULES = (
    Module("M1", "초대형", (1000, 5000), (40, 120), "이중깃형",  "대둔근, 대퇴사두근"),
    Module("M2", "대형",   (400, 1200),  (30,  80), "깃형",      "삼각근, 대흉근"),
    Module("M3", "중형",   (150,  450),  (25,  60), "방추형",    "이두근, 삼두근"),
    Module("M4", "소형",   (40,   200),  (35,  90), "평행형",    "지굴근, 봉공근"),
    Module("M5", "초소형", (2,     50),  ( 5,  20), "복합형",    "충양근, 골간근"),
    Module("M6", "평판형", (50,  1500),  (20,  80), "판상형",    "복사근, 기립근"),
)

BY_NAME = {m.name: m for m in MODULES}


# ─────────────────────────────────────────────────────────────────────────────
# 힘 모델
# ─────────────────────────────────────────────────────────────────────────────


def stroke_from_length(l0: float, eps_max: float = EPS_MAX) -> float:
    """활성 길이 `l0`(m)인 모듈의 스트로크(m). 힘이 아니라 **길이의 비율**이다."""
    return eps_max * l0


def active_range(f_max: float):
    """능동 힘 범위. 당기기만 한다 → 길항 배치가 필수다."""
    return 0.0, f_max


def passive_force(length, free_length, l0, f_max, *, eps_passive=EPS_PASSIVE,
                  cap_ratio=PASSIVE_CAP_RATIO, slack=0.0):
    """이완 상태의 수동 힘 (N). 식은 근육은 스프링이다.

    `slack` 은 변형 B의 텐던 여유 길이 `s_i`. 0이면 변형 A(항상 연결)다.
    상한은 `cap_ratio · f_max` — 그 이유는 파일 상단 `PASSIVE_CAP_RATIO` 참고.
    """
    stretch = length - free_length - slack
    if stretch <= 0.0 or l0 <= 0.0:
        return 0.0
    force = f_max * stretch / (eps_passive * l0)
    return min(force, cap_ratio * f_max) if cap_ratio is not None else force


def force_bounds(length, free_length, l0, f_max, **kw):
    """근육 i가 낼 수 있는 힘 범위 `[f_passive, f_passive + F_max]`.

    **하한이 0이 아니라는 것이 Hill 근육과의 차이**이고, 실현가능성 판정에
    그대로 들어간다.
    """
    fp = passive_force(length, free_length, l0, f_max, **kw)
    return fp, fp + f_max


def module_mass_g(f_max: float, rho_f: float = RHO_F_TARGET) -> float:
    return f_max / rho_f


def required_length_for_travel(delta_l: float, eps_max: float = EPS_MAX) -> float:
    """이동량 `delta_l`(m)을 내려면 필요한 활성 길이 `l0 ≥ Δl / ε_max`."""
    return delta_l / eps_max


def positive_spanning_floor(n_dof: int) -> int:
    """당김 전용 액추에이터로 n차원 토크 공간을 덮는 데 필요한 최소 개수.

    당기기만 하는 액추에이터의 모멘트암 열들이 R^n 을 **양의 결합**으로 덮으려면
    최소 n+1 개가 필요하다. 하지 6 DOF → 다리당 7개, 양다리 14개가 이론적 바닥이다.
    """
    return n_dof + 1
