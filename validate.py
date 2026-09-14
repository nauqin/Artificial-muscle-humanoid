"""검산 밴드와 리포트.

목적은 **자릿수**가 맞는지 보는 것이다. 계산은 입력이 틀려도 불평 없이 숫자를
내놓기 때문에, 검산이 없으면 조용히 틀린 결과로 몇 달을 갈 수 있다.

  PASS  정상범위 안
  WARN  밴드를 벗어났지만 자릿수는 맞음 — 피험자·속도 차이일 수 있다
  FAIL  상한의 2배 초과 또는 하한의 절반 미만 — 파이프라인 오류를 의심
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# 밴드
# ─────────────────────────────────────────────────────────────────────────────

#: 항목 -> (하한, 상한, 단위, 벗어나면 의심할 것)
BANDS = {
    "ankle":    (1.2, 1.6, "Nm/kg", "GRF-발 매핑 / 좌표계 / 단위"),
    "knee":     (0.4, 0.8, "Nm/kg", "GRF-발 매핑 / 좌표계 / 단위"),
    "hip":      (0.7, 1.2, "Nm/kg", "GRF-발 매핑 / 좌표계 / 단위"),
    "grf_peak": (0.9, 1.3, "BW",    "GRF 단위(N vs BW) / 체중"),
    "cop_dist": (None, 0.30, "m",   "COP 단위(mm vs m) / 좌우 반전"),
    "residual": (None, 5.0, "%",    "모델-모션 부정합 / 동역학 불일치"),
    "recon":    (None, 5.0, "%",    "SO 미수렴 / 시간축 불일치 / 근육 누락"),
    "reserve":  (None, 10.0, "%",   "근육만으로 모멘트 부족 — 모델·구간·GRF 확인"),
    "moment_arm": (0.02, 0.08, "m",  "모델 스케일 / 길이 단위(mm vs m)"),
}

#: 한 걸음의 대략적 길이. 분석 구간이 이보다 짧으면 push-off가 잘렸을 수 있다.
ONE_STEP_S = 0.8

VERDICT_ORDER = {"PASS": 0, "WARN": 1, "FAIL": 2}


def verdict(value, lo, hi) -> str:
    """FAIL은 상한의 2배 초과 또는 하한의 절반 미만일 때만."""
    if value is None or not np.isfinite(value):
        return "FAIL"
    if hi is not None and value > 2.0 * hi:
        return "FAIL"
    if lo is not None and value < 0.5 * lo:
        return "FAIL"
    if hi is not None and value > hi:
        return "WARN"
    if lo is not None and value < lo:
        return "WARN"
    return "PASS"


def worst(verdicts) -> str:
    return max(verdicts, key=lambda v: VERDICT_ORDER[v]) if verdicts else "PASS"


@dataclass
class Check:
    name: str
    value: float
    band: str          # BANDS 의 key
    verdict: str = "PASS"
    note: str = ""


@dataclass
class Report:
    checks: list = field(default_factory=list)
    hints: list = field(default_factory=list)

    def add(self, name, value, band, note=""):
        lo, hi, _unit, _suspect = BANDS[band]
        self.checks.append(Check(name, value, band, verdict(value, lo, hi), note))
        return self.checks[-1]

    @property
    def verdict(self) -> str:
        return worst([c.verdict for c in self.checks])

    def get(self, name):
        for c in self.checks:
            if c.name == name:
                return c
        return None

    def render(self, indent="  ") -> str:
        lines = []
        for c in self.checks:
            lo, hi, unit, suspect = BANDS[c.band]
            if lo is None:
                rng = f"< {hi:g}"
            elif hi is None:
                rng = f"> {lo:g}"
            else:
                rng = f"{lo:g} ~ {hi:g}"
            line = f"{indent}{c.name:<16s} {c.value:8.3f} {unit:<6s} [{rng:>11s}]  {c.verdict}"
            if c.verdict != "PASS":
                line += f"   의심: {suspect}"
            if c.note:
                line += f"   {c.note}"
            lines.append(line)
        for h in self.hints:
            lines.append(f"{indent}※ {h}")
        return "\n".join(lines)


# ─────────────────────────────────────────────────────────────────────────────
# ID 검산
# ─────────────────────────────────────────────────────────────────────────────

#: 관절 -> ID 출력 컬럼 접두사
JOINTS = {
    "ankle": "ankle_angle",
    "knee": "knee_angle",
    "hip": "hip_flexion",
}


def peak_moment(id_mot, coord: str):
    """`<coord>_moment` 컬럼의 피크 절대값 (Nm). 없으면 None."""
    for suffix in ("_moment", "_force"):
        name = coord + suffix
        if name in id_mot.names:
            return float(np.max(np.abs(id_mot.column(name))))
    return None


def residual_percent(id_mot, grf_peak_n: float):
    """골반 잔차력 피크 / 피크 GRF (%)."""
    cols = [f"pelvis_t{a}_force" for a in "xyz"]
    if not all(c in id_mot.names for c in cols):
        return None
    mag = np.linalg.norm(np.column_stack([id_mot.column(c) for c in cols]), axis=1)
    if grf_peak_n <= 0:
        return None
    return float(mag.max() / grf_peak_n * 100.0)


def validate_id(id_mot, *, mass_kg, grf_peak_bw, grf_peak_n, cop_dist,
                window=None) -> Report:
    """ID 결과 한 건을 검산한다."""
    rep = Report()
    near_zero = []

    for joint, prefix in JOINTS.items():
        for side in ("r", "l"):
            coord = f"{prefix}_{side}"
            peak = peak_moment(id_mot, coord)
            if peak is None:
                continue
            per_kg = peak / mass_kg
            rep.add(f"{joint}_{side}", per_kg, joint)
            lo = BANDS[joint][0]
            if lo is not None and per_kg < 0.5 * lo:
                near_zero.append(f"{joint}_{side}")

    rep.add("grf_peak", grf_peak_bw, "grf_peak")
    if cop_dist is not None:
        rep.add("cop_dist", cop_dist, "cop_dist")
    resid = residual_percent(id_mot, grf_peak_n)
    if resid is not None:
        rep.add("residual", resid, "residual")

    # 한쪽 관절만 0에 가까우면 매핑보다 먼저 구간 길이를 본다.
    # 구간이 한 걸음보다 짧으면 그쪽 다리의 push-off가 통째로 잘렸을 수 있다.
    if near_zero and window is not None:
        span = window[1] - window[0]
        if span < ONE_STEP_S:
            rep.hints.append(
                f"{', '.join(near_zero)} 이(가) 0에 가깝다. 분석 구간이 {span:.2f} s로 "
                f"한 걸음({ONE_STEP_S:g} s)보다 짧다 — GRF-발 매핑이 아니라 midstance 골이 "
                f"지지 구간을 쪼갠 경우일 수 있다.")
        else:
            rep.hints.append(
                f"{', '.join(near_zero)} 이(가) 0에 가깝다. 구간({span:.2f} s)은 충분하므로 "
                f"GRF-발 매핑·좌표계를 본다. 좌우를 따로 볼 것 — 한쪽만 틀리면 매핑 문제다.")
    return rep


# ─────────────────────────────────────────────────────────────────────────────
# SO 검산
# ─────────────────────────────────────────────────────────────────────────────


def validate_so(id_mot, force_mot, ma_by_coord, coords, *, reserve_prefix="reserve_") -> Report:
    """`Σ(근육힘 × 모멘트암) + 보조액추에이터 =? ID 모멘트` 를 좌표마다.

    `recon`(재현 오차) 하나만 보면 안 된다. 보조 액추에이터는 어떤 모멘트든
    만들어낼 수 있으므로 재현 오차는 언제나 0에 가깝게 나온다. 근육이 전혀 일을
    안 해도 통과한다. **`reserve`가 실제 품질 지표다.**
    """
    rep = Report()
    t = id_mot.time

    for coord in coords:
        ma = ma_by_coord.get(coord)
        if ma is None:
            continue
        id_col = f"{coord}_moment" if f"{coord}_moment" in id_mot.names else f"{coord}_force"
        if id_col not in id_mot.names:
            continue
        id_moment = id_mot.column(id_col)

        # 근육 힘 × 모멘트암 — 모멘트암은 ID 시간격자로 보간한다.
        muscle_moment = np.zeros_like(t)
        for muscle in ma.names:
            if muscle not in force_mot.names:
                continue
            arm = np.interp(t, ma.time, ma.column(muscle))
            force = np.interp(t, force_mot.time, force_mot.column(muscle))
            muscle_moment += arm * force

        # 보조 액추에이터. 좌표 이름은 **정확히 일치**시켜 찾는다 —
        # 부분일치를 쓰면 knee_angle_r 이 knee_angle_r_beta 것까지 끌어온다.
        reserve = np.zeros_like(t)
        for cand in (f"{reserve_prefix}{coord}", coord):
            if cand in force_mot.names:
                reserve = np.interp(t, force_mot.time, force_mot.column(cand))
                break

        peak = float(np.max(np.abs(id_moment)))
        if peak <= 0:
            continue
        rms = float(np.sqrt(np.mean((muscle_moment + reserve - id_moment) ** 2)))
        rep.add(f"recon/{coord}", rms / peak * 100.0, "recon")
        rep.add(f"reserve/{coord}", float(np.max(np.abs(reserve))) / peak * 100.0, "reserve")

    return rep


def max_of(rep: Report, prefix: str):
    vals = [c.value for c in rep.checks if c.name.startswith(prefix)]
    return max(vals) if vals else None


# ─────────────────────────────────────────────────────────────────────────────
# 그림
# ─────────────────────────────────────────────────────────────────────────────


def plot_id(path, id_mot, grf_mot, plates, mass_kg, weight_n, window=None):
    """관절 모멘트·GRF 그래프. 정상범위를 회색 띠로 깐다."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(4, 1, figsize=(9, 10), sharex=True)
    t = id_mot.time

    for ax, (joint, prefix) in zip(axes[:3], JOINTS.items()):
        lo, hi, unit, _ = BANDS[joint]
        if lo is not None:
            ax.axhspan(lo, hi, color="0.85", zorder=0)
            ax.axhspan(-hi, -lo, color="0.85", zorder=0)
        for side, color in (("r", "tab:red"), ("l", "tab:blue")):
            col = f"{prefix}_{side}_moment"
            if col in id_mot.names:
                ax.plot(t, id_mot.column(col) / mass_kg, color=color, label=side)
        ax.set_ylabel(f"{joint}\n[{unit}]")
        ax.axhline(0, color="0.6", lw=0.6)
        ax.legend(loc="upper right", fontsize=8)

    ax = axes[3]
    lo, hi, _, _ = BANDS["grf_peak"]
    ax.axhspan(lo, hi, color="0.85", zorder=0)
    for key, plate in plates.items():
        if "y" in plate.force:
            ax.plot(grf_mot.time, grf_mot.column(plate.force["y"]) / weight_n, label=key)
    ax.set_ylabel("vGRF\n[BW]")
    ax.set_xlabel("time [s]")
    ax.legend(loc="upper right", fontsize=8)

    if window:
        for ax in axes:
            ax.axvline(window[0], color="tab:green", ls="--", lw=0.8)
            ax.axvline(window[1], color="tab:green", ls="--", lw=0.8)
    axes[0].set_title(os.path.basename(path).replace(".png", ""))
    fig.tight_layout()
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    fig.savefig(path, dpi=110)
    plt.close(fig)
