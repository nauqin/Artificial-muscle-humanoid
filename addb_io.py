"""파일 탐색·분류, .mot 파싱, GRF-발 자동 매핑, 지지 구간 탐지.

이 모듈의 방침은 하나다 — **이름표를 믿지 않는다.**

- GRF 파일인지는 컬럼 이름에 'force'가 들어있는지가 아니라 힘 x/y/z 와
  작용점 x/y/z 가 모두 있는지로 판정한다.
- 힘판이 어느 발에 붙는지는 컬럼 이름(`..._r_...`)이 아니라 COP와 발 사이의
  실제 수평 거리로 정한다.

둘 다 틀려도 계산은 에러 없이 그럴듯한 그래프를 내놓기 때문이다.
"""

from __future__ import annotations

import os
import re
import sys
from dataclasses import dataclass, field

import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# 상수
# ─────────────────────────────────────────────────────────────────────────────

#: 지지 구간 판정 문턱 (총 수직 GRF / 체중)
SUPPORT_THRESHOLD = 0.70
#: midstance 골로 끊긴 두 구간을 이어붙일 최대 간격 (초)
SUPPORT_GAP_MAX = 0.25
#: 그 간격 동안 총 수직 GRF가 이보다 높게 유지되어야 이어붙인다 (체중 배수)
SUPPORT_GAP_FLOOR = 0.50
#: 힘판-발 매핑에 쓸 프레임의 수직력 하한 (체중 배수)
CONTACT_THRESHOLD = 0.10

#: 파일 이름에서 trial을 구분하는 데 도움이 안 되는 토큰
GENERIC_TOKENS = {
    "ik", "id", "grf", "mot", "sto", "results", "result", "output", "outputs",
    "motion", "motions", "kinematics", "coordinates", "coords", "forces",
    "force", "data", "subject", "scaled", "model", "states", "inverse",
    "kinematic", "dynamics", "analyze", "processed", "filtered", "final",
}

#: 모델 자동 선택에서 뒤로 미룰 템플릿 이름 조각
TEMPLATE_HINTS = (
    "unscaled", "generic", "template", "rajagopal2016", "gait2392", "gait2354",
    "rajagopallaiuhlrich2023", "preScale".lower(),
)

_AXES = ("x", "y", "z")


# ─────────────────────────────────────────────────────────────────────────────
# .mot / .sto 읽고 쓰기 (OpenSim 없이 동작한다)
# ─────────────────────────────────────────────────────────────────────────────


@dataclass
class Mot:
    """`.mot`/`.sto` 한 개."""

    path: str
    header: dict = field(default_factory=dict)
    names: list = field(default_factory=list)   # time 을 제외한 컬럼 이름
    time: np.ndarray = None
    data: np.ndarray = None                     # (T, len(names))

    @property
    def in_degrees(self) -> bool:
        return str(self.header.get("indegrees", "no")).strip().lower() in ("yes", "true", "1")

    @property
    def t0(self) -> float:
        return float(self.time[0])

    @property
    def t1(self) -> float:
        return float(self.time[-1])

    @property
    def duration(self) -> float:
        return self.t1 - self.t0

    def column(self, name: str) -> np.ndarray:
        return self.data[:, self.names.index(name)]

    def at(self, t) -> np.ndarray:
        """시간 t(스칼라 또는 배열)에서 전 컬럼을 선형보간."""
        t = np.atleast_1d(np.asarray(t, dtype=float))
        out = np.empty((t.size, len(self.names)))
        for j in range(len(self.names)):
            out[:, j] = np.interp(t, self.time, self.data[:, j])
        return out


def read_mot(path: str) -> Mot:
    """`.mot`/`.sto` 파싱. endheader 앞은 `key=value` 헤더로 읽는다."""
    header, labels, rows = {}, None, []
    with open(path, "r", errors="replace") as fh:
        in_header = True
        for raw in fh:
            line = raw.strip()
            if in_header:
                if not line:
                    continue
                if line.lower() == "endheader":
                    in_header = False
                    continue
                if "=" in line:
                    k, v = line.split("=", 1)
                    header[k.strip().lower()] = v.strip()
                continue
            if labels is None:
                labels = line.split()
                continue
            if not line:
                continue
            rows.append([float(x) for x in line.split()])

    if labels is None or not rows:
        raise ValueError(f"데이터가 없다: {path}")

    arr = np.asarray(rows, dtype=float)
    if arr.shape[1] != len(labels):
        raise ValueError(f"컬럼 수 불일치 ({arr.shape[1]} vs {len(labels)}): {path}")

    return Mot(path=path, header=header, names=labels[1:], time=arr[:, 0], data=arr[:, 1:])


def write_mot(path: str, names, time, data, name=None, in_degrees=False):
    """OpenSim이 읽는 `.mot`으로 저장."""
    time = np.asarray(time, dtype=float)
    data = np.asarray(data, dtype=float)
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write(f"{name or os.path.splitext(os.path.basename(path))[0]}\n")
        fh.write("version=1\n")
        fh.write(f"nRows={data.shape[0]}\n")
        fh.write(f"nColumns={data.shape[1] + 1}\n")
        fh.write(f"inDegrees={'yes' if in_degrees else 'no'}\n")
        fh.write("endheader\n")
        fh.write("time\t" + "\t".join(names) + "\n")
        for i in range(data.shape[0]):
            fh.write(f"{time[i]:.8f}\t" + "\t".join(f"{v:.8f}" for v in data[i]) + "\n")


# ─────────────────────────────────────────────────────────────────────────────
# GRF 컬럼 해석
# ─────────────────────────────────────────────────────────────────────────────

#: plate key를 만들 때 지울 토큰. 'm'(자유모멘트)을 여기 넣어야 v/p/m 세 항목이
#: 같은 plate로 묶인다. 이걸 빠뜨리면 자유모멘트가 힘으로 분류되어 나중에
#: 지면반력을 덮어쓴다 — 에러는 나지 않는다.
_KIND_TOKENS = {
    "v": "force", "f": "force", "force": "force",
    "p": "point", "point": "point", "cop": "point",
    "m": "torque", "torque": "torque", "moment": "torque",
}
_DROP_TOKENS = {"ground", "grf", "vec", "vector"} | set(_KIND_TOKENS)


def parse_grf_column(label: str):
    """GRF 컬럼 이름 하나를 (plate_key, kind, axis)로. GRF가 아니면 None.

    지원하는 규약:
      ground_force_vx / ground_force_px / ground_torque_x       (OpenSim 기본)
      ground_force_r_vx / ground_torque_r_x                     (좌우 접미)
      1_ground_force_vx / 1_ground_torque_x                     (번호 접두)
      ground_force_calcn_r_vx / _px / _mx                       (AddBiomechanics)
    """
    low = label.strip().lower()
    if not low or low[-1] not in _AXES:
        return None
    axis = low[-1]
    stem = low[:-1].rstrip("_")
    if not stem:
        return None

    tokens = [t for t in re.split(r"[_\s.]+", stem) if t]
    if not tokens:
        return None

    # kind는 뒤에서부터 찾는다 — 마지막 토큰(v/p/m)이 가장 구체적이다.
    kind = None
    for tok in reversed(tokens):
        if tok in _KIND_TOKENS:
            kind = _KIND_TOKENS[tok]
            break
    if kind is None:
        return None

    key = "_".join(t for t in tokens if t not in _DROP_TOKENS)
    return key, kind, axis


@dataclass
class Plate:
    """힘판 하나 — 힘/작용점/자유모멘트 컬럼 묶음."""

    key: str
    force: dict = field(default_factory=dict)    # axis -> 컬럼 이름
    point: dict = field(default_factory=dict)
    torque: dict = field(default_factory=dict)

    @property
    def complete(self) -> bool:
        return all(a in self.force for a in _AXES) and all(a in self.point for a in _AXES)

    def identifier(self, kind: str):
        """ExternalForce가 쓰는 접두사. 컬럼 이름에서 축 문자만 뗀 것."""
        cols = getattr(self, kind)
        if not all(a in cols for a in _AXES):
            return None
        prefixes = {cols[a][: -1] for a in _AXES}
        return prefixes.pop() if len(prefixes) == 1 else None


def find_plates(mot: Mot) -> dict:
    """`.mot` 안의 힘판들을 찾는다. key -> Plate."""
    plates = {}
    for label in mot.names:
        parsed = parse_grf_column(label)
        if parsed is None:
            continue
        key, kind, axis = parsed
        plate = plates.setdefault(key, Plate(key=key))
        getattr(plate, kind)[axis] = label
    return plates


def classify_mot(path: str):
    """`.mot`을 'grf' / 'kinematics' / 'other'로 분류.

    문자열 'force' 검색은 쓰지 않는다 — states 파일의
    `/forceset/soleus_r/activation` 같은 라벨이 GRF로 오분류된다.
    힘 x/y/z 와 작용점 x/y/z 가 **모두** 있어야 GRF다. IK 파일의
    pelvis_tx/ty/tz 는 x/y/z 삼총사를 이루지만 작용점이 없어서 걸러진다.
    """
    try:
        mot = read_mot(path)
    except Exception:
        return "other", None

    plates = find_plates(mot)
    if any(p.complete for p in plates.values()):
        return "grf", mot

    lower = [n.lower() for n in mot.names]
    pelvis = sum(1 for n in lower if n.startswith("pelvis_"))
    joints = sum(1 for n in lower if any(k in n for k in ("hip_", "knee_", "ankle_", "lumbar_")))
    if pelvis >= 3 and joints >= 4:
        return "kinematics", mot
    return "other", mot


# ─────────────────────────────────────────────────────────────────────────────
# 지지 구간
# ─────────────────────────────────────────────────────────────────────────────


def total_vertical(mot: Mot, plates: dict) -> np.ndarray:
    """모든 힘판의 수직력 합."""
    total = np.zeros_like(mot.time)
    for plate in plates.values():
        if "y" in plate.force:
            total = total + mot.column(plate.force["y"])
    return total


def _runs(mask: np.ndarray):
    """True 구간들의 (시작, 끝) 인덱스(포함)."""
    out, start = [], None
    for i, v in enumerate(mask):
        if v and start is None:
            start = i
        elif not v and start is not None:
            out.append((start, i - 1))
            start = None
    if start is not None:
        out.append((start, len(mask) - 1))
    return out


def grf_support_window(time, total_vy, weight_n, *, threshold=SUPPORT_THRESHOLD,
                       gap_max=SUPPORT_GAP_MAX, gap_floor=SUPPORT_GAP_FLOOR):
    """총 수직 GRF가 체중의 `threshold` 이상인 최장 연속 구간 (t0, t1).

    보행의 총 수직 GRF는 쌍봉이고 midstance에 골이 있다. 그 골이 문턱을 스치면
    지지 구간이 반으로 쪼개져 한쪽 다리의 push-off가 통째로 잘린다. 그래서
    **짧고(≤gap_max) 얕은(≥gap_floor) 골은 이어붙인다.**

    히스테리시스(들어갈 때 0.7, 나올 때 0.3)로는 안 된다 — 골은 고쳐지지만 양
    끝까지 같이 늘어나서, 피험자가 힘판을 빠져나가는 구간이 딸려 들어온다.
    """
    time = np.asarray(time, dtype=float)
    total_vy = np.asarray(total_vy, dtype=float)
    if weight_n <= 0:
        return float(time[0]), float(time[-1])

    runs = _runs(total_vy >= threshold * weight_n)
    if not runs:
        return float(time[0]), float(time[-1])

    merged = [list(runs[0])]
    for start, end in runs[1:]:
        prev_end = merged[-1][1]
        gap = total_vy[prev_end + 1: start]
        gap_dt = time[start] - time[prev_end]
        if gap_dt <= gap_max and (gap.size == 0 or gap.min() >= gap_floor * weight_n):
            merged[-1][1] = end          # midstance 골 — 이어붙인다
        else:
            merged.append([start, end])

    start, end = max(merged, key=lambda r: time[r[1]] - time[r[0]])
    return float(time[start]), float(time[end])


# ─────────────────────────────────────────────────────────────────────────────
# 모델 자세 (OpenSim 필요)
# ─────────────────────────────────────────────────────────────────────────────


class Pose:
    """IK 결과로 모델을 세워 바디 위치를 재는 도구."""

    def __init__(self, model, ik_path: str):
        import opensim as osim

        self.osim = osim
        self.model = model
        # initSystem() 을 먼저 불러야 한다. 그 전에 getSimbodyEngine() 을 건드리면
        # 트레이스백도 없이 세그폴트로 죽는다.
        self.state = model.initSystem()

        sto = osim.Storage(ik_path)
        if sto.isInDegrees():
            model.getSimbodyEngine().convertDegreesToRadians(sto)

        labels = sto.getColumnLabels()
        names = [labels.get(i) for i in range(labels.getSize())][1:]
        n = sto.getSize()
        time = np.empty(n)
        data = np.empty((n, len(names)))
        for i in range(n):
            sv = sto.getStateVector(i)
            time[i] = sv.getTime()
            vec = sv.getData()
            for j in range(min(vec.getSize(), len(names))):
                data[i, j] = vec.get(j)

        self.time, self.data, self.names = time, data, names
        coord_set = model.getCoordinateSet()
        model_coords = {coord_set.get(i).getName(): i for i in range(coord_set.getSize())}
        self.pairs = [(model_coords[nm], j) for j, nm in enumerate(names) if nm in model_coords]

    def set(self, t: float):
        coord_set = self.model.getCoordinateSet()
        for ci, j in self.pairs:
            value = float(np.interp(t, self.time, self.data[:, j]))
            coord_set.get(ci).setValue(self.state, value, False)
        self.model.assemble(self.state)
        self.model.realizePosition(self.state)
        return self.state

    def body_position(self, body_name: str) -> np.ndarray:
        body = self.model.getBodySet().get(body_name)
        p = body.getPositionInGround(self.state)
        return np.array([p.get(0), p.get(1), p.get(2)])

    def station_position(self, body_name: str, offset) -> np.ndarray:
        """바디 좌표계의 한 점을 지면 좌표계로. 몸통 위쪽처럼 원점이 아닌 곳을 잡을 때."""
        body = self.model.getBodySet().get(body_name)
        local = self.osim.Vec3(float(offset[0]), float(offset[1]), float(offset[2]))
        p = body.findStationLocationInGround(self.state, local)
        return np.array([p.get(0), p.get(1), p.get(2)])


# ─────────────────────────────────────────────────────────────────────────────
# GRF-발 매핑
# ─────────────────────────────────────────────────────────────────────────────


def assign_plates_to_feet(model, ik_path, grf: Mot, plates: dict, feet=("calcn_r", "calcn_l"),
                          weight_n=None, max_frames=40, log=print):
    """힘판을 발에 배정한다. 컬럼 이름은 **보지 않는다.**

    수직력이 체중의 10%를 넘는 프레임에서 모델을 IK 자세로 세우고, COP와 각 발의
    수평 거리를 재서 가까운 쪽에 배정한다. 택배 주소 라벨 대신 GPS를 보는 셈.

    반환: {plate_key: {"body", "dist", "margin", "frames"}}
    """
    if weight_n is None:
        weight_n = model_weight_n(model)

    usable = {k: p for k, p in plates.items() if p.complete}
    if not usable:
        return {}

    pose = Pose(model, ik_path)
    result = {}

    for key, plate in usable.items():
        fy = grf.column(plate.force["y"])
        cop = np.column_stack([grf.column(plate.point[a]) for a in _AXES])
        idx = np.flatnonzero(fy > CONTACT_THRESHOLD * weight_n)
        if idx.size == 0:
            log(f"  [WARN] plate '{key}': 접지 프레임이 없다 — 건너뜀")
            continue
        if idx.size > max_frames:
            idx = idx[np.linspace(0, idx.size - 1, max_frames).astype(int)]

        sums = {f: 0.0 for f in feet}
        for i in idx:
            pose.set(float(grf.time[i]))
            for foot in feet:
                fp = pose.body_position(foot)
                sums[foot] += float(np.hypot(cop[i, 0] - fp[0], cop[i, 2] - fp[2]))
        means = {f: sums[f] / idx.size for f in feet}
        order = sorted(means, key=means.get)
        best, second = order[0], order[1] if len(order) > 1 else order[0]
        result[key] = {
            "body": best,
            "dist": means[best],
            "margin": means[second] - means[best],
            "means": means,
            "frames": int(idx.size),
        }

    _resolve_conflicts(result, feet, log)

    for key, info in result.items():
        detail = "  ".join(f"{f}={info['means'][f]:.3f}" for f in feet)
        log(f"  {key:>10s} -> {info['body']}   {detail}   +{info['margin']:.3f}")
    return result


def _resolve_conflicts(result, feet, log):
    """두 힘판이 같은 발을 고르면, 총 거리가 작은 조합으로 바꾼다."""
    if len(result) != 2 or len(feet) != 2:
        return
    (k1, a), (k2, b) = result.items()
    if a["body"] != b["body"]:
        return
    f0, f1 = feet
    cost_same = a["means"][f0] + b["means"][f1]
    cost_swap = a["means"][f1] + b["means"][f0]
    pick = (f0, f1) if cost_same <= cost_swap else (f1, f0)
    log(f"  [WARN] 두 힘판이 같은 발을 골랐다 — 총 거리로 재배정 ({pick[0]}, {pick[1]})")
    for info, foot in ((a, pick[0]), (b, pick[1])):
        other = f1 if foot == f0 else f0
        info["body"] = foot
        info["dist"] = info["means"][foot]
        info["margin"] = info["means"][other] - info["means"][foot]


def model_weight_n(model, gravity=9.80665) -> float:
    """모델 체중 (N)."""
    return float(total_mass(model) * gravity)


def total_mass(model) -> float:
    """모델 총 질량 (kg)."""
    bodies = model.getBodySet()
    return float(sum(bodies.get(i).getMass() for i in range(bodies.getSize())))


# ─────────────────────────────────────────────────────────────────────────────
# 폴더 탐색
# ─────────────────────────────────────────────────────────────────────────────


def _tokens(stem: str) -> set:
    raw = {t for t in re.split(r"[^a-z0-9]+", stem.lower()) if t}
    return {t for t in raw if t not in GENERIC_TOKENS}


def find_models(root: str) -> list:
    """폴더 아래 `.osim`을 전부. 템플릿으로 보이는 것은 뒤로 민다."""
    found = []
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in filenames:
            if fn.lower().endswith(".osim"):
                found.append(os.path.join(dirpath, fn))

    def rank(path):
        low = os.path.basename(path).lower()
        template = any(h in low for h in TEMPLATE_HINTS)
        scaled = "scaled" in low or "subject" in low
        return (1 if template else 0, 0 if scaled else 1, len(path))

    return sorted(found, key=rank)


def pick_model(root: str, log=print):
    """모델 자동 선택. 고른 모델의 체중을 반드시 찍는다.

    한 번 배포판 기본 모형(Rajagopal2016.osim)을 집어 잔차가 84 %까지 뛴 적이
    있다. 체중이 예상과 다르면 모델을 잘못 고른 것이다.
    """
    models = find_models(root)
    if not models:
        return None
    if len(models) > 1:
        log(f"  ({len(models)}개 중 선택 — 나머지: "
            + ", ".join(os.path.basename(m) for m in models[1:4]) + ")")
    return models[0]


def find_mots(root: str) -> list:
    out = []
    for dirpath, _dirnames, filenames in os.walk(root):
        for fn in filenames:
            if fn.lower().endswith((".mot", ".sto")):
                out.append(os.path.join(dirpath, fn))
    return sorted(out)


@dataclass
class Trial:
    name: str
    ik: str
    grf: str


def find_trials(root: str, min_duration=0.5, log=print) -> list:
    """subject 폴더에서 (IK, GRF) 짝을 찾는다.

    짝을 확정할 수 없으면 **짝짓지 않고 건너뛴다.** 엉뚱한 조합으로 계산해
    조용히 틀리는 것보다 빠뜨리는 쪽이 낫다.
    """
    kinematics, grfs = [], []
    for path in find_mots(root):
        kind, mot = classify_mot(path)
        if mot is None or mot.duration < min_duration:
            continue
        if kind == "kinematics":
            kinematics.append(path)
        elif kind == "grf":
            grfs.append(path)

    trials = []
    # IK도 GRF도 하나뿐일 때만 이름을 안 보고 묶는다. GRF가 여럿이면 반드시
    # 토큰을 확인한다 — grf_run.mot 이 walk IK 와 짝지어진 적이 있다.
    unambiguous = len(kinematics) == 1 and len(grfs) == 1
    for grf_path in grfs:
        grf_tokens = _tokens(_stem(grf_path))
        if unambiguous:
            best = kinematics[0]
        else:
            scored = [(len(grf_tokens & _tokens(_stem(k))), k) for k in kinematics]
            scored.sort(key=lambda s: (-s[0], s[1]))
            if not scored:
                log(f"  [WARN] IK 후보가 없다 — 건너뜀: {os.path.basename(grf_path)}")
                continue
            if scored[0][0] == 0:
                log(f"  [WARN] 짝을 확정할 수 없어 건너뜀: {os.path.basename(grf_path)}")
                continue
            if len(scored) > 1 and scored[0][0] == scored[1][0]:
                log(f"  [WARN] 짝이 모호해 건너뜀: {os.path.basename(grf_path)}")
                continue
            best = scored[0][1]
        trials.append(Trial(name=_trial_name(best, grf_path), ik=best, grf=grf_path))

    trials.sort(key=lambda t: t.name)
    return trials


def _stem(path: str) -> str:
    return os.path.splitext(os.path.basename(path))[0]


def _trial_name(ik_path: str, grf_path: str) -> str:
    """두 파일 이름의 공통 앞부분. `walking1_segment_0_ik` + `..._grf` → `walking1_segment_0`.

    공통 접두가 없으면(`ik_output_walk` + `grf_walk`) 공통 토큰으로 물러난다.
    토큰만 쓰면 순서가 뒤집히고(`0_segment_walking1`) 대소문자도 잃으므로
    (`walkingTS1` → `walkingts1`) 접두를 먼저 본다.
    """
    a, b = _stem(ik_path), _stem(grf_path)
    n = 0
    while n < min(len(a), len(b)) and a[n] == b[n]:
        n += 1
    common = a[:n].strip("_-. ")
    if len(common) >= 3:
        return common

    shared = _tokens(a) & _tokens(b)
    for stem in (b, a):
        ordered = [t for t in re.split(r"[^A-Za-z0-9]+", stem) if t.lower() in shared]
        if ordered:
            return "_".join(ordered)
    return a


def find_subjects(data_root: str) -> list:
    """`data/` 아래 subject 폴더들. 없으면 data_root 자체를 subject로 본다."""
    subs = [os.path.join(data_root, d) for d in sorted(os.listdir(data_root))
            if os.path.isdir(os.path.join(data_root, d)) and not d.startswith(".")]
    subs = [s for s in subs if find_models(s)]
    return subs or ([data_root] if find_models(data_root) else [])


def log_stderr(*args):
    print(*args, file=sys.stderr)
