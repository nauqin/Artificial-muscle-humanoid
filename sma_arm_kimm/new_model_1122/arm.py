"""기계연 1자유도 팔 + SMA 근육 두 개 (이두 · 삼두) — 가설 바꿔 보기용 (11/20, 김은서)
팀 코드 sma_kimm/kimm_arm.py 와 같은 팔 · 근육 · 시간표에, 아래 네 가지를 옵션으로 바꿀 수 있게 한 것

  --mode ours   지금 팀 코드: 각 근육이 반대쪽 끝 각도에서 딱 0.40 (다 식은 신전근 저항 없음)  [기본]
  --mode prof   교수님 elbow_rig.py 방식: 중립 각도 --qmid 에서 두 근육 모두 0.40 → 굽히면 삼두가 0.40 넘게 늘어나 다 식어도 버팀
  --qmid 0      (prof) 중립 각도 [°]
  --qt 82.5     (ours) 삼두가 다 펴지는 각도 [°] — 작게 하면 그 각도 넘을 때만 다 식은 삼두가 버팀
  --hx 1        팔 실험 열 배율 (전류 RMS 를 잰 구간이 다르면 ×2 · ×4)
  --cx 1        팔 실험 팬 냉각 배율 (강제 · 반냉각만, 켜는 동안 냉각 0.72 는 그대로)
  --damp 0.02   관절 감쇠 [N·m·s/rad] (가정) · --fric 0  관절 마찰 [N·m] (가정, 교수님 elbow_rig 은 감쇠 0.3 · 마찰 0.2)

돌리기 예:  python arm.py                          → 지금 팀 코드 (평균 오차 15.9°)
            python arm.py --mode prof               → 교수님 방식 (12.7°)
            python arm.py --mode prof --cx 0.6      → 교수님 방식 + 팬 냉각 0.6 (4.7°)
필요: pip install mujoco numpy  (같은 폴더에 sma_model.py)

[실제 값] 근육 0.18 m · 최대 수축 0.40 · T3.1 힘 식, 가열 267 W @ 36 A · 강제 냉각 4.62 W/K (근육 하나 시험에 맞춤),
          팔 실험 전류 25 / 29 / 29 / 29.5 A, 모멘트 팔 50 mm, 하중 0.5 kg @ 320 mm, 블레이드 냉각 시간표
[가정 값] 켜는 동안 냉각 0.72 W/K · 전완 0.2 kg · 감쇠 0.02 · 삼두 모멘트 팔 = 이두 · 삼두 전류 = 이두"""
import argparse, math, numpy as np, mujoco, sma_model as sm

CUR = {0.125: 25.0, 0.25: 29.0, 0.5: 29.0, 1.0: 29.5}           # 기계연 팔 실험 전류 [A]
MEAS = {0.125: 59.8, 0.25: 50.4, 0.5: 35.6, 1.0: 14.1}          # 기계연 팔 실험 굽힘 [°]
L0, EPSL = sm.M3["L0"], sm.M3["epsL"]
R, LOAD, ARM, MF, DAMP = 0.05, 0.5, 0.32, 0.2, 0.02             # 모멘트 팔 · 하중 · 하중 거리 [m, kg, m] · 전완 [kg] · 감쇠 [가정]
QFULL = L0 * EPSL / R                                           # 근육 하나가 다 줄 때 각도 = 82.5°

class Arm:
    def __init__(s, mode="ours", qmid=0.0, qt=None, hx=1.0, cx=1.0, damp=DAMP, fric=0.0):
        s.mode, s.hx, s.cx, s.damp, s.fric = mode, hx, cx, damp, fric
        s.qmid = math.radians(qmid); s.qt = math.radians(qt) if qt is not None else QFULL
        s.qmax = s.qmid + QFULL if mode == "prof" else QFULL
        s.HF = sm.KIMM_HEAT["h_forced"] * cx; s.H_HALF = s.HF / 2; s.H_ON = 0.72
        xml = f"""<mujoco><compiler angle="radian"/><option timestep="0.002" gravity="0 0 -9.81"/>
<worldbody><body name="fore" pos="0 0 0">
<joint name="elbow" type="hinge" axis="0 -1 0" range="0 {s.qmax}" limited="true" damping="{damp}" frictionloss="{fric}"/>
<geom type="capsule" fromto="0 0 0 0 0 -{ARM}" size="0.01" mass="{MF}"/>
<geom type="sphere" pos="0 0 -{ARM}" size="0.03" mass="{LOAD}"/>
</body></worldbody>
<tendon><fixed name="bic"><joint joint="elbow" coef="-{R}"/></fixed><fixed name="tri"><joint joint="elbow" coef="{R}"/></fixed></tendon>
<actuator><general tendon="bic" ctrlrange="0 1" ctrllimited="true" gainprm="-1000" biastype="none"/>
<general tendon="tri" ctrlrange="0 1" ctrllimited="true" gainprm="-1000" biastype="none"/></actuator></mujoco>"""
        s.m = mujoco.MjModel.from_xml_string(xml)

    def strains(s, q):
        """팔 각도 → 이두 · 삼두 변형률 (0.40 = 다 펴짐, 넘으면 다 식어도 버팀)"""
        if s.mode == "prof":
            return max(EPSL - R * (q - s.qmid) / L0, 0.0), max(EPSL + R * (q - s.qmid) / L0, 0.0)
        return max(EPSL - R * q / L0, 0.0), max(EPSL - R * (s.qt - q) / L0, 0.0)

    def phase(s, x):
        """한 주기 안 위치 x (0~1) → (켜기, 냉각 세기) 이두 · 삼두. 블레이드: 데우기 1/4 → 반냉각 1/4 → 강제 1/2"""
        b = (True, s.H_ON) if x < .25 else (False, s.H_HALF) if x < .5 else (False, s.H_HALF if x >= .75 else s.HF)
        t = (False, s.H_HALF if .25 <= x < .5 else s.HF) if x < .5 else (True, s.H_ON) if x < .75 else (False, s.H_HALF)
        return b, t

    def run(s, f, n=None):
        n = n or max(6, int(round(30 * f))); per = 1 / f; P = sm.heat_from_current(CUR[f]) * s.hx
        mus = [sm.SMA(), sm.SMA()]; d = mujoco.MjData(s.m); mujoco.mj_forward(s.m, d); t = 0.0; L = []; DT = s.m.opt.timestep
        while t < n * per:
            cmd = s.phase((t % per) / per); eps = s.strains(d.qpos[0])
            for j in range(2):
                on, h = cmd[j]; F = mus[j].step(P if on else 0.0, DT, eps=eps[j], hS=h); d.ctrl[j] = min(max(F, 0.0) / 1000.0, 1.0)
            mujoco.mj_step(s.m, d); t += DT; L.append((t, math.degrees(d.qpos[0]), mus[0].T, mus[1].T))
        L = np.array(L); r = []
        for k in range(n):
            sl = (L[:, 0] >= k * per) & (L[:, 0] < (k + 1) * per); r.append((L[sl, 1].max() - L[sl, 1].min(), L[sl, 1].max(), L[sl, 1].min()))
        r = np.array(r)                                         # 처음 두 주기는 빼고 평균
        return dict(P=P, bend=r[2:, 0].mean(), hi=r[2:, 1].mean(), lo=r[2:, 2].mean(), Tmax=L[:, 2:].max(), log=L, first=r[:5, 0])

    def label(s):
        g = f"교수님 방식 중립 {math.degrees(s.qmid):.0f}°" if s.mode == "prof" else f"지금 방식 (삼두 {math.degrees(s.qt):.1f}° 에서 다 펴짐)"
        return f"{g} · 열 ×{s.hx:g} · 팬 냉각 ×{s.cx:g} · 감쇠 {s.damp:g} · 마찰 {s.fric:g}"

def err(bends):
    return float(np.sqrt(np.mean([(b - MEAS[f]) ** 2 for b, f in zip(bends, MEAS)])))

def args():
    ap = argparse.ArgumentParser(); ap.add_argument("--mode", choices=["ours", "prof"], default="ours")
    ap.add_argument("--qmid", type=float, default=0.0); ap.add_argument("--qt", type=float, default=None)
    ap.add_argument("--hx", type=float, default=1.0); ap.add_argument("--cx", type=float, default=1.0)
    ap.add_argument("--damp", type=float, default=DAMP); ap.add_argument("--fric", type=float, default=0.0); return ap.parse_args()

if __name__ == "__main__":
    a = args(); A = Arm(a.mode, a.qmid, a.qt, a.hx, a.cx, a.damp, a.fric); out = []
    for f in MEAS:
        r = A.run(f); out.append(r["bend"])
        print(f"{f:5.3f} Hz | 전류 {CUR[f]} A → {r['P']:4.0f} W | 굽힘 {r['bend']:5.1f}° (실측 {MEAS[f]}°) | 위/아래 {r['hi']:5.1f}/{r['lo']:5.1f}° | 최고 {r['Tmax']:5.1f} °C | 처음 5회 {np.round(r['first'], 1).tolist()}", flush=True)
    print(f"{A.label()} | 평균 오차 {err(out):.1f}°")
