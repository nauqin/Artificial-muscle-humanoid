"""기계연 1자유도 팔 + SMA 근육 두 개 (이두 · 삼두) — 새 설정 (11/20, 김은서)
예전 설정 kimm_arm.py 에서 두 가지만 바꿈 (둘 다 가설, 기계연 확인 안 됨):
  ① 근육 다는 방식 = 교수님 elbow_rig.py: 팔을 쭉 폈을 때 (0°) 두 근육 모두 최대 길이 → 굽히면 삼두가 최대 길이를 넘어 다 식어도 버팀
  ② 팔 실험 팬 냉각 = 근육 하나 시험 냉각의 0.6 배 (이두를 끈 뒤 팔이 내려오는 모양으로 정함, 굽힘 크기는 안 봄)
돌리기:  python kimm_arm_new.py      (필요: mujoco, numpy, 같은 폴더의 sma_model.py)
출력:    팔 실험 빠르기 4개 (0.125 / 0.25 / 0.5 / 1 Hz) 의 굽힘 각도 vs 기계연 실측

[실제 값 — 기계연 회신 · 원본 데이터]
  근육: 길이 0.18 m, 최대 수축 0.40, 힘 식 T3.1  (sma_model.M3)
  열:   근육 하나 실험 (5 kg, 0.5 Hz, 36 A) 에 맞춘 값 — 가열 267 W @ 36 A, 강제 냉각 4.62 W/K  (sma_model.KIMM_HEAT)
        팔 실험 전류 25 / 29 / 29 / 29.5 A → 열 = 267 × (전류 / 36)²
        ※ 팔 실험은 맞추는 데 안 씀 (확인용)
  팔:   모멘트 팔 (팔꿈치 ~ 와이어 부착점) 50 mm, 팔꿈치 ~ 하중 320 mm, 하중 0.5 kg
  냉각: 블레이드 시간표 — 한 주기 = 데우기 1/4 → 반냉각 1/4 → 강제 냉각 1/2, 삼두는 반 주기 늦게
        데우는 근육은 자연 냉각, 반대쪽은 강제, 반냉각 구간은 두 근육 모두 강제의 절반 (블레이드가 가운데)

[가정 값 — 회신에 없음]
  켜는 동안 냉각 0.72 W/K · 전완 0.2 kg (막대 가운데) · 관절 감쇠 0.02 (베어링)
  삼두 모멘트 팔 = 이두와 같은 50 mm · 모멘트 팔이 각도와 상관없이 일정 (도르래처럼)
  근육 다는 방식: 0° 에서 이두 · 삼두 모두 최대 길이 (교수님 elbow_rig.py 와 같음) · 팔 팬 냉각 0.6 배
  삼두 전류 = 이두와 같음 (회신 전류는 이두만 있음)

결과 (이 설정 그대로): 굽힘 57.7 / 48.2 / 27.9 / 9.7°  vs  실측 59.8 / 50.4 / 35.6 / 14.1°  (평균 오차 4.7°)
  (예전 설정 kimm_arm.py 는 77.9 / 73.2 / 32.8 / 1.4°, 15.9°)
가정 값 바꿔 보기 (환경 변수): MF=전완 질량, DAMP=감쇠, CX=팬 냉각 배율 (기본 0.6)
"""
import os, math, numpy as np, mujoco, sma_model as sm

CX = float(os.environ.get('CX', 0.6))                          # ② 팔 실험 팬 냉각 배율 [가정, 내려오는 모양으로 정함]
P36, HF = sm.KIMM_HEAT["P"], sm.KIMM_HEAT["h_forced"] * CX     # 36 A 일 때 열 [W], 팔 실험 강제 냉각 [W/K]
H_ON, H_HALF = 0.72, HF / 2                                    # 켜는 동안 [가정], 반냉각 = 절반
L0, EPSL = sm.M3["L0"], sm.M3["epsL"]
R, LOAD, ARM = 0.05, 0.5, 0.32                                 # 모멘트 팔, 하중, 팔꿈치 ~ 하중 [m, kg, m]
MF = float(os.environ.get('MF', 0.2)); DAMP = float(os.environ.get('DAMP', 0.02))   # [가정]
QMAX = L0 * EPSL / R                                           # 각도 범위 위쪽 = 이두가 다 줄어드는 각도 (82.5°)
CUR = {0.125: 25.0, 0.25: 29.0, 0.5: 29.0, 1.0: 29.5}           # 기계연 팔 실험 전류 [A]
MEAS = {0.125: 59.8, 0.25: 50.4, 0.5: 35.6, 1.0: 14.1}          # 기계연 팔 실험 굽힘 [°]

XML = f"""<mujoco><compiler angle="radian"/><option timestep="0.002" gravity="0 0 -9.81"/>
<worldbody><body name="fore" pos="0 0 0">
<joint name="elbow" type="hinge" axis="0 -1 0" range="0 {QMAX}" limited="true" damping="{DAMP}"/>
<geom type="capsule" fromto="0 0 0 0 0 -{ARM}" size="0.01" mass="{MF}"/>
<geom type="sphere" pos="0 0 -{ARM}" size="0.03" mass="{LOAD}"/>
</body></worldbody>
<tendon><fixed name="bic"><joint joint="elbow" coef="-{R}"/></fixed><fixed name="tri"><joint joint="elbow" coef="{R}"/></fixed></tendon>
<actuator><general tendon="bic" ctrlrange="0 1" ctrllimited="true" gainprm="-1000" biastype="none"/>
<general tendon="tri" ctrlrange="0 1" ctrllimited="true" gainprm="-1000" biastype="none"/></actuator></mujoco>"""
m = mujoco.MjModel.from_xml_string(XML); DT = m.opt.timestep; GAIN = 1000.0

def phase(q):
    """한 주기 안 위치 q (0~1) → (켜기, 냉각 세기) 이두 · 삼두"""
    b = (True, H_ON) if q < .25 else (False, H_HALF) if q < .5 else (False, H_HALF if q >= .75 else HF)
    t = (False, H_HALF if .25 <= q < .5 else HF) if q < .5 else (True, H_ON) if q < .75 else (False, H_HALF)
    return b, t

def run(f, n=None):
    n = n or max(6, int(round(30 * f))); per = 1 / f; P = P36 * (CUR[f] / 36.0) ** 2
    mus = [sm.SMA(), sm.SMA()]                                                # 이두, 삼두
    d = mujoco.MjData(m); mujoco.mj_forward(m, d); t = 0.0; L = []
    while t < n * per:
        cmd = phase((t % per) / per); q = d.qpos[0]
        eps = (max(EPSL - R * q / L0, 0.0), EPSL + R * q / L0)       # ① 이두 · 삼두 변형률: 0° 에서 둘 다 0.40, 굽히면 삼두가 0.40 넘음
        for j in range(2):
            on, h = cmd[j]
            F = mus[j].step(P if on else 0.0, DT, eps=eps[j], hS=h)
            d.ctrl[j] = min(max(F, 0.0) / GAIN, 1.0)
        mujoco.mj_step(m, d); t += DT; L.append((t, math.degrees(d.qpos[0]), mus[0].T, mus[1].T))
    L = np.array(L); r = []
    for k in range(n):
        s = (L[:, 0] >= k * per) & (L[:, 0] < (k + 1) * per); r.append((L[s, 1].max() - L[s, 1].min(), L[s, 1].max(), L[s, 1].min()))
    r = np.array(r)                                            # 처음 두 주기는 빼고 평균
    return dict(P=P, bend=r[2:, 0].mean(), hi=r[2:, 1].mean(), lo=r[2:, 2].mean(), Tmax=L[:, 2:].max(), log=L)

if __name__ == "__main__":
    out = []
    for f in MEAS:
        r = run(f); out.append(r['bend'])
        print(f"{f:5.3f} Hz | 전류 {CUR[f]} A → {r['P']:4.0f} W | 굽힘 {r['bend']:5.1f}° (실측 {MEAS[f]}°) | 위/아래 {r['hi']:5.1f}/{r['lo']:5.1f}° | 최고 {r['Tmax']:5.1f} °C", flush=True)
    e = np.sqrt(np.mean([(a - MEAS[f]) ** 2 for a, f in zip(out, MEAS)]))
    print(f"새 설정 (교수님 방식 · 팬 냉각 ×{CX}) · 전완 {MF} kg · 감쇠 {DAMP} | 평균 오차 {e:.1f}°")
