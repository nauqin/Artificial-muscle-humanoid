# -*- coding: utf-8 -*-
"""②-3 유빈 최종 A안 오른팔로 물병 들기 동작 따라가기 (SMA 근육, MuJoCo)

목표 동작: 유빈 Moco 결과 gui_upper_robot40h160_finalA_lift_s6_PFWJ_2.mot (6.4 초, 물병 3.79 kg)
  - 팔꿈치 목표 각도 = 그 파일의 elbow_flex_r (6° ~ 143°)
  - 어깨·몸통이 움직이는 효과 = 매 순간 위팔 방향을 계산해서 "중력 방향"을 바꿔 넣음 (어깨가 빨리 흔드는 관성은 무시)
  - 비교용 "필요한 힘" = 같은 파일의 Moco 활성도 × 설계 힘 (근육이 이상적일 때 필요한 케이블 힘)

제어 (두 단계, 기계연 장치의 "전류제어 = 온도제어" 방식을 흉내냄)
  ① 바깥: 목표 각도와의 차이(PD) → 명령 −1 ~ 1. + 면 굽힘 근육, − 면 폄 근육을 데움
  ② 안쪽: 명령 → 목표 온도 (대기 52 °C ~ 최대 92 °C) → 온도 차이만큼 가열 (상한 UMAX = 32번째 줄, 발표 숫자는 0.35 · 보정값은 0.32)
     목표보다 뜨거우면 강제 냉각 4.29 W/K, 아니면 켜는 동안 냉각 0.72 W/K (기계연 13장 값)
  관절 마찰(감쇠) 5 N·m·s/rad [가정 — 실물 측정 필요. 0.3 이면 물병 때문에 크게 출렁임]
근육: 은서 모델 (εL 0.33, 잠열, 히스테리시스). M3(450 N, 50 mm) 와 M1(설계 힘, 86 mm) 비교. ⚠ M1 의 열 값은 M3 그대로 [가정]

실행 (이 폴더에서, conda activate osim — elbow_finalA.py build 를 먼저 해 둘 것)
  python lift_finalA.py run              → 표 + 그림 out\\fig_lift_finalA.png
  python lift_finalA.py view M1          → MuJoCo 창 (M3 / M1)
필요: 이 폴더에 finalA_elbow_r.xml, sma_model.py, 유빈 osim, 유빈 .mot
     (.mot 이 없으면 C:\\work\\addb_opensim\\report\\moco 에서 찾음)
"""
import sys, os, math, copy, json, time, logging
import xml.etree.ElementTree as ET
import numpy as np, mujoco
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'sma_arm_kimm'))   # 김은서 sma_model.py v3 (수정 없이 사용)
import sma_model as sm
logging.getLogger('matplotlib').setLevel(logging.ERROR)

NAME = 'gui_upper_robot40h160_finalA_lift_s6_PFWJ_2'
XML = os.path.join(HERE, 'finalA_elbow_r.xml')
EPSL = 0.40
UMAX = 1.0                                   # 은서 v3: power_duty(1) = 267 W (36 A)
KP, KD = float(os.environ.get('KP', '0.03')), float(os.environ.get('KD', '0.04'))                  # 명령 = KP × 각도 오차(°) + KD × 속도 오차(°/s), −1 ~ 1
T_STANDBY, T_SPAN, KT = 52.0, float(os.environ.get('TSPAN', 40.0)), 0.1   # 목표 온도 = 52 + 40 × 명령, 가열 = 0.1 × (목표 − 지금 온도)
DAMP = 5.0                           # 관절 감쇠 [가정]
H_FORCED = sm.KIMM_HEAT["h_forced"]; H_HALF = H_FORCED / 2; H_ON = sm.KIMM_HEAT["h_on"]   # 은서 v3 값
DESIGN_F = {'elbow_M0_r': 2327.05, 'elbow_M1_r': 1834.86}


def find(fname):
    for c in [os.path.join(HERE, fname), os.path.join(HERE, '..', 'report', 'moco', fname), os.path.join('report', 'moco', fname), os.path.join(r'C:\work\addb_opensim\report\moco', fname)]:
        if os.path.exists(c): return c
    sys.exit(f'{fname} 을 못 찾았어요 → 이 폴더로 복사해 주세요:\n  copy C:\\work\\addb_opensim\\report\\moco\\{fname} .')


# ── 위팔 방향 (땅 → 골반 → 몸통 → 위팔, 회전만) ──
def rx(a): c, s = math.cos(a), math.sin(a); return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
def ry(a): c, s = math.cos(a), math.sin(a); return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
def rz(a): c, s = math.cos(a), math.sin(a); return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
def euler_xyz(v): return rx(v[0]) @ ry(v[1]) @ rz(v[2])


def axis_rot(axis, q):
    a = np.asarray(axis, float); a = a / np.linalg.norm(a); K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + math.sin(q) * K + (1 - math.cos(q)) * K @ K


def chain(osim_path, names=('ground_pelvis', 'back', 'acromial_r')):
    model = ET.parse(osim_path).getroot().find('Model'); out = []
    for j in model.find('JointSet/objects'):
        if j.get('name') not in names: continue
        fr = {f.get('name'): f for f in j.find('frames')}
        pf, cf = fr[j.findtext('socket_parent_frame').strip()], fr[j.findtext('socket_child_frame').strip()]
        axes = []
        for ta in j.iter('TransformAxis'):
            if not ta.get('name').startswith('rotation'): continue
            co = (ta.findtext('coordinates') or '').split(); lf = ta.find('LinearFunction')
            ab = [float(x) for x in lf.findtext('coefficients').split()] if lf is not None else [1.0, 0.0]
            axes.append((np.array([float(x) for x in ta.findtext('axis').split()]), co[0] if co else None, ab))
        out.append((names.index(j.get('name')), euler_xyz([float(x) for x in pf.findtext('orientation').split()]),
                    euler_xyz([float(x) for x in cf.findtext('orientation').split()]), axes))
    return [c[1:] for c in sorted(out, key=lambda c: c[0])]


def humerus_R(CH, q):
    R = np.eye(3)
    for Rp, Rc, axes in CH:
        Rj = np.eye(3)
        for ax, co, (a, b) in axes: Rj = Rj @ axis_rot(ax, a * q.get(co, 0.0) + b)
        R = R @ Rp @ Rj @ Rc.T
    return R


def read_mot(path):
    L = open(path, encoding='utf-8').read().splitlines(); i = [k for k, l in enumerate(L) if l.strip() == 'endheader'][0]
    hdr = L[i + 1].split('\t'); D = np.array([[float(x) for x in l.split('\t')] for l in L[i + 2:] if l.strip()])
    deg = any('inDegrees=yes' in l for l in L[:i])
    return hdr, D, deg


# ── 근육 ──
def params(module, cable):
    p = copy.deepcopy(sm.M3); p['epsL'] = EPSL; p['force_model'] = 't31'
    if module == 'M1':
        p['FMAX'] = DESIGN_F[cable]; p['L0'] = 0.26; p['A_wire'] = p['FMAX'] / 400e6
    p['K_mod'] = p['FMAX'] / EPSL
    return p


def prepare():
    if not os.path.exists(XML): sys.exit('finalA_elbow_r.xml 이 없어요 → 먼저  python elbow_finalA.py build')
    hdr, D, deg = read_mot(find(NAME + '.mot')); CH = chain(find(NAME + '_model.osim'))
    t = D[:, 0] - D[0, 0]; col = lambda c: D[:, hdr.index(c)] * (math.pi / 180 if deg and not c.startswith('/') else 1.0)
    qref = col('elbow_flex_r')
    rot = ['pelvis_tilt', 'pelvis_list', 'pelvis_rotation', 'lumbar_extension', 'lumbar_bending', 'lumbar_rotation', 'arm_flex_r', 'arm_add_r', 'arm_rot_r']
    Rx90 = rx(math.pi / 2)                                                   # MuJoCo 화면용으로 세운 만큼
    grav = np.array([Rx90 @ (humerus_R(CH, {c: col(c)[k] for c in rot}).T @ np.array([0, -9.81, 0])) for k in range(len(t))])
    Freq = {n: col(f'/forceset/{n}/activation') * DESIGN_F[n] for n in DESIGN_F}
    return dict(t=t, qref=qref, grav=grav, Freq=Freq)


def setup(module, REF):
    m = mujoco.MjModel.from_xml_string(open(XML, encoding='utf-8').read()); m.opt.timestep = 0.0005
    m.dof_damping[:] = DAMP
    d = mujoco.MjData(m); qa = m.jnt_qposadr[m.joint('elbow_flex_r').id]; va = m.jnt_dofadr[m.joint('elbow_flex_r').id]
    ib, it = m.tendon('elbow_M0_r').id, m.tendon('elbow_M1_r').id
    gb, gt = -m.actuator_gainprm[m.actuator('elbow_M0_r').id, 0], -m.actuator_gainprm[m.actuator('elbow_M1_r').id, 0]
    d.qpos[qa] = 0.0; mujoco.mj_forward(m, d); Lb_ref = d.ten_length[ib]
    d.qpos[qa] = math.radians(150); mujoco.mj_forward(m, d); Lt_ref = d.ten_length[it]
    d.qpos[qa] = REF['qref'][0]; mujoco.mj_forward(m, d)
    pb, pt = params(module, 'elbow_M0_r'), params(module, 'elbow_M1_r')
    return dict(m=m, d=d, qa=qa, va=va, ib=ib, it=it, gb=gb, gt=gt, Lb_ref=Lb_ref, Lt_ref=Lt_ref, pb=pb, pt=pt,
                b=sm.SMA(pb, latent=True), t=sm.SMA(pt, latent=True))


def step(S, REF):
    m, d = S['m'], S['d']; tt = min(d.time, REF['t'][-1])
    qr = np.interp(tt, REF['t'], REF['qref']); dqr = np.interp(tt, REF['t'], np.gradient(REF['qref'], REF['t']))
    m.opt.gravity[:] = [np.interp(tt, REF['t'], REF['grav'][:, k]) for k in range(3)]
    e = math.degrees(qr - d.qpos[S['qa']]); de = math.degrees(dqr - d.qvel[S['va']])
    cmd = float(np.clip(KP * e + KD * de, -1, 1))
    b, t = S['b'], S['t']

    def temp_ctl(sma, T_des):                      # 안쪽 온도 제어
        return float(np.clip(KT * (T_des - sma.T), 0, UMAX)), (H_FORCED if sma.T > T_des + 1 else H_ON)
    ub, hb = temp_ctl(b, T_STANDBY + T_SPAN * max(cmd, 0.0))
    ut, ht = temp_ctl(t, T_STANDBY + T_SPAN * max(-cmd, 0.0))
    eb = max(EPSL - (S['Lb_ref'] - d.ten_length[S['ib']]) / S['pb']['L0'], 0.0)
    et = max(EPSL - (S['Lt_ref'] - d.ten_length[S['it']]) / S['pt']['L0'], 0.0)
    fb = b.step(b.power_duty(min(ub, b.u_cap())), m.opt.timestep, eps=eb, hS=hb)
    ft = t.step(t.power_duty(min(ut, t.u_cap())), m.opt.timestep, eps=et, hS=ht)
    d.ctrl[0] = min(fb / S['gb'], 1.0); d.ctrl[1] = min(ft / S['gt'], 1.0)
    mujoco.mj_step(m, d)
    return qr, fb, ft, ub, ut


def simulate(module, REF):
    S = setup(module, REF); R = []
    while S['d'].time < REF['t'][-1]:
        qr, fb, ft, ub, ut = step(S, REF)
        R.append((S['d'].time, math.degrees(qr), math.degrees(S['d'].qpos[S['qa']]), S['b'].T, S['t'].T, fb, ft, ub, ut))
    return np.array(R)


mode = sys.argv[1] if len(sys.argv) > 1 else 'run'
if mode == 'run':
    REF = prepare(); res = {}
    print(f"목표: 유빈 물병 들기 {REF['t'][-1]:.1f} 초, 팔꿈치 {math.degrees(REF['qref'].min()):.0f}° ~ {math.degrees(REF['qref'].max()):.0f}°")
    print(f"필요한 힘 (Moco): 굽힘 최대 {REF['Freq']['elbow_M0_r'].max():.0f} N, 폄 최대 {REF['Freq']['elbow_M1_r'].max():.0f} N\n")
    print(f"{'근육':4s} {'평균 오차':>8s} {'최대 오차':>8s} {'10° 안 시간':>10s} {'최고 온도':>8s} {'최대 힘 굽힘/폄':>16s} {'가열 에너지':>9s}")
    for module in ['M3', 'M1']:
        R = simulate(module, REF); res[module] = R
        err = R[:, 2] - R[:, 1]; dt = np.diff(R[:, 0], prepend=0)
        E = float(np.sum((R[:, 7] + R[:, 8]) * sm.M3['Qon'] * dt))
        print(f"{module:4s} {np.sqrt(np.mean(err ** 2)):7.1f}° {np.abs(err).max():7.1f}° {np.mean(np.abs(err) < 10) * 100:9.0f} % "
              f"{R[:, 3:5].max():7.1f}°C {R[:, 5].max():7.0f} / {R[:, 6].max():5.0f} N {E:7.0f} J")
    print('\n평균 오차 = RMS (°). 10° 안 시간 = 목표와 10° 이내였던 시간 비율. 가열 에너지 = 6.4 초 동안 넣은 전기 열 (근육 2개 합)')
    os.makedirs(os.path.join(HERE, 'out'), exist_ok=True)
    json.dump({k: {'rms_deg': float(np.sqrt(np.mean((v[:, 2] - v[:, 1]) ** 2))), 'max_err_deg': float(np.abs(v[:, 2] - v[:, 1]).max()),
                   'Tmax': float(v[:, 3:5].max()), 'Fmax_flex': float(v[:, 5].max()), 'Fmax_ext': float(v[:, 6].max())} for k, v in res.items()},
              open(os.path.join(HERE, 'out', 'lift_finalA.json'), 'w', encoding='utf-8'), indent=1)
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['Malgun Gothic', 'NanumGothic', 'DejaVu Sans']; plt.rcParams['axes.unicode_minus'] = False
    fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
    R3, R1 = res['M3'], res['M1']
    ax[0].plot(R1[:, 0], R1[:, 1], 'k-', lw=3, label='목표 (유빈 들기 동작)')
    ax[0].plot(R3[:, 0], R3[:, 2], '-', c='#9AA5B5', lw=2, label='M3 (450 N, 50 mm)'); ax[0].plot(R1[:, 0], R1[:, 2], '-', c='#B8501E', lw=2, label='M1 크기 (설계 힘, 86 mm)')
    ax[0].set(ylabel='팔꿈치 각도 (°)', title='물병 3.79 kg 들기 — 팔꿈치 따라가기'); ax[0].legend(fontsize=9)
    ax[1].plot(REF['t'], REF['Freq']['elbow_M0_r'], 'k--', lw=1.5, label='필요한 굽힘 힘 (Moco)')
    ax[1].plot(R1[:, 0], R1[:, 5], c='#C0392B', lw=1.8, label='M1 굽힘 힘'); ax[1].plot(R1[:, 0], R1[:, 6], c='#2B6CB0', lw=1.8, label='M1 폄 힘')
    ax[1].plot(R3[:, 0], R3[:, 5], c='#C0392B', lw=1, alpha=.5, label='M3 굽힘 힘'); ax[1].set(ylabel='케이블 힘 (N)'); ax[1].legend(fontsize=8, ncol=2)
    ax[2].plot(R1[:, 0], R1[:, 3], c='#C0392B', lw=1.8, label='M1 굽힘 근육'); ax[2].plot(R1[:, 0], R1[:, 4], c='#2B6CB0', lw=1.8, label='M1 폄 근육')
    ax[2].axhline(72, c='#888', ls='--', lw=1); ax[2].axhline(52, c='#888', ls=':', lw=1)
    ax[2].set(ylabel='온도 (°C)', xlabel='시간 (초)'); ax[2].legend(fontsize=8)
    for a in ax: a.grid(alpha=.3)
    plt.tight_layout(); out = os.path.join(HERE, 'out', 'fig_lift_finalA.png'); plt.savefig(out, dpi=140); print('그림:', out)

elif mode == 'view':
    import mujoco.viewer
    module = sys.argv[2] if len(sys.argv) > 2 else 'M1'
    REF = prepare(); S = setup(module, REF); m, d = S['m'], S['d']; nxt = 0.0
    print(f'{module} · 물병 들기 따라가기 — 끝나면 처음부터 다시. 창을 닫으면 끝')
    with mujoco.viewer.launch_passive(m, d) as v:
        while v.is_running():
            t0 = time.time()
            for _ in range(int(0.01 / m.opt.timestep)): qr, fb, ft, ub, ut = step(S, REF)
            if d.time >= nxt:
                print(f"t={d.time:4.1f}s 목표 {math.degrees(qr):6.1f}° 지금 {math.degrees(d.qpos[S['qa']]):6.1f}°  굽힘 {S['b'].T:5.1f}°C {fb:6.0f} N  폄 {S['t'].T:5.1f}°C {ft:6.0f} N"); nxt += 0.5
            if d.time >= REF['t'][-1]:                                   # 끝나면 처음부터 다시
                mujoco.mj_resetData(m, d); d.qpos[S['qa']] = REF['qref'][0]; mujoco.mj_forward(m, d)
                S['b'], S['t'] = sm.SMA(S['pb'], latent=True), sm.SMA(S['pt'], latent=True); nxt = 0.0; print('— 처음부터 다시 —')
            v.sync(); time.sleep(max(0.0, 0.01 - (time.time() - t0)))
else:
    print(__doc__)
