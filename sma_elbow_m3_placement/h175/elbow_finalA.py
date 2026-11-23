# -*- coding: utf-8 -*-
"""②-1 유빈 최종 A안 팔꿈치 → MuJoCo 팔 (오른팔, 팔꿈치 케이블 2개)

하는 일
  build : 유빈 들기 모델(.osim)에서 위팔·아래팔(+요골·손·물병) 뼈대와 elbow_M0_r(굽힘)·elbow_M1_r(폄) 부착점을 읽어
          MuJoCo 팔(finalA_elbow_r.xml)을 만든다. 같은 osim 을 파이썬으로 직접 계산한 값과 케이블 길이·모멘트암을 비교한다.
  check : OpenSim(osim 환경)으로 같은 각도에서 케이블 길이·모멘트암을 재서 MuJoCo 와 비교한다 → 표 + 그림. 통과 기준
          길이 0.5 % 이내, 모멘트암 3 mm 이내 (유빈 README "변환 검증" 기준).
  view  : MuJoCo 창에서 팔꿈치를 0° → 150° → 0° 로 천천히 움직이며 케이블 길이를 보여 준다 (근육 힘 없음, 눈으로 확인용).

실행 (C:\\work\\addb_opensim 에서, conda activate osim)
  python elbow_finalA\\elbow_finalA.py build
  python elbow_finalA\\elbow_finalA.py check
  python elbow_finalA\\elbow_finalA.py view
  osim 파일이 다른 곳에 있으면 마지막에 경로를 붙인다:  ... build D:\\경로\\파일.osim

규약: 각도는 OpenSim 그대로 (elbow_flex_r 0 = 쭉 폄, + = 굽힘, 0~150°). 위팔·아래팔 안의 좌표(부착점)도 OpenSim 그대로이고,
     화면에서 똑바로 서 보이게 팔 전체를 x축으로 90° 돌려 놓았다 (y 위 → z 위). 길이·모멘트암은 돌림과 상관없다.
     어깨는 0° 로 고정(팔을 아래로 내린 자세), 아래팔 회전(pro_sup)·손목은 0 으로 고정. 손 질량에 물병 3.79 kg 이 들어 있다.
"""
import sys, os, json, math
import xml.etree.ElementTree as ET
import numpy as np
import logging; logging.getLogger('matplotlib').setLevel(logging.ERROR)

HERE = os.path.dirname(os.path.abspath(__file__))
OSIM_DEFAULT = os.path.join('report', 'moco', 'gui_upper_robot40h160_finalA_lift_s6_PFWJ_2_model.osim')
XML_OUT = os.path.join(HERE, 'finalA_elbow_r.xml')
CABLES = ['elbow_M0_r', 'elbow_M1_r']                  # M0 = 굽힘(이두 역할), M1 = 폄(삼두 역할)
NAMES = {'elbow_M0_r': 'M0 굽힘', 'elbow_M1_r': 'M1 폄'}
Q_DEG = np.arange(0, 151, 10)                          # 비교할 팔꿈치 각도
TOL_LEN, TOL_MA = 0.5, 3.0                             # 통과 기준: 길이 %, 모멘트암 mm


# ───────────────────────── osim 읽기 (OpenSim 없이) ─────────────────────────
def vec(t): return np.array([float(x) for x in t.split()])


def rx(a): c, s = math.cos(a), math.sin(a); return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])
def ry(a): c, s = math.cos(a), math.sin(a); return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])
def rz(a): c, s = math.cos(a), math.sin(a); return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
def euler_xyz(v): return rx(v[0]) @ ry(v[1]) @ rz(v[2])     # OpenSim: 몸체 고정 x-y-z 순서


def axis_rot(axis, q):
    a = np.asarray(axis, float); a = a / np.linalg.norm(a); K = np.array([[0, -a[2], a[1]], [a[2], 0, -a[0]], [-a[1], a[0], 0]])
    return np.eye(3) + math.sin(q) * K + (1 - math.cos(q)) * K @ K


def quat(R):
    w = math.sqrt(max(0.0, 1 + R[0, 0] + R[1, 1] + R[2, 2])) / 2
    x = math.copysign(math.sqrt(max(0.0, 1 + R[0, 0] - R[1, 1] - R[2, 2])) / 2, R[2, 1] - R[1, 2])
    y = math.copysign(math.sqrt(max(0.0, 1 - R[0, 0] + R[1, 1] - R[2, 2])) / 2, R[0, 2] - R[2, 0])
    z = math.copysign(math.sqrt(max(0.0, 1 - R[0, 0] - R[1, 1] + R[2, 2])) / 2, R[1, 0] - R[0, 1])
    return np.array([w, x, y, z])


def read_osim(path):
    model = ET.parse(path).getroot().find('Model')
    bodies = {}
    for b in model.find('BodySet/objects'):
        bodies[b.get('name')] = dict(mass=float(b.findtext('mass')), com=vec(b.findtext('mass_center')), inertia=vec(b.findtext('inertia')))
    joints = {}
    for j in model.find('JointSet/objects'):
        fr = {f.get('name'): f for f in j.find('frames')}
        pf, cf = fr[j.findtext('socket_parent_frame').strip()], fr[j.findtext('socket_child_frame').strip()]
        axes = []
        if j.tag == 'PinJoint': axes = [np.array([0, 0, 1.0])]
        elif j.tag == 'UniversalJoint': axes = [np.array([1.0, 0, 0]), np.array([0, 1.0, 0])]
        joints[j.get('name')] = dict(kind=j.tag, parent=pf.findtext('socket_parent').strip().split('/')[-1],
                                     child=cf.findtext('socket_parent').strip().split('/')[-1],
                                     t_p=vec(pf.findtext('translation')), R_p=euler_xyz(vec(pf.findtext('orientation'))),
                                     t_c=vec(cf.findtext('translation')), R_c=euler_xyz(vec(cf.findtext('orientation'))), axes=axes,
                                     coords=[c.get('name') for c in j.findall('coordinates/Coordinate')],
                                     ranges=[vec(c.findtext('range')) for c in j.findall('coordinates/Coordinate')])
    cables = {}
    for tag in ['Millard2012EquilibriumMuscle', 'Thelen2003Muscle', 'DeGrooteFregly2016Muscle', 'PathActuator']:
        for m in model.iter(tag):
            if m.get('name') in CABLES:
                pts = [(p.findtext('socket_parent_frame').strip().split('/')[-1], vec(p.findtext('location'))) for p in m.iter('PathPoint')]
                wraps = [w.get('name') for w in m.iter('PathWrap')]
                F = float(m.findtext('max_isometric_force') or m.findtext('optimal_force'))
                cables[m.get('name')] = dict(points=pts, wraps=wraps, F=F)
    return bodies, joints, cables


def x_child(J, qs):
    """부모 몸체 좌표 → 자식 몸체 좌표 변환 (R, t): X_p · 관절회전(q) · X_c⁻¹"""
    Rj = np.eye(3)
    for ax, q in zip(J['axes'], qs): Rj = Rj @ axis_rot(ax, q)
    R = J['R_p'] @ Rj @ J['R_c'].T
    return R, J['t_p'] - R @ J['t_c']


def py_length(joints, cables, name, q):
    """파이썬 직접 계산 (위팔 좌표계): 직선 2점 케이블 길이"""
    R, t = x_child(joints['elbow_r'], [q])
    (b0, p0), (b1, p1) = cables[name]['points']
    P = {'humerus_r': lambda p: p, 'ulna_r': lambda p: t + R @ p}
    return float(np.linalg.norm(P[b1](p1) - P[b0](p0)))


def find_osim():
    for a in sys.argv[2:]:
        if a.lower().endswith('.osim'): return a
    for c in [OSIM_DEFAULT, os.path.join(HERE, '..', OSIM_DEFAULT), os.path.join(HERE, os.path.basename(OSIM_DEFAULT))]:
        if os.path.exists(c): return c
    sys.exit(f'osim 파일을 못 찾았어요: {OSIM_DEFAULT}\n  C:\\work\\addb_opensim 에서 실행하거나, 명령 끝에 osim 경로를 붙여 주세요.')


# ───────────────────────── MuJoCo xml 만들기 ─────────────────────────
def f(v): return ' '.join(f'{x:.9g}' for x in np.ravel(v))


def inertial(b):
    I = b['inertia']
    if np.all(np.abs(I[:3]) < 1e-9): return f'<inertial pos="{f(b["com"])}" mass="{b["mass"]:.6g}" diaginertia="1e-6 1e-6 1e-6"/>'
    return f'<inertial pos="{f(b["com"])}" mass="{b["mass"]:.6g}" fullinertia="{f(I)}"/>'


def build_xml(bodies, joints, cables):
    Je, Jr, Jh = joints['elbow_r'], joints['radioulnar_r'], joints['radius_hand_r']
    R0, t0 = x_child(Je, [0.0])                                  # 위팔 → 아래팔(척골), 팔꿈치 0°
    Rr, tr = x_child(Jr, [0.0])                                  # 척골 → 요골 (회내·회외 0 고정)
    Rh, th = x_child(Jh, [0.0, 0.0])                             # 요골 → 손 (손목 0 고정)
    ax_c = Je['R_c'] @ np.array([0, 0, 1.0])                     # 팔꿈치 축 (척골 좌표)
    lo, hi = Je['ranges'][0]
    hand_c = bodies['hand_r']['com']
    sites_h = ''.join(f'<site name="{n}_p0" pos="{f(c["points"][0][1])}" size="0.006" rgba="1 0.2 0.2 1"/>' for n, c in cables.items())
    sites_u = ''.join(f'<site name="{n}_p1" pos="{f(c["points"][1][1])}" size="0.006" rgba="1 0.2 0.2 1"/>' for n, c in cables.items())
    tendons = ''.join(f'<spatial name="{n}" width="0.004" rgba="{"0.9 0.1 0.1" if "M0" in n else "0.1 0.3 0.9"} 1"><site site="{n}_p0"/><site site="{n}_p1"/></spatial>' for n in cables)
    acts = ''.join(f'<general name="{n}" tendon="{n}" ctrlrange="0 1" ctrllimited="true" gainprm="-{c["F"]:.6g}" biastype="none"/>' for n, c in cables.items())
    elbow_world = Je['t_p']
    return f"""<mujoco model="finalA_elbow_r">
  <!-- 자동 생성: elbow_finalA.py build — 유빈 최종 A안 들기 모델의 오른팔 (40 kg · 160 cm). 각도는 OpenSim 그대로, 위팔 안 좌표도 OpenSim 그대로 (화면용으로 통째로 z 위로 세움) -->
  <compiler angle="radian"/>
  <option timestep="0.0005" gravity="0 0 -9.81" integrator="implicitfast"/>
  <visual><global offwidth="1280" offheight="720"/></visual>
  <worldbody>
    <light pos="0.5 1 1" dir="-0.5 -1 -1"/>
    <geom type="plane" pos="0 0 -0.75" size="0.8 0.8 0.05" rgba="0.8 0.8 0.8 1"/>
    <light pos="0 -1 1.5" dir="0 0.6 -1"/>
    <body name="humerus_r" pos="0 0 0.3" quat="0.70710678 0.70710678 0 0">  <!-- OpenSim(y 위) → MuJoCo(z 위): x축으로 90° 돌림 -->
      {inertial(bodies['humerus_r'])}
      <geom type="capsule" fromto="0 0 0 {f(elbow_world)}" size="0.022" rgba="0.7 0.7 0.72 1" contype="0" conaffinity="0"/>
      {sites_h}
      <body name="ulna_r" pos="{f(t0)}" quat="{f(quat(R0))}">
        <joint name="elbow_flex_r" type="hinge" pos="{f(Je['t_c'])}" axis="{f(ax_c)}" range="{lo:.6g} {hi:.6g}" limited="true" damping="0.3"/>
        {inertial(bodies['ulna_r'])}
        <geom type="capsule" fromto="{f(Je['t_c'])} 0 -0.24 0" size="0.018" rgba="0.7 0.7 0.72 1" contype="0" conaffinity="0"/>
        {sites_u}
        <body name="radius_r" pos="{f(tr)}" quat="{f(quat(Rr))}">
          {inertial(bodies['radius_r'])}
          <body name="hand_r" pos="{f(th)}" quat="{f(quat(Rh))}">
            {inertial(bodies['hand_r'])}
            <geom type="sphere" pos="{f(hand_c)}" size="0.045" rgba="0.55 0.1 0.15 1" contype="0" conaffinity="0"/>
          </body>
        </body>
      </body>
    </body>
  </worldbody>
  <tendon>{tendons}</tendon>
  <actuator>{acts}</actuator>
</mujoco>
"""


def load_mj(path):
    """MuJoCo 는 윈도우 한글 경로의 파일을 못 열어서, 파이썬으로 읽어 문자열로 넘긴다 (메시 파일이 없어서 가능)"""
    import mujoco
    return mujoco.MjModel.from_xml_string(open(path, encoding='utf-8').read())


def mj_table(xml_path, names, q_deg):
    import mujoco
    m = load_mj(xml_path); d = mujoco.MjData(m)
    jid = m.joint('elbow_flex_r').id; dof = m.jnt_dofadr[jid]; qa = m.jnt_qposadr[jid]
    L = {n: [] for n in names}; MA = {n: [] for n in names}

    def lengths(q):
        d.qpos[:] = 0; d.qpos[qa] = q; mujoco.mj_forward(m, d)
        return {n: float(d.ten_length[m.tendon(n).id]) for n in names}
    h = 1e-6
    for qd in q_deg:
        q = math.radians(qd); l0, lp, lm = lengths(q), lengths(q + h), lengths(q - h)
        for n in names:
            L[n].append(l0[n]); MA[n].append(-(lp[n] - lm[n]) / (2 * h))   # 모멘트암 = −dL/dq (OpenSim 과 같은 부호)
    return {n: np.array(v) for n, v in L.items()}, {n: np.array(v) for n, v in MA.items()}


def compare(title, Lref, MAref, Lmj, MAmj, names):
    ok = True
    print(f'\n[{title}]  통과 기준: 길이 {TOL_LEN} % · 모멘트암 {TOL_MA} mm 이내')
    print(f"{'케이블':10s} {'각도':>5s} {'길이 기준':>9s} {'MuJoCo':>8s} {'차이%':>7s} | {'모멘트암 기준':>11s} {'MuJoCo':>8s} {'차이mm':>7s}")
    for n in names:
        for k, qd in enumerate(Q_DEG):
            dl = (Lmj[n][k] - Lref[n][k]) / Lref[n][k] * 100; dm = (MAmj[n][k] - MAref[n][k]) * 1000
            if k % 3 == 0 or k == len(Q_DEG) - 1:
                print(f"{NAMES[n]:10s} {qd:5.0f} {Lref[n][k] * 1000:8.1f}mm {Lmj[n][k] * 1000:7.1f}mm {dl:+7.3f} | {MAref[n][k] * 1000:9.1f}mm {MAmj[n][k] * 1000:7.1f}mm {dm:+7.2f}")
        el = np.max(np.abs((Lmj[n] - Lref[n]) / Lref[n] * 100)); em = np.max(np.abs(MAmj[n] - MAref[n]) * 1000)
        good = el <= TOL_LEN and em <= TOL_MA; ok &= good
        print(f"  → {NAMES[n]}: 길이 최대 차이 {el:.3f} %, 모멘트암 최대 차이 {em:.2f} mm  {'통과' if good else '미달'}")
        print(f"     길이 변화(0→150°) {abs(Lref[n][-1] - Lref[n][0]) * 1000:.1f} mm · 모멘트암 {MAref[n].min() * 1000:.1f} ~ {MAref[n].max() * 1000:.1f} mm")
    print('\n결과:', '모두 통과 — MuJoCo 팔이 유빈 설계와 같다' if ok else '미달 있음 — 화면을 캡처해서 보내 주세요')
    return ok


def plot(Lref, MAref, Lmj, MAmj, names, ref_label, out):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    plt.rcParams['font.family'] = ['Malgun Gothic', 'DejaVu Sans']; plt.rcParams['axes.unicode_minus'] = False
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for n, c in zip(names, ['#C0392B', '#2B6CB0']):
        ax[0].plot(Q_DEG, Lref[n] * 1000, '-', c=c, lw=2.5, label=f'{NAMES[n]} · {ref_label}'); ax[0].plot(Q_DEG, Lmj[n] * 1000, 'o', c=c, mfc='white', label=f'{NAMES[n]} · MuJoCo')
        ax[1].plot(Q_DEG, MAref[n] * 1000, '-', c=c, lw=2.5, label=f'{NAMES[n]} · {ref_label}'); ax[1].plot(Q_DEG, MAmj[n] * 1000, 'o', c=c, mfc='white', label=f'{NAMES[n]} · MuJoCo')
    ax[0].set(xlabel='팔꿈치 각도 (°) · 0 = 쭉 폄', ylabel='케이블 길이 (mm)', title='케이블 길이')
    ax[1].set(xlabel='팔꿈치 각도 (°) · 0 = 쭉 폄', ylabel='모멘트암 (mm) · + = 굽히는 쪽', title='모멘트암 (축 ~ 케이블 거리)')
    ax[1].axhline(0, c='#999', lw=1)
    for a in ax: a.grid(alpha=.3); a.legend(fontsize=8)
    fig.suptitle(f'유빈 최종 A안 오른팔 팔꿈치 케이블 · {ref_label} vs MuJoCo', fontsize=13, fontweight='bold')
    plt.tight_layout(); plt.savefig(out, dpi=140); print('그림:', out)


# ───────────────────────── 실행 ─────────────────────────
mode = sys.argv[1] if len(sys.argv) > 1 else 'build'
osim_path = find_osim()
bodies, joints, cables = read_osim(osim_path)
names = [n for n in CABLES if n in cables]

if mode == 'build':
    print('osim:', osim_path)
    for n in names:
        c = cables[n]
        print(f"  {n}: {c['points'][0][0]} → {c['points'][1][0]}, 힘 {c['F']:.0f} N, 감기는 면 {c['wraps'] or '없음 (직선)'}")
    print(f"  손 질량 {bodies['hand_r']['mass']:.2f} kg (물병 포함), 아래팔 {bodies['ulna_r']['mass'] + bodies['radius_r']['mass']:.2f} kg")
    open(XML_OUT, 'w', encoding='utf-8').write(build_xml(bodies, joints, cables))
    print('MuJoCo 팔 저장:', XML_OUT)
    Lpy = {n: np.array([py_length(joints, cables, n, math.radians(q)) for q in Q_DEG]) for n in names}
    h = 1e-6
    MApy = {n: np.array([-(py_length(joints, cables, n, math.radians(q) + h) - py_length(joints, cables, n, math.radians(q) - h)) / (2 * h) for q in Q_DEG]) for n in names}
    Lmj, MAmj = mj_table(XML_OUT, names, Q_DEG)
    compare('osim 직접 계산 vs MuJoCo', Lpy, MApy, Lmj, MAmj, names)
    print('\n다음: python elbow_finalA\\elbow_finalA.py check   (OpenSim 으로 한 번 더 확인)')

elif mode == 'check':
    try:
        import opensim as osim
    except ImportError:
        sys.exit('opensim 이 없어요 → conda activate osim 한 뒤 다시 실행하세요.')
    if not os.path.exists(XML_OUT): sys.exit('먼저 build 를 실행하세요.')
    M = osim.Model(osim_path); s = M.initSystem()
    coord = M.getCoordinateSet().get('elbow_flex_r'); mus = M.getMuscles()
    Los = {n: [] for n in names}; MAos = {n: [] for n in names}
    for qd in Q_DEG:
        for cn in [c.getName() for c in M.getCoordinateSet()]:
            C = M.getCoordinateSet().get(cn)
            if not C.getLocked(s): C.setValue(s, 0.0, False)
        coord.setValue(s, math.radians(qd)); M.realizePosition(s)
        for n in names:
            mm = mus.get(n); Los[n].append(mm.getLength(s)); MAos[n].append(mm.computeMomentArm(s, coord))
    Los = {n: np.array(v) for n, v in Los.items()}; MAos = {n: np.array(v) for n, v in MAos.items()}
    Lmj, MAmj = mj_table(XML_OUT, names, Q_DEG)
    ok = compare('OpenSim vs MuJoCo', Los, MAos, Lmj, MAmj, names)
    os.makedirs(os.path.join(HERE, 'out'), exist_ok=True)
    plot(Los, MAos, Lmj, MAmj, names, 'OpenSim', os.path.join(HERE, 'out', 'fig_elbow_check.png'))
    json.dump({'q_deg': Q_DEG.tolist(), **{f'{n}_{k}': (v[n] * 1000).round(3).tolist() for k, v in [('L_osim_mm', Los), ('L_mj_mm', Lmj), ('MA_osim_mm', MAos), ('MA_mj_mm', MAmj)] for n in names}, 'pass': bool(ok)},
              open(os.path.join(HERE, 'out', 'elbow_check.json'), 'w', encoding='utf-8'), indent=1)

elif mode == 'view':
    import time, mujoco, mujoco.viewer
    if not os.path.exists(XML_OUT): sys.exit('먼저 build 를 실행하세요.')
    m = load_mj(XML_OUT); d = mujoco.MjData(m); qa = m.jnt_qposadr[m.joint('elbow_flex_r').id]
    print('팔꿈치 0° → 150° → 0° 반복 (근육 힘 없이 자세만). 창을 닫으면 끝'); t0 = time.time(); nxt = 0.0
    with mujoco.viewer.launch_passive(m, d) as v:
        while v.is_running():
            tt = time.time() - t0; qd = 75 * (1 - math.cos(2 * math.pi * tt / 8))     # 8초에 한 번
            d.qpos[qa] = math.radians(qd); mujoco.mj_forward(m, d)
            if tt >= nxt:
                print(f"각도 {qd:6.1f}°  " + '  '.join(f"{NAMES[n]} 길이 {d.ten_length[m.tendon(n).id] * 1000:6.1f} mm" for n in names)); nxt += 0.5
            v.sync(); time.sleep(0.01)
else:
    print(__doc__)
