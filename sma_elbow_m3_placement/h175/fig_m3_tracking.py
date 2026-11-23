import os, sys, json, numpy as np
os.environ.setdefault('NM3', '5'); sys.argv = ['x']
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import opt_m3 as o
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from matplotlib import font_manager
for f in font_manager.findSystemFonts():
    if 'NotoSansCJK' in f or 'Nanum' in f: font_manager.fontManager.addfont(f)
plt.rcParams['font.family'] = [n.name for n in font_manager.fontManager.ttflist if 'CJK' in n.name or 'Nanum' in n.name][:1] or ['DejaVu Sans']
o._init(); g = o._W['g']; ns = o._W['ns']; G, yub = g['task']('path')
st = json.load(open('out/static_m3.json'))['opt']; top = json.load(open(f'out/optm3_{o.NM3}_top.json'))[0]
C = [('유빈 175 배치', yub['elbow_M0_r'], yub['elbow_M1_r'], '#9AA5B5'),
     ('정적 최적화', np.array(st['elbow_M0_r']['x']), np.array(st['elbow_M1_r']['x']), '#D69E2E'),
     ('MuJoCo 최적화', np.array(top['x']['elbow_M0_r']), np.array(top['x']['elbow_M1_r']), '#2E7D4F')]
S = {'elbow_M0_r': o.S_M3, 'elbow_M1_r': o.S_M3}
fig, ax = plt.subplots(figsize=(8, 3.6)); first = True
for name, xf, xe, c in C:
    o.set_muscle(g['xml_from']({'elbow_M0_r': xf, 'elbow_M1_r': xe}, o._W['tmp']), S); R = ns['simulate']('custom', o._W['REF'])
    if first: ax.plot(R[:, 0], R[:, 1], 'k--', lw=1.5, label='목표 (유빈 물병 들기)'); first = False
    e = np.sqrt(np.mean((R[:, 2] - R[:, 1]) ** 2)); ax.plot(R[:, 0], R[:, 2], color=c, lw=2, label=f'{name} (RMS {e:.1f}°)')
ax.set_xlabel('시간 (s)'); ax.set_ylabel('팔꿈치 각도 (°)'); ax.set_title(f'175 cm 팔 · 기계연 M3×{o.NM3} — 물병 3.79 kg 들기 추종')
ax.legend(fontsize=8, frameon=False, loc='upper right'); ax.grid(alpha=.3); fig.tight_layout(); fig.savefig(f'out/fig175_M3x{o.NM3}.png', dpi=160)
print('ok')
