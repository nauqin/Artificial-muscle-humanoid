# -*- coding: utf-8 -*-
"""MuJoCo 물병 들기 오차를 직접 줄이는 배치 최적화 (CMA-ES)
변수 12개: 굽힘 케이블 6 좌표 + 폄 케이블 6 좌표 (유빈 elbow_r 상자 안)
근육 길이: 케이블마다 "들어갈 수 있는 가장 긴 근육" = 케이블 최소 길이 (최대 261 mm, 스트로크 86 mm 와 같음) → 스트로크 = 0.33 × 근육 길이
조건(기하): 스트로크 · 모멘트암 ≥ 5 mm · 뼈 거리 ≥ 20 mm  → 어기면 다시 뽑기 (시뮬레이션 안 함)
목적: 부착점을 0.01 mm 흔들어 3번 돌린 평균 오차 (근육 힘 = M1 설계 힘, 유빈과 같음)
"""
import sys, os, json, time, math
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
_W = {}


def _init():
    sys.argv = ['x']
    g = {'__file__': os.path.join(HERE, 'place_elbow.py'), '__name__': 'pe'}
    exec(open(os.path.join(HERE, 'place_elbow.py'), encoding='utf-8').read().split("if __name__ == '__main__':")[0], g)
    G, yubin = g['task']('path'); ns = g['lift_ns'](); REF = ns['prepare']()
    _W.update(g=g, G=G, ns=ns, REF=REF, tmp=os.path.join(HERE, 'out', f'_opt_{os.getpid()}.xml'))


S_M3 = 0.40 * 0.18   # 기계연 M3: 근육 0.18 m × 0.40 = 72 mm (은서 v3)
NM3 = int(os.environ.get('NM3', '4'))   # 케이블마다 M3 묶음 수
JIT = float(os.environ.get('JIT', '1e-5'))


def decode(z):
    g = _W['g']; z = np.clip(np.asarray(z, float), 0, 1)
    lo, hi = g['LO'], g['HI']
    xf = lo + z[0:6] * (hi - lo); xe = lo + z[6:12] * (hi - lo)
    S = {}
    for n, x in (('elbow_M0_r', xf), ('elbow_M1_r', xe)):
        S[n] = S_M3
    return xf, xe, S


def geom_violation(xf, xe, S):
    g = _W['g']; v = 0.0
    for n, x in (('elbow_M0_r', xf), ('elbow_M1_r', xe)):
        obj, F, *_ = g['evaluate'](_W['G'], x[None], n, max(S[n], 1e-3), True)
        if not np.isfinite(F[0]): v += (obj[0] / 1e5 - 1) / 10
    return v


def set_muscle(xml, S):
    g, ns = _W['g'], _W['ns']; sm = ns['sm']; import copy

    def params(module, cable):
        p = copy.deepcopy(sm.M3); p['epsL'] = 0.40; p['force_model'] = 't31'
        p['FMAX'] = float(os.environ.get('FM3', 450.0)) * NM3; p['L0'] = 0.18; p['A_wire'] = p['FMAX'] / 400e6; p['K_mod'] = p['FMAX'] / 0.40
        return p
    ns['params'] = params; ns['XML'] = xml


def run_once(xf, xe, S, seed):
    g, ns = _W['g'], _W['ns']; rng = np.random.default_rng(seed)
    X = {'elbow_M0_r': xf + rng.normal(0, JIT, 6), 'elbow_M1_r': xe + rng.normal(0, JIT, 6)}
    set_muscle(g['xml_from'](X, _W['tmp']), S); R = ns['simulate']('custom', _W['REF'])
    e = R[:, 2] - R[:, 1]; return float(np.sqrt(np.mean(e ** 2)))


def objective(args):
    z, k = args
    if not _W: _init()
    xf, xe, S = decode(z)
    v = geom_violation(xf, xe, S)
    if v > 0: return 100.0 + 100.0 * v, None
    r = [run_once(xf, xe, S, 1000 * k + j) for j in range(int(os.environ.get('OPT_NJ', '3')))]
    return float(np.mean(np.minimum(r, 90.0))), r


def encode(xf, xe):
    g = _W['g']; lo, hi = g['LO'], g['HI']
    return np.concatenate([(xf - lo) / (hi - lo), (xe - lo) / (hi - lo)])


if __name__ == '__main__':
    import cma
    from multiprocessing import Pool
    gens = int(sys.argv[1]) if len(sys.argv) > 1 else 18
    sig0 = float(sys.argv[2]) if len(sys.argv) > 2 else 0.05
    _init()
    start = os.environ.get('OPT_START')
    if start:
        best0 = json.load(open(start, encoding='utf-8'))[0]
    else:
        st = json.load(open(os.path.join(HERE, 'out', 'static_m3.json'), encoding='utf-8'))['opt']; best0 = dict(x={n: st[n]['x'] for n in st})
    z0 = encode(np.array(best0['x']['elbow_M0_r']), np.array(best0['x']['elbow_M1_r']))
    es = cma.CMAEvolutionStrategy(z0, sig0, {'bounds': [0, 1], 'popsize': 10, 'seed': 3, 'verbose': -9})
    log = []; t0 = time.time(); k = 0
    with Pool(2, initializer=_init) as pool:
        for gen in range(gens):
            Z = es.ask(); tries = 0
            for i in range(len(Z)):                       # 기하 조건을 어기는 후보는 다시 뽑기 (계산이 빨라서 시뮬레이션 전에 거름)
                for _ in range(300):
                    if geom_violation(*decode(Z[i])) == 0: break
                    Z[i] = es.ask(1)[0]; tries += 1
            res = pool.map(objective, [(z, k + i + int(os.environ.get('OPT_SEED0', '0'))) for i, z in enumerate(Z)]); k += len(Z)
            f = [r[0] for r in res]; es.tell(Z, f)
            for z, (fv, rr) in zip(Z, res): log.append(dict(z=list(map(float, z)), f=fv, runs=rr))
            b = min(log, key=lambda d: d['f'])
            print(f'gen {gen + 1:2d}: 이번 세대 최고 {min(f):5.1f}° · 지금까지 최고 {b["f"]:5.1f}° · 다시 뽑기 {tries} · {time.time() - t0:4.0f} s', flush=True)
            json.dump(log, open(os.path.join(HERE, 'out', f'optm3_{NM3}_log.json'), 'w'), indent=0)
    top = sorted([d for d in log if d['f'] < 100], key=lambda d: d['f'])[:5]
    out = []
    for d in top:
        xf, xe, S = decode(d['z']); out.append(dict(f=d['f'], runs=d['runs'], S=S, x={'elbow_M0_r': xf.tolist(), 'elbow_M1_r': xe.tolist()}))
    json.dump(out, open(os.path.join(HERE, 'out', f'optm3_{NM3}_top.json'), 'w'), indent=1)
    print('top:', [round(d['f'], 1) for d in out], 'S(mm):', [{k[6:8]: round(v * 1000) for k, v in d['S'].items()} for d in out])
