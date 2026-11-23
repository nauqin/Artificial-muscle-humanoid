"""기계연 M3 (근육 0.18 m, 스트로크 0.40 × 0.18 = 72 mm) 조건 정적 최적화 — place_elbow 의 evaluate/optimize 그대로, S 만 72 mm"""
import sys, os, json, math, numpy as np
D = os.path.dirname(os.path.abspath(__file__)); os.chdir(D); sys.path.insert(0, D); sys.argv = ['x']
g = {'__file__': os.path.join(D, 'place_elbow.py'), '__name__': 'pe'}
exec(open('place_elbow.py', encoding='utf-8').read().split("if __name__ == '__main__':")[0], g)
S = 0.40 * 0.18
G, yubin = g['task']('path'); out = {'S': S, 'opt': {}, 'yubin': {}}
for n in g['CABLES']:
    o, F, use, rmin, Lmin = g['evaluate'](G, yubin[n][None], n, S, True)
    _, F2, *_ = g['evaluate'](G, yubin[n][None], n, S, False)
    out['yubin'][n] = dict(F_fit=float(F[0]), F_nofit=float(F2[0]), stroke=float(use[0]), r_min=float(rmin[0]), L_min=float(Lmin[0]))
    r = g['optimize'](G, n, S, True); out['opt'][n] = r
    print(n, '유빈:', {k: round(v, 4) for k, v in out['yubin'][n].items()}, '| 정적 최적:', {k: (round(v, 4) if not isinstance(v, list) else '') for k, v in r.items()}, 'M3 개수', math.ceil(r['F_need'] / 450), flush=True)
json.dump(out, open(os.environ.get('STATIC_OUT', 'out/static_m3.json'), 'w'), indent=1)
