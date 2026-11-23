"""배치들을 기계연 M3×n 으로 N 번 (부착점 JIT 흔들기) 돌려 비교"""
import sys, os, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import opt_m3 as o
def job(a):
    tag, xf, xe, seed = a
    if not o._W: o._init()
    S = {'elbow_M0_r': o.S_M3, 'elbow_M1_r': o.S_M3}
    return tag, o.run_once(np.array(xf), np.array(xe), S, seed)
if __name__ == '__main__':
    from multiprocessing import Pool
    N = int(sys.argv[1]); o._init(); g = o._W['g']; G, yubin = g['task']('path')
    st = json.load(open(os.path.join(HERE, 'out', 'static_m3.json')))['opt']
    cand = [('유빈 160 배치', yubin['elbow_M0_r'].tolist(), yubin['elbow_M1_r'].tolist()),
            ('정적 최적화 (M3)', st['elbow_M0_r']['x'], st['elbow_M1_r']['x'])]
    tp = os.path.join(HERE, 'out', f'optm3_{o.NM3}_top.json')
    if os.path.exists(tp):
        d = json.load(open(tp))[0]; cand.append(('MuJoCo 최적화 (M3)', d['x']['elbow_M0_r'], d['x']['elbow_M1_r']))
    yp = os.path.join(HERE, 'out', 'static_yubin.json')
    if os.path.exists(yp):
        sy = json.load(open(yp))['opt']; cand.append(('유빈식 정적 (β 0.25)', sy['elbow_M0_r']['x'], sy['elbow_M1_r']['x']))
    mp = os.path.join(HERE, 'out', 'static_mine2.json')
    if os.path.exists(mp):
        s2 = json.load(open(mp))['opt']; cand.append(('내 정적 2 (스트로크 끝까지)', s2['elbow_M0_r']['x'], s2['elbow_M1_r']['x']))
    jobs = [(c[0], c[1], c[2], 70000 + 100 * i + j) for i, c in enumerate(cand) for j in range(N)]
    with Pool(2, initializer=o._init) as p: res = p.map(job, jobs)
    out = {}
    for t, r in res: out.setdefault(t, []).append(r)
    summ = {}
    for c in cand:
        v = np.array(out[c[0]]); summ[c[0]] = dict(runs=v.tolist(), median=float(np.median(v)), fail=int((v > 25).sum()), x=[c[1], c[2]])
        print(f"M3×{o.NM3} (M3 {os.environ.get('FM3', '450')} N) 흔들기 {o.JIT*1000:g} mm | {c[0]:16s} 중앙값 {np.median(v):5.1f}°  범위 {v.min():5.1f}~{v.max():5.1f}°  실패(>25°) {int((v > 25).sum())}/{N}", flush=True)
    json.dump(summ, open(os.path.join(HERE, 'out', f"evalm3_M3x{o.NM3}_F{os.environ.get('FM3', '450')}_jit{o.JIT}.json"), 'w', encoding='utf-8'), ensure_ascii=False, indent=1)
