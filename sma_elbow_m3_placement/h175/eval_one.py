import sys, os, json, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import opt_m3 as o, eval_m3 as E
if __name__ == '__main__':
    from multiprocessing import Pool
    f, N, idx = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]); d = json.load(open(f))['opt']
    jobs = [('x', d['elbow_M0_r']['x'], d['elbow_M1_r']['x'], 70000 + 100 * idx + j) for j in range(N)]
    with Pool(2, initializer=o._init) as p: v = np.array([r for _, r in p.map(E.job, jobs)])
    print(f"M3×{o.NM3} {os.path.basename(f)} 중앙값 {np.median(v):.1f}°  범위 {v.min():.1f}~{v.max():.1f}  실패(>25°) {(v > 25).sum()}/{N}", flush=True)
