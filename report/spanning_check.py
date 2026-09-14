"""시너지 k개의 토크 방향벡터가 하지 토크 공간을 양의 결합으로 덮는가.

"시너지당 새 근육 하나씩 달아서 걷게 할 수 있나"의 첫 관문. 당김 전용 액추에이터
n+1 개가 R^n 을 양으로 스패닝해야 임의 방향 토크를 낼 수 있다. n+1 개일 때 판정은
rank = n 이고 영공간 벡터의 부호가 전부 같은지다.

    conda activate osim
    python report/spanning_check.py                 # 기본: 4 DOF, k=5, 6 subject 표본
    python report/spanning_check.py --k 6 --coords hip_flexion,hip_adduction,hip_rotation,knee_angle,ankle_angle
"""
import argparse, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import synergy as sy, addb_io as io, run_muscle_analysis as rma  # noqa: E402

DEFAULT_COORDS = ("hip_flexion", "hip_adduction", "knee_angle", "ankle_angle")
DEFAULT_TRIALS = (("subject3", "walking1_segment_0"), ("subject3", "walkingTS2_segment_0"),
                  ("subject4", "walking1_segment_0"), ("subject5", "walking1_segment_0"),
                  ("subject10", "walking1_segment_0"), ("subject11", "walking2_segment_0"))


def torque_vectors(sub, trial, side, k, coords, out_root):
    """시너지별 단위 활성당 평균 토크 벡터 (k × n_coords) [Nm]."""
    base = os.path.join(out_root, sub, trial)
    frc = io.read_mot(f"{base}/so/{trial}_StaticOptimization_force.sto")
    mus = sy.muscle_names(os.path.join(ROOT, "data", sub, f"{sub}_scaled.osim"))
    cs = [f"{c}_{side}" for c in coords]
    ma = rma.load_moment_arms(f"{base}/ma", trial, cs)
    A, names, _, _ = sy.load_activations(f"{base}/so/{trial}_StaticOptimization_activation.sto",
                                         side, muscles=mus, log=lambda *a: None)
    C, W = sy.nmf(A, k, seed=0); C, W, _ = sy.sort_synergies(C, W)
    share = C[:, :, None] * W[None, :, :]
    share /= np.maximum(share.sum(axis=1, keepdims=True), 1e-12)   # 근육 힘을 시너지에 귀속
    F = np.column_stack([frc.column(n) for n in names]); t = frc.time
    R = {c: np.column_stack([np.interp(t, ma[c].time, ma[c].column(n)) for n in names]) for c in cs}
    T = np.stack([np.stack([(share[:, i, :] * F * R[c]).sum(axis=1) for c in cs], axis=1) for i in range(k)])
    return np.stack([(T[i] * C[:, i:i + 1]).sum(axis=0) / max(C[:, i].sum(), 1e-9) for i in range(k)])


def positive_spanning(V):
    """(rank, 성립 여부, |영공간 λ|). k = n+1 일 때만 엄밀하다; k > n+1 이면 LP 가 필요하다."""
    M = V.T; n = M.shape[0]
    rank = np.linalg.matrix_rank(M)
    if rank < n or M.shape[1] != n + 1:
        return rank, None if M.shape[1] != n + 1 else False, None
    _, _, vt = np.linalg.svd(M); null = vt[-1] / np.abs(vt[-1]).max()
    ok = bool(np.all(null > 1e-6) or np.all(null < -1e-6))
    return rank, ok, np.abs(null)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=5)
    ap.add_argument("--coords", default=",".join(DEFAULT_COORDS))
    ap.add_argument("--out-root", default=os.path.join(ROOT, "out", "real"))
    args = ap.parse_args()
    coords = tuple(c for c in args.coords.split(",") if c)
    print(f"k={args.k}  DOF={len(coords)} {coords}   (n+1={len(coords)+1})")
    print(f"{'subject/trial/side':28s} rank  양의스패닝   |영공간 λ|")
    hits, total = 0, 0
    for sub, trial in DEFAULT_TRIALS:
        for side in ("r", "l"):
            try:
                V = torque_vectors(sub, trial, side, args.k, coords, args.out_root)
            except Exception as exc:
                print(f"{sub}/{trial}/{side}: {exc}"); continue
            rank, ok, lam = positive_spanning(V); total += 1; hits += bool(ok)
            lam_s = "  ".join(f"{x:.2f}" for x in lam) if lam is not None else "-"
            verdict = {True: "예", False: "아니오", None: "k≠n+1"}[ok]
            print(f"{sub + '/' + trial.replace('_segment_0', '') + '/' + side:28s} {rank:4d}   {verdict:>5s}    {lam_s}")
    print(f"\n양의 스패닝 성립: {hits}/{total}")


if __name__ == "__main__":
    main()
