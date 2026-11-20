"""`.trc` (마커 3D 궤적) 읽고 쓰기 — OpenSim 없이 numpy 만."""
import os
import numpy as np


def read_trc(path):
    """반환 (names, time (T,), X (T, M, 3) m). 단위가 mm 면 m 로 바꾼다. 빈 칸은 NaN."""
    with open(path, encoding="utf-8", errors="replace") as fh: lines = fh.read().splitlines()
    keys = lines[1].split("\t"); vals = lines[2].split("\t"); meta = dict(zip([k.strip() for k in keys], [v.strip() for v in vals]))
    names = [n for n in lines[3].split("\t")[2:] if n.strip()]
    rows = []
    for ln in lines[5:]:
        if not ln.strip(): continue
        c = ln.split("\t")
        rows.append([float(x) if x.strip() else np.nan for x in c[:2 + 3 * len(names)]] + [np.nan] * max(0, 2 + 3 * len(names) - len(c)))
    A = np.asarray(rows, float); X = A[:, 2:2 + 3 * len(names)].reshape(len(A), len(names), 3)
    if meta.get("Units", "m").lower() == "mm": X = X / 1000.0
    return names, A[:, 1], X


def write_trc(path, names, time, X, rate=None):
    """X (T, M, 3) m. NaN 은 빈 칸으로 쓴다."""
    rate = rate or (1.0 / float(np.median(np.diff(time))) if len(time) > 1 else 60.0)
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(f"PathFileType\t4\t(X/Y/Z)\t{os.path.basename(path)}\n")
        fh.write("DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n")
        fh.write(f"{rate:g}\t{rate:g}\t{len(time)}\t{len(names)}\tm\t{rate:g}\t1\t{len(time)}\n")
        fh.write("Frame#\tTime\t" + "\t".join(f"{n}\t\t" for n in names) + "\n")
        fh.write("\t\t" + "\t".join(f"X{k+1}\tY{k+1}\tZ{k+1}" for k in range(len(names))) + "\n\n")
        for i in range(len(time)):
            vals = ["\t".join("" if np.isnan(a) else f"{a:.5f}" for a in X[i, j]) for j in range(len(names))]
            fh.write(f"{i+1}\t{time[i]:.5f}\t" + "\t".join(vals) + "\n")
    return path
