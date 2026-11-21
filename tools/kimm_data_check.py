"""기계연 원자료(팔 실험 각도 4개 + 근육 하나 5 kg 변위)와 김은서 님 SMA 모델(kimm_arm_new.py)을 시간 파형으로 비교 (2026-11-21).
목적: 모델이 근육 하나 0.5 Hz 시험에 맞춘 것이라 빠른 구간(1 Hz ≈ 걷기 주기)에서 얼마나 어긋나는지, 열 속도를 몇 배로 보면 실측에 맞는지 → 걷기 검증의 범위.
  원자료: '기계연구원 실험 파일 공유 건/' (공개 저장소에 올리지 않는다). 모델: git origin/main 의 sma_arm_kimm/ (복사하지 않고 읽는다).
  열 속도 배율 k = 가열·냉각(강제·반냉각·켜는 동안)을 모두 k 배 → 온도가 변하는 속도만 k 배, 도달 온도는 같다.
    /c/Users/USER/miniforge3/envs/mj/python.exe tools/kimm_data_check.py [--k 1,1.5,2,3]
출력: report/mujoco160/kimm_check.json, fig_kimm_check.png
"""
import argparse, json, os, subprocess, sys, types, zipfile, re
os.environ.setdefault("MKL_THREADING_LAYER", "SEQUENTIAL")
import numpy as np                                           # noqa: E402
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); DD = os.path.join(ROOT, "기계연구원 실험 파일 공유 건"); OUT = os.path.join(ROOT, "report", "mujoco160")
FILES = {0.125: "0.5kg_0.125hz_Angle.xlsx", 0.25: "0.5kg_0.25hz_Angle.xlsx.xlsx", 0.5: "0.5kg_0.5hz_Angle.xlsx.xlsx", 1.0: "0.5kg_1hz_Angle.xlsx.xlsx"}


def read_xlsx(path, sheet=1):
    """openpyxl 없이 숫자 시트 읽기 (mj 환경에 openpyxl 이 없다)."""
    z = zipfile.ZipFile(path); xml = z.read(f"xl/worksheets/sheet{sheet}.xml").decode("utf-8"); rows = []
    for r in re.findall(r"<row[^>]*>(.*?)</row>", xml, flags=re.S):
        vals = [v for c in re.findall(r"<c\b[^>]*?>(.*?)</c>", r, flags=re.S) for v in re.findall(r"<v>([^<]*)</v>", c)]
        try: rows.append([float(v) for v in vals if v != ""])
        except ValueError: pass
    n = max(len(r) for r in rows); return np.array([r for r in rows if len(r) == n])


def cycle_avg(t, y, f, skip=2, nb=100):
    """주기 1/f 로 접어 평균 파형 (처음 skip 주기 제외). 위상은 최솟값이 0 에 오게 돌린다."""
    per = 1.0 / f; m = t >= t[0] + skip * per; ph = ((t[m] - t[0]) % per) / per; b = np.minimum((ph * nb).astype(int), nb - 1)
    w = np.array([y[m][b == i].mean() if np.any(b == i) else np.nan for i in range(nb)]); ok = ~np.isnan(w); w = np.interp(np.arange(nb), np.arange(nb)[ok], w[ok], period=nb)
    return np.roll(w, -int(np.argmin(w)))


def pp_per_cycle(t, y, f, skip=2):
    per = 1.0 / f; out = []; k = skip
    while t[0] + (k + 1) * per <= t[-1]:
        m = (t >= t[0] + k * per) & (t < t[0] + (k + 1) * per); out.append(y[m].max() - y[m].min()); k += 1
    return np.array(out)


def load_model():
    g = lambda p: subprocess.run(["git", "-C", ROOT, "show", f"origin/main:sma_arm_kimm/{p}"], capture_output=True, check=True).stdout.decode("utf-8")
    sm = types.ModuleType("sma_model"); exec(compile(g("sma_model.py"), "sma_model.py", "exec"), sm.__dict__); sys.modules["sma_model"] = sm
    arm = types.ModuleType("kimm_arm_new"); arm.__dict__["__name__"] = "kimm_arm_new"; exec(compile(g("kimm_arm_new.py"), "kimm_arm_new.py", "exec"), arm.__dict__)
    return sm, arm


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--k", default="1,1.5,2,3"); a = ap.parse_args(); ks = [float(x) for x in a.k.split(",")]
    meas = {}
    for f, fn in FILES.items():
        d = read_xlsx(os.path.join(DD, fn)); t, y = d[:, 0], d[:, -1]; on = np.argmax(y > 1.0); t0 = max(t[on] - 0.5 / f * 0.0, 0.0); m = t >= t0
        pp = pp_per_cycle(t[m], y[m], f); meas[f] = dict(t=t[m] - t0, y=y[m], pp=pp, wave=cycle_avg(t[m], y[m], f))
        print(f"실측 {f:5.3f} Hz: 기록 {t[-1]:.0f} s, 주기별 굽힘 평균 {pp.mean():.1f}° (처음 3회 {np.round(pp[:3], 1)}, 마지막 3회 {np.round(pp[-3:], 1)}), 최대 {y.max():.1f}°, 주기 안 최솟값 평균 {np.mean([meas[f]['wave'].min()]):.1f}°", flush=True)
    sm, arm = load_model(); P0, HF0, HON0 = arm.P36, arm.HF, arm.H_ON; res = {}; waves = {}
    for k in ks:
        arm.P36, arm.HF, arm.H_ON, arm.H_HALF = P0 * k, HF0 * k, HON0 * k, HF0 * k / 2; row = {}
        for f in FILES:
            r = arm.run(f); L = r["log"]; row[f] = float(r["bend"]); waves[(k, f)] = cycle_avg(L[:, 0], L[:, 1], f)
        e = float(np.sqrt(np.mean([(row[f] - meas[f]["pp"].mean()) ** 2 for f in FILES]))); res[k] = dict(bend=row, rmse_deg=e)
        print(f"모델 열 속도 ×{k:g}: 굽힘 " + " / ".join(f"{row[f]:.1f}" for f in FILES) + f"°  (실측 " + " / ".join(f"{meas[f]['pp'].mean():.1f}" for f in FILES) + f"°)  오차 RMS {e:.1f}°", flush=True)
    json.dump(dict(measured_bend_deg={str(f): float(meas[f]["pp"].mean()) for f in FILES}, model={str(k): v for k, v in ((k, dict(bend={str(f): b for f, b in v["bend"].items()}, rmse_deg=v["rmse_deg"])) for k, v in res.items())}),
              open(os.path.join(OUT, "kimm_check.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from matplotlib import font_manager
    for fnt in ("Malgun Gothic", "NanumGothic"):
        if any(x.name == fnt for x in font_manager.fontManager.ttflist): plt.rcParams["font.family"] = fnt; break
    plt.rcParams["axes.unicode_minus"] = False
    fig, ax = plt.subplots(1, 5, figsize=(24, 4.8)); x = np.arange(100); cols = ["#c0392b", "#e67e22", "#27ae60", "#2980b9", "#8e44ad"]
    for i, f in enumerate(FILES):
        ax[i].plot(x, meas[f]["wave"], "k-", lw=2.5, label="기계연 실측 (주기 평균)")
        for j, k in enumerate(ks): ax[i].plot(x, waves[(k, f)], color=cols[j % 5], lw=1.5, label=f"모델, 열 속도 ×{k:g}")
        ax[i].set_title(f"팔 실험 {f:g} Hz (주기 {1 / f:g} s)", weight="bold"); ax[i].set_xlabel("주기 [%] (가장 편 순간 = 0)"); ax[i].set_ylabel("팔꿈치 굽힘 [°]"); ax[i].grid(alpha=0.3)
        if i == 0: ax[i].legend(fontsize=8)
    fs = list(FILES); ax[4].semilogx(fs, [meas[f]["pp"].mean() for f in fs], "ko-", lw=2.5, label="기계연 실측")
    for j, k in enumerate(ks): ax[4].semilogx(fs, [res[k]["bend"][f] for f in fs], "o-", color=cols[j % 5], label=f"모델 ×{k:g} (오차 {res[k]['rmse_deg']:.1f}°)")
    ax[4].set_title("빠르기별 굽힘 크기", weight="bold"); ax[4].set_xlabel("빠르기 [Hz]"); ax[4].set_ylabel("굽힘 [°]"); ax[4].grid(alpha=0.3, which="both"); ax[4].legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig_kimm_check.png"), dpi=110); print("그림 → report/mujoco160/fig_kimm_check.png")


if __name__ == "__main__":
    main()
