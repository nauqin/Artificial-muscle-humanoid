"""관절각만으로 지면반력을 추정한다 (힘판 없는 데이터: OpenCap·영상). **`osim` 환경.**

    python tools/grf_estimate.py --model <osim> --ik <ik.mot> --out <grf.mot>       # 한 건
    python tools/grf_estimate.py --validate                                         # 힘판이 있는 우리 데이터로 정확도 재기

원리 (뉴턴–오일러, 몸 전체):
    발이 받는 힘의 합      F   = M · (무게중심 가속도 − 중력)
    무게중심 둘레 모멘트   M_c = d(전신 각운동량)/dt
이 둘은 관절각 q(t) 와 모델(질량·관성)만 있으면 정해진다. 남는 문제는 "두 발에 어떻게 나누나" 뿐이다.
  - 합력의 작용점 p (바닥 위, 수평 모멘트가 0 이 되는 점)를 구한다.
  - 발마다 발바닥 선분(뒤꿈치–발끝)에서 p 에 가장 가까운 점 q_r, q_l 을 찾는다.
  - p 가 한 발 위에 있으면 그 발이 전부 받는다(한 발 지지). 두 발 사이에 있으면 q_l→q_r 선 위의 위치 비율로 나눈다(지렛대 비).
  - 발이 바닥에서 4 cm 넘게 떠 있으면 그 발 몫은 0 (부드럽게 줄인다).
  - p 가 발 밖으로 3 cm 넘게 나가면 발 가장자리로 당긴다(clamp). 이때 남는 모멘트는 골반 잔차로 간다 — 숨기지 않고 보고한다.
이렇게 만든 지면반력으로 같은 모델·같은 관절각의 ID 를 돌리면 골반 잔차가 (clamp·공중 구간을 빼면) 0 에 가깝다.
즉 "같은 모델 짝" 원칙을 지킨다. 대신 **관절각의 잡음이 그대로 힘이 된다** — 필터 컷오프(기본 6 Hz)가 중요하다.

정확도 (2026-11-19, `--validate`, 40 kg·1.60 m 로봇 뼈대, 힘판 값과 비교): 결과는 docs 의 표. 요약하면 스쿼트·서서 하는 동작은
수직력·관절 토크가 잘 맞고, 걷기는 두 발 지지 구간의 배분과 앞뒤 힘이 덜 맞는다.

한계: (1) 두 발 지지에서 배분은 가정이다(정역학적으로 정해지지 않는다). (2) 손·엉덩이 등 발 말고 다른 곳이 닿는 동작(의자, 난간,
벽 짚기)은 못 다룬다 — 그 힘이 전부 발로 간 것처럼 나온다. (3) 발이 둘 다 뜬 구간은 힘 0. (4) 미끄러지는 바닥·경사는 가정 밖.
"""
import argparse, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import addb_io as io                                   # noqa: E402

G = 9.80665
FEET = ("calcn_r", "calcn_l")
LIFT_H0, LIFT_W = 0.04, 0.008                          # 발이 이 높이(m)를 넘으면 뜬 것으로 본다, 전환 폭
COP_MARGIN = 0.04                                      # 합력 작용점이 발 선분에서 이만큼(m)까지 벗어나는 것은 허용 (발 폭)
MIN_FY = 0.03                                          # 수직력이 체중의 이 비율보다 작으면 공중으로 본다
PIVOT = 1.4                                            # 두 발 지지에서 각 발의 힘이 향하는 점의 높이 = 무게중심 높이 × 이 값 (0 이면 끔). --validate 로 고름


def smooth_resample(mot, fs=100.0, lowpass=6.0):
    """.mot → 균일 시간축 (fs Hz), 라디안, 영위상 4차 Butterworth 저역통과. 반환 (names, t, Q, Qd)."""
    from scipy.signal import butter, filtfilt
    from scipy.interpolate import CubicSpline
    D = mot.data.copy()
    if mot.in_degrees:
        for j, n in enumerate(mot.names):
            if n not in ("pelvis_tx", "pelvis_ty", "pelvis_tz"): D[:, j] = np.radians(D[:, j])
    t = np.arange(mot.time[0], mot.time[-1] + 1e-9, 1.0 / fs)
    Q = CubicSpline(mot.time, D, axis=0)(t)
    if lowpass and lowpass > 0 and len(t) > 30:
        b, a = butter(4, lowpass / (fs / 2.0)); Q = filtfilt(b, a, Q, axis=0, padlen=min(len(t) - 1, int(fs)))
    Qd = CubicSpline(t, Q, axis=0)(t, 1)
    return list(mot.names), t, Q, Qd


class Body:
    """모델을 세워 프레임마다 무게중심·운동량·발바닥 점을 잰다."""
    def __init__(self, model_path, feet=FEET):
        import opensim as osim
        osim.Logger.setLevelString("Error"); self.osim = osim
        proc = osim.ModelProcessor(os.path.abspath(model_path)); proc.append(osim.ModOpRemoveMuscles())
        self.model = proc.process(); self.state = self.model.initSystem()
        self.mass = self.model.getTotalMass(self.state); self.cs = self.model.getCoordinateSet()
        self.bodies = [self.model.getBodySet().get(i) for i in range(self.model.getBodySet().getSize())]
        self.m = np.array([b.getMass() for b in self.bodies])
        self.com_local = [b.getMassCenter() for b in self.bodies]
        self.I = []
        for b in self.bodies:
            mo = b.getInertia().getMoments(); pr = b.getInertia().getProducts()
            self.I.append(np.array([[mo.get(0), pr.get(0), pr.get(1)], [pr.get(0), mo.get(1), pr.get(2)], [pr.get(1), pr.get(2), mo.get(2)]]))
        # 발바닥 점 (calcn 좌표계): 뒤꿈치, 중족골두(mtp 관절), 발끝. calcn 원점이 뒤꿈치 바닥이다.
        self.sole = {}
        for f in feet:
            side = f[-1]; js = self.model.getJointSet(); mtp = js.get(f"mtp_{side}").get_frames(0).get_translation()
            mx = mtp.get(0); self.sole[f] = [osim.Vec3(0.0, 0.0, 0.0), osim.Vec3(mx, 0.0, 0.0), osim.Vec3(mx * 1.38, 0.0, 0.0)]
        self.feet = feet

    def set(self, names, q, qd):
        cs = self.cs
        for n, v, vd in zip(names, q, qd):
            if not cs.contains(n): continue
            c = cs.get(n)
            if c.get_locked(): continue
            c.setValue(self.state, float(v), False); c.setSpeedValue(self.state, float(vd))
            if cs.contains(n + "_beta"):                 # 슬개골 종속 좌표 (무릎각과 1:1)
                cb = cs.get(n + "_beta"); cb.setValue(self.state, float(v), False); cb.setSpeedValue(self.state, float(vd))
        self.model.realizeVelocity(self.state)

    def momentum(self):
        """(무게중심 위치, 선운동량/질량 = 무게중심 속도, 무게중심 둘레 각운동량)"""
        s = self.state; P = np.zeros((len(self.bodies), 3)); V = np.zeros_like(P); H = np.zeros(3); spin = []
        for i, b in enumerate(self.bodies):
            p = b.findStationLocationInGround(s, self.com_local[i]); v = b.findStationVelocityInGround(s, self.com_local[i])
            P[i] = [p.get(k) for k in range(3)]; V[i] = [v.get(k) for k in range(3)]
            w = b.getVelocityInGround(s).get(0); R = b.getTransformInGround(s).R()
            Rm = np.array([[R.get(r, c) for c in range(3)] for r in range(3)]); wv = np.array([w.get(k) for k in range(3)])
            spin.append(Rm @ self.I[i] @ Rm.T @ wv)
        c = (self.m[:, None] * P).sum(0) / self.mass; vc = (self.m[:, None] * V).sum(0) / self.mass
        H = (np.cross(P - c, self.m[:, None] * V)).sum(0) + np.sum(spin, axis=0)
        return c, vc, H

    def sole_points(self):
        out = {}
        for f in self.feet:
            b = self.model.getBodySet().get(f); pts = []
            for loc in self.sole[f]:
                p = b.findStationLocationInGround(self.state, loc); pts.append([p.get(k) for k in range(3)])
            out[f] = np.array(pts)
        return out


def closest_on_segment(p, a, b):
    """수평면(x, z)에서 선분 a–b 위 p 에 가장 가까운 점."""
    ab = b - a; L2 = float(ab @ ab)
    s = 0.0 if L2 < 1e-12 else float(np.clip((p - a) @ ab / L2, 0.0, 1.0))
    return a + s * ab


def estimate(model_path, ik_path, out_path=None, *, lowpass=6.0, fs=100.0, floor="auto", feet=FEET, log=print):
    """지면반력 추정. 반환 dict(time, names, Q(필터된 관절각), grf(열 이름→배열), diag). out_path 를 주면 .mot 저장."""
    mot = io.read_mot(ik_path); names, t, Q, Qd = smooth_resample(mot, fs, lowpass)
    body = Body(model_path, feet); T = len(t)
    C = np.zeros((T, 3)); VC = np.zeros((T, 3)); H = np.zeros((T, 3)); SP = {f: np.zeros((T, 3, 3)) for f in feet}
    for i in range(T):
        body.set(names, Q[i], Qd[i]); C[i], VC[i], H[i] = body.momentum()
        sp = body.sole_points()
        for f in feet: SP[f][i] = sp[f]
    dt = 1.0 / fs; A = np.gradient(VC, dt, axis=0); Hd = np.gradient(H, dt, axis=0)
    M = body.mass; W = M * G
    F = M * (A + np.array([0.0, G, 0.0]))                       # 발이 받아야 하는 합력
    Mo = np.cross(C, F) + Hd                                    # 지면 원점 둘레 모멘트
    low = np.min([SP[f][:, :, 1].min(axis=1) for f in feet], axis=0)
    y0 = float(np.percentile(low, 5)) if floor == "auto" else float(floor)

    cols = {f"ground_force_{f}_{k}{a}": np.zeros(T) for f in feet for k in "vpm" for a in "xyz"}
    clamp = np.zeros(T); flight = np.zeros(T, bool); share = np.zeros((T, 2)); cop = np.full((T, 3), np.nan)
    for i in range(T):
        h = np.array([SP[f][i][:, 1].min() - y0 for f in feet])
        g = 1.0 / (1.0 + np.exp((h - LIFT_H0) / LIFT_W))        # 접촉 정도 0~1
        if F[i, 1] < MIN_FY * W or g.sum() < 0.05: flight[i] = True; continue
        # 합력 작용점 p (y = y0): Mo = p × F + T_y ŷ  →  수평 성분으로 p_x, p_z
        px = (Mo[i, 2] + y0 * F[i, 0]) / F[i, 1]; pz = (y0 * F[i, 2] - Mo[i, 0]) / F[i, 1]
        p = np.array([px, pz]); cop[i] = [px, y0, pz]
        seg = [(SP[f][i][0][[0, 2]], SP[f][i][2][[0, 2]]) for f in feet]
        q = [closest_on_segment(p, a, b) for a, b in seg]; delta = np.zeros(2)
        for _ in range(4):                                      # 발 중심선 위 두 점 q 와 배분 w 를 w0·q0 + w1·q1 = p 가 되게 맞춰 간다
            d = q[1] - q[0]; L2 = float(d @ d)
            s = 0.5 if L2 < 1e-10 else float(np.clip((p - q[0]) @ d / L2, 0.0, 1.0))  # 0 = 첫째 발, 1 = 둘째 발
            w = np.array([(1 - s) * g[0], s * g[1]])
            if w.sum() < 1e-6: w = g.copy()                     # p 쪽 발이 떠 있으면 닿은 발이 받는다
            w = w / w.sum(); delta = p - (w[0] * q[0] + w[1] * q[1])
            q = [closest_on_segment(q[k] + delta, *seg[k]) for k in range(2)]
        share[i] = w; delta = p - (w[0] * q[0] + w[1] * q[1]); nd = float(np.linalg.norm(delta))
        if nd > COP_MARGIN: clamp[i] = nd - COP_MARGIN; delta = delta * COP_MARGIN / nd
        Ty = Mo[i, 1] - (pz * F[i, 0] - px * F[i, 2])            # 수직축 자유모멘트
        pk = [q[k] + delta for k in range(2)]; Fk = [F[i] * w[k] for k in range(2)]
        # 두 발이 서로 밀어내는 힘 (합력·합모멘트는 그대로): 각 발의 힘이 무게중심 위 한 점(pivot)을 향하도록.
        # 두 발 지지에서 힘의 합만으로는 이 성분을 알 수 없다. 사람은 벌려 선 발로 바깥을, 걷기에서 뒷발은 앞으로·앞발은 뒤로 민다.
        dd = pk[1] - pk[0]; Ld = float(np.linalg.norm(dd))
        if PIVOT and Ld > 0.05:
            dh = dd / Ld; Hp = PIVOT * (C[i, 1] - y0); P = C[i, [0, 2]]
            c = [float(Fk[k][[0, 2]] @ dh) for k in range(2)]; hdes = [Fk[k][1] * float((P - pk[k]) @ dh) / Hp for k in range(2)]
            f_int = 4 * w[0] * w[1] * ((c[0] - hdes[0]) - (c[1] - hdes[1])) / 2.0
            Fk[0] = Fk[0] - f_int * np.array([dh[0], 0.0, dh[1]]); Fk[1] = Fk[1] + f_int * np.array([dh[0], 0.0, dh[1]])
        for k, f in enumerate(feet):
            for a, v in zip("xyz", Fk[k]): cols[f"ground_force_{f}_v{a}"][i] = v
            cols[f"ground_force_{f}_px"][i] = pk[k][0]; cols[f"ground_force_{f}_py"][i] = y0; cols[f"ground_force_{f}_pz"][i] = pk[k][1]
            cols[f"ground_force_{f}_my"][i] = Ty * w[k]
    contact = ~flight
    diag = dict(mass_kg=round(M, 2), floor_y=round(y0, 4), frames=T, duration_s=round(float(t[-1] - t[0]), 2),
                fy_peak_bw=round(float(F[contact, 1].max() / W), 3) if contact.any() else 0.0,
                fy_mean_bw=round(float(F[contact, 1].mean() / W), 3) if contact.any() else 0.0,
                flight_pct=round(100.0 * flight.mean(), 1), clamp_pct=round(100.0 * (clamp > 0).mean(), 1),
                clamp_max_cm=round(100.0 * float(clamp.max()), 1), fx_peak_bw=round(float(np.abs(F[contact][:, [0, 2]]).max() / W), 3) if contact.any() else 0.0)
    log(f"  지면반력 추정: 질량 {M:.1f} kg, 바닥 y={y0:+.3f} m, 수직력 평균 {diag['fy_mean_bw']:.2f}·최대 {diag['fy_peak_bw']:.2f} BW, 수평력 최대 {diag['fx_peak_bw']:.2f} BW, "
        f"공중 {diag['flight_pct']:.0f} %, 작용점이 발 밖 {diag['clamp_pct']:.0f} % (최대 {diag['clamp_max_cm']:.1f} cm)")
    if abs(diag["fy_mean_bw"] - 1.0) > 0.08 and diag["flight_pct"] < 5:
        log(f"  [주의] 수직력 평균이 체중의 {diag['fy_mean_bw']:.2f}배 — 구간 앞뒤로 몸이 오르내리는 중이거나 골반 높이에 추세(드리프트)가 있다")
    if diag["clamp_pct"] > 20: log("  [주의] 합력 작용점이 발 밖으로 나가는 프레임이 20 % 를 넘는다 — 관절각 잡음이 크거나 발 말고 다른 곳이 닿는 동작")
    cn = list(cols); Gm = np.column_stack([cols[c] for c in cn])
    if out_path: io.write_mot(out_path, cn, t, Gm)
    return dict(time=t, names=names, Q=Q, grf_names=cn, grf=Gm, diag=diag, com=C, share=share, flight=flight, clamp=clamp, cop=cop, F=F)


# ─────────────────────────────────────────────────────────────────────────────
# 검증: 힘판이 있는 우리 데이터에서 추정값과 비교
# ─────────────────────────────────────────────────────────────────────────────

VALID = [  # (이름, 모델, ik, 힘판 grf, 힘판 ID, 컷오프)
    ("스쿼트 (LabValidation subject3)", "data/robot40_h160/robot40_h160.osim", "data/robot40_h160/lower/squats1_segment_0_ik.mot",
     "data/robot40_h160/lower/squats1_segment_0_grf.mot", "out/lower40_h160/subject3/squats1_segment_0_id.sto", 6.0),
    ("걷기 (LabValidation subject3)", "data/robot40_h160/robot40_h160.osim", "data/robot40_h160/lower/walking1_segment_0_ik.mot",
     "data/robot40_h160/lower/walking1_segment_0_grf.mot", "out/lower40_h160/subject3/walking1_segment_0_id.sto", 6.0),
    ("걷기 트레드밀 (Moore T047_12)", "data/robot40_h160/robot40_h160.osim", "data/robot40_h160/moore/T047_segment_12_ik.mot",
     "data/robot40_h160/moore/T047_segment_12_grf.mot", "out/lower40_h160/moore/T047_segment_12_id.sto", 6.0),
]
ID_COORDS = ("hip_flexion", "hip_adduction", "knee_angle", "ankle_angle")


def validate(out_csv=None, log=print):
    """힘판 값·힘판 ID 와 비교. 추정 지면반력으로 ID 를 다시 돌려 관절 토크 차이를 잰다."""
    import csv, tempfile
    import run_trial as rt
    rows = []
    for label, model, ik, grf, idp, lp in VALID:
        model, ik, grf, idp = [os.path.join(ROOT, x) for x in (model, ik, grf, idp)]
        if not all(os.path.exists(x) for x in (model, ik, grf, idp)): log(f"[건너뜀] {label}: 파일 없음"); continue
        log(f"\n=== {label} ===")
        tmp = tempfile.mkdtemp(prefix="grfval_"); est_path = os.path.join(tmp, "est_grf.mot"); ikf = os.path.join(tmp, "ikf.mot")
        r = estimate(model, ik, est_path, lowpass=lp, log=log)
        io.write_mot(ikf, r["names"], r["time"], r["Q"])
        ref = io.read_mot(grf); idr = io.read_mot(idp); t0, t1 = idr.t0, idr.t1
        W = r["diag"]["mass_kg"] * G; sel = (r["time"] >= t0) & (r["time"] <= t1); tt = r["time"][sel]
        est = {n: r["grf"][sel, j] for j, n in enumerate(r["grf_names"])}
        refi = {n: np.interp(tt, ref.time, ref.column(n)) for n in ref.names}
        tot_e = sum(est[f"ground_force_{f}_vy"] for f in FEET); tot_r = sum(refi[f"ground_force_{f}_vy"] for f in FEET)
        hx_e = sum(est[f"ground_force_{f}_vx"] for f in FEET); hx_r = sum(refi[f"ground_force_{f}_vx"] for f in FEET)
        per_foot = np.sqrt(np.mean([(est[f"ground_force_{f}_vy"] - refi[f"ground_force_{f}_vy"]) ** 2 for f in FEET])) / W * 100
        # 작용점 차이: 힘판 수직력이 체중 20 % 를 넘는 프레임에서
        cop_err = []
        for f in FEET:
            on = refi[f"ground_force_{f}_vy"] > 0.2 * W
            if on.any(): cop_err.append(np.hypot(est[f"ground_force_{f}_px"][on] - refi[f"ground_force_{f}_px"][on], est[f"ground_force_{f}_pz"][on] - refi[f"ground_force_{f}_pz"][on]).mean())
        # 추정 지면반력으로 ID (필터는 이미 했으므로 ID 툴 필터는 끈다)
        gm = io.read_mot(est_path); plates = {k: p for k, p in io.find_plates(gm).items() if p.complete}
        mapping = {k: {"body": k} for k in plates}
        ext = rt.write_external_loads(os.path.join(tmp, "ext.xml"), est_path, plates, mapping)
        id_est, _ = rt.run_inverse_dynamics(model, ikf, ext, tmp, "est", t0, t1, lowpass=-1)
        ide = io.read_mot(id_est); row = dict(trial=label, t0=round(t0, 2), t1=round(t1, 2), fy_total_rmse_pctBW=round(float(np.sqrt(np.mean((tot_e - tot_r) ** 2)) / W * 100), 1),
                                              fx_total_rmse_pctBW=round(float(np.sqrt(np.mean((hx_e - hx_r) ** 2)) / W * 100), 1), fy_per_foot_rmse_pctBW=round(float(per_foot), 1),
                                              cop_err_cm=round(float(np.mean(cop_err)) * 100, 1) if cop_err else None)
        log(f"  수직력 합 오차 {row['fy_total_rmse_pctBW']} %BW, 앞뒤 힘 합 {row['fx_total_rmse_pctBW']} %BW, 발별 수직력 {row['fy_per_foot_rmse_pctBW']} %BW, 작용점 {row['cop_err_cm']} cm")
        res = np.sqrt(sum(ide.column(f"pelvis_t{a}_force") ** 2 for a in "xyz")); row["residual_pct_peakF"] = round(float(res.max() / max(tot_e.max(), 1e-9) * 100), 1)
        resm = np.sqrt(sum(ide.column(f"pelvis_{a}_moment") ** 2 for a in ("tilt", "list", "rotation"))); row["residual_moment_max_Nm"] = round(float(resm.max()), 1)
        log(f"  추정 지면반력 ID 의 골반 잔차: 힘 최대 {row['residual_pct_peakF']} % (수직력 피크 대비), 모멘트 최대 {row['residual_moment_max_Nm']} N·m")
        for c in ID_COORDS:
            for side in "rl":
                n = f"{c}_{side}_moment"; a = np.interp(tt, ide.time, ide.column(n)); b = np.interp(tt, idr.time, idr.column(n)); pk = np.abs(b).max()
                row[f"{c}_{side}_rmse_Nm"] = round(float(np.sqrt(np.mean((a - b) ** 2))), 2); row[f"{c}_{side}_peak_ref_Nm"] = round(float(pk), 1); row[f"{c}_{side}_peak_est_Nm"] = round(float(np.abs(a).max()), 1)
            log(f"  {c:14s} 토크 RMSE 오른쪽 {row[f'{c}_r_rmse_Nm']:5.1f} / 왼쪽 {row[f'{c}_l_rmse_Nm']:5.1f} N·m   (힘판 ID 피크 {row[f'{c}_r_peak_ref_Nm']:.0f} / {row[f'{c}_l_peak_ref_Nm']:.0f}, 추정 피크 {row[f'{c}_r_peak_est_Nm']:.0f} / {row[f'{c}_l_peak_est_Nm']:.0f})")
        rows.append(row)
    if out_csv and rows:
        os.makedirs(os.path.dirname(os.path.abspath(out_csv)), exist_ok=True)
        with open(out_csv, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
        log(f"\n저장 {out_csv}")
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model"); ap.add_argument("--ik"); ap.add_argument("--out")
    ap.add_argument("--lowpass", type=float, default=6.0); ap.add_argument("--floor", default="auto")
    ap.add_argument("--validate", action="store_true", help="힘판이 있는 우리 데이터 3건으로 정확도를 잰다")
    a = ap.parse_args()
    if a.validate: validate(os.path.join(ROOT, "report", "custom", "grf_estimate_validation.csv")); return
    if not (a.model and a.ik and a.out): ap.error("--model --ik --out 이 필요하다 (또는 --validate)")
    estimate(a.model, a.ik, a.out, lowpass=a.lowpass, floor=a.floor)


if __name__ == "__main__":
    main()
