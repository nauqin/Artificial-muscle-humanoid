"""아래팔 회내외 전용 케이블 2개 (회내·회외) — 노뼈 드럼(원통 감김) 설계와 Moco 확인 (2026-11-19, README 13-10).

왜: 13-9 에서 회내외는 직선 케이블로 안 됐다 (회내외 토크가 작아 LP 가 전용 케이블을 안 두고, 노뼈 직선 케이블의 회내외 모멘트암 1~2 mm).
분리: 회내외 관절은 자뼈–노뼈 사이라 팔꿈치 케이블(위팔→자뼈)의 회내외 모멘트는 정확히 0, 자뼈→노뼈 케이블의 팔꿈치 모멘트도 정확히 0.
  → 팔꿈치 2개는 175 cm 최종 세트(td_br50_rb, report/td/best_set.md)를 그대로 두고 회내외 2개만 설계한다.
드럼: 노뼈에 회내외 축과 같은 축의 원통 (반지름 R, OpenSim WrapCylinder, 케이블마다 하나 — 홈 두 개짜리 드럼).
  케이블 = 자뼈 가이드 A → 드럼을 감아 → 노뼈 끝점 B (드럼 표면 + b_off). 회내 케이블은 엎칠수록 감김이 풀려 짧아지고 (모멘트암 ≈ +R),
  회외 케이블은 반대 (≈ −R). 감김 각 θ 가 모든 자세에서 [θmin, θmax] 안이어야 접촉이 유지된다 (θ 는 회내외 각과 1:1 로 변한다).
변수 10개: 공통 R, zc (드럼 중심의 축 방향 위치, 노뼈 원점(노뼈머리) 기준, − = 먼쪽) + 케이블마다 A (축 거리 ρ, 둘레 각 φ, 홈 기준 축 방향 Δz),
  B 의 둘레 각 φB. A 는 회내외 0 (뒤침) 자세의 노뼈 축 좌표계로 정하고 자뼈 좌표로 옮긴다. 홈(B 의 축 방향)은 zc ± groove.
설계 쪽 목적 = 두 케이블 모듈 일 용량 (timedomain.size_cable; 장력 = 1.1 × 회내외 토크 ÷ 그 순간 모멘트암, 길항 없음 — 분리라 정확한 최소 장력)
  + 감김 각 범위 벌점 + 자뼈 가이드 브래킷 범위 벌점. 프레임 = 과제 시행 (dt 0.02) + 강건성 (5 프레임마다 회내·회외 ±X, × 1.1, path_design --robust-frac 과 같은 꼴).
Moco: 팔꿈치 2 + 회내외 2 (좌우) 모델로 MocoInverse (run_moco.sh, 175 cm). --outer N: Moco 에서 회내외 보조가 남은 순간의
  (자세, 케이블 토크 + 보조) 를 설계 프레임으로 더해 다시 설계. --moco-nm N: Moco 회내외 보조 95퍼센타일을 목적으로 Nelder–Mead (평가 N 회, 1 회 = Moco 2 건).

    python report/td/prosup_wrap.py --robust-nm 0.25 --out-tag _prosup_wrap [--outer 3] [--moco-nm 0] [--design-only] [--final]
→ report/newmus_paths_forearm_r_k4_40kg_h175sk<tag>.json, report/moco/newmus_paths_forearm_r_k4_40kg_h175sk<tag>{,_load3.79}_model.osim,
  report/td/prosup_wrap<tag>/ (OpenSim 길이 비교 csv·그림, 반복 Moco 지표·로그)
"""
import argparse, csv, json, os, shutil, subprocess, sys, time
import numpy as np
from scipy.optimize import differential_evolution, minimize
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import path_design as pd                                     # noqa: E402
import timedomain as td                                      # noqa: E402
import matplotlib.pyplot as plt                              # noqa: E402  (글꼴은 path_design 이 정한다)

MODEL = "data/robot40_h175sk/robot40_h175sk.osim"; TASKS = "report/upper/tasks_40kg_h175sk.csv"
ELBOW = "report/newmus_paths_elbow_r_k2_40kg_h175sk_td_br50_rb.json"; STEM = "newmus_paths_forearm_r_k4_40kg_h175sk"
MOCO_ENV = dict(B="1.0802", UP="out/upper40_h175sk", LO="out/lower40_h175sk", DR="data/robot40_h175sk")   # run_forearm_175.sh 와 같다
CAB = (("M2", "회내", -1.0, +1.0), ("M3", "회외", +1.0, -1.0))   # 이름, 역할, 감김 방향 s (A→B 축 둘레, +1 반시계), 회내외 모멘트암 부호
PN = ["R", "zc", "rhoA_pro", "phiA_pro", "dzA_pro", "phiB_pro", "rhoA_sup", "phiA_sup", "dzA_sup", "phiB_sup"]
J_PS = 1                                                     # forearm_r COORDS 안의 pro_sup_r 번호 → Frames 섭동 열 1+2·1 (+δ), 2+2·1 (−δ)


class Geo:
    """회내외 축 좌표계 (노뼈 좌표, 회내외 0 자세): 원점 O (관절 중심), z = 축 a (위쪽 = 몸쪽), x0 = 노뼈 x 를 축에 수직으로, y0 = a × x0. 자뼈 ← 노뼈 변환 (회내외 0)."""
    def __init__(self, model_path):
        import opensim as osim
        osim.Logger.setLevelString("Error"); m = osim.Model(model_path); s = m.initSystem()
        T = m.getJointSet().get("radioulnar_r").getChildFrame().findTransformInBaseFrame()
        Rc = np.array([[T.R().get(i, k) for k in range(3)] for i in range(3)]); self.O = np.array([T.p().get(i) for i in range(3)])
        self.a = Rc[:, 2] / np.linalg.norm(Rc[:, 2]); x = np.eye(3)[0] - self.a[0] * self.a; self.x0 = x / np.linalg.norm(x); self.y0 = np.cross(self.a, self.x0)
        m.getCoordinateSet().get("pro_sup_r").setValue(s, 0.0); m.realizePosition(s); B = m.getBodySet(); rad, uln = B.get("radius_r"), B.get("ulna_r")
        self.t_ur = rad.findStationLocationInAnotherFrame(s, osim.Vec3(0, 0, 0), uln).to_numpy()
        self.R_ur = np.column_stack([rad.findStationLocationInAnotherFrame(s, osim.Vec3(*e), uln).to_numpy() - self.t_ur for e in np.eye(3)])

    def rpt(self, z, rho, phi):
        return self.O + z * self.a + rho * (np.cos(phi) * self.x0 + np.sin(phi) * self.y0)

    def to_ulna(self, p):
        return self.t_ur + self.R_ur @ p


def bounds(a):
    pi = np.pi; b = [(a.r_min, a.r_max), (a.zc_min, a.zc_max)]
    for _ in CAB: b += [(0.016, 0.05), (-pi, pi), (-a.dz_max, a.dz_max), (-pi, pi)]
    return b


def cables(x, geo, a, theta_max=None):
    """변수 x → 회내외 케이블 2개 (json 꼴). theta_max (케이블별, rad) 를 주면 원통 좌표계 x 를 그 케이블이 감는 호 전체의 가운데로 (quadrant +x)."""
    R, zc = float(x[0]), float(x[1]); out = []
    for c, (nm, role, s, sg) in enumerate(CAB):
        rhoA, phiA, dzA, phiB = (float(v) for v in x[2 + 4 * c: 6 + 4 * c]); zg = zc + a.groove * sg
        A_r = geo.rpt(zg + dzA, rhoA, phiA); B_r = geo.rpt(zg, R + a.b_off, phiB)
        aB = np.arccos(R / (R + a.b_off)); mid = phiB - s * aB - s * (0.5 * theta_max[c] if theta_max is not None else np.pi / 2)
        xd = np.cos(mid) * geo.x0 + np.sin(mid) * geo.y0; frame = np.column_stack([xd, np.cross(geo.a, xd), geo.a])
        wc = dict(body="radius_r", center=geo.rpt(zc, 0.0, 0.0).tolist(), axis=geo.a.tolist(), radius=R, sign=s, length=a.cyl_length, quadrant="+x", frame=frame.tolist())
        out.append(dict(name=nm, label=f"{role} (노뼈 드럼)", chain=["ulna_r", "radius_r"], points=[geo.to_ulna(A_r).tolist(), B_r.tolist()], wrap_cyl=wc))
    return out


class Design:
    """설계 쪽 계산: 자세 Frames (고유 자세) + 장력 프레임 (자세 번호, 회내외 토크 × 1.1)."""
    def __init__(self, a, geo):
        self.a, self.geo = a, geo; pd.configure("forearm_r")
        Q, T, meta = pd.poses_time_from_tasks(os.path.join(ROOT, a.tasks), dt=a.dt); self.meta = meta
        self.Q = Q; self.idx = np.arange(len(Q)); self.tau = pd.MARGIN * T[:, J_PS]; self.kind = np.zeros(len(Q), int)      # 0 데이터, 1 강건성, 2 Moco 추가
        iv = np.arange(0, len(Q), a.robust_every)
        for sg in (1.0, -1.0):
            self.idx = np.r_[self.idx, iv]; self.tau = np.r_[self.tau, np.full(len(iv), sg * pd.MARGIN * a.robust_nm)]; self.kind = np.r_[self.kind, np.ones(len(iv), int)]
        print(f"자세 {len(Q)}개 (시행 {len(meta['names'])}), 회내외 토크 {T[:, J_PS].min():+.2f} ~ {T[:, J_PS].max():+.2f} N·m, 각 {np.degrees(Q[:, J_PS].min()):.1f} ~ {np.degrees(Q[:, J_PS].max()):.1f}°, "
              f"강건성 ±{a.robust_nm} N·m (× 1.1) 자세 {len(iv)}개", flush=True)
        tic = time.time(); self.fr = pd.Frames(os.path.join(ROOT, a.model), Q); print(f"뼈 변환 저장 {time.time() - tic:.0f}s", flush=True)
        self.curves = td.load_curves(a.curves, beta=a.beta); self.grades = [g.strip() for g in a.grades.split(",") if g.strip()]
        br = pd.setup_bracket(a.bracket, a.box_scale, dict(ulna_r=tuple(a.ulna_scale)))
        self.tube = br["bones"][pd.BODIES.index("ulna_r")]; self.chain = [pd.BODIES.index("ulna_r"), pd.BODIES.index("radius_r")]
        self.cad = None
        if a.cad: self.setup_cad(a)

    def setup_cad(self, a):
        """--cad (2026-11-19, 연구실 팔 시제품): 가이드 A 높이 = CAD 자뼈 표면 기준 (tools/cad_sdf.py 거리장),
        드럼 = 회내외 축 둘레 원통이라 자뼈와의 간격은 자세와 무관 → 축 방향 위치마다 CAD 자뼈 표면의 축 거리 최솟값 (ρ_u) 과
        노뼈 부품의 축 거리 최댓값 (ρ_r) 을 미리 재 두고, 드럼 폭 안에서 R ≤ ρ_u − 간격, R ≥ ρ_r + 벽 을 벌점으로."""
        sys.path.insert(0, os.path.join(ROOT, "tools")); import cad_sdf as cs
        pd.BRACKET = dict(H=a.bracket, bones={}, clear={}); pd.BRACKET["cad"] = pd.setup_bracket_cad(a.cad); fj = json.load(open(os.path.join(ROOT, a.cad), encoding="utf-8"))
        Ar, tr, E = np.array(fj["A"]["radius_r"]), np.array(fj["t"]["radius_r"]), np.array(fj["E_C"])
        def prof(parts, fn):
            pts = np.vstack([cs.surface_samples(cs.read_stl_mm(os.path.join(ROOT, cs.STL, q + ".stl")), spacing=0.5) for q in parts])
            pr = tr + (pts - E) @ Ar.T                                                     # CAD → 노뼈 좌표 (회내외 0 = CAD 조립 자세)
            z = (pr - self.geo.O) @ self.geo.a; rho = np.linalg.norm(pr - self.geo.O - np.outer(z, self.geo.a), axis=1)
            edges = np.arange(-0.30, 0.03, 0.0005); k = np.digitize(z, edges); out = np.full(len(edges) + 1, np.nan)
            for i in np.unique(k): out[i] = fn(rho[k == i])
            return edges, out
        self.cad = dict(edges=prof(fj["parts"]["ulna"], np.min)[0], rho_u=prof(fj["parts"]["ulna"], np.min)[1], rho_r=prof(["Elbow__03_Radius_Proximal"], np.max)[1],
                        clear=a.drum_clear, wall=a.drum_wall, half=a.groove + a.drum_margin)
        rp = tr + (cs.surface_samples(cs.read_stl_mm(os.path.join(ROOT, cs.STL, "Elbow__03_Radius_Proximal.stl")), spacing=0.5) - E) @ Ar.T   # 노뼈 돌기 각도 보고용
        zr = (rp - self.geo.O) @ self.geo.a; v = rp - self.geo.O - np.outer(zr, self.geo.a)
        self.cad["rad_pts"] = (zr, np.linalg.norm(v, axis=1), np.degrees(np.arctan2(v @ self.geo.y0, v @ self.geo.x0)))
        z0 = np.arange(a.zc_min, a.zc_max + 1e-9, 0.005)
        print("CAD 드럼 자리 (축 방향 위치 → 자뼈까지 / 노뼈 바깥까지 축 거리, mm): " + ", ".join(f"{z*1000:+.0f}: {self.drum_room(z)[0]*1000:.1f}/{self.drum_room(z)[1]*1000:.1f}" for z in z0), flush=True)

    def drum_room(self, zc):
        """드럼 폭 [zc − half, zc + half] 안의 (자뼈 표면 축 거리 최소, 노뼈 부품 축 거리 최대) m."""
        c = self.cad; k0, k1 = np.digitize([zc - c["half"], zc + c["half"]], c["edges"]); ru = c["rho_u"][k0:k1 + 1]; rr = c["rho_r"][k0:k1 + 1]
        return (float(np.nanmin(ru)) if np.isfinite(ru).any() else 1.0), (float(np.nanmax(rr)) if np.isfinite(rr).any() else 0.0)

    def add_frames(self, Qn, taun):
        """Moco 에서 보조가 남은 순간 (자세, 필요 토크) 를 더한다 (바깥 반복). 토크는 이미 × 1.1 한 값."""
        fr2 = pd.Frames(os.path.join(ROOT, self.a.model), Qn); n0 = len(self.Q)
        self.fr.R = np.concatenate([self.fr.R, fr2.R]); self.fr.p = np.concatenate([self.fr.p, fr2.p]); self.Q = np.vstack([self.Q, Qn])
        self.idx = np.r_[self.idx, n0 + np.arange(len(Qn))]; self.tau = np.r_[self.tau, taun]; self.kind = np.r_[self.kind, np.full(len(Qn), 2)]

    def tube_h(self, p_u):
        if self.cad: return float(pd.cad_height(pd.BODIES.index("ulna_r"), p_u))
        p0, p1, Rt = self.tube; ax = p1 - p0; t = np.clip((p_u - p0) @ ax / (ax @ ax), 0.0, 1.0); return float(np.linalg.norm(p_u - (p0 + t * ax)) - Rt)

    def eval(self, x, detail=False, cabs=None):
        a = self.a; cabs = cabs or cables(x, self.geo, a); work, pen, det = 0.0, 0.0, []
        th_lo, th_hi = np.radians(a.theta_min), np.radians(a.theta_max)
        for (nm, role, s, sg), cab in zip(CAB, cabs):
            P = np.array(cab["points"]); L, inside, th = pd.path_lengths(self.fr, self.chain, P, with_inside=True, wrap_cyl=cab["wrap_cyl"], with_theta=True)
            r = -(L[:, 1 + 2 * J_PS] - L[:, 2 + 2 * J_PS]) / (2 * pd.DELTA); r_el = -(L[:, 1] - L[:, 2]) / (2 * pd.DELTA); th0 = th[:, 0]
            pen += 1e3 * float(inside.mean()) + a.w_theta * float(np.mean(np.maximum(0.0, th_lo - th0) ** 2 + np.maximum(0.0, th0 - th_hi) ** 2)) * 100
            h = self.tube_h(P[0]); pen += pd.W_BRACKET_WALL * (max(0.0, -h) ** 2 + max(0.0, h - a.bracket) ** 2)
            arm = sg * r[self.idx]; need = np.maximum(sg * self.tau, 0.0); act = need > 1e-9
            pen += 1e4 * float(np.sum(np.maximum(0.0, 0.002 - arm[act]) ** 2))
            T = np.where(act, need / np.maximum(arm, 0.002), 0.0); sL = (L[:, 0].max() - L[:, 0])[self.idx]
            z = td.size_cable(T, sL, self.curves, self.grades); work += z["obj"]
            if detail: det.append(dict(name=nm, role=role, r=r, r_el=r_el, L=L[:, 0], theta=th0, h=h, T=T, s=sL, size=z, inside=bool(inside.any())))
        if self.cad:
            ru, rr = self.drum_room(float(x[1])); R = float(x[0])
            pen += pd.W_BRACKET_WALL * (max(0.0, R + self.cad["clear"] - ru) ** 2 + (max(0.0, rr + self.cad["wall"] - R) ** 2 if self.cad["wall"] is not None else 0.0))
        if detail: return dict(total=work + pen, work=work, pen=pen, cables=det, cabs=cabs)
        return work + pen

    def solve(self, x0=None, seed=0):
        a = self.a; b = bounds(a); tic = time.time(); gen = [0]
        def cb(intermediate_result):
            gen[0] += 1
            if gen[0] % 25 == 0: print(f"  DE 세대 {gen[0]}: 최선 {intermediate_result.fun:.4f}  {time.time() - tic:.0f}s", flush=True)
        res = differential_evolution(self.eval, b, seed=seed, maxiter=a.de_iter if x0 is None else max(50, a.de_iter // 3), popsize=a.de_pop, tol=1e-8, polish=False, x0=x0, callback=cb)   # 반복 2회차부터는 앞 해에서 출발, 세대 1/3
        r2 = minimize(self.eval, res.x, method="Nelder-Mead", bounds=b, options=dict(maxfev=4000, xatol=1e-6, fatol=1e-7))
        x = r2.x if r2.fun < res.fun else res.x
        print(f"설계 쪽 최적화: DE {res.fun:.4f} → NM {r2.fun:.4f}  ({res.nfev + r2.nfev} 회, {time.time() - tic:.0f}s)", flush=True)
        return x

    def report(self, x):
        d = self.eval(x, detail=True); thmax = [float(c["theta"].max()) for c in d["cables"]]
        d["cabs"] = cables(x, self.geo, self.a, theta_max=thmax)                                   # 원통 좌표계 x = 감는 호 전체의 가운데 (quadrant +x)
        print(f"  R {x[0]*1000:.1f} mm, 드럼 중심 축 방향 {x[1]*1000:+.1f} mm (노뼈머리 기준)   목적 {d['total']:.3f} (일 용량 {d['work']:.3f} J, 벌점 {d['pen']:.4f})", flush=True)
        if self.cad:
            ru, rr = self.drum_room(float(x[1]))
            print(f"  CAD: 드럼 폭 ±{self.cad['half']*1000:.1f} mm 안에서 자뼈 표면까지 축 거리 {ru*1000:.1f} mm → 드럼과 틈 {(ru - x[0])*1000:+.1f} mm, "
                  f"노뼈 부품 바깥 {rr*1000:.1f} mm → 드럼 벽 {(x[0] - rr)*1000:+.1f} mm" + (" (벽 조건 끔: 드럼 = 노뼈에 고정하는 링, 돌기는 아래 각도)" if self.cad["wall"] is None else ""), flush=True)
            zr, rho, ph = self.cad["rad_pts"]; sel = (np.abs(zr - x[1]) <= self.cad["half"]) & (rho > x[0])
            if sel.any():
                hist = np.histogram(ph[sel], bins=np.arange(-180, 181, 10))[0] > 0; spans = [int(b) for b in np.arange(-180, 180, 10)[hist]]
                print(f"  노뼈 부품이 드럼 반지름 밖으로 나온 각도 (노뼈 x0 기준, 10° 칸): {spans} — 와이어가 감기는 호와 겹치면 안 된다 (검사 단계)", flush=True)
        dist_ax = float(np.linalg.norm(np.cross(self.geo.to_ulna(self.geo.rpt(x[1], 0, 0)) - self.tube[0], self.tube[1] - self.tube[0])) / np.linalg.norm(self.tube[1] - self.tube[0]))
        print(f"  드럼 높이에서 회내외 축 ↔ 자뼈 관 축 거리 {dist_ax*1000:.1f} mm (자뼈 관 반지름 {self.tube[2]*1000:.0f} + 드럼 {x[0]*1000:.1f} mm → 틈 {(dist_ax - self.tube[2] - x[0])*1000:+.1f} mm)", flush=True)
        for c, (nm, role, s, sg) in zip(d["cables"], CAB):
            data = self.kind == 0; arm = sg * c["r"][self.idx]
            print(f"  {nm} {role}: 회내외 모멘트암 {c['r'].min()*1000:+.1f}~{c['r'].max()*1000:+.1f} mm "
                  f"(그 방향 유효 {arm.min()*1000:.1f}~{arm.max()*1000:.1f}), 팔꿈치 모멘트암 최대 |{np.abs(c['r_el']).max()*1000:.3f}| mm, 감김 각 {np.degrees(c['theta'].min()):.0f}~{np.degrees(c['theta'].max()):.0f}°, "
                  f"가이드 높이 {c['h']*1000:.1f} mm  | {td.fmt_size(c['size'])}  (데이터 프레임 장력 최대 {c['T'][data].max():.0f} N)", flush=True)
        return d


# ─────────────────────────────── OpenSim 길이 비교 ───────────────────────────────
def check_opensim(js, model_path, Q, out_csv):
    """OpenSim WrapCylinder 경로 길이·모멘트암 vs 설계 계산 (cylinder_path). 오른쪽은 같은 자세, 왼쪽은 같은 각 (거울상이면 같아야 한다)."""
    import opensim as osim
    osim.Logger.setLevelString("Error"); m = osim.Model(model_path); s = m.initSystem(); cs = m.getCoordinateSet(); fs = m.getForceSet()
    fr = pd.Frames(model_path, Q); chain = [pd.BODIES.index("ulna_r"), pd.BODIES.index("radius_r")]; rows = []; summ = {}
    for mm in js["muscles"]:
        if not mm.get("wrap_cyl"): continue
        L = pd.path_lengths(fr, chain, np.array(mm["points"]), wrap_cyl=mm["wrap_cyl"]); r_np = -(L[:, 1 + 2 * J_PS] - L[:, 2 + 2 * J_PS]) / (2 * pd.DELTA)
        for side in ("r", "l"):
            gp = osim.PathActuator.safeDownCast(fs.get(f"{mm['name']}_{side}")).getGeometryPath(); Lo, ro, re = [], [], []
            for q in Q:
                cs.get(f"elbow_flex_{side}").setValue(s, float(q[0]), False); cs.get(f"pro_sup_{side}").setValue(s, float(q[1]), False); m.realizePosition(s)
                Lo.append(gp.getLength(s)); ro.append(gp.computeMomentArm(s, cs.get(f"pro_sup_{side}"))); re.append(gp.computeMomentArm(s, cs.get(f"elbow_flex_{side}")))
            Lo, ro, re = np.array(Lo), np.array(ro), np.array(re); dL = (Lo - L[:, 0]) * 1000; dr = (ro - r_np) * 1000
            summ[f"{mm['name']}_{side}"] = dict(dL_max_mm=float(np.abs(dL).max()), dr_max_mm=float(np.abs(dr).max()), r_os_mm=(float(ro.min() * 1000), float(ro.max() * 1000)), r_el_max_mm=float(np.abs(re).max() * 1000))
            print(f"  OpenSim {mm['name']}_{side}: 길이 차 최대 {np.abs(dL).max():.3f} mm, 회내외 모멘트암 차 최대 {np.abs(dr).max():.3f} mm "
                  f"(OpenSim {ro.min()*1000:+.2f}~{ro.max()*1000:+.2f}, 설계 {r_np.min()*1000:+.2f}~{r_np.max()*1000:+.2f} mm), 팔꿈치 모멘트암 최대 |{np.abs(re).max()*1000:.3f}| mm", flush=True)
            for i, q in enumerate(Q):
                rows.append(dict(cable=mm["name"], side=side, elbow_deg=np.degrees(q[0]), prosup_deg=np.degrees(q[1]), L_osim_mm=Lo[i] * 1000, L_design_mm=L[i, 0] * 1000,
                                 r_osim_mm=ro[i] * 1000, r_design_mm=r_np[i] * 1000, r_elbow_osim_mm=re[i] * 1000))
    cs.get("elbow_flex_r").setValue(s, 0.0, False); cs.get("pro_sup_r").setValue(s, 0.0, False)
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    return summ, rows


def check_poses(design, n=40):
    """비교 자세: 회내외 격자 (0 ~ 최대) × 팔꿈치 3 + 데이터 자세 n 개."""
    ps = np.linspace(0.0, max(1.6, float(design.Q[:, J_PS].max())), 33); grid = np.array([[e, p] for e in (0.0, 1.2, 2.4) for p in ps])
    pick = design.Q[np.linspace(0, len(design.meta["trial"]) - 1, n).astype(int)]
    return np.vstack([grid, pick])


def plot(design, d, rows, png):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2)); q = np.degrees(design.Q[:, J_PS]); o = np.argsort(q)
    cols = ("#c0392b", "#2471a3")
    for c, col in zip(d["cables"], cols):
        ax[0].plot(q[o], c["r"][o] * 1000, color=col, lw=1.5, label=f"{c['name']} {c['role']} (설계 계산)")
        rr = [r for r in rows if r["cable"] == c["name"] and r["side"] == "r"]
        ax[0].plot([r["prosup_deg"] for r in rr], [r["r_osim_mm"] for r in rr], "o", ms=3, mfc="none", color=col, label=f"{c['name']} OpenSim (오른쪽)")
        rl = [r for r in rows if r["cable"] == c["name"] and r["side"] == "l"]
        ax[0].plot([r["prosup_deg"] for r in rl], [r["r_osim_mm"] for r in rl], "x", ms=3, color=col, label=f"{c['name']} OpenSim (왼쪽)")
        ax[1].plot(q[o], np.degrees(c["theta"][o]), color=col, label=c["name"])
        data = design.kind[: len(design.idx)] == 0
        ax[2].plot(np.degrees(design.Q[design.idx[data], J_PS]), c["T"][data], ".", ms=2, color=col, label=f"{c['name']} 데이터")
        rb = design.kind == 1; ax[2].plot(np.degrees(design.Q[design.idx[rb], J_PS]), c["T"][rb], "_", ms=4, color=col, alpha=0.5, label=f"{c['name']} 강건성")
    ax[0].axhline(0, color="k", lw=0.6); ax[0].set_xlabel("회내외 각 (°, + 엎침)"); ax[0].set_ylabel("회내외 모멘트암 (mm)"); ax[0].legend(fontsize=7); ax[0].set_title("모멘트암: 설계 계산 vs OpenSim WrapCylinder")
    ax[1].axhspan(design.a.theta_min, design.a.theta_max, color="0.9"); ax[1].set_xlabel("회내외 각 (°)"); ax[1].set_ylabel("드럼 감김 각 (°)"); ax[1].legend(fontsize=7); ax[1].set_title("감김 각 (회색 = 허용 범위)")
    ax[2].set_xlabel("회내외 각 (°)"); ax[2].set_ylabel("장력 (N)"); ax[2].legend(fontsize=7); ax[2].set_title("필요 장력 (1.1 × 토크 ÷ 모멘트암)")
    fig.tight_layout(); fig.savefig(png, dpi=140); plt.close(fig)


# ─────────────────────────────── Moco ───────────────────────────────
def write_json(design, d, elbow, path, extra):
    out = dict(space="forearm_r", coords=list(pd.COORDS), bodies=list(pd.BODIES), model=design.a.model, elbow_from=design.a.elbow, muscles=[])
    for m in elbow["muscles"]: out["muscles"].append(m)                                   # 팔꿈치 M0·M1 그대로
    for c, cab in zip(d["cables"], d["cabs"]):
        z = c["size"]; cab = dict(cab, force_N=float(z["T_peak"]), T_peak_N=float(z["T_peak"]), force_mode="time", s_max_m=float(z["s_max"]),
                                  L_min_m=float(c["L"].min()), L_max_m=float(c["L"].max()), arm_min=float(np.abs(c["r"]).min()), arm_mean=float(np.abs(c["r"]).mean()),
                                  theta_deg=[float(np.degrees(c["theta"].min())), float(np.degrees(c["theta"].max()))], guide_h_m=c["h"],
                                  module=dict(name=z["grade"], bundles=z["n"], n_cont=z["n_cont"], work_J=z["work_J"], mass_park_g=z["mass_park_g"], feasible=z["feasible"]))
        out["muscles"].append(cab)
    out.update(extra); json.dump(out, open(os.path.join(ROOT, path), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    return out


def run_moco(js_rel, model_rel, load_rel, tasks, prefix, outdir, official=False):
    """run_moco.sh 를 과제마다 병렬로. official 이면 results_paths.csv 에 줄을 더하고 report/moco 에 둔다 (최종 확인). 반환 {과제: 지표 dict}."""
    bash = shutil.which("bash"); env = dict(os.environ, PYTHONUTF8="1", MKL_THREADING_LAYER="SEQUENTIAL", LOGDIR=os.path.relpath(outdir, ROOT).replace("\\", "/"), **MOCO_ENV)
    procs = []; tic = time.time()
    for t in tasks:
        nm = f"{prefix}_{t}"; mj = os.path.join(outdir, f"{nm}_metrics.json")
        extra = (["--prefilter"] if t.startswith("walk") else []) + ["--metrics-json", os.path.relpath(mj, ROOT).replace("\\", "/")]
        if not official: extra += ["--out-dir", os.path.relpath(outdir, ROOT).replace("\\", "/"), "--no-csv"]
        cmd = [bash, "report/modulefirst/run_moco.sh", "forearm_r", js_rel, load_rel if t.startswith("lift") else model_rel, t, nm] + extra
        procs.append((t, nm, mj, subprocess.Popen(cmd, cwd=ROOT, env=env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8", errors="replace")))
    res = {}
    for t, nm, mj, p in procs:
        txt, _ = p.communicate()
        if not os.path.exists(mj): print(f"  Moco {nm} 실패 (지표 없음): {txt.strip()[-300:]}", flush=True); res[t] = None; continue
        res[t] = json.load(open(mj, encoding="utf-8")); cr = res[t]["coords"]
        print(f"  Moco {nm}: 성공={res[t]['success']}  회내외 보조 95퍼센타일 오른쪽 {cr['pro_sup_r']['assist_p95_pct']:.1f} % / 왼쪽 {cr['pro_sup_l']['assist_p95_pct']:.1f} % "
              f"(피크 {cr['pro_sup_r']['assist_peak_pct']:.1f} / {cr['pro_sup_l']['assist_peak_pct']:.1f} %), 팔꿈치 95퍼센타일 {cr['elbow_flex_r']['assist_p95_pct']:.1f} / {cr['elbow_flex_l']['assist_p95_pct']:.1f} %, "
              f"케이블 토크 vs ID 95퍼센타일 팔꿈치 {cr['elbow_flex_r']['err_p95_pct']:.1f} %·회내외 {cr['pro_sup_r']['err_p95_pct']:.1f} %  | 전체 최악 95퍼센타일 {res[t]['worst_p95_pct']:.1f} %", flush=True)
    print(f"  Moco {len(tasks)} 건 {time.time() - tic:.0f}s", flush=True)
    return res


def moco_score(res):
    """판정: 모든 좌표·양쪽 보조 95퍼센타일 < 10 %, 설계 쪽 케이블 토크 vs ID 95퍼센타일 < 100 %. 목적값 = 회내외 오른쪽 보조 95퍼센타일 최악 (%)."""
    if any(v is None for v in res.values()): return 1e3, False
    ps = max(v["coords"]["pro_sup_r"]["assist_p95_pct"] for v in res.values())
    ok = all(v["worst_p95_pct"] < 10 and v["err_p95_pct"] < 100 for v in res.values())
    return ps, ok


def moco_frames(res, thr_frac):
    """바깥 반복: 회내외 보조 |reserve| > thr_frac × ID 피크 인 순간 → (자세, 1.1 × (케이블 토크 + 보조))."""
    Qn, tn = [], []
    for t, v in res.items():
        if v is None: continue
        sr = v["series"]["pro_sup_r"]; rsv = np.array(sr["reserve_Nm"]); cab = np.array(sr["cable_Nm"]); Q = np.array(v["Q"]); pk = v["coords"]["pro_sup_r"]["id_peak_Nm"]
        sel = np.abs(rsv) > thr_frac * pk
        if sel.any(): Qn.append(Q[sel]); tn.append(pd.MARGIN * (cab[sel] + rsv[sel])); print(f"  {t}: 보조가 남은 순간 {int(sel.sum())}개 → 설계 프레임 추가 (필요 토크 {np.abs(cab[sel] + rsv[sel]).max():.2f} N·m 까지)", flush=True)
    return (np.vstack(Qn), np.concatenate(tn)) if Qn else (None, None)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=MODEL); ap.add_argument("--tasks", default=TASKS); ap.add_argument("--elbow", default=ELBOW, help="팔꿈치 케이블 json (그대로 쓴다)")
    ap.add_argument("--out-tag", default="_prosup_wrap"); ap.add_argument("--dt", type=float, default=0.02)
    ap.add_argument("--robust-nm", type=float, default=0.25, help="강건성: 모든 자세(5 프레임마다)에서 회내·회외 양쪽으로 이 토크 N·m (× 1.1 여유, path_design 과 같은 꼴)")
    ap.add_argument("--robust-every", type=int, default=5)
    ap.add_argument("--grades", default="M3,M4,M5"); ap.add_argument("--curves", default=None); ap.add_argument("--beta", type=float, default=None)
    ap.add_argument("--r-min", type=float, default=0.006); ap.add_argument("--r-max", type=float, default=0.015, help="드럼 반지름 범위 m")
    ap.add_argument("--zc-min", type=float, default=-0.08); ap.add_argument("--zc-max", type=float, default=-0.012, help="드럼 중심 축 방향 범위 m (노뼈머리 기준, − 먼쪽)")
    ap.add_argument("--dz-max", type=float, default=0.02, help="가이드 A 의 축 방향 거리 (홈 기준) 상한 m"); ap.add_argument("--groove", type=float, default=0.004, help="두 홈의 축 방향 간격의 절반 m")
    ap.add_argument("--b-off", type=float, default=0.001, help="노뼈 끝점 B 를 드럼 표면에서 띄우는 거리 m"); ap.add_argument("--cyl-length", type=float, default=0.06, help="OpenSim 원통 길이 m")
    ap.add_argument("--theta-min", type=float, default=15.0); ap.add_argument("--theta-max", type=float, default=165.0, help="감김 각 허용 범위 ° (접촉 유지 / OpenSim 한쪽 감김)")
    ap.add_argument("--w-theta", type=float, default=10.0); ap.add_argument("--bracket", type=float, default=0.02, help="자뼈 가이드 A 높이 상한 m (자뼈 관 표면 기준)")
    ap.add_argument("--box-scale", type=float, default=1.0802); ap.add_argument("--ulna-scale", type=float, nargs=3, default=(1.0, 1.0495, 1.0), help="자뼈 관 축 추가 배율 (run_forearm_175.sh 의 --body-scale)")
    ap.add_argument("--cad", default=None, help="CAD 맞춤 json (tools/cad_arm_model.py): 가이드 높이·드럼 간격을 연구실 팔 시제품 부품 기준으로 (2026-11-19). 없으면 지금과 같다")
    ap.add_argument("--drum-clear", type=float, default=0.002, help="--cad: 드럼과 자뼈 표면 사이 틈 하한 m"); ap.add_argument("--drum-wall", type=float, default=None, help="--cad: 드럼이 노뼈 부품을 통째로 감싸게 할 벽 두께 하한 m (기본 끔 — CAD 에서는 자뼈와 겹쳐 불가능, 드럼은 노뼈에 고정하는 링)")
    ap.add_argument("--drum-margin", type=float, default=0.003, help="--cad: 드럼 폭 = 두 홈 ± 이 값 m")
    ap.add_argument("--moco-up", default=None, help="Moco 의 상체 ID 폴더 (기본 out/upper40_h175sk). CAD 팔 모델이면 out/upper40_h175sk_cadarm")
    ap.add_argument("--de-iter", type=int, default=300); ap.add_argument("--de-pop", type=int, default=20); ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--design-only", action="store_true", help="설계·OpenSim 비교·모델까지만 (Moco 안 함)")
    ap.add_argument("--outer", type=int, default=3, help="바깥 반복 최대 횟수 (Moco lift + walk 통과하면 멈춤)"); ap.add_argument("--outer-thr", type=float, default=0.02, help="보조 |reserve| > 이 비율 × ID 피크 인 순간을 추가")
    ap.add_argument("--moco-nm", type=int, default=0, help="바깥 반복 뒤 Nelder–Mead 평가 횟수 (평가 1회 = Moco 2 건, 약 1.5 분). 0 이면 안 함")
    ap.add_argument("--reuse", action="store_true", help="첫 설계를 다시 풀지 않고 저장된 json (같은 --out-tag) 의 prosup_wrap.x 에서 시작 (설계만 돌린 뒤 Moco 로 이어 갈 때)")
    ap.add_argument("--final", action="store_true", help="마지막 설계로 walk·walks9·lift·lifts6 공식 Moco (results_paths.csv 에 td175_forearm<tag>_* 줄)")
    a = ap.parse_args(); a.elbow_abs = os.path.join(ROOT, a.elbow)
    if a.moco_up: MOCO_ENV["UP"] = a.moco_up
    outdir = os.path.join(ROOT, f"report/td/prosup_wrap{a.out_tag}"); os.makedirs(outdir, exist_ok=True)
    js_rel = f"report/{STEM}{a.out_tag}.json"; mdl_rel = f"report/moco/{STEM}{a.out_tag}_model.osim"; load_rel = f"report/moco/{STEM}{a.out_tag}_load3.79_model.osim"
    x_saved = None
    if a.reuse:
        x_saved = np.array([json.load(open(os.path.join(ROOT, js_rel), encoding="utf-8"))["prosup_wrap"]["x"][k] for k in PN]); os.environ["PROSUP_OVERWRITE"] = "1"
    for p in (js_rel, mdl_rel, load_rel):
        if os.path.exists(os.path.join(ROOT, p)) and not os.environ.get("PROSUP_OVERWRITE"): raise SystemExit(f"{p} 가 이미 있다 — 다른 --out-tag 를 쓰거나 PROSUP_OVERWRITE=1")
    elbow = json.load(open(a.elbow_abs, encoding="utf-8")); assert all(m["chain"] == ["humerus_r", "ulna_r"] for m in elbow["muscles"]), "팔꿈치 json 은 위팔→자뼈 케이블만 (회내외 모멘트 0) 이어야 한다"
    geo = Geo(os.path.join(ROOT, a.model)); design = Design(a, geo); pd.configure("forearm_r")
    print(f"회내외 축 (노뼈 좌표) {np.round(geo.a, 4)}, 원점 {np.round(geo.O, 4)};  팔꿈치 케이블 {a.elbow} ({[m['module']['name'] + '×' + str(m['module']['bundles']) for m in elbow['muscles']]})", flush=True)
    hist = []; x = None; base = os.path.join(ROOT, a.model); res = None

    def build(x, it_tag):
        d = design.report(x); extra = dict(prosup_wrap=dict(x=dict(zip(PN, map(float, x))), robust_nm=a.robust_nm, grades=a.grades, theta_deg=[a.theta_min, a.theta_max], bracket_m=a.bracket, cad=a.cad,
                                                             n_moco_frames=int((design.kind == 2).sum()), iteration=it_tag, history=hist))
        out = write_json(design, d, elbow, js_rel, extra); os.environ["PROSUP_OVERWRITE"] = "1"
        pd.build_model(out, base, os.path.join(ROOT, mdl_rel)); pd.build_model(out, base.replace(".osim", "_load3.79.osim"), os.path.join(ROOT, load_rel))
        return d, out

    for it in range(max(1, a.outer if not a.design_only else 1)):
        print(f"\n── 설계 {it + 1}회 (Moco 추가 프레임 {int((design.kind == 2).sum())}개) ──", flush=True)
        x = x_saved if (it == 0 and x_saved is not None) else design.solve(x0=x, seed=a.seed + it); d, out = build(x, it + 1)
        if it == 0:
            summ, rows = check_opensim(out, os.path.join(ROOT, mdl_rel), check_poses(design), os.path.join(outdir, "osim_vs_design.csv"))
            plot(design, d, rows, os.path.join(outdir, "fig_prosup_wrap.png")); hist.append(dict(iteration=0, osim_check=summ))
            bad = [k for k, v in summ.items() if v["dL_max_mm"] > 0.5 or v["dr_max_mm"] > 0.5]
            if bad: print(f"  ※ OpenSim 과 설계 계산이 0.5 mm 넘게 다르다: {bad} — Moco 전에 확인 필요", flush=True)
        if a.design_only: break
        res = run_moco(js_rel, mdl_rel, load_rel, ["lift", "walk"], f"pw{a.out_tag}_it{it + 1}", outdir)
        sc, ok = moco_score(res); hist.append(dict(iteration=it + 1, x=list(map(float, x)), prosup_p95_pct=sc, passed=ok, work_J=d["work"],
                                                   per_task={t: (None if v is None else dict(worst_p95_pct=v["worst_p95_pct"], err_p95_pct=v["err_p95_pct"], prosup=v["coords"]["pro_sup_r"])) for t, v in res.items()}))
        print(f"  → 반복 {it + 1}: 회내외 보조 95퍼센타일 최악 {sc:.1f} %, 판정 {'통과' if ok else '미달'}", flush=True)
        if ok: break
        Qn, tn = moco_frames(res, a.outer_thr)
        if Qn is None: print("  보조가 남은 순간이 없다 — 회내외 밖(팔꿈치 등)의 문제", flush=True); break
        design.add_frames(Qn, tn)

    if a.moco_nm > 0 and not a.design_only:
        print(f"\n── Nelder–Mead (Moco 평가 {a.moco_nm}회) ──", flush=True); cache = {}; n = [0]; b = bounds(a); lo = np.array([v[0] for v in b]); hi = np.array([v[1] for v in b])
        def f(xx):
            xx = np.clip(xx, lo, hi); key = tuple(np.round(xx, 5))
            if key in cache: return cache[key]
            dd = design.eval(xx, detail=True)
            if dd["pen"] > 1e-3: cache[key] = 100.0 + dd["pen"]; return cache[key]                      # 설계 쪽에서 이미 안 되는 점 (감김 범위·브래킷) 은 Moco 를 안 돌린다
            n[0] += 1; dd2, out2 = build(xx, f"nm{n[0]}"); r = run_moco(js_rel, mdl_rel, load_rel, ["lift", "walk"], f"pw{a.out_tag}_nm{n[0]}", outdir)
            sc2, ok2 = moco_score(r); val = sc2 + 0.01 * dd["work"]; cache[key] = val
            hist.append(dict(iteration=f"nm{n[0]}", x=list(map(float, xx)), prosup_p95_pct=sc2, passed=ok2, work_J=dd["work"])); print(f"  NM {n[0]}: 목적 {val:.3f}", flush=True)
            return val
        step = np.array([0.001, 0.005] + [0.003, 0.2, 0.003, 0.2] * 2); sim = [x] + [np.clip(x + np.eye(len(x))[i] * step[i], lo, hi) for i in range(len(x))]
        r = minimize(f, x, method="Nelder-Mead", options=dict(initial_simplex=np.array(sim), maxfev=a.moco_nm, xatol=1e-4, fatol=0.1))
        best = min((h for h in hist if "x" in h and h.get("prosup_p95_pct") is not None), key=lambda h: (not h["passed"], h["prosup_p95_pct"] + 0.01 * h["work_J"]))
        x = np.array(best["x"]); print(f"  NM 끝: 가장 좋은 점 회내외 보조 95퍼센타일 {best['prosup_p95_pct']:.1f} % ({best['iteration']})", flush=True); build(x, f"nm_best")

    if a.final and not a.design_only:
        print("\n── 최종 Moco (walk·walks9·lift·lifts6, results_paths.csv) ──", flush=True)
        fin = run_moco(js_rel, mdl_rel, load_rel, ["walk", "walks9", "lift", "lifts6"], f"td175_forearm{a.out_tag}", outdir, official=True)
        sc, ok = moco_score(fin); hist.append(dict(iteration="final", prosup_p95_pct=sc, passed=ok,
                                                   per_task={t: (None if v is None else dict(worst_p95_pct=v["worst_p95_pct"], err_p95_pct=v["err_p95_pct"], coords=v["coords"])) for t, v in fin.items()}))
        print(f"  → 최종: 회내외 보조 95퍼센타일 최악 {sc:.1f} %, 판정 {'통과' if ok else '미달'}", flush=True)
    js = json.load(open(os.path.join(ROOT, js_rel), encoding="utf-8")); js["prosup_wrap"]["history"] = hist
    json.dump(js, open(os.path.join(ROOT, js_rel), "w", encoding="utf-8"), indent=1, ensure_ascii=False); print("저장", js_rel, mdl_rel, load_rel, os.path.relpath(outdir, ROOT))


if __name__ == "__main__":
    main()
