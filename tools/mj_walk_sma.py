"""160 cm 하체 설계의 MuJoCo 걷기 검증 — 인공근육(SMA) 열 특성 반영 (2026-11-21).

교수님 연구실 저장소(artificial_muscle_humanoid/mujoco_track)의 순동역학 걷기 모의실험(walk_sim.WalkSim: 기준 궤적 추종 PD + ID 피드포워드
→ 케이블 힘 배분, 골반 하네스)에서 **액추에이터 부분만** 바꿔 가며 비교한다.
  ideal : 지연 없음 (명령한 힘이 바로 나옴)
  lag   : 교수님 코드의 1차 열지연 (가열 0.12 s / 냉각 0.13 s)
  sma   : 김은서 님 SMA 열 모델 (sma_arm_kimm/sma_model.py — 온도·상변태·잠열·하중 문턱 이동, 힘 식 T3.1), 기계연 M3 실측에 맞춘 값
          (36 A = 267 W, 강제 냉각 4.62 W/K, 길이 0.18 m, 최대 수축 0.40, 힘 450 N)
  sma --cool X : 냉각 세기만 X 배 (수냉 등 가정)

등급·묶음 환산 (가정, 사용자 결정 2026-11-21 "비율로 늘리는 걸로"):
  등급 g 모듈 하나 = M3 모듈의 힘·열용량·가열·냉각을 (g 의 규격 최대 힘 ÷ 450 N) 배 → 가열·냉각 시간은 M3 와 같다.
  길이 = 0.18 m × (g 의 규격 최대 수축 ÷ M3 의 60 mm) → 최대 수축 40 % 가 규격 수축 × 1.2 (M3 는 72 mm = 기계연 실측).
  묶음 n 개 = 같은 모듈 n 개 병렬 (힘 n 배, 온도는 같다).
다는 방식: 케이블이 가장 길 때(설계 json 의 L_max) 근육 변형률 = 최대 수축(0.40). 관절이 움직여 케이블이 짧아지면 변형률이 준다 (김은서 kimm_arm_new.py 와 같음).
가열 제어 (단순 비례): 듀티 = clip(KU × (원하는 힘 − 지금 힘) ÷ 모듈 최대 힘, 0, 1). 듀티 > 0 이면 '켜는 동안 냉각'(0.72 W/K × 배율), 0 이면 강제 냉각.

    MJ=/c/Users/USER/miniforge3/envs/mj/python.exe
    $MJ tools/mj_walk_sma.py --mode ideal ; $MJ tools/mj_walk_sma.py --mode lag ; $MJ tools/mj_walk_sma.py --mode sma [--cool 3]
출력: report/mujoco160/walk_<tag>.json / .npz
"""
import argparse, json, os, subprocess, sys, types
os.environ.setdefault("MKL_THREADING_LAYER", "SEQUENTIAL")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.environ.setdefault("AMH_DATA", ROOT)
LAB = os.path.join(ROOT, "..", "artificial_muscle_humanoid", "mujoco_track", "src"); sys.path.insert(0, LAB); sys.path.insert(0, ROOT)
import numpy as np                                          # noqa: E402
import mujoco                                               # noqa: E402
from common import Leg16                                    # noqa: E402
import walk_sim as ws                                       # noqa: E402
from actuator_model import BY_NAME                          # noqa: E402


def load_sma_model():
    """김은서 님 sma_model.py: 작업 폴더에 있으면 그것, 없으면 git 의 origin/main 에서 읽는다 (복사본을 만들지 않는다)."""
    p = os.path.join(ROOT, "sma_arm_kimm", "sma_model.py")
    src = open(p, encoding="utf-8").read() if os.path.exists(p) else subprocess.run(["git", "-C", ROOT, "show", "origin/main:sma_arm_kimm/sma_model.py"], capture_output=True, check=True).stdout.decode("utf-8")
    mod = types.ModuleType("sma_model"); exec(compile(src, "sma_model.py", "exec"), mod.__dict__); return mod


class SMABank:
    """케이블 16개의 SMA 열 모델. walk_sim.SMA 와 같은 자리(.a, .step, .invert)에 끼운다.
    .a = 케이블 힘 ÷ 설계 힘(MuJoCo 액추에이터 gain) — 모듈 용량이 설계 힘보다 크면 1 을 넘을 수 있다."""

    def __init__(self, L, P, sm, cool=1.0, ku=8.0, speed=1.0):
        self.L, self.m, self.d = L, L.m, L.d; self.sm = sm; self.ku = ku; self.cool = cool; n = self.m.nu
        self.gain = L.Fmax.copy(); self.mus, self.nb, self.L0, self.Lmax, self.grade = [], np.zeros(n), np.zeros(n), np.zeros(n), []
        for i in range(n):
            mu = P["muscles"][i // 2]; g, k = mu["module"]["name"], int(mu["module"]["bundles"]); G = BY_NAME[g]
            r = G.force_n[1] / sm.M3["FMAX"]; p = dict(sm.M3)
            p.update(name=f"{g} (M3 × {r:.2f})", m=sm.M3["m"] * r, hS=sm.M3["hS"] * r * cool * speed, Qon=sm.M3["Qon"] * r * speed, FMAX=G.force_n[1], L0=sm.M3["L0"] * G.stroke_mm[1] / BY_NAME["M3"].stroke_mm[1])
            p["A_wire"] = p["FMAX"] / 400e6; p["K_mod"] = p["FMAX"] / p["epsL"]
            self.mus.append(sm.SMA(p)); self.nb[i] = k; self.L0[i] = p["L0"]; self.Lmax[i] = mu["L_max_m"]; self.grade.append(f"{g}×{k}")
        self.epsL = sm.M3["epsL"]; self.h_on = sm.KIMM_HEAT["h_on"] * speed;   # speed: 가열·냉각을 모두 이 배수로 (온도가 변하는 속도만 빨라진다)
        self.cap = np.array([mu.p["FMAX"] for mu in self.mus]) * self.nb
        self.tid = self.m.actuator_trnid[:, 0].copy(); self.a = np.zeros(n); self.F = np.zeros(n); self.u = np.zeros(n); self.f_des = np.zeros(n)
        self.m.actuator_ctrlrange[:, 1] = np.maximum(1.0, self.cap / self.gain)          # 모듈이 설계 힘보다 세게 당길 수 있게
        self.T = np.array([mu.T for mu in self.mus]); self.eps_now = np.full(n, self.epsL)

    def eps(self):
        return np.clip(self.epsL - (self.Lmax - self.d.ten_length[self.tid]) / self.L0, 0.0, None)

    def f_hot(self):
        """지금 자세에서 다 데웠을 때 낼 수 있는 힘 (T3.1: F = K_A · ε) — 배분 상한으로 쓴다."""
        return self.cap * np.minimum(self.eps() / self.epsL, 1.0)

    def invert(self, a_des, dt_c, gamma=None):            # 원하는 힘을 기억만 한다 (가열 듀티는 step 에서 힘 오차로 정함)
        self.f_des = a_des * self.gain; return a_des

    def step(self, u_cmd, dt):
        e = self.eps(); self.eps_now = e
        for i, mu in enumerate(self.mus):
            f1 = self.F[i] / self.nb[i]; fd = self.f_des[i] / self.nb[i]
            u = min(max(self.ku * (fd - f1) / mu.p["FMAX"], 0.0), 1.0) * mu.u_cap()
            h = self.h_on * (mu.p["FMAX"] / self.sm.M3["FMAX"]) if u > 0 else mu.p["hS"]      # 데우는 근육에는 바람을 안 보낸다 (기계연 블레이드) — 냉각 배율은 끌 때만
            self.F[i] = self.nb[i] * mu.step(mu.power_duty(u), dt, eps=e[i], hS=h); self.u[i] = u; self.T[i] = mu.T
        self.a = self.F / self.gain
        return self.a


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--mode", choices=("ideal", "lag", "sma"), default="sma"); ap.add_argument("--cool", type=float, default=1.0)
    ap.add_argument("--cycles", type=float, default=1.0); ap.add_argument("--harness", default="yz", choices=("yz", "full", "none")); ap.add_argument("--ku", type=float, default=8.0)
    ap.add_argument("--dir", default="report/mujoco160"); ap.add_argument("--xml", default="leg16_h160_dyn.xml"); ap.add_argument("--osim", default="data/robot40_h160/robot40_h160.osim")
    ap.add_argument("--paths", default="report/newmus_paths_lower_k8_40kg_h160_td_br20.json"); ap.add_argument("--data", default="data/robot40_h160/lower"); ap.add_argument("--trial", default="walking1")
    ap.add_argument("--tag", default=""); ap.add_argument("--warm", type=float, default=3.0, help="sma: 정해진 동작으로 미리 도는 주기 수 (마지막 주기로 통계)"); ap.add_argument("--kin-only", action="store_true"); ap.add_argument("--dump-table", action="store_true"); ap.add_argument("--speed", type=float, default=1.0, help="sma: 열 속도 배율 (가열·냉각 모두). 모델이 실제보다 느릴 가능성을 보는 민감도"); ap.add_argument("--slow", type=float, default=1.0, help="sma 정해진 동작: 같은 동작·같은 토크를 이 배수만큼 천천히 (근육 속도 한계를 보려는 것 — 느린 걷기의 실제 토크는 조금 다르다)"); ap.add_argument("--lead", type=float, default=0.0, help="sma 정해진 동작: 원하는 힘을 이 초만큼 앞당겨 가열 (동작을 미리 아는 제어)"); a = ap.parse_args()
    O = os.path.join(ROOT, a.dir); pj = os.path.join(ROOT, a.paths); P = json.load(open(pj, encoding="utf-8"))
    L = Leg16(os.path.join(O, a.xml), os.path.join(ROOT, a.osim), pj)
    ref = ws.Reference(L, os.path.join(ROOT, a.data), a.trial, id_npz=os.path.join(O, "full", f"id_{a.trial}.npz"), n_cycles=a.cycles)
    if a.dump_table:                                       # 걷기 한 주기의 자세별 케이블 길이·모멘트암·ID 토크 표 → tools/walk_kin_fast.py (env fit) 가 근육만 빠르게 계산
        dtk = 0.001; nc = int(round(ref.T / dtk)); tid = L.m.actuator_trnid[:, 0]; TL, RR, TD = np.zeros((nc, L.m.nu)), np.zeros((nc, 2, 4, 8)), np.zeros((nc, 8))
        for k in range(nc):
            L.set_pose(ref.coords_at(k * dtk)); TL[k] = L.d.ten_length[tid]
            for si, s_ in enumerate("rl"): RR[k, si] = L.moment_arms(s_); TD[k, 4 * si:4 * si + 4] = ref.ff_at(k * dtk, s_)
        np.savez(os.path.join(O, "walk_table.npz"), TL=TL, RR=RR, TD=TD, dt=dtk, T=ref.T, Lmax=np.array([P["muscles"][i // 2]["L_max_m"] for i in range(L.m.nu)]), gain=L.Fmax,
                 grade=np.array([P["muscles"][i // 2]["module"]["name"] for i in range(L.m.nu)]), bundles=np.array([P["muscles"][i // 2]["module"]["bundles"] for i in range(L.m.nu)]),
                 act_r=np.array(L.act_idx["r"]), act_l=np.array(L.act_idx["l"]), paths=a.paths)
        print("표 →", os.path.join(a.dir, "walk_table.npz"), TL.shape); return
    sim = ws.WalkSim(L, ref, harness=(False if a.harness == "none" else (True if a.harness == "full" else "yz")), ff=True)
    bank = None
    if a.mode == "ideal": sim.sma = ws.SMA(L.m.nu, 1e-3, 1e-3); sim.gamma = 0.0
    elif a.mode == "sma":
        bank = SMABank(L, P, load_sma_model(), cool=a.cool, ku=a.ku, speed=a.speed); sim.sma = bank
        print("모듈:", " ".join(bank.grade[::2]), "| 길이 mm:", np.round(bank.L0[::2] * 1000).astype(int).tolist(), "| 용량 N:", np.round(bank.cap[::2]).astype(int).tolist(), "| 설계 힘 N:", np.round(bank.gain[::2]).astype(int).tolist(), flush=True)
    kin = None
    if bank is not None:
        # ① 정해진 동작(기준 걷기)을 그대로 따라가며 근육만 계산 — 균형과 무관하게 "근육이 걷기에 필요한 힘을 제때 내는가".
        #    매 순간 ID 토크를 케이블에 배분(상한 = 지금 자세에서 다 데웠을 때의 힘)하고 가열 제어로 따라간다. 마지막 주기로 통계. 끝 상태가 ②의 시작 온도.
        orig = ws.allocate; dtk = 0.001; nc = int(round(ref.T / dtk)); K = {k: [] for k in ("t", "tau_des", "tau_act", "tau_cap", "Tmp", "u", "F", "fdes", "eps")}
        # 한 주기 표를 먼저 만든다: 자세별 케이블 길이·모멘트암·원하는 힘 (동작이 정해져 있으므로 가열을 lead 초 앞당겨 걸 수 있다 — --lead)
        TL, RR, FD, TD = np.zeros((nc, L.m.nu)), np.zeros((nc, 2, 4, 8)), np.zeros((nc, L.m.nu)), np.zeros((nc, 8))
        for k in range(nc):
            L.set_pose(ref.coords_at(k * dtk)); TL[k] = L.d.ten_length[bank.tid]; fh = bank.f_hot()
            for si, s in enumerate("rl"):
                idx = L.act_idx[s]; R = L.moment_arms(s); tau = ref.ff_at(k * dtk, s); RR[k, si] = R; TD[k, 4 * si:4 * si + 4] = tau; FD[k, idx] = orig(R, tau, np.maximum(fh[idx], 1e-3))
        nl = int(round(a.lead / (dtk * a.slow))); eps_tab = np.clip(bank.epsL - (bank.Lmax - TL) / bank.L0, 0.0, None); bank.eps = lambda: eps_tab[bank._k]
        for k in range(int(round(a.warm * nc))):
            kc = k % nc; bank._k = kc; bank.f_des = FD[(kc + nl) % nc].copy(); bank.step(None, dtk * a.slow)
            ta = np.concatenate([RR[kc, si] @ bank.F[L.act_idx[s]] for si, s in enumerate("rl")]); tc = np.concatenate([RR[kc, si] @ FD[kc, L.act_idx[s]] for si, s in enumerate("rl")])
            for key, v in (("t", k * dtk), ("tau_des", TD[kc]), ("tau_act", ta), ("tau_cap", tc), ("Tmp", bank.T.copy()), ("u", bank.u.copy()), ("F", bank.F.copy()), ("fdes", FD[kc]), ("eps", bank.eps_now.copy())): K[key].append(v)
        del bank.eps                                          # 순동역학에서는 다시 MuJoCo 의 케이블 길이로
        K = {k: np.array(v) for k, v in K.items()}; last = K["t"] >= (a.warm - 1) * ref.T; pk = np.abs(K["tau_des"][last]).max(axis=0)
        e_act = K["tau_act"][last] - K["tau_des"][last]; e_cap = K["tau_cap"][last] - K["tau_des"][last]
        kin = dict(warm_cycles=a.warm, speed=a.speed, lead_s=a.lead, slow=a.slow, period_s=float(ref.T * a.slow), id_peak_Nm=pk.round(1).tolist(),
                   torque_err_pct_of_peak_rms=(100 * np.sqrt(np.mean(e_act ** 2, axis=0)) / pk).round(1).tolist(), torque_err_pct_of_peak_p95=(100 * np.percentile(np.abs(e_act), 95, axis=0) / pk).round(1).tolist(),
                   static_limit_err_pct_of_peak_p95=(100 * np.percentile(np.abs(e_cap), 95, axis=0) / pk).round(1).tolist(),
                   T_max_per_cable_r=K["Tmp"][last].max(axis=0)[::2].round(0).tolist(), T_min_per_cable_r=K["Tmp"][last].min(axis=0)[::2].round(0).tolist(),
                   duty_mean_per_cable_r=K["u"][last].mean(axis=0)[::2].round(2).tolist(), eps_min_per_cable_r=K["eps"][last].min(axis=0)[::2].round(3).tolist(),
                   heat_W_mean_total=float(sum(K["u"][last][:, i].mean() * bank.mus[i].p["Qon"] * bank.nb[i] for i in range(L.m.nu))))
        np.savez(os.path.join(O, f"kin_{a.tag or ('sma' + (f'_cool{a.cool:g}' if a.cool != 1 else '') + (f'_speed{a.speed:g}' if a.speed != 1 else '') + (f'_lead{a.lead:g}' if a.lead else '') + (f'_slow{a.slow:g}' if a.slow != 1 else ''))}.npz"), **K); print("① 정해진 동작:", json.dumps(kin, ensure_ascii=False), flush=True)
        if a.kin_only:
            json.dump(dict(mode="sma", cool=a.cool, modules=bank.grade[::2], kinematic=kin), open(os.path.join(O, f"kin_{a.tag or ('sma' + (f'_cool{a.cool:g}' if a.cool != 1 else '') + (f'_speed{a.speed:g}' if a.speed != 1 else '') + (f'_lead{a.lead:g}' if a.lead else '') + (f'_slow{a.slow:g}' if a.slow != 1 else ''))}.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False); return
        # ② 순동역학: 배분 상한 = 지금 자세에서 다 데웠을 때의 힘 (설계 힘이 아니라 근육이 실제 낼 수 있는 힘)
        def alloc(R, tau, Fm, reg=1e-4, w=ws.W_ROWS, _s=[0]):
            idx = L.act_idx["rl"[_s[0] % 2]]; _s[0] += 1
            return orig(R, tau, np.maximum(bank.f_hot()[idx], 1e-3), reg, w)
        ws.allocate = alloc
    sim.reset(0.0)
    if bank is not None: sim.sma.a = bank.F / bank.gain
    T_end = a.cycles * ref.T; log = {k: [] for k in ("t", "q", "qref", "a", "Fh", "pelvis", "grf", "tau_des", "tau_act", "Tmp", "u", "eps", "fdes", "F")}; fell = False
    while sim.t < T_end:
        u, a_des, tau_des, qref = sim.nominal_control(); Fh = sim.step(u, qref); d = sim.d
        tau_act = np.concatenate([L.moment_arms(s) @ (sim.sma.a[L.act_idx[s]] * L.Fmax[L.act_idx[s]]) for s in "rl"])
        log["t"].append(sim.t); log["q"].append(d.qpos[sim.jq].copy()); log["qref"].append(ref.leg_ref(sim.t)[1]); log["a"].append(np.array(sim.sma.a).copy()); log["Fh"].append(np.linalg.norm(Fh))
        log["pelvis"].append(d.xpos[sim.pelvis].copy()); log["grf"].append([sim.contact_force("r"), sim.contact_force("l")]); log["tau_des"].append(tau_des.copy()); log["tau_act"].append(tau_act)
        if bank is not None: log["Tmp"].append(bank.T.copy()); log["u"].append(bank.u.copy()); log["eps"].append(bank.eps_now.copy()); log["fdes"].append(bank.f_des.copy()); log["F"].append(bank.F.copy())
        if sim.fallen(): fell = True; break
    for k in log: log[k] = np.array(log[k])
    err = np.degrees(log["q"] - log["qref"]); W = L.mass * 9.81; te = log["tau_act"] - log["tau_des"]
    res = dict(mode=a.mode, cool=a.cool, harness=a.harness, cycles=a.cycles, period_s=float(ref.T), sim_time_s=float(log["t"][-1]), fell=bool(fell),
               rmse_deg=float(np.sqrt(np.mean(err ** 2))), rmse_deg_per_joint=np.sqrt(np.mean(err ** 2, axis=0)).round(2).tolist(),
               harness_BW_rms=float(np.sqrt(np.mean(log["Fh"] ** 2)) / W), harness_BW_peak=float(log["Fh"].max() / W), pelvis_drop_m=float(log["pelvis"][0, 1] - log["pelvis"][:, 1].min()),
               torque_rmse_Nm_per_joint=np.sqrt(np.mean(te ** 2, axis=0)).round(2).tolist(), torque_des_rms_Nm_per_joint=np.sqrt(np.mean(log["tau_des"] ** 2, axis=0)).round(2).tolist(),
               dist_m=float(log["pelvis"][-1, 0] - log["pelvis"][0, 0]), grf_peak_BW=(log["grf"].max(axis=0) / W).round(2).tolist())
    if bank is not None:
        res.update(kinematic=kin, T_max_C=float(log["Tmp"].max()), T_max_per_cable_r=log["Tmp"].max(axis=0)[::2].round(0).tolist(), duty_mean=float(log["u"].mean()),
                   force_short_frac=float(np.mean(log["F"] < 0.8 * log["fdes"] - 5.0)), modules=bank.grade[::2])
    tag = a.tag or (a.mode + (f"_cool{a.cool:g}" if a.mode == "sma" and a.cool != 1 else "") + (f"_speed{a.speed:g}" if a.mode == "sma" and a.speed != 1 else "") + (f"_lead{a.lead:g}" if a.mode == "sma" and a.lead else "") + (f"_c{a.cycles:g}" if a.cycles != 1 else "") + ("" if a.harness == "yz" else f"_h{a.harness}"))
    json.dump(res, open(os.path.join(O, f"walk_{tag}.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False); np.savez(os.path.join(O, f"walk_{tag}.npz"), **log)
    print(json.dumps(res, ensure_ascii=False))


if __name__ == "__main__":
    main()
