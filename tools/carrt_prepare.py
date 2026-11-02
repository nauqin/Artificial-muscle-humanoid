"""CARRT(Zenodo 8034000) .c3d → 마커 이름 매핑 → .trc → IK → (손 하중 모델) 까지. **`osim` 환경.**

    python tools/carrt_prepare.py --raw ../raw_carrt/MocapDataPaperData --subject 1 --trials PB1 PFWJ2 PEWJ2 OCH2 --out data/carrt
→ data/carrt/subject1/{subject1_scaled.osim, subject1_load<kg>.osim, <trial>_ik.mot, <trial>_markers.trc}, ik_log.csv

뼈대: Moore2015 subject3 모델(Rajagopal 팔 포함, 마커 45개)을 generic 으로 삼아 OpenSim ScaleTool 로 분절별
스케일(마커 쌍 길이 비) + 마커 재배치(정지 자세 IK 뒤 마커를 실측 위치로) 를 한다. 스케일 없이 쓰면 마커 오차
7 cm(골반 폭 29 vs 23 cm). 정지 자세는 `--static` trial 의 처음 0.3 s (서 있는 시작 자세). 마커는 CARRT 43개 중
35개를 모델 마커 이름으로 바꿔 쓰고, 모델에 없는 8개(CLAV, LBACK, RSHOP, LSHOP, RFIN, LFIN, RMFC, LMFC)는 버린다.
손 하중: 물체 질량을 hand_r (양손 과제는 hand_r·hand_l 에 반씩) 질량에 더한 모델을 따로 쓴다 — 외력 대신 질량.
좌표계: Vicon z-up → OpenSim y-up (x 축 둘레 −90° 회전).
"""
import argparse, csv, os, re, sys
import numpy as np
import opensim as osim

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
osim.Logger.setLevelString("Error")

#: CARRT 마커 → Moore2015(AddB Rajagopal 팔 포함) 모델 마커
MAP = {"RPSI": "RPSIS", "LPSI": "LPSIS", "LASI": "LASIS", "RASI": "RASIS", "LGT": "LGTRO", "RGT": "RGTRO",
       "C7": "C7", "T10": "T10", "STRN": "STRN", "RSHO": "RSHO", "LSHOA": "LSHO", "RUPA": "RDELT", "LUPA": "LDELT",
       "RELB": "RLEE", "LELB": "LLEE", "RELBM": "RMEE", "LELBM": "LMEE", "RFRA": "RFRM", "LFRA": "LFRM",
       "RWRA": "RLW", "RWRB": "RMW", "LWRA": "LLW", "LWRB": "LMW",
       "RTHI": "FRTHI", "LTHI": "FLTHI", "RLRC": "RLEK", "LLFC": "LLEK", "RTBI": "RATI", "LTBI": "LATI",
       "RLMAL": "RLM", "LLMAL": "LLM", "RCAL": "RHEE", "LCAL": "LHEE", "RTOE": "RTOE", "LTOE": "LTOE"}
#: IK 가중치. 골반·몸통·어깨·팔꿈치 1, 손목 0.5(라디얼/울나 쪽 구분이 불확실), 허벅지/정강이 클러스터 0.5
WEIGHT = {"RLW": 0.5, "RMW": 0.5, "LLW": 0.5, "LMW": 0.5, "FRTHI": 0.5, "FLTHI": 0.5, "RATI": 0.5, "LATI": 0.5, "RDELT": 0.5, "LDELT": 0.5, "RSHO": 0.5, "LSHO": 0.5}
#: 모델에서 몸통에 붙어 있지만 CARRT 에서는 상완 마커(RUPA/LUPA)인 것 — 상완으로 옮긴다 (안 옮기면 팔 들 때 22 cm 오차)
REPARENT = {"RDELT": "/bodyset/humerus_r", "LDELT": "/bodyset/humerus_l"}
DEMO = {1: (60, 1.626), 2: (65, 1.651), 3: (58, 1.600), 4: (86, 1.727), 5: (58, 1.524), 6: (68, 1.626), 7: (51, 1.651), 8: (79, 1.626), 9: (60, 1.727), 10: (76, 1.676)}
#: 과제 코드 → (설명, 손 하중 kg, 양손?). 하중은 논문 Table 3 공칭값 — 확인 전 None
TASKS = {"BH": ("빗질", 0.0, False), "FPC": ("컵으로 마시기", 0.0, False), "FDC": ("컵으로 마시기", 0.0, False), "OCL": ("아래 찬장 열기", 0.0, False), "OCH": ("위 찬장 열기", 0.0, False), "OC": ("찬장 열기", 0.0, False), "OCF": ("찬장 열기(F)", 0.0, False),
         "PB": ("상자 들기", 0.52, True), "PDC": ("먼지떨이 들고 청소", 0.18, False), "PEWJ": ("빈 물병 들기", 0.90, False),
         "PFWJ": ("찬 물병 들기", 3.79, False), "PWJ": ("물병 들고 따르기", 3.79, False)}
#: 하중은 논문(Sensors 2023, 23:8354) Table 3 공칭값: 상자 0.52, 먼지떨이 0.18, 빈 1갤런 물병 0.90, 찬 물병 3.79, 컵 0.02 kg.
#: 논문은 양손 과제가 하나라고만 적었고 어느 것인지 안 써서 상자를 양손으로 본다.


def code_of(tr):
    """'PFWJ! 2' → 'PFWJ', 'PB4' → 'PB' (피험자마다 파일 이름 꼴이 조금씩 다르다)."""
    return re.sub(r"[^A-Za-z].*$", "", tr)


def c3d_to_trc(c3d_path, trc_path):
    ad = osim.C3DFileAdapter(); ad.setLocationForForceExpression(osim.C3DFileAdapter.ForceLocation_CenterOfPressure)
    tables = ad.read(c3d_path); mk = ad.getMarkersTable(tables)
    labels = list(mk.getColumnLabels()); t = np.array(mk.getIndependentColumn())
    X = np.zeros((mk.getNumRows(), len(labels), 3))
    for i in range(mk.getNumRows()):
        row = mk.getRowAtIndex(i)
        for j in range(len(labels)): X[i, j] = [row[j].get(a) for a in range(3)]
    X = X / 1000.0                                                # mm → m
    # z-up 이면 y-up 으로: (x, y, z) → (x, z, −y)
    keep = [j for j, l in enumerate(labels) if l in MAP]
    if np.nanstd(X[:, keep, 2]) > np.nanstd(X[:, keep, 1]): X = np.stack([X[:, :, 0], X[:, :, 2], -X[:, :, 1]], axis=2); rotated = True
    else: rotated = False
    names = [MAP[labels[j]] for j in keep]
    rate = float(mk.getTableMetaDataAsString("DataRate"))
    with open(trc_path, "w") as fh:
        fh.write(f"PathFileType\t4\t(X/Y/Z)\t{os.path.basename(trc_path)}\n")
        fh.write("DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\tOrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n")
        fh.write(f"{rate:g}\t{rate:g}\t{len(t)}\t{len(names)}\tm\t{rate:g}\t1\t{len(t)}\n")
        fh.write("Frame#\tTime\t" + "\t".join(f"{n}\t\t" for n in names) + "\n")
        fh.write("\t\t" + "\t".join(f"X{k+1}\tY{k+1}\tZ{k+1}" for k in range(len(names))) + "\n\n")
        for i in range(len(t)):
            vals = []
            for j in keep:
                v = X[i, j]; vals.append("\t".join("" if np.isnan(a) else f"{a:.5f}" for a in v))
            fh.write(f"{i+1}\t{t[i]:.5f}\t" + "\t".join(vals) + "\n")
    return names, t, rotated, [l for l in labels if l not in MAP]


#: 분절별 스케일에 쓸 마커 쌍: (바디들, 축, [(마커, 마커), ...]). 여러 쌍이면 평균 비.
#: 어깨(RSHO)가 견봉이 아니라 삼각근 위라 상완 길이를 못 재므로(측정 20.6 cm) 상완은 **전완 비**를 대신 쓴다.
#: 골반은 폭(ASIS–ASIS, Z)과 깊이(ASIS–PSIS, X·Y)를 따로. 몸통은 C7–PSIS 로 균일.
MEASURE = [(("pelvis",), "Z", [("RASIS", "LASIS")]), (("pelvis",), "XY", [("RASIS", "RPSIS"), ("LASIS", "LPSIS")]),
           (("torso",), "XYZ", [("C7", "RPSIS"), ("C7", "LPSIS")]),
           (("humerus_r", "ulna_r", "radius_r", "hand_r"), "XYZ", [("RLEE", "RLW"), ("RMEE", "RMW")]),
           (("humerus_l", "ulna_l", "radius_l", "hand_l"), "XYZ", [("LLEE", "LLW"), ("LMEE", "LMW")]),
           (("femur_r",), "XYZ", [("RGTRO", "RLEK")]), (("femur_l",), "XYZ", [("LGTRO", "LLEK")]),
           (("tibia_r",), "XYZ", [("RLEK", "RLM")]), (("tibia_l",), "XYZ", [("LLEK", "LLM")]),
           (("talus_r", "calcn_r", "toes_r"), "XYZ", [("RHEE", "RTOE")]), (("talus_l", "calcn_l", "toes_l"), "XYZ", [("LHEE", "LTOE")])]


def place_markers(model_path, static_trc, t_static, out_path, log=print):
    """정지 자세에서 IK 를 풀고, 모델 마커를 실측 위치(부모 뼈 좌표계, 프레임 평균)로 옮긴다.
    ScaleTool 의 MarkerPlacer 가 파이썬에서 마커를 안 옮기는 것을 확인(2026-11-02)해 직접 한다."""
    model = osim.Model(model_path)
    for n, parent in REPARENT.items():
        if model.getMarkerSet().contains(n): model.getMarkerSet().get(n).setParentFrameName(parent)
    model.finalizeConnections(); state = model.initSystem()
    ik = osim.InverseKinematicsTool(); ik.setModel(model); ik.setMarkerDataFileName(static_trc)
    ik.setStartTime(t_static[0]); ik.setEndTime(t_static[1]); tmp = os.path.join(os.path.dirname(out_path), "_static_ik.mot"); ik.setOutputMotionFileName(tmp)
    md = osim.MarkerData(static_trc); names = [md.getMarkerNames().get(i) for i in range(md.getMarkerNames().getSize())]
    tasks = ik.getIKTaskSet()
    for n in names:
        if model.getMarkerSet().contains(n):
            tk = osim.IKMarkerTask(); tk.setName(n); tk.setWeight(WEIGHT.get(n, 1.0)); tk.setApply(True); tasks.cloneAndAppend(tk)
    ik.run()
    import addb_io as io
    q = io.read_mot(tmp); cs = model.getCoordinateSet(); acc = {n: [] for n in names}
    for fi in range(md.getNumFrames()):
        fr = md.getFrame(fi); tt = fr.getFrameTime()
        if tt < t_static[0] or tt > t_static[1]: continue
        row = q.at(tt)[0]
        for j, cn in enumerate(q.names):
            c = cs.get(cn); v = row[j]
            if q.in_degrees and c.getMotionType() == 1: v = np.deg2rad(v)
            if not c.get_locked(): c.setValue(state, float(v), False)
        model.realizePosition(state)
        for j, n in enumerate(names):
            if not model.getMarkerSet().contains(n): continue
            g = fr.getMarker(j)
            if any(np.isnan(g.get(k)) for k in range(3)): continue
            mk = model.getMarkerSet().get(n); body = mk.getParentFrame()
            loc = model.getGround().findStationLocationInAnotherFrame(state, g, body); acc[n].append([loc.get(k) for k in range(3)])
    moved = 0
    for n, L in acc.items():
        if L:
            model.getMarkerSet().get(n).set_location(osim.Vec3(*np.mean(L, axis=0))); moved += 1
    model.finalizeConnections(); model.printToXML(out_path); log(f"  마커 재배치 {moved}개 → {out_path}")
    return out_path


def scale_model(base_model, static_trc, out_dir, subj, mass_kg, t_static=(0.0, 0.3), log=print):
    """ScaleTool: 측정 기반 분절 스케일 → 정지 자세 마커 재배치. 결과 <subj>_scaled.osim."""
    st = osim.ScaleTool(); st.setName(subj); st.setSubjectMass(mass_kg)
    st.getGenericModelMaker().setModelFileName(os.path.abspath(base_model))
    ms = st.getModelScaler(); ms.setApply(True)
    order = osim.ArrayStr(); order.append("measurements"); ms.setScalingOrder(order)
    for bodies, axes, pairs in MEASURE:
        meas = osim.Measurement(); meas.setName("+".join(bodies) + "_" + axes); meas.setApply(True)
        for a, b in pairs:
            mp = osim.MarkerPair(); mp.setMarkerName(0, a); mp.setMarkerName(1, b); meas.getMarkerPairSet().cloneAndAppend(mp)
        for bn in bodies:
            bs = osim.BodyScale(); bs.setName(bn); ax = osim.ArrayStr(); [ax.append(x) for x in axes]; bs.setAxisNames(ax); meas.getBodyScaleSet().cloneAndAppend(bs)
        ms.getMeasurementSet().cloneAndAppend(meas)
    ms.setMarkerFileName(os.path.abspath(static_trc)); tr = osim.ArrayDouble(); tr.append(t_static[0]); tr.append(t_static[1]); ms.setTimeRange(tr)
    ms.setPreserveMassDist(True); ms.setOutputModelFileName(os.path.abspath(os.path.join(out_dir, f"{subj}_scaled_nomarker.osim")))
    pl = st.getMarkerPlacer(); pl.setApply(True); pl.setStaticPoseFileName(os.path.abspath(static_trc)); pl.setTimeRange(tr)
    tasks = osim.IKTaskSet(); md = osim.MarkerData(static_trc); names = {md.getMarkerNames().get(i) for i in range(md.getMarkerNames().getSize())}
    for n in sorted(names):
        tk = osim.IKMarkerTask(); tk.setName(n); tk.setWeight(WEIGHT.get(n, 1.0)); tk.setApply(True); tasks.cloneAndAppend(tk)
    ts = pl.getIKTaskSet()
    for i in range(tasks.getSize()): ts.cloneAndAppend(tasks.get(i))
    out = os.path.abspath(os.path.join(out_dir, f"{subj}_scaled.osim")); pl.setOutputModelFileName(out)
    st.printToXML(os.path.join(out_dir, f"{subj}_scale_setup.xml"))
    if not st.run(): raise RuntimeError("ScaleTool 실패")
    nomk = os.path.join(out_dir, f"{subj}_scaled_nomarker.osim")
    place_markers(nomk, static_trc, t_static, out, log=log)
    m = osim.Model(out); m.initSystem(); tot = sum(m.getBodySet().get(i).getMass() for i in range(m.getBodySet().getSize()))
    log(f"  스케일 완료: 총질량 {tot:.2f} kg (목표 {mass_kg}), 모델 {out}")
    return out


def with_hand_load(model_path, out_path, load_kg, bimanual):
    """물체 질량을 손 바디에 더한다 (양손이면 반씩). 관성은 점질량으로 손 무게중심에."""
    m = osim.Model(model_path); m.initSystem()
    for hand, share in (("hand_r", 0.5 if bimanual else 1.0), ("hand_l", 0.5 if bimanual else 0.0)):
        if share == 0: continue
        b = m.getBodySet().get(hand); b.setMass(b.getMass() + load_kg * share)
    m.printToXML(out_path); return out_path


def run_ik(model_path, trc_path, out_mot, log):
    ik = osim.InverseKinematicsTool(); ik.setName(os.path.basename(out_mot)[:-4])
    model = osim.Model(model_path); model.initSystem(); ik.setModel(model)
    ik.setMarkerDataFileName(trc_path); ik.setOutputMotionFileName(out_mot)
    md = osim.MarkerData(trc_path); ik.setStartTime(md.getStartFrameTime()); ik.setEndTime(md.getLastFrameTime())
    tasks = osim.IKTaskSet(); names = set(md.getMarkerNames().get(i) for i in range(md.getMarkerNames().getSize()))
    for i in range(model.getMarkerSet().getSize()):
        mk = model.getMarkerSet().get(i).getName()
        if mk in names:
            tk = osim.IKMarkerTask(); tk.setName(mk); tk.setWeight(WEIGHT.get(mk, 1.0)); tk.setApply(True); tasks.cloneAndAppend(tk)
    ik.set_IKTaskSet(tasks); ik.set_report_errors(True); ik.setResultsDir(os.path.dirname(out_mot))
    ik.run()
    err = os.path.join(os.path.dirname(out_mot), f"{ik.getName()}_ik_marker_errors.sto")
    rms = mx = None
    if os.path.exists(err):
        import addb_io as io
        e = io.read_mot(err); rms = float(np.mean(e.column("marker_error_RMS"))); mx = float(np.max(e.column("marker_error_max")))
    log(f"  IK {os.path.basename(out_mot)}: RMS {rms*100 if rms else float('nan'):.2f} cm, 최대 {mx*100 if mx else float('nan'):.2f} cm")
    return rms, mx


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", required=True); ap.add_argument("--subject", type=int, required=True); ap.add_argument("--trials", nargs="+", required=True)
    ap.add_argument("--out", default=os.path.join(ROOT, "data/carrt")); ap.add_argument("--base-model", default=os.path.join(ROOT, "data/moore2015/subject3/subject3_scaled.osim"))
    ap.add_argument("--load", type=float, default=None, help="손 하중 kg (TASKS 표 대신 강제)")
    ap.add_argument("--static", default=None, help="정지 자세로 쓸 trial (기본 첫 trial 의 0~0.3 s)")
    a = ap.parse_args(); subj = f"subject{a.subject}"; out_dir = os.path.join(a.out, subj); os.makedirs(out_dir, exist_ok=True)
    mass, height = DEMO[a.subject]

    import glob as _glob
    def find_c3d(tr):
        """전신(Full Body) 폴더만 — 'Half Body Right Hand' 는 오른팔 20마커뿐이라 쓰지 않는다.
        tr 이 코드만(PB)이면 그 코드의 첫 번째 파일(PB1 …)을 고른다."""
        code = code_of(tr)
        for d in sorted(os.listdir(os.path.join(a.raw, "Data", f"Subject {a.subject}"))):
            if "full" not in d.lower(): continue
            folder = os.path.join(a.raw, "Data", f"Subject {a.subject}", d, "ADL" if code in TASKS else "ROM")
            cand = os.path.join(folder, tr + ".c3d")
            if os.path.exists(cand): return cand
            if tr == code:
                hits = sorted(_glob.glob(os.path.join(folder, code + "[0-9! ]*.c3d")))
                if hits: return hits[0]
    def to_trc(tr):
        trc = os.path.join(out_dir, f"{tr}_markers.trc")
        if not os.path.exists(trc): c3d_to_trc(find_c3d(tr), trc)
        return trc
    a.trials = [os.path.basename(find_c3d(tr))[:-4] if find_c3d(tr) else tr for tr in a.trials]
    static = a.static or a.trials[0]
    model = os.path.join(out_dir, f"{subj}_scaled.osim")
    if not os.path.exists(model):
        model = scale_model(a.base_model, to_trc(static), out_dir, subj, mass, log=print)
    print(f"[{subj}] 모델 {model} (질량 {mass} kg, 키 {height} m, 정지 자세 {static} 0~0.3 s)")
    rows = []
    for tr in a.trials:
        code = code_of(tr); desc, load, bimanual = TASKS.get(code, (tr, 0.0, False))
        if a.load is not None: load = a.load
        sub_dir = find_c3d(tr)
        if not sub_dir: print(f"  [없음] {tr}"); continue
        trc = os.path.join(out_dir, f"{tr}_markers.trc"); names, t, rotated, dropped = c3d_to_trc(sub_dir, trc)
        print(f"[{subj}/{tr}] {desc}: 마커 {len(names)}개 (버림 {dropped}), {t[-1]:.2f} s, z-up 회전 {rotated}")
        rms, mx = run_ik(model, trc, os.path.join(out_dir, f"{tr}_ik.mot"), print)
        mp = model
        if load: mp = with_hand_load(model, os.path.join(out_dir, f"{subj}_load{load:g}{'_bi' if bimanual else ''}.osim"), load, bimanual)
        rows.append(dict(trial=tr, task=desc, load_kg=load, bimanual=bimanual, ik_rms_cm=round((rms or 0) * 100, 2), ik_max_cm=round((mx or 0) * 100, 2), model=os.path.relpath(mp, ROOT)))
    with open(os.path.join(out_dir, "ik_log.csv"), "a", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); fh.seek(0, 2)
        if fh.tell() == 0: w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
