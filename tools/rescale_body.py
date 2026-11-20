"""로봇 크기(예: 40 kg, 1.35 m)로 뼈대를 줄이고, 그 뼈대로 상체 ID 를 다시 돌려 tasks.csv 를 만든다.

왜 다시 돌리나: 몸 자체의 토크는 질량×키에 비례해 그냥 줄지만, 손 하중(물병 3.79 kg)은 안 줄어든다. 그래서 상체는
줄인 뼈대 + 그대로인 하중으로 ID 를 다시 해야 한다. 하체는 손 하중이 없으니 기존 토크를 질량×키 비율로 환산한다
(`design_directions.py --ref-mass --ref-height`). 다만 하체 Moco 검증용으로 지면반력·골반 위치를 같은 비율로 줄인 파일도 만든다.

    python tools/rescale_body.py --mass 40 --height 1.35 --out data/robot40
→ data/robot40/robot40.osim (+ 손 하중 모델), out/upper40/<dataset>/<subject>/<trial>/, report/upper/tasks_40kg.csv,
   data/robot40/lower/<trial>_ik.mot·_grf.mot (하체 검증용, 골반 위치 ×키비, 힘 ×질량비, 모멘트 ×질량비×키비)
"""
import argparse, csv, os, shutil, sys
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "tools"))
import addb_io as io                                   # noqa: E402
BASE = "data/moore2015/subject3/subject3_scaled.osim"; BASE_HEIGHT = 1.62   # 팔 포함 Rajagopal, Moore subject3 (63.68 kg)


def scale_model(base, out, length_ratio, mass, extra=None):
    """모든 뼈를 length_ratio 로. extra = {뼈 이름: (x, y, z) 추가 배율} 이면 그 뼈만 그 배율을 더 곱한다 (2026-11-19, CAD 마디 길이 맞춤)."""
    import opensim as osim
    m = osim.Model(base); st = m.initSystem(); ss = osim.ScaleSet(); extra = extra or {}
    for i in range(m.getBodySet().getSize()):
        bn = m.getBodySet().get(i).getName(); f = np.array(extra.get(bn, (1.0, 1.0, 1.0))) * length_ratio
        sc = osim.Scale(); sc.setSegmentName(bn); sc.setScaleFactors(osim.Vec3(*f)); sc.setApply(True); ss.cloneAndAppend(sc)
    ok = m.scale(st, ss, True, float(mass)); st = m.initSystem()
    assert ok, "scale 실패"
    m.setName(os.path.splitext(os.path.basename(out))[0]); m.printToXML(out)
    print(f"모델 {out}: 질량 {m.getTotalMass(st):.1f} kg, 길이 ×{length_ratio:.4f}" + (f", 추가 {extra}" if extra else ""))
    return out


# 하체 검증 입력의 출처 사람 (2026-11-19 --per-source): 길이비는 로봇 키 ÷ 그 사람 키, 질량비는 로봇 질량 ÷ 그 사람 질량.
# 지금까지(기본)는 둘 다 키 1.62 m·질량 63.5 kg 로 나눴다 — LabValidation 피험자3(1.69 m)은 길이가 4 % 크게, Moore 피험자3(63.68 kg)은 질량이 0.3 % 크게 환산됐다.
SOURCES = {"data/subject3": (1.69, 63.5), "data/moore2015/subject3": (1.62, 63.68)}


def scale_lower_inputs(src_dir, trials, out_dir, length_ratio, mass_ratio):
    """하체 검증용 IK·GRF: 골반 위치 ×길이비, 힘 ×질량비, COP ×길이비, 자유모멘트 ×질량비×길이비."""
    os.makedirs(out_dir, exist_ok=True)
    for tr in trials:
        ik = io.read_mot(os.path.join(src_dir, f"{tr}_ik.mot")); D = ik.data.copy()
        for j, n in enumerate(ik.names):
            if n in ("pelvis_tx", "pelvis_ty", "pelvis_tz"): D[:, j] *= length_ratio
        io.write_mot(os.path.join(out_dir, f"{tr}_ik.mot"), ik.names, ik.time, D, in_degrees=ik.in_degrees)
        g = io.read_mot(os.path.join(src_dir, f"{tr}_grf.mot")); G = g.data.copy()
        for j, n in enumerate(g.names):
            if "_v" in n and n[-3:-1] == "_v": G[:, j] *= mass_ratio                      # 힘
            elif n.endswith(("_px", "_py", "_pz")): G[:, j] *= length_ratio               # COP
            elif n.endswith(("_mx", "_my", "_mz")) or "torque" in n: G[:, j] *= mass_ratio * length_ratio   # 자유모멘트
        io.write_mot(os.path.join(out_dir, f"{tr}_grf.mot"), g.names, g.time, G)
        print(f"  하체 입력 {tr}: IK·GRF 저장")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mass", type=float, required=True); ap.add_argument("--height", type=float, required=True)
    ap.add_argument("--out", default=None, help="모델 폴더 (기본 data/robot<mass>)"); ap.add_argument("--tasks", default="report/upper/tasks.csv")
    ap.add_argument("--lowpass", type=float, default=6.0); ap.add_argument("--skip-id", action="store_true")
    ap.add_argument("--tag", default="", help="출력 이름 뒤에 붙일 말 (예: _h160 → data/robot40_h160, out/upper40_h160, tasks_40kg_h160.csv). 같은 질량·다른 키를 나란히 둘 때")
    ap.add_argument("--segment-scale", action="append", default=[], metavar="BODY=X,Y,Z",
                    help="그 뼈에만 추가 배율 (키 배율 위에 곱함, 여러 번 가능). 예: humerus_r=1,1.0967,1 (CAD 위팔 길이). 2026-11-19")
    ap.add_argument("--per-source", action="store_true", help="하체 검증 입력의 길이·질량비를 출처 사람별로 (LabValidation 1.69 m·63.5 kg, Moore 1.62 m·63.68 kg). 기본은 지금까지처럼 1.62 m·63.5 kg")
    a = ap.parse_args()
    import opensim as osim; osim.Logger.setLevelString("Error")
    import run_trial as rt; from carrt_prepare import with_hand_load
    name = f"robot{a.mass:g}{a.tag}"; out_dir = os.path.join(ROOT, a.out or f"data/{name}"); os.makedirs(out_dir, exist_ok=True)
    lr = a.height / BASE_HEIGHT; base = os.path.join(ROOT, BASE)
    extra = {s.split("=")[0]: tuple(float(v) for v in s.split("=")[1].split(",")) for s in a.segment_scale}
    model = scale_model(base, os.path.join(out_dir, f"{name}.osim"), lr, a.mass, extra)
    if extra or a.per_source:
        import json
        json.dump(dict(height_m=a.height, mass_kg=a.mass, base=BASE, base_height_m=BASE_HEIGHT, length_ratio=lr, segment_extra=extra, per_source=a.per_source,
                       note="segment_extra: CAD(step/) 실측 마디 길이 ÷ 균일 축소 모델 길이. 출처는 README 12절"), open(os.path.join(out_dir, "segment_scales.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    # 상체: 손 하중별 모델 + ID
    rows_out = []; tasks_in = a.tasks if os.path.isabs(a.tasks) else os.path.join(ROOT, a.tasks)
    for row in csv.DictReader(open(tasks_in, encoding="utf-8")):
        load = float(row.get("hand_load_kg") or 0); bimanual = "_bi" in (row.get("note") or "") or row["task"] == "lift_box"
        if load > 0:
            mp = os.path.join(out_dir, f"{name}_load{load:g}{'_bi' if bimanual else ''}.osim")
            if not os.path.exists(mp): with_hand_load(model, mp, load, bimanual)
        else: mp = model
        od = os.path.join(ROOT, f"out/upper{a.mass:g}{a.tag}", row["dataset"].lower(), row["subject"], row["trial"]); idp = os.path.join(od, f"{row['trial']}_id.sto")
        if not a.skip_id and not os.path.exists(idp):
            print(f"[ID] {row['dataset']} {row['subject']} {row['trial']} 하중 {load} kg 모델 {os.path.basename(mp)}", flush=True)
            rt.run_trial_no_grf(mp, os.path.join(ROOT, row["ik"]), od, row["trial"], lowpass=a.lowpass, log=lambda *x: None)
        r = dict(row); r.update(id=os.path.relpath(idp, ROOT), mass_kg=f"{a.mass:g}", height_m=f"{a.height:g}", note=f"{name} 뼈대로 ID 재계산, 손 하중 {load} kg 그대로"); rows_out.append(r)
    tasks_out = os.path.join(ROOT, f"report/upper/tasks_{a.mass:g}kg{a.tag}.csv")
    with open(tasks_out, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys())); w.writeheader(); w.writerows(rows_out)
    print("tasks", tasks_out, len(rows_out), "행")
    # 하체 검증 입력 (data/subject3 걷기·스쿼트)
    for src, trials, sub in (("data/subject3", ("walking1_segment_0", "squats1_segment_0"), "lower"), ("data/moore2015/subject3", ("T047_segment_12",), "moore")):   # moore = 전신 걷기 검증용
        h_src, m_src = SOURCES[src] if a.per_source else (BASE_HEIGHT, 63.5)
        scale_lower_inputs(os.path.join(ROOT, src), trials, os.path.join(out_dir, sub), a.height / h_src, a.mass / m_src)


if __name__ == "__main__":
    main()
