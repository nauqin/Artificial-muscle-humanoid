"""모듈 먼저 실험(2026-11-05) 후처리: (1) 어깨 힘 ×1.5 json, (2) Moco 검증용 모델(손 하중 없는 것·찬 물병 3.79 kg 든 것).

맥에서 한 방식 그대로: _f15 json 은 경로 설계 결과의 force_N 만 1.5배(다른 값은 그대로), 들기 모델은 같은 경로를
data/robot40_h160/robot40_h160_load3.79.osim 위에 path_design.build_model 로 다시 만든 것이다(감는 구는 box-scale 배율).

    python report/modulefirst/prepare_models.py <space> <경로 json> [--force-scale 1.5] [--box-scale 0.9877]
→ (힘 배율이 있으면) <json>_f15.json, report/moco/newmus_paths_<이름>_model.osim, ..._load3.79_model.osim
"""
import argparse, json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "report"))
import path_design as pd                                     # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("space"); ap.add_argument("paths"); ap.add_argument("--force-scale", type=float, default=1.0); ap.add_argument("--box-scale", type=float, default=0.9877)
ap.add_argument("--base", default="data/robot40_h160/robot40_h160.osim"); ap.add_argument("--load-base", default="data/robot40_h160/robot40_h160_load3.79.osim")
a = ap.parse_args(); pd.configure(a.space)
if pd.WRAP and a.box_scale != 1.0: pd.WRAP = dict(pd.WRAP, center=tuple(v * a.box_scale for v in pd.WRAP["center"]), radius=pd.WRAP["radius"] * a.box_scale)
pj = a.paths if os.path.isabs(a.paths) else os.path.join(ROOT, a.paths); out = json.load(open(pj, encoding="utf-8"))
if a.force_scale != 1.0:
    tag = f"_f{int(round(a.force_scale * 10))}"
    for m in out["muscles"]: m["force_N"] = m["force_N"] * a.force_scale
    out["note"] = f"힘 ×{a.force_scale:g} (경로 설계 뒤, 맥 _f15 와 같은 방식)"; pj = pj.replace(".json", f"{tag}.json")
    json.dump(out, open(pj, "w", encoding="utf-8"), indent=1, ensure_ascii=False); print("저장", os.path.relpath(pj, ROOT))
stem = os.path.basename(pj).replace(".json", "")                      # newmus_paths_<space>_k<k><tag>
if a.force_scale != 1.0: pd.build_model(out, os.path.join(ROOT, a.base), os.path.join(ROOT, f"report/moco/{stem}_model.osim"))   # 힘이 바뀐 모델 (배율 1 이면 path_design 이 이미 만들었다)
if a.space != "lower": pd.build_model(out, os.path.join(ROOT, a.load_base), os.path.join(ROOT, f"report/moco/{stem}_load3.79_model.osim"))
