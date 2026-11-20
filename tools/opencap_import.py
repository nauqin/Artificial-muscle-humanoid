"""OpenCap 세션 → 우리 형식 (`data/opencap/<이름>/<trial>_ik.mot` + 모델 + 몸 치수). **`osim` 환경.**

OpenCap(https://www.opencap.ai)은 아이폰 2대 이상으로 찍어 클라우드에서 3D 관절각까지 계산해 준다. 결과 모델이
LaiUhlrich2022 인데, 좌표 이름·관절 축이 우리 Rajagopal(팔 포함) 모델과 같아서 관절각을 그대로 옮길 수 있다
(2026-11-19 확인: 좌표 35개 이름 일치, hip/knee/ankle/back/acromial/elbow 관절 축 방향 일치).

    # 1) 세션 번호나 주소로 바로 (공개 세션은 로그인 없이, 비공개는 토큰 필요 — 아래)
    python tools/opencap_import.py 4d5c3eb1-1a59-4ea1-9178-d3634610561c --name opencap_sample
    python tools/opencap_import.py https://app.opencap.ai/session/<세션 번호> --trials squat STS
    # 2) OpenCap 웹에서 "Download" 로 받은 zip 이나 푼 폴더
    python tools/opencap_import.py ~/Downloads/OpenCapData_xxxx.zip --name lift_box
→ data/opencap/<이름>/
     <이름>_scaled.osim        OpenCap 이 스케일한 사람 모델 (참고·영상용)
     <trial>_ik.mot            관절각, 라디안, 우리 좌표 이름 (knee_angle_*_beta 는 뺌)
     session.json              질량·키·trial 목록·출처
     raw/                      받은 원본 (.mot 도 단위, .trc, sessionMetadata.yaml, --videos 면 영상)

그 다음은 `tools/add_motion.py --opencap data/opencap/<이름> ...` 가 로봇 뼈대로 옮기고 지면반력·토크를 구한다.

비공개 세션: OpenCap 계정 토큰을 환경변수 `OPENCAP_TOKEN` 에 넣거나 `~/.opencap_token` 파일에 한 줄로 둔다.
토큰은 `python tools/opencap_import.py --login` 으로 받는다(아이디·비밀번호를 물어보고 `~/.opencap_token` 에 저장).
※ 로그인·비공개 세션 경로는 계정이 없어 **실행해 보지 못했다**(공개 세션과 폴더·zip 경로만 시험). 안 되면 웹에서 zip 으로 받아 2) 로 한다.
"""
import argparse, glob, json, os, re, shutil, sys, tempfile, urllib.request, zipfile
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0, ROOT)
import addb_io as io                                   # noqa: E402

API = "https://api.opencap.ai/"
UUID = re.compile(r"[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", re.I)
SKIP_TRIALS = ("calibration", "neutral")              # 관절각이 없는 trial
TOKEN_FILE = os.path.join(os.path.expanduser("~"), ".opencap_token")
TRANSLATIONS = ("pelvis_tx", "pelvis_ty", "pelvis_tz")


def token():
    t = os.environ.get("OPENCAP_TOKEN")
    if not t and os.path.exists(TOKEN_FILE): t = open(TOKEN_FILE).read().strip()
    return t or None


def api_get(path):
    req = urllib.request.Request(API + path, headers={"Accept": "application/json", "User-Agent": "addb-opensim-import"})
    if token(): req.add_header("Authorization", f"Token {token()}")
    try:
        with urllib.request.urlopen(req, timeout=60) as r: return json.load(r)
    except urllib.error.HTTPError as e:
        if e.code in (401, 403, 404):
            raise SystemExit(f"OpenCap API {e.code}: 세션이 비공개이거나 번호가 틀리다. 토큰을 넣거나(--login), 웹에서 zip 으로 받아 그 파일을 넘긴다.")
        raise


def login():
    import getpass, urllib.parse
    user = input("OpenCap 아이디: "); pw = getpass.getpass("비밀번호: ")
    req = urllib.request.Request(API + "login/", data=urllib.parse.urlencode({"username": user, "password": pw}).encode())
    with urllib.request.urlopen(req, timeout=60) as r: tok = json.load(r).get("token")
    if not tok: raise SystemExit("토큰을 못 받았다 (이메일 인증 단계가 있으면 웹에서 zip 으로 받는다)")
    open(TOKEN_FILE, "w").write(tok); print(f"토큰 저장: {TOKEN_FILE}")


def download_session(sid, raw, trials=None, videos=False, log=print):
    """API 로 세션을 받아 raw/ 에 OpenCap zip 과 비슷한 모양으로 둔다. 반환 (meta dict, 세션 이름)."""
    s = api_get(f"sessions/{sid}/"); os.makedirs(raw, exist_ok=True)
    json.dump(s, open(os.path.join(raw, "session_api.json"), "w"), indent=1)
    def get(url, rel):
        dst = os.path.join(raw, rel); os.makedirs(os.path.dirname(dst), exist_ok=True)
        if not os.path.exists(dst): urllib.request.urlretrieve(url, dst)
        return dst
    n = 0
    for t in s["trials"]:
        name = t["name"]; tags = {}
        for r in t["results"]: tags.setdefault(r["tag"], []).append(r["media"])
        if name == "neutral":
            if "opensim_model" in tags: get(tags["opensim_model"][0], "OpenSimData/Model/LaiUhlrich2022_scaled.osim")
            if "session_metadata" in tags: get(tags["session_metadata"][0], "sessionMetadata.yaml")
            if "marker_data" in tags: get(tags["marker_data"][0], "MarkerData/neutral.trc")
        if name in SKIP_TRIALS or (trials and name not in trials): continue
        if "ik_results" not in tags: log(f"  [건너뜀] {name}: 관절각 결과 없음"); continue
        get(tags["ik_results"][0], f"OpenSimData/Kinematics/{name}.mot")
        if "marker_data" in tags: get(tags["marker_data"][0], f"MarkerData/{name}.trc")
        if videos:
            for ci, u in enumerate(tags.get("video-sync", [])): get(u, f"Videos/Cam{ci}/InputMedia/{name}/{name}_sync.mp4")
        n += 1; log(f"  받음 {name}")
    if not n: raise SystemExit("받을 trial 이 없다 (--trials 이름을 확인)")
    return s.get("name") or s.get("sessionName") or sid[:8]


def read_metadata(folder):
    """sessionMetadata.yaml 에서 질량·키. 없으면 None."""
    hits = glob.glob(os.path.join(folder, "**", "sessionMetadata.yaml"), recursive=True)
    if not hits: return None, None
    mass = height = None
    for line in open(hits[0], encoding="utf-8", errors="replace"):
        if line.startswith("mass_kg"): mass = float(line.split(":")[1])
        if line.startswith("height_m"): height = float(line.split(":")[1])
    return mass, height


def convert_mot(src, dst, log=print):
    """OpenCap .mot(도) → 우리 .mot(라디안). 이동 좌표(pelvis_t*)는 그대로, knee_angle_*_beta(슬개골 종속 좌표)는 뺀다."""
    m = io.read_mot(src); keep = [j for j, n in enumerate(m.names) if not n.endswith("_beta")]
    D = m.data[:, keep].copy(); names = [m.names[j] for j in keep]
    if m.in_degrees:
        for j, n in enumerate(names):
            if n not in TRANSLATIONS: D[:, j] = np.radians(D[:, j])
    io.write_mot(dst, names, m.time, D, in_degrees=False)
    return names, m.time


def import_folder(folder, out_dir, name, trials=None, log=print):
    """OpenCap 폴더(zip 을 푼 것 또는 download_session 결과)에서 모델·관절각을 찾아 우리 형식으로."""
    mots = sorted(glob.glob(os.path.join(folder, "**", "Kinematics", "*.mot"), recursive=True)) or \
        [f for f in sorted(glob.glob(os.path.join(folder, "**", "*.mot"), recursive=True)) if "Dynamics" not in f]
    models = sorted(glob.glob(os.path.join(folder, "**", "Model", "*.osim"), recursive=True)) or sorted(glob.glob(os.path.join(folder, "**", "*.osim"), recursive=True))
    if not mots: raise SystemExit(f"관절각(.mot)을 못 찾았다: {folder}")
    mass, height = read_metadata(folder); os.makedirs(out_dir, exist_ok=True)
    if models: shutil.copyfile(models[0], os.path.join(out_dir, f"{name}_scaled.osim"))
    done = []
    for f in mots:
        tr = os.path.splitext(os.path.basename(f))[0]
        if tr in SKIP_TRIALS or (trials and tr not in trials): continue
        tr_safe = re.sub(r"[^0-9A-Za-z_]+", "_", tr)
        names, t = convert_mot(f, os.path.join(out_dir, f"{tr_safe}_ik.mot")); done.append(tr_safe)
        log(f"  {tr_safe}: {len(t)} 프레임, {t[0]:.2f}~{t[-1]:.2f} s, 좌표 {len(names)}개")
    meta = dict(name=name, source="opencap", mass_kg=mass, height_m=height, trials=done, model=(f"{name}_scaled.osim" if models else None),
                note="관절각은 OpenCap(LaiUhlrich2022) 결과를 라디안으로 바꾼 것. 지면반력 없음.")
    json.dump(meta, open(os.path.join(out_dir, "session.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    if mass is None or height is None: log("  [주의] sessionMetadata.yaml 에서 질량·키를 못 읽었다 — add_motion.py 에 --mass --height 를 직접 준다")
    log(f"저장 {out_dir}  (질량 {mass} kg, 키 {height} m, trial {len(done)}개)")
    return meta


def run(src, name=None, out=os.path.join(ROOT, "data", "opencap"), trials=None, videos=False, log=print):
    """src: 세션 번호·주소·폴더·zip. 반환 (out_dir, meta)."""
    m = UUID.search(src) if not os.path.exists(src) else None
    if m:
        sid = m.group(0).lower(); tmp_raw = tempfile.mkdtemp(prefix="opencap_")
        log(f"OpenCap 세션 {sid} 받는 중 …"); sname = download_session(sid, tmp_raw, trials, videos, log)
        name = name or re.sub(r"[^0-9A-Za-z_]+", "_", sname); out_dir = os.path.join(out, name); raw = os.path.join(out_dir, "raw")
        if os.path.exists(raw): shutil.rmtree(raw)
        os.makedirs(out_dir, exist_ok=True); shutil.move(tmp_raw, raw)
    elif os.path.isdir(src) or zipfile.is_zipfile(src):
        name = name or re.sub(r"[^0-9A-Za-z_]+", "_", os.path.splitext(os.path.basename(os.path.normpath(src)))[0]); out_dir = os.path.join(out, name); raw = os.path.join(out_dir, "raw")
        if os.path.abspath(src) != os.path.abspath(raw):
            if os.path.exists(raw): shutil.rmtree(raw)
            if os.path.isdir(src): shutil.copytree(src, raw)
            else:
                os.makedirs(raw, exist_ok=True)
                with zipfile.ZipFile(src) as z:
                    for info in z.infolist():                # 경로가 raw/ 밖으로 나가는 항목은 풀지 않는다
                        dst = os.path.realpath(os.path.join(raw, info.filename))
                        if dst.startswith(os.path.realpath(raw) + os.sep): z.extract(info, raw)
    else:
        raise SystemExit(f"세션 번호(8-4-4-4-12 꼴)도, 있는 폴더·zip 도 아니다: {src}")
    return out_dir, import_folder(raw, out_dir, name, trials, log)


def main():
    ap = argparse.ArgumentParser(description="OpenCap 세션을 우리 형식으로")
    ap.add_argument("src", nargs="?", help="세션 번호 / app.opencap.ai 주소 / 받은 zip / 푼 폴더")
    ap.add_argument("--name", default=None, help="저장할 이름 (기본: 세션 이름)"); ap.add_argument("--out", default=os.path.join(ROOT, "data", "opencap"))
    ap.add_argument("--trials", nargs="+", default=None, help="이 trial 만"); ap.add_argument("--videos", action="store_true", help="동기화된 영상도 받는다 (API 경로만)")
    ap.add_argument("--login", action="store_true", help="OpenCap 토큰을 받아 ~/.opencap_token 에 저장")
    a = ap.parse_args()
    if a.login: return login()
    if not a.src: ap.error("src 가 필요하다")
    run(a.src, a.name, a.out, a.trials, a.videos)


if __name__ == "__main__":
    main()
