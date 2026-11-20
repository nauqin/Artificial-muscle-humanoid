"""영상 한 개 → 프레임마다 사람 몸 표면 3D 점 (카메라 좌표, m). **`pose` 환경에서 돈다** (torch + CUDA).

    python tools/video_pose.py --video my.mp4 --out data/custom/my/my_pose.npz

보통은 직접 부르지 않는다 — `tools/add_motion.py --video ...` 가 `pose` 환경의 파이썬으로 이 파일을 대신 돌린다.

모델: NLF (Neural Localizer Fields, Sárándi & Pons-Moll, NeurIPS 2024) `nlf_l_multi_0.3.2.torchscript`.
  - 영상 한 장에서 SMPL 몸 표면 6,890점과 관절 24개의 **절대 3D 위치**(카메라 좌표, 미터 단위 추정)를 낸다.
  - 가중치는 처음 한 번 GitHub release 에서 받아 `~/.cache/nlf/` 에 둔다 (493 MB).
  - **비상업 연구용** 라이선스다. 논문·발표에 쓰면 NLF 를 인용한다.
카메라는 **고정**이라고 가정한다(삼각대). 화각을 모르면 55° 로 둔다 — 화각이 틀리면 깊이(카메라에서 떨어진 거리)가
비례해서 틀리지만 몸의 모양·관절각에는 거의 영향이 없다. 폰 기종을 알면 `--fov` 로 준다(아이폰 기본 카메라 세로 약 69°·가로 약 58°).

출력 npz: verts (T, 6890, 3) m · joints (T, 24, 3) m · vert_unc (T, 1024, 모델이 내는 대로) · time (T,) s · fps · size (W, H) · fov
          사람이 안 잡힌 프레임은 NaN.
"""
import argparse, os, sys, time, urllib.request
import numpy as np

URL = "https://github.com/isarandi/nlf/releases/download/v0.3.2/nlf_l_multi_0.3.2.torchscript"
DEFAULT_MODEL = os.path.join(os.path.expanduser("~"), ".cache", "nlf", "nlf_l_multi_0.3.2.torchscript")


def fetch_model(path):
    if os.path.exists(path) and os.path.getsize(path) > 400e6: return path
    os.makedirs(os.path.dirname(path), exist_ok=True)
    print(f"NLF 가중치를 받는다 (493 MB, 처음 한 번): {URL}", flush=True)
    urllib.request.urlretrieve(URL, path + ".part"); os.replace(path + ".part", path)
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--video", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--model", default=DEFAULT_MODEL)
    ap.add_argument("--fov", type=float, default=55.0, help="화면 긴 변 기준 화각(도). 모르면 55")
    ap.add_argument("--start", type=float, default=None); ap.add_argument("--end", type=float, default=None, help="쓸 구간 (초)")
    ap.add_argument("--batch", type=int, default=8)
    ap.add_argument("--person", default="largest", help="여러 명이 찍혔을 때: largest(첫 프레임에서 가장 큰 사람) 또는 center(화면 가운데)")
    a = ap.parse_args()

    import cv2, torch
    import torchvision  # noqa: F401  — 모델 안의 torchvision::nms 연산을 등록한다 (없으면 jit.load 가 죽는다)
    if not torch.cuda.is_available(): print("[경고] CUDA 가 안 잡힌다 — CPU 로 돌면 프레임당 수 초가 걸린다", flush=True)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model = torch.jit.load(fetch_model(a.model), map_location=dev).eval()

    cap = cv2.VideoCapture(a.video)
    if not cap.isOpened(): raise SystemExit(f"영상을 못 연다: {a.video}")
    fps = cap.get(cv2.CAP_PROP_FPS) or 30.0; n_total = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    W, H = int(cap.get(cv2.CAP_PROP_FRAME_WIDTH)), int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT))
    f0 = int(round((a.start or 0.0) * fps)); f1 = int(round(a.end * fps)) if a.end is not None else n_total
    if f0: cap.set(cv2.CAP_PROP_POS_FRAMES, f0)
    print(f"영상 {os.path.basename(a.video)}: {W}×{H}, {fps:.2f} fps, 프레임 {f0}~{f1} / {n_total}", flush=True)

    def infer(frames):
        x = torch.from_numpy(np.stack(frames)[..., ::-1].copy()).permute(0, 3, 1, 2).to(dev)   # BGR → RGB, N×3×H×W uint8
        with torch.inference_mode():
            try: return model.detect_smpl_batched(x, default_fov_degrees=a.fov)
            except Exception as e:                       # 인자 이름이 다른 판이면 기본 화각으로
                if not getattr(infer, "warned", False): print(f"[경고] 화각 인자를 못 넘겼다({type(e).__name__}) — 모델 기본 화각 사용", flush=True); infer.warned = True
                return model.detect_smpl_batched(x)

    V, J, U, T = [], [], [], []; prev = None; tic = time.time(); fi = f0; buf = []; missed = 0

    def flush():
        nonlocal prev, missed
        if not buf: return
        p = infer(buf)
        for b in range(len(buf)):
            v = p["vertices3d"][b].float().cpu().numpy() / 1000.0; j = p["joints3d"][b].float().cpu().numpy() / 1000.0   # mm → m, (사람 수, 점, 3)
            u = p["vertex_uncertainties"][b].float().cpu().numpy() if "vertex_uncertainties" in p else np.zeros(v.shape[:2])
            if len(v) == 0:                                                                  # 사람을 못 잡은 프레임 (불확실도는 점 1024개짜리라 크기를 그대로 따른다)
                V.append(np.full(v.shape[1:], np.nan)); J.append(np.full(j.shape[1:], np.nan)); U.append(np.full(u.shape[1:], np.nan)); missed += 1; continue
            if prev is None:
                box = p["boxes"][b].float().cpu().numpy()                                    # (사람 수, 5) x, y, w, h, 점수
                if a.person == "center": k = int(np.argmin(np.hypot(box[:, 0] + box[:, 2] / 2 - W / 2, box[:, 1] + box[:, 3] / 2 - H / 2)))
                else: k = int(np.argmax(box[:, 2] * box[:, 3]))
            else: k = int(np.argmin(np.linalg.norm(j[:, 0] - prev, axis=1)))                 # 앞 프레임 골반에 가장 가까운 사람
            prev = j[k, 0]; V.append(v[k]); J.append(j[k]); U.append(u[k])
        buf.clear()

    while fi < f1:
        ok, fr = cap.read()
        if not ok: break
        buf.append(fr); T.append(fi / fps); fi += 1
        if len(buf) == a.batch:
            flush()
            if len(V) % (a.batch * 10) == 0: print(f"  {len(V)} 프레임  {len(V) / (time.time() - tic):.1f} fps", flush=True)
    flush(); cap.release()
    if not V: raise SystemExit("프레임을 하나도 못 읽었다")
    V = np.asarray(V, np.float32); J = np.asarray(J, np.float32); U = np.asarray(U, np.float32)
    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    np.savez_compressed(a.out, verts=V, joints=J, vert_unc=U, time=np.asarray(T[:len(V)]), fps=fps, size=np.array([W, H]), fov=a.fov, video=os.path.abspath(a.video))
    print(f"저장 {a.out}: {len(V)} 프레임 ({time.time() - tic:.0f} s), 사람 못 잡은 프레임 {missed}개", flush=True)
    if missed > 0.2 * len(V): print("[경고] 사람을 못 잡은 프레임이 20 % 를 넘는다 — 전신이 화면에 다 들어오게 다시 찍는 게 좋다", flush=True)


if __name__ == "__main__":
    main()
