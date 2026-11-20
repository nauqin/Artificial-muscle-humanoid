# 새 동작 넣기 — OpenCap·영상 → 로봇 관절 토크 → 최종 A안 검사·재설계 (2026-11-19)

10/28 계획 H9·U3("OpenCap 으로 원하는 동작을 찍어 6-1 부터 재설계")를 지금 파이프라인에 붙였다.
**명령 하나**로 동작을 로봇 뼈대(40 kg · 1.60 m)의 관절 토크로 바꾸고, 설계 입력표(`tasks.csv`)에 넣고, 지금 설계가 그 동작을 낼 수 있는지 본다.

```
OpenCap 세션 (아이폰 2대+) ── opencap_import ─────────────────────┐
영상 한 개 (폰 1대, 삼각대) ── video_pose (NLF) ── video_ik ──────┤→ 사람 관절각 → 로봇 뼈대로 → 지면반력 추정 → ID → tasks_<묶음>.csv
관절각 .mot ──────────────────────────────────────────────────────┘                                                 │
                                                       fit_check (최종 A안 가능?)  ←────────────────────────────────┤
                                                       design_directions → path_design → moco (재설계) ←────────────┘
```

## 1. 쓰는 법

```bash
conda activate osim
export PYTHONUTF8=1 MKL_THREADING_LAYER=SEQUENTIAL     # 윈도 (TUTORIAL 9-9)

# (가) OpenCap 세션 — 세션 번호·app.opencap.ai 주소·웹에서 받은 zip·푼 폴더 아무거나
python tools/add_motion.py --opencap <세션 번호 또는 zip> --trials lift1 lift2 --set lift --task lift_box --hand-load 3 --bimanual --check

# (나) 영상 한 개 — 찍힌 사람 키·몸무게를 꼭 준다
python tools/add_motion.py --video ~/Videos/reach.mp4 --height 1.70 --mass 62 --set reach --task reach_high --check

# (다) 이미 있는 관절각 (우리 좌표 이름의 .mot)
python tools/add_motion.py --ik some_ik.mot --src-height 1.69 --name trial1 --set misc --task carry

# 그 묶음을 기존 과제에 더해 다시 설계 (부위마다 수 분). 경로 설계 명령은 끝에 찍어 준다
python tools/add_motion.py --opencap ... --set lift --redesign lower,trunk,shoulder_r,elbow_r
```

결과
| 파일 | 내용 |
| --- | --- |
| `report/custom/tasks_<묶음>.csv` | 설계 입력표 (`report/upper/tasks_40kg_h160.csv` 와 같은 열). `design_directions.py`·`path_design.py` 의 `--tasks` 로 쓴다 |
| `out/custom40_h160/<묶음>/<이름>/<이름>_id.sto`, `_id.png` | 로봇 관절 토크, 그래프 (**먼저 그림을 본다** — 튀는 곳이 있으면 6절) |
| `data/custom/<묶음>/<이름>_ik.mot`, `_grf.mot` | 로봇 뼈대 관절각(라디안, 100 Hz, 필터 후), 추정 지면반력 |
| `report/custom/fit_check_<묶음>_last.csv` | `--check` 결과 |
| `data/opencap/<세션>/` | OpenCap 에서 받은 원본(`raw/`)과 사람 관절각 |
| `data/custom/<묶음>/video/` | 영상: 3D 점(`_pose.npz`), 마커(`.trc`), 사람 모델, IK, 진단(`_video.json`) |

주요 옵션 (`python tools/add_motion.py -h`)
| 옵션 | 뜻 |
| --- | --- |
| `--set` | 묶음 이름. 같은 묶음끼리 한 과제표에 쌓인다 (같은 trial 이름이면 바꿔 쓴다) |
| `--task` | 과제 이름 (설계 결과의 "과제별 필요 용량"이 이 이름으로 묶인다) |
| `--hand-load KG`, `--bimanual` | 손에 든 물건. 구간 **내내** 들고 있다고 본다 → 들기 전·내려놓은 뒤는 `--t0 --t1` 로 자른다 |
| `--t0 --t1` | 쓸 구간 (초) |
| `--lowpass HZ` | 관절각 컷오프. 기본 OpenCap·영상 **4 Hz**, `--ik` 6 Hz. 점프·달리기는 8 |
| `--upper-only` | 지면반력 없이 상체 토크만 (의자·난간·벽처럼 발 말고 다른 곳이 닿는 동작, 발이 안 보이는 영상) |
| `--check` | 최종 A안으로 낼 수 있나 (`report/fit_check.py`, 2절) |
| `--redesign 부위,…` | 기존 과제 + 이 묶음으로 관절 배분비·토크 용량 재설계 |
| `--robot MODEL MASS HEIGHT` | 다른 로봇 몸 (예: `data/robot40_h175sk/robot40_h175sk.osim 40 1.75`). 이때 `--check`/`--redesign` 은 꺼진다(최종 A안이 160 cm 기준) |
| `--fov`, `--footlock auto/on/off` | 영상: 화각(모르면 55°), 디딘 발 기준 위치 재구성 |

## 2. `--check` 판정 (`report/fit_check.py`, Moco 없이 수 초)

trial × 부위(하체·몸통·어깨·팔꿈치)마다
1. **토크형**: 관절 배분비 V·토크 용량 c 만으로(0 ≤ u ≤ 1) 그 토크를 내나 — 근육 개수·용량이 모자란가
2. **케이블형**: 실제 케이블 경로의 그 자세 모멘트암 × 케이블 힘으로 내나 — 그 자세에서 지렛대가 모자란가
3. **길이**: 그 자세들이 케이블 길이 범위를 설계 때보다 얼마나 넓히나, 그 결과 모듈 스트로크를 넘나

케이블형 남는 토크 95퍼센타일 ≤ 피크의 10 % 이고 길이 초과가 없으면 **가능**. 정적 검사라 통과해도 확정은 `report/modulefirst/run_moco.sh`·`moco_path_muscles.py` 로 한다.
(설계 범위 자체가 이미 상한을 몇 mm 넘는 케이블이 있다 — 팔꿈치 M0 83.7 mm vs 80. 그래서 "1 mm 넘게 **넓히고** 그 결과 상한도 넘을 때"만 초과로 본다.)

예시 (`report/custom/tasks_demo.csv`, OpenCap 공개 세션의 키 1.86 m 남성):
| trial | 하체 | 몸통 | 어깨 | 팔꿈치 |
| --- | --- | --- | --- | --- |
| 스쿼트 (OpenCap, 4 Hz) | 가능 (피크 57 N·m) | 가능 (89) | 가능 | 가능 |
| 앉았다 일어나기 (`--upper-only`) | — | 가능 (101, 케이블형 7.9 %) | **길이 초과** 11 mm (팔을 앞으로 뻗어 범위 20 mm 넓힘) | 가능 |
| 점프 (`--ik`, 8 Hz) | **길이 초과** 21 mm (착지 발목) | **토크 부족** 22 % (154 N·m) | 길이 초과 3 mm | 가능 |

점프 결과는 9/15 결정("동작 범위 = 걷기 + 스쿼트, 점프 제외")과 맞는다.

## 3. 무엇을 하나

### 3-1. OpenCap → 사람 관절각 (`tools/opencap_import.py`)
OpenCap 결과 모델 LaiUhlrich2022 는 우리 Rajagopal(팔 포함) 모델과 **좌표 35개 이름이 같고 관절 좌표계 방향도 같다**(둘 다 Rajagopal2015 에서 나왔다. hip·knee·ankle·subtalar·back·acromial·elbow·radioulnar·손목 관절 프레임을 두 모델에서 직접 비교). 그래서 관절각을 그대로 옮긴다(도 → 라디안, 슬개골 종속 좌표 `*_beta` 는 뺌).
공개 세션은 로그인 없이 API 로 받는다. 비공개는 `--login`(토큰을 `~/.opencap_token` 에) — **계정이 없어 이 경로는 실행해 보지 못했다.** 안 되면 웹에서 zip 으로 받아 넘긴다.

### 3-2. 영상 한 개 → 사람 관절각 (`tools/video_pose.py`, `tools/video_ik.py`)
1. **NLF**(Neural Localizer Fields, Sárándi & Pons-Moll, NeurIPS 2024)가 프레임마다 SMPL 몸 표면 6,890점의 3D 위치(카메라 좌표, m)를 낸다. `pose` 환경(4절), RTX 5080 에서 약 15 fps.
2. 표면 점 84개를 뼈 마커로 쓴다. 어느 점을 어느 뼈의 어디에 붙일지는 `tools/smpl_calibrate.py` 가 한 번 정해 `tools/smpl_markers.json` 에 둔다 — OpenCap 공개 세션의 영상(아이폰 3대)과 OpenCap 관절각이 같은 순간이라, 영상을 NLF 로 돌린 뒤 표면 점마다 "어느 뼈 좌표계에서 가장 안 움직이나"로 고른다(스쿼트·걷기·드롭점프로 보정, 뼈 기준 흔들림 2.4~5.4 cm). 관절 바로 옆 피부는 두 뼈 사이에서 늘어나므로 뺀다(뺀 전/후 앉았다 일어나기 고관절 굴곡, 카메라 0: 12.1 → 9.8°).
3. 키: SMPL 관절 사슬 길이로 키를 추정하고 `--height` 에 맞춰 전체 크기를 맞춘다(카메라 한 대는 절대 크기를 모른다).
4. 위쪽: 카메라 옆 기울기 0(삼각대), 앞뒤 기울기는 "무게중심이 평균적으로 발 위"로, 여러 걸음이면 발 디딘 자리들이 한 평면인 것도 쓴다.
5. **깊이 평활**: 카메라 시선 방향의 몸 전체 이동만 1 Hz 로 거른다. 카메라 한 대는 거리가 프레임마다 2~6 cm 출렁이고, 이것이 무게중심 가속도 → 지면반력 → 하체 토크를 3~5배로 부풀렸다(고치기 전 달리기 하체 피크 OpenCap 100 → 영상 378~446 N·m, 고친 뒤 115~142).
6. 스케일(뼈마다 SMPL 뼈 길이 비) → OpenSim IK.
7. **디딘 발 기준 위치**: 디딘 발이 미끄러지지 않게 골반 위치를 관절각에서 다시 만든다(디딤 판정은 깊이 성분을 뺀 발 속도 + 높이, 히스테리시스). 트레드밀처럼 디딘 발이 움직이면 자동으로 꺼진다.

### 3-3. 로봇 뼈대로 (`add_motion.retarget`)
관절각은 그대로, 골반 위치는 키 비(1.60 / 사람 키)로, 로봇에서 잠긴 좌표(subtalar·mtp·손목)는 0, 발바닥이 y=0 에 닿게 높이를 맞춘다. 상체 토크가 사람 몸 크기로 안 줄어드는 손 하중 문제(11/3)는 로봇 뼈대 + 하중 그대로로 ID 를 해서 피한다(`rescale_body.py` 와 같은 원칙).

### 3-4. 지면반력 추정 (`tools/grf_estimate.py`)
힘판이 없으니 몸 전체 뉴턴–오일러로 구한다: 합력 = M·(무게중심 가속도 − 중력), 무게중심 둘레 모멘트 = 각운동량 변화율. 남는 것은 두 발 배분뿐이다.
- 합력 작용점이 발 중심선 위 어디에 있는지로 두 발에 나눈다(지렛대 비). 발이 4 cm 넘게 뜨면 그 발 몫은 0.
- 두 발이 서로 밀어내는 힘(합력에는 안 나타나는 성분)은 각 발의 힘이 무게중심 위 1.4배 높이의 한 점을 향한다고 둔다(`PIVOT`). 0·1.0·1.4·1.8 중 힘판 데이터에 가장 맞는 값.
- 이 지면반력으로 같은 모델·같은 관절각의 ID 를 돌리면 골반 잔차가 0 에 가깝다("같은 모델 짝" 원칙). 작용점이 발 밖으로 4 cm 넘게 나가는 프레임만 잔차가 남고, 로그에 비율이 찍힌다.

## 4. 환경

`osim` 은 그대로. 영상(나)만 `pose` 환경이 더 필요하다 (torch 와 opensim 의 numpy·MKL 이 부딪혀 따로 둔다).
```bash
conda create -n pose python=3.11 -y
~/miniforge3/envs/pose/python.exe -m pip install torch torchvision --index-url https://download.pytorch.org/whl/cu128   # RTX 50xx 는 cu128 이상
~/miniforge3/envs/pose/python.exe -m pip install -r tools/requirements_pose.txt
```
`add_motion.py` 는 `~/miniforge3/envs/pose/python.exe` 를 찾는다. 다른 곳이면 `POSE_PYTHON` 환경변수로 알려 준다.
NLF 가중치(493 MB)는 처음 한 번 GitHub release 에서 `~/.cache/nlf/` 로 받는다. **NLF 는 비상업 연구용 라이선스**다 — 논문·발표에 쓰면 인용한다.
`tools/opencap_core/LaiUhlrich2022.osim` 은 opencap-core(Apache-2.0, `LICENSE` 같이 둠)의 일반 모델이다(영상 경로의 사람 모델 바탕).

## 5. 믿을 수 있는 정도

### 5-1. 지면반력 추정 vs 힘판 (`python tools/grf_estimate.py --validate` → `grf_estimate_validation.csv`)
40 kg · 1.60 m 로봇 뼈대, 힘판 지면반력을 질량·키 비로 줄인 것과 비교. 토크 RMSE 는 오른쪽 / 왼쪽, 괄호는 힘판 ID 피크.
| | 스쿼트 (LabValidation) | 걷기 (LabValidation) | 트레드밀 걷기 (Moore)* |
| --- | --- | --- | --- |
| 수직력 합 RMSE | 4.8 %BW | 4.4 %BW | 7.0 %BW |
| 발별 수직력 RMSE | 3.6 %BW | 6.2 %BW | 6.3 %BW |
| 작용점 오차 | 2.6 cm | 2.8 cm | 3.8 cm |
| 고관절 굴곡 토크 RMSE | 1.4 / 4.3 N·m (45 / 52) | 5.9 / 6.2 (33 / 32) | 9.2 / 12.7 (76 / 93) |
| 고관절 내전 | 7.5 / 6.4 (15 / 16) | 3.4 / 5.0 (30 / 35) | 5.5 / 6.5 (43 / 45) |
| 무릎 | 3.9 / 6.5 (48 / 35) | 5.5 / 7.6 (16 / 16) | 7.6 / 10.5 (80 / 58) |
| 발목 | 2.4 / 6.4 (24 / 32) | 7.3 / 10.5 (60 / 56) | 8.1 / 11.7 (73 / 56) |

\* Moore 힘판 ID 는 15 Hz(TUTORIAL 9-6), 추정은 6 Hz 라 피크 비교에 필터 차이가 섞여 있다(추정 고관절 피크 44 / 39 가 작게 나온다).
처음 판에서 작용점 식의 부호가 틀려 토크 RMSE 가 50~115 N·m 였다 — 힘판 비교가 없었으면 몰랐을 것이다.

### 5-2. 영상 한 개 vs OpenCap (`python tools/video_validate.py` → `video_validation.csv`, `.log`)
보정에 **안 쓴** jump·run_2_5ms(트레드밀)·STS 를 카메라 3대 각각 **한 대씩만** 써서(9건) OpenCap(3대) 결과와 비교했다. 기준 OpenCap 도 영상 추정이다(실험실 마커 대비 관절각 평균 4~5°, Uhlrich 2023).

관절각 RMSE 중앙값 (범위): 고관절 굴곡 7.1° (2.8~19.2), 내전 3.1°, 회전 3.9°, 무릎 5.0° (4.0~8.3), 발목 8.8° (5.5~10.0), 요추 굴곡 4.1°, 옆굽힘 2.3°, 회전 2.8°, 어깨 굴곡 6.3°, 내전 3.1°, 팔꿈치 5.9°.
나쁜 쪽(고관절 19°)은 STS 의 정면에 가까운 카메라 — 굴곡이 깊이 방향으로 숨는다. → **옆에서 45° 앞쪽**에서 찍는다.

로봇 토크 (95퍼센타일, 영상 ÷ OpenCap):
| 부위 | 중앙값 | 범위 | 판단 |
| --- | --- | --- | --- |
| 팔꿈치 | 0.97 | 0.88~1.04 | 영상 한 개로 충분 |
| 몸통 | 0.99 | 0.84~2.21 (트레드밀 달리기만 1.7~2.2) | 대체로 충분 |
| 어깨 | 1.23 | 0.81~1.52 (절댓값 3 → 5 N·m) | 충분 (작은 값의 비) |
| **하체** | **1.54** | 1.22~1.76 | **크게 나온다(보수적)**. 순간 최대는 2~4배까지 튄다 |

→ **상체 설계는 영상 한 개로, 하체 설계는 OpenCap(아이폰 2대+)으로.** 영상 하체 토크로 재설계하면 용량을 과하게 잡는다.

## 6. 주의

- **설계 LP 는 모든 프레임을 다 내야 한다** — 튀는 프레임 하나가 토크 용량을 정한다. 넣기 전에 `_id.png` 를 보고, 튀면 `--t0 --t1` 로 자르거나 `--lowpass` 를 낮춘다. 예: 점프(8 Hz)를 넣으니 몸통 총 토크 용량이 155.5 → 289.5 N·m(+86 %).
- 의자·난간·벽·무릎 짚기: 그 힘이 전부 발로 간 것처럼 나와 하체 토크가 틀리다 → `--upper-only`. 로그의 "작용점이 발 밖 %"가 20 % 를 넘으면 의심한다.
- 체공(점프·달리기): 체공 중 몸이 포물선을 안 그리면(관절각 잡음) 그만큼 골반 잔차로 남는다. 점프 8 Hz 에서 잔차 95퍼센타일 17 %BW.
- 트레드밀: 디딘 발 기준 위치가 꺼진다(발이 벨트와 같이 움직임). 영상 하체 토크가 더 나빠진다.
- 손 하중은 구간 내내 손에 붙인 질량이다(물체를 놓는 순간을 모른다).
- 두 발 지지의 배분은 가정이다(정역학적으로 정해지지 않는다). 걷기의 두 발 지지 구간이 가장 덜 맞는다(5-1).
- `smpl_markers.json` 은 한 사람(남, 1.86 m, 83 kg)으로 만들었다. 체형·옷이 많이 다르면 표면 점과 뼈 사이 거리가 달라진다. 우리 사람을 OpenCap 으로 한 번 찍으면 `tools/smpl_calibrate.py --session data/opencap/<그 세션>` 로 다시 만들 수 있다(그 세션 영상을 `video_pose.py` 로 먼저 돌린다).

## 7. 찍는 법

**OpenCap** (권장, 하체까지): app.opencap.ai 안내대로 아이폰·아이패드 2대 이상, 체크보드 보정, 중립 자세 → 동작. 끝나면 세션 번호를 `--opencap` 에 준다.
**영상 한 개** (상체, 빠른 확인):
- 삼각대에 고정. 손에 들고 찍으면 안 된다(카메라가 움직이면 몸이 움직인 것으로 나온다).
- 머리부터 **발끝까지** 처음부터 끝까지 화면 안에. 한 사람만.
- 옆에서 45° 앞쪽. 정면·정옆은 굽힘이 깊이 방향으로 숨는다.
- 2~4 m 거리, 가능하면 60 fps, 밝게. 몸에 붙는 옷.
- 평평한 바닥. 동작 앞뒤로 1초씩 가만히 서 있기(필터 가장자리 0.2 s 를 버린다).

## 8. 파일

| 파일 | 역할 |
| --- | --- |
| `tools/add_motion.py` | 전체 흐름 (입력 → 로봇 뼈대 → 지면반력 → ID → 과제표 → 검사·재설계) |
| `tools/opencap_import.py` | OpenCap 세션(번호·zip·폴더) → `data/opencap/<이름>/` |
| `tools/video_pose.py` | 영상 → SMPL 표면 3D 점 (`pose` 환경, NLF) |
| `tools/video_ik.py` | 3D 점 → 세계 좌표·마커·스케일·IK·디딘 발 기준 위치 |
| `tools/smpl_calibrate.py`, `tools/smpl_markers.json` | 표면 점 → 뼈 마커 보정 (한 번) |
| `tools/grf_estimate.py` | 관절각 → 지면반력 (`--validate` 로 힘판 비교) |
| `tools/video_validate.py` | 영상 한 개 vs OpenCap 비교 |
| `tools/trc_io.py` | `.trc` 읽기·쓰기 |
| `tools/opencap_core/` | OpenCap 일반 모델 LaiUhlrich2022 + Apache-2.0 라이선스 |
| `tools/requirements_pose.txt` | `pose` 환경 |
| `report/fit_check.py` | 최종 A안 정적 검사 (`FINAL_A` 표에 설계 파일) |
| `report/custom/` | 과제표·검사·검증 결과 |

## 9. 다시 돌리기 (검증)

```bash
python tools/opencap_import.py 4d5c3eb1-1a59-4ea1-9178-d3634610561c --name opencap_sample --videos   # 공개 세션 (키 1.86 m)
for t in squat walk_1_25ms DJ STS jump run_2_5ms; do for c in 0 1 2; do
  ~/miniforge3/envs/pose/python.exe -X utf8 tools/video_pose.py --video data/opencap/opencap_sample/raw/Videos/Cam$c/InputMedia/$t/${t}_sync.mp4 \
      --out data/opencap/opencap_sample/pose/${t}_cam$c.npz; done; done                                # 18개, 약 12분
python tools/smpl_calibrate.py --trials squat walk_1_25ms DJ        # → tools/smpl_markers.json
python tools/video_validate.py --trials jump run_2_5ms STS          # 5-2 표 (약 15분)
python tools/grf_estimate.py --validate                             # 5-1 표
```
