# 후배용 따라하기 — 사람 걸음 데이터에서 근육 시너지까지 (그리고 그 다음)

이 문서는 처음 받은 사람이 **명령어를 그대로 쳐서** 시너지 분석까지 끝내고, 상체 같은 다른 부위로 바꿀 때 어디를 고치면 되는지 알기 위한 것이다. 왜 그렇게 했는지는 `NOTES.md`, 전체 현황은 `README.md`, 설계 단계 설명은 `report/인수인계_보행시너지_새근육_2026-09-14.md`에 있다.

맥(Apple Silicon) 기준이다. 한 단계가 끝날 때마다 "확인" 줄의 결과가 나와야 다음으로 간다.

---

## 0. 준비 (최초 1회, 30분)

### 0-1. 코드 받기
```bash
cd ~/Desktop
git clone https://github.com/nauqin/Artificial-muscle-humanoid.git addb_opensim
cd addb_opensim
```
저장소는 비공개라 GitHub 계정을 초대받아야 한다. 연습용 데이터(`data/subject3`, `dryrun`)는 저장소에 들어 있다.

### 0-2. conda 환경 `osim` 만들기
```bash
conda create -n osim -c opensim-org -c conda-forge python=3.11 opensim=4.6 numpy pandas matplotlib ffmpeg -y
conda activate osim
python -m pip install scipy imageio-ffmpeg pyyaml python-pptx pypdf
```
pip는 반드시 `python -m pip`로 친다(다른 환경의 pip가 잡히는 사고 방지). `.b3d` 원본을 직접 변환할 때만 두 번째 환경 `b3d`가 필요한데(README "환경 만들기"), 이 문서에서는 쓰지 않는다.

### 0-3. 환경 확인
```bash
python check_env.py
```
확인: 마지막 줄에 `RESULT: PASS`. 20초쯤 걸린다. FAIL이면 그 위 `[FAIL]` 줄에 원인이 적혀 있다.

---

## 1. 데이터가 어떻게 생겼나 (5분)

```bash
ls data/subject3
```
확인: `subject3_scaled.osim`(피험자 몸에 맞춘 근골격 모델, 근육 80개), `walking1_segment_0_ik.mot`(관절 각도), `walking1_segment_0_grf.mot`(지면반력)이 보인다. 걷기 6시행(walking1~3, walkingTS2~4)과 스쿼트, 앉았다 서기가 있다.

- `_ik.mot`: 시간 × 관절 좌표(라디안). AddBiomechanics가 마커에서 계산한 것.
- `_grf.mot`: 시간 × 발 두 개의 힘·작용점·모멘트.
- 출처: OpenCap 논문(Uhlrich 2023) 데이터를 AddBiomechanics가 처리한 것. 다른 피험자나 다른 동작은 공유 드라이브 `인공근육휴머노이드/Uhlrich2023_Formatted_With_Arm_test`에서 받아 `data/<피험자>/`에 같은 형식으로 둔다.

파일을 직접 열어 보려면:
```bash
head -12 data/subject3/walking1_segment_0_ik.mot | cut -c1-120
```

---

## 2. 한 번에 돌리기: 관절 토크 → 근육 힘 → 시너지 (15분)

```bash
python batch.py --data data --out out/tutorial --subject subject3 --so --synergy --kmax 12
```
무엇이 돌아가나(시행마다):
1. 역동역학(ID): 관절 각도 + 지면반력 → 관절 토크. `run_trial.py`
2. 모멘트암 계산(Muscle Analysis). `run_muscle_analysis.py`
3. 정적 최적화(SO): 관절 토크를 근육 80개에 나눠 근육별 힘·활성도. `run_static_optimization.py`
4. NMF 시너지: 활성도 행렬을 k개 패턴으로 분해, k는 VAF 90 % 기준. `synergy.py`

확인: 끝나면 `out/tutorial/summary.csv`가 생기고, 시행마다 아래 폴더가 생긴다.
```bash
ls out/tutorial/subject3/walking1_segment_0
ls out/tutorial/subject3/walking1_segment_0/so
```
- `walking1_segment_0_id.sto` 관절 토크, `..._id.png` 그래프
- `so/..._StaticOptimization_activation.sto` 근육 활성도, `..._force.sto` 근육 힘
- `so/..._synergy_r.npz` 오른다리 시너지(W: 근육 가중치, C: 시간 패턴), `..._synergy_r.png` 그림

---

## 3. 결과 읽기 (10분)

```bash
column -s, -t < out/tutorial/summary.csv | cut -c1-200
```
한 줄이 시행 하나다. 볼 열:
- `ankle_r`, `knee_r`, `hip_r`: 관절 토크 피크(Nm/kg). 걷기면 발목 1.4~1.9, 무릎 0.5~1.0, 고관절 0.8~1.2가 정상 범위.
- `residual`: 골반 잔차력(체중 대비). 0.1 아래면 지면반력과 운동학이 잘 맞는 것.
- `so_reserve_max`: 근육이 못 내서 보조 액추에이터가 떠맡은 토크. 작을수록 좋다.
- `syn_k_r`, `syn_vaf_r`: 오른다리 시너지 개수와 그때 VAF.
- `id_verdict`, `so_verdict`, `syn_verdict`: PASS / WARN / FAIL. WARN은 `validate.py`의 정상 범위 밖이라는 뜻이고 README "검산 기준" 절에 범위가 적혀 있다.

그림을 열어 본다.
```bash
open out/tutorial/subject3/walking1_segment_0/walking1_segment_0_id.png
open out/tutorial/subject3/walking1_segment_0/so/walking1_segment_0_synergy_r.png
```
시너지 그림은 위가 W(어느 근육이 어느 패턴에 묶이나), 아래가 C(걸음 중 언제 켜지나)다.

파이썬에서 숫자로 꺼내기:
```bash
python - <<'EOF2'
import numpy as np
d = np.load("out/tutorial/subject3/walking1_segment_0/so/walking1_segment_0_synergy_r.npz", allow_pickle=True)
print(list(d.keys()))
print("W", d["W"].shape, "C", d["C"].shape)
EOF2
```

---

## 4. 시너지 하나만 따로 돌리기 (옵션 바꿔 볼 때)

batch가 만든 활성도 파일로 시너지만 다시 뽑는다. k 상한, 좌우, 조용한 근육 제외 문턱을 바꿔 가며 결과가 어떻게 변하는지 본다.
```bash
python synergy.py --activation out/tutorial/subject3/walking1_segment_0/so/walking1_segment_0_StaticOptimization_activation.sto \
    --model data/subject3/subject3_scaled.osim --out out/tutorial/syn_test --name w1 --side r --kmax 10 --min-peak 0.05
```
`--model`은 꼭 준다. 활성도 파일에 팔 토크 액추에이터 열이 섞여 있어서 모델로 근육만 걸러야 한다.

---

## 5. 시너지 수가 곧 근육 수가 아닌 이유 (20분, 선택)

시너지 5개를 당기기만 하는 인공근육 5개로 바꾸면 안 되는 것을 직접 확인한다.
```bash
python report/spanning_check.py          # 시너지 토크 방향이 4자유도를 양의 결합으로 덮나 → 12개 중 1개만
python report/lumped_actuator_test.py    # 새 근육 5개 vs 8개로 측정 토크를 낼 수 있나 (준정적)
```
확인: 5개는 무릎 토크가 모자라고 8개는 낸다. 결과 CSV와 그림은 `report/`에 생긴다.

실측 근전도로 SO를 검산하려면 OpenCap 원자료 `LabValidation_withoutVideos`가 `addb_opensim` 옆 폴더에 있어야 한다(공유 드라이브). 있으면:
```bash
python report/emg_check.py
```

---

## 6. 그 다음 단계 (설계) — 시너지에서 인공근육 8개와 케이블 경로까지

여기서부터는 "사람은 어떻게 걷나"가 아니라 "인공근육을 몇 개, 어떤 비율로, 어디에 붙이나"를 정하는 단계다. 원리는 `report/인수인계_보행시너지_새근육_2026-09-14.md` 10장에 있고, 여기서는 돌리는 순서와 확인할 숫자만 적는다. 전부 `addb_opensim` 폴더에서 `conda activate osim` 한 뒤에 친다.

### 6-0. 먼저 있어야 하는 것
- `out/real/<피험자>/<시행>/<시행>_id.sto`: 7명 걷기 12시행 + subject3 스쿼트의 관절 토크. 2절을 `--subject` 없이 `--out out/real`로 돌리면 생긴다(SO·시너지는 설계에 안 쓰이니 `--no-muscle`로 ID만 돌려도 된다, 10분).
  ```bash
  python batch.py --data data --out out/real --no-muscle
  ```
- `../LabValidation_withoutVideos/<피험자>/sessionMetadata.yaml`: 피험자 체중(정규화에 씀). 공유 드라이브의 OpenCap 원자료를 `addb_opensim` 옆에 둔다.
- 다른 피험자 데이터(`data/subject4` 등)는 공유 드라이브 `Uhlrich2023_Formatted_With_Arm_test`에서 받는다. 연습만 하려면 subject3 하나로도 돌아가지만 숫자는 문서와 달라진다.

### 6-1. 근육 8개의 토크 벡터와 용량 정하기 (10~20분)
```bash
python -u report/design_directions.py --frames 800 --maxfev 1200 --starts 3
```
- 하는 일: 모든 프레임의 관절 토크를 신호 0~1로 낼 수 있는 근육 8개의 토크 방향을, 총 용량이 최소가 되게 찾는다(안쪽 선형계획, 바깥 Powell). `--frames`는 탐색에 쓸 프레임 표본 수, `--maxfev`는 Powell 평가 횟수 상한, `--starts`는 무작위 시작점 수(축 정렬·시너지 시작점은 항상 포함).
- 결과: `report/newmus_V_design_k8.npz`(키 `V` 단위 방향, `c` 용량, `V_nm` = 방향 × 용량), `report/design_directions.csv`, `report/fig_design_directions.png`.
- 확인: 로그 끝의 총용량. 현재 설계는 **394 Nm**이고 축 정렬 522, 시너지 방향 521이 같이 찍힌다. 그림에서 8개 막대가 네 관절에 어떤 비율인지 본다.
- 100평가마다 진행 줄이 찍힌다. 안 찍히면 `-u`를 빼먹은 것.

### 6-2. 토크형 근육으로 걷기 재현 (검증 1, 1분/시행)
케이블로 만들기 전에 "그 토크 벡터 8개로 정말 걷나"를 본다. 사람 근육 80개를 떼고 좌표 액추에이터 4개 × 시너지 제어기로 근육 하나를 흉내 낸 모델이다.
```bash
python report/moco_new_muscles.py --subject subject3 --trial walking1_segment_0 --k 8 --umax 1.0 \
    --reserve-weight 100 --v-from design report/newmus_V_design_k8.npz --tag _tutorial
```
- `--v-from design <npz>`: 6-1의 설계를 쓴다. `--reserve-weight 100`은 꼭 준다(1이면 보조 액추에이터가 싸서 결과가 착시).
- 결과: 로그의 "reserve 최악 N %"와 `report/moco/results_newmus.csv`에 한 줄 추가. 해 파일 `report/moco/newmus_<이름>_solution.sto`.
- 확인: 보조 최악이 1 % 아래면 재현. subject5·10처럼 100 %가 넘게 나오는 시행은 Moco가 발 접촉 순간에 토크를 잘못 내는 노드 때문이다(NOTES "Moco 튐 노드"). 영상은 `python report/newmus_video.py --subject subject3 --trial walking1 --k 8 --tag _tutorial --umax 1.0`.

### 6-3. 케이블 부착점 정하기 (1단계 + 2단계, 5분)
```bash
python -u report/path_design.py
```
- 하는 일: 1단계는 근육 하나씩, 자세 121개에서 케이블 모멘트암 방향이 6-1의 토크 벡터와 나란하도록 부착점을 차분진화 + Powell로 찾는다(근육당 20초). 2단계는 8개 부착점과 힘 배율을 함께 조정해 모든 자세에서 측정 토크 × 1.1을 재현하는 잔차를 줄인다(90초).
- 결과: `report/newmus_paths_k8.json`(근육별 부착점·힘·모멘트암·스트로크), 검증용 모델 `report/moco/newmus_paths_subject3_model.osim`, 그림 `report/fig_path_design.png`.
- 확인: 로그에서 "2단계 끝: 최악 잔차 N %"가 5 % 아래, 각 근육의 필요 힘이 5,000 N 아래, 스트로크가 0.40 아래. 현재 설계는 잔차 4.1 %, 힘 590~3,460 N, 스트로크 0.02~0.19.
- 1단계 결과는 두고 2단계만 다시 돌리려면 `python -u report/path_design.py --refine-from-json`.
- 뼈별 부착 허용 범위는 `report/path_design.py` 맨 위 `BOX`(단위 m). 로봇 뼈대 치수가 오면 여기를 바꾼다.

그림으로 꼭 확인한다(숫자가 좋아도 케이블이 뼈를 뚫을 수 있다).
```bash
python report/plot_path_layout.py          # report/fig_path_layout.png: 서기·걷기·깊은 스쿼트 세 자세
python report/static_check.py              # Moco 없이, 케이블 실제 모멘트암으로 걷기·스쿼트 전 프레임 재현 검사 (0 % 가 정상)
```

### 6-4. 케이블 근육으로 걷기·스쿼트 재현 (검증 2, 1~4분/시행)
```bash
python report/moco_path_muscles.py --subject subject3 --trial walking1_segment_0 --tag _tutorial
python report/moco_path_muscles.py --subject subject3 --trial squats1_segment_0 --t0 2.6 --t1 4.4 --tag _tutorial
python report/moco_path_muscles.py --subject subject3 --trial walking1_segment_0 --t0 0.03 --t1 1.27 --name paths_subject3_walking1_full_tutorial
```
- 하는 일: 6-3의 케이블 모델(좌우 16개)로 MocoInverse. 관절 각도와 지면반력은 사람 것으로 고정하고 케이블 신호만 푼다.
- 결과: 로그에 좌표별 "보조 피크 / ID 피크 / 보조 %"와 오른다리 "새 근육 토크 vs ID 최악 %", `report/moco/results_paths.csv`에 한 줄, 해 `report/moco/paths_<이름>_solution.sto`.
- 확인: "새 근육 토크 vs ID 최악"이 10 % 아래면 재현. 현재 설계는 걷기 2.8~6.3 %, 스쿼트 9.5 %. walkingTS2처럼 100 %가 넘으면 Moco 튐 노드다. 튐 노드를 빼고 다시 계산하려면
  ```bash
  python report/artifact_metric.py paths_subject3_walkingTS2_tutorial report/newmus_paths_k8.json walkingTS2
  ```
- 세 번째 명령은 지지 구간이 아니라 시행 전체(0.03~1.27 s)를 돌린 것으로, GUI 영상용이다. 이때 "토크 vs ID" 지표는 ID 창(0.25~1.12 s) 밖에서는 의미가 없다.
- 입력을 바꿀 때(`--ik --ext --id`) 세 파일이 **같은 모델에서 나온 짝**인지 꼭 확인한다. 다른 출처를 섞으면 발목 오차가 수백 %로 튄다(NOTES "입력은 반드시 같은 모델 짝").

### 6-5. 영상과 GUI 파일 (3분)
```bash
python report/paths_video.py --name paths_subject3_walking1_tutorial --grf data/subject3/walking1_segment_0_grf.mot \
    --id out/real/subject3/walking1_segment_0/walking1_segment_0_id.sto --title "subject3 걷기" --slow 3 --loops 3
python report/paths_video.py --name paths_subject3_squats1_tutorial --grf data/subject3/squats1_segment_0_grf.mot \
    --id out/real/subject3/squats1_segment_0/squats1_segment_0_id.sto --title "subject3 스쿼트" --slow 1.5 --loops 2
python report/gui_export.py paths_subject3_walking1_full_tutorial paths_subject3_squats1_tutorial
```
- `paths_video.py`: `report/video_<이름>.mp4`. 화면은 뼈대 + 케이블(신호가 클수록 진함) + 지면반력 + 신호 막대 + 케이블 토크 vs 측정 토크 + 보조 토크. `--name`은 6-4에서 쓴 Moco 이름(`paths_` 포함)이다.
- `gui_export.py`: OpenSim GUI는 케이블 구동기의 색을 활성도로 안 바꾸므로, 같은 경로를 근육 모델로 바꾼 `report/moco/gui_muscles_model.osim`과 관절 각도 + 활성도 열을 넣은 `report/moco/gui_<이름>.mot`을 만든다. GUI에서 File → Open Model로 .osim을 열고 File → Load Motion으로 .mot을 열면 케이블이 활성도 색으로 움직인다. 뼈 메시는 `report/moco/Geometry`.
- GUI가 안 열리면 NOTES의 "OpenSim 4.6 GUI 맥" 절(Java 지정).

### 6-6. 한 번에 다시 돌리기
6-1부터 6-5까지를 순서대로 치면 된다. 설계가 바뀌면(동작 세트, 허용 상자, 근육 수) 6-1부터, 부착 범위만 바뀌면 6-3부터, 검증 시행만 추가하면 6-4부터 돌린다. 결과 CSV(`results_newmus.csv`, `results_paths.csv`)는 덮어쓰지 않고 줄이 쌓이니 `--tag`로 구분한다.

## 7. 상체로 바꾸려면 어디를 고치나

이 코드는 하지 기준으로 짜였지만 근육·좌표 이름을 본문에 박아 두지 않았다. 바꿀 곳은 정해져 있다.

| 무엇 | 어디 | 어떻게 |
| --- | --- | --- |
| 관절 좌표 목록 | `models.py`의 `LOWER_LIMB_COORDS`, `DIRECTIONS` | 상지 좌표(`arm_flex`, `arm_add`, `arm_rot`, `elbow_flex`, `pro_sup` 등)와 부호 뜻을 추가. 좌표 이름은 모델 .osim 안의 `<Coordinate name=...>`에서 확인 |
| 모멘트암·SO 좌표 | `batch.py --coords a,b,c` | 명령줄에서 상지 좌표를 넘긴다. 기본값은 하지 6개 |
| 시너지 좌우 | `synergy.py --side` | 팔은 `r`/`l`로 한쪽씩 |
| 정상 범위 검산 | `validate.py`의 `BANDS` | 하지 보행 기준이라 상지에는 전부 WARN이 뜬다. 문헌값으로 새 밴드를 넣거나 검산을 끄고 본다 |
| 모델 | `data/<피험자>/*.osim` | **지금 데이터의 팔은 근육이 아니라 토크 액추에이터다.** SO·시너지를 상지에 하려면 팔 근육이 있는 모델(예: MoBL-ARMS, Rajagopal 팔 근육판)로 스케일한 데이터가 필요하다. ID(관절 토크)까지는 지금 모델로도 된다 |
| 설계 단계 | `report/design_directions.py`의 `COORDS`, `report/path_design.py`의 `COORDS`·`BODIES`·`BOX` | 좌표 4개와 뼈 4개, 부착 허용 상자를 상지용으로 바꾼다 |

순서 추천: ① 팔 동작 데이터로 ID만 먼저 돌려 관절 토크 범위를 본다(`batch.py --no-muscle`) → ② 팔 근육 모델을 구한다 → ③ `models.py`와 `--coords`를 바꿔 SO·시너지 → ④ `spanning_check.py`로 근육 수 하한 → ⑤ 설계.

### 7-1. 팔 포함 보행 데이터(Moore2015)로 ID 까지 해 보기

`data/moore2015/` 에 Moore2015(트레드밀 보행, 팔 포함) subject3·subject9 가 변환돼 있다. 기존
`data/subject3` 과 다른 사람이다. 트레드밀 데이터라 옵션이 다르다 — **GRF 도 좌표와 같은 컷오프로
필터**해야 잔차가 5 % 대로 떨어진다(안 하면 16 %).

```bash
conda activate osim
# 한 조각만 (1.6 m/s, 20초)
python run_trial.py --model data/moore2015/subject3/subject3_scaled.osim \
  --ik data/moore2015/subject3/T048_segment_26_ik.mot \
  --grf data/moore2015/subject3/T048_segment_26_grf.mot \
  --out out/moore2015 --name s3_T048_26 --lowpass 15 --grf-lowpass 15
# 전부 (189조각, 20~30분)
python batch.py --data data/moore2015 --out out/moore2015 --no-muscle --lowpass 15 --grf-lowpass 15
```

`out/moore2015/summary.csv` 의 `residual` 이 8 % 이하이고 `t1-t0` 가 5초 이상인 조각만 쓴다.
trial 이름 뜻: subject3 T046/T047/T048, subject9 T025/T026/T027 = 0.8/1.2/1.6 m/s. 0.8 m/s 는
쓸 만한 조각이 없다. 보행 방향은 +z 다. ID 결과(`*_id.sto`)에 `lumbar_*_moment`, `arm_flex_r_moment`,
`elbow_flex_r_moment` 등 상체 토크가 들어 있으므로 상체 설계의 첫 입력으로 쓴다.

새 `.b3d` 를 더 받으면(맥 `b3d` 환경, `nimblephysics` 는 윈도 휠이 없다):

```bash
conda activate b3d
python b3d_to_opensim.py /Volumes/MRLAB/인공근육휴머노이드/subjectN.b3d --out data/moore2015
```

`--walking-only` 는 쓰지 않는다(이름에 walk 가 없어 다 걸러진다). 로그의 `[GRF 보정]` 줄은 떠 있는
발의 잡음 힘을 0 으로 만들고 발 밖에 찍힌 COP 를 발 아래로 옮긴 프레임 수다(NOTES 2026-10-30 참고).

---

## 8. 막혔을 때

- `[환경 오류]`로 시작하는 메시지: `conda activate osim`을 안 한 것. 메시지에 적힌 대로.
- `batch.py`가 어떤 시행을 건너뛴다: 이미 결과가 있어서다. 다시 하려면 `--force`.
- ID 잔차가 크다(`residual` > 0.2): 지면반력이 없는 구간이 섞인 것. `run_trial.py`가 자동으로 자르지만, 로그의 "trial … 프레임 (…/… 사용)" 줄을 본다.
- 시너지 k가 시행마다 다르다: 정상이다. `--criterion global`로 바꾸면 보수적으로 나온다. NOTES "시너지 개수" 절.
- 그 밖의 모든 것: `NOTES.md`를 먼저 검색한다. 같은 실수를 이미 한 번 했을 확률이 높다.
