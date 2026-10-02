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

## 6. 그 다음 단계 (설계) — 명령만

시너지 분석이 끝난 뒤의 설계 단계는 `report/인수인계_보행시너지_새근육_2026-09-14.md` 10~11장을 읽고 돌린다. 입력은 `out/real/`(7명 걷기 ID)이 필요해서, 2절을 `--out out/real --subject` 없이 7명 전부로 돌린 뒤에 가능하다.
```bash
python report/design_directions.py --frames 800 --maxfev 1200 --starts 3     # 근육 8개 토크 벡터·용량 (10~20분)
python report/path_design.py                                                # 케이블 부착점 1단계 + 2단계 (5분)
python report/plot_path_layout.py                                           # 경로 그림
python report/moco_path_muscles.py --subject subject3 --trial walking1_segment_0 --tag _test   # 케이블 근육으로 걷기 재현 (1분)
python report/paths_video.py --name paths_subject3_walking1_test --grf data/subject3/walking1_segment_0_grf.mot \
    --id out/real/subject3/walking1_segment_0/walking1_segment_0_id.sto --title "subject3 걷기" --slow 3 --loops 3
```

---

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

---

## 8. 막혔을 때

- `[환경 오류]`로 시작하는 메시지: `conda activate osim`을 안 한 것. 메시지에 적힌 대로.
- `batch.py`가 어떤 시행을 건너뛴다: 이미 결과가 있어서다. 다시 하려면 `--force`.
- ID 잔차가 크다(`residual` > 0.2): 지면반력이 없는 구간이 섞인 것. `run_trial.py`가 자동으로 자르지만, 로그의 "trial … 프레임 (…/… 사용)" 줄을 본다.
- 시너지 k가 시행마다 다르다: 정상이다. `--criterion global`로 바꾸면 보수적으로 나온다. NOTES "시너지 개수" 절.
- 그 밖의 모든 것: `NOTES.md`를 먼저 검색한다. 같은 실수를 이미 한 번 했을 확률이 높다.
