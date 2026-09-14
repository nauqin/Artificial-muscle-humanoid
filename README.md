# AddBiomechanics → OpenSim 보행 분석 파이프라인

마커·힘판 데이터에서 **관절 모멘트 → 근육별 힘 → 근육 시너지**까지 자동으로
계산한다. AddBiomechanics는 스케일된 `.osim`과 IK 결과 `.mot`을 함께 배포하므로
Scale/IK는 직접 하지 않는다.

각 단계마다 결과가 물리적으로 말이 되는지 자동으로 검산한다. 계산은 입력이 틀려도
불평 없이 숫자를 내놓기 때문에, 검산이 없으면 조용히 틀린 결과로 몇 달을 갈 수 있다.

**현재 상태** — macOS(Apple Silicon)에서 재구축했고 드라이런이 예전 실측값을 재현한다.
실제 데이터 Uhlrich2023은 **7명 42 trial**로 확정했다 — subject6·9는 AddBiomechanics
자체 잔차가 나빠 제외, subject10은 `manualReview` 플래그를 열어 복구. 시너지는
`global` 기준 **5개**(4~7)이고 보행 수정 조건으로 개수는 바뀌지 않는다 → 결과 절.
다음은 개수가 아니라 구성·타이밍 비교, 그리고 시너지당 새 근육 하나로 걷기다.

예전 윈도우 환경에서는 `.b3d`를 읽는 nimblephysics에 Windows 휠이 없어 변환만
리눅스에서 했지만, **맥에서는 그럴 필요가 없다.** 변환도 분석도 이 기계 한 대에서 끝난다.

```
마커·힘판  →  IK(모션)  →  ID(관절 모멘트)  →  SO(근육별 힘)  →  NMF(근육 시너지)
                              ↑                     ↑                  ↑
                         MuscleAnalysis(모멘트암) ───┴── 검증에 사용 ───┘
```

설계 근거, 실측 수치, 밟은 함정은 NOTES.md에 있다. 값이 이상하면 그쪽을 본다.

---

## 빠른 시작

```bash
cd ~/Desktop/scHooL/의료로봇연구실/addb_opensim
conda activate osim
python check_env.py          # RESULT: PASS 확인 (20초)
python batch.py --data dryrun --out out/batchtest --force --so --synergy
```

---

## 실행 방법

### 환경

환경이 둘이다. 서로 numpy 요구가 달라 분리했다.

| 환경 | 무엇 | 언제 |
| --- | --- | --- |
| `osim` | python 3.11 + OpenSim 4.6 + numpy 2.4 + pandas + matplotlib | 분석 전부 |
| `b3d` | python 3.11 + nimblephysics 0.10.52.1 + numpy 1.x | `.b3d` 변환만 |

```bash
conda activate osim        # 여러 명령을 이어서 칠 때
./run.sh check_env.py      # 어느 폴더에서든 한 방만 돌릴 때
```

`run.sh`는 `conda activate osim` 후 이 폴더로 옮겨 스크립트를 돌린다.

### 배치 — 평소 쓰는 것

`data/` 아래 각 하위 폴더를 subject 하나로 보고 전수 처리한다. 폴더 구조는 탐색으로
처리하므로 AddBiomechanics가 어떤 레이아웃으로 주든 신경 쓰지 않아도 된다.

```bash
python batch.py --data data --out out --so --synergy
```

| 옵션 | 기본값 | 뜻 |
| --- | --- | --- |
| `--data <폴더>` | `data` | subject 폴더들의 상위 폴더 |
| `--out <폴더>` | `out` | 결과를 쌓을 곳 |
| `--subject <이름>` | 전체 | 특정 subject만. 여러 번 지정 가능 |
| `--model <파일>` | 자동 선택 | 모델을 강제 지정 |
| `--force` | 끔 | 기존 결과도 다시 계산 (기본은 건너뜀) |
| `--no-muscle` | 끔 | MuscleAnalysis 생략, ID만 |
| `--so` | 끔 | Static Optimization으로 근육별 힘까지 |
| `--synergy` | 끔 | NMF 시너지 분석 (`--so` 포함) |
| `--concat` | 끔 | 같은 조건의 trial을 이어붙여 시너지를 한 번 더 분해 |
| `--side r\|l\|both` | `both` | 시너지를 분석할 다리 |
| `--kmax <N>` | `10` | 시너지 개수 탐색 상한 |
| `--criterion dual\|global` | `dual` | 시너지 개수 선택 기준 (아래 검산 절 참고) |
| `--seed <N>` | `0` | NMF 초기값 난수 시드 |
| `--coords a,b,c` | 발목·무릎·고관절 좌우 6개 | 모멘트암을 뽑을 좌표 |
| `--start` / `--end` | 자동 | 시간 구간. 보통 건드릴 필요 없다 |
| `--lowpass <Hz>` | `6.0` | 좌표 저역통과 컷오프 |

### 개별 실행

```bash
python run_trial.py               --model <osim> --ik <ik.mot> --grf <grf.mot> --out <폴더> --name <이름>
python run_muscle_analysis.py     --model <osim> --ik <ik.mot> --out <폴더> --name <이름>
python run_static_optimization.py --model <osim> --ik <ik.mot> --grf <grf.mot> --id <id.sto> --ma-dir <ma폴더> --out <폴더> --name <이름>
python synergy.py                 --activation <*_StaticOptimization_activation.sto> --model <osim> --out <폴더> --name <이름> --side both
python animate.py                 --model <osim> --ik <ik.mot>
python gait_video.py              --model <osim> --ik <ik.mot> --grf <grf.mot> --out walk.mp4
```

`synergy.py`의 `--model`은 선택이 아니다 — 활성도 파일에는 근육 말고 **팔 토크
액추에이터**(`shoulder_flex_r` 등)도 들어있어서, 모델로 걸러내지 않으면 근육이 아닌
시너지 성분이 만들어진다. 그 밖에 `--min-peak`(기본 0.05, 이보다 조용한 근육은
분석에서 제외)을 받는다. `run_static_optimization.py`는 `--ma-dir`을 생략하면
모멘트암을 직접 계산한다.

---

## 결과 읽는 법

### 폴더

```
out/
  summary.csv                          ← 전체 요약. 여기부터 본다
  <subject>/<trial>/
    <trial>_id.sto                     관절 모멘트 (본 결과)
    <trial>_id.png                     모멘트·GRF 그래프, 정상범위 회색 띠
    <trial>_id_setup.xml               재현용 설정. GUI에서 열어 확인 가능
    <trial>_external_loads.xml         GRF를 어느 발에 붙였는지
    ma/            <trial>_MuscleAnalysis_MomentArm_<좌표>.sto 등
    so/            <trial>_StaticOptimization_force.sto      근육별 힘 (N)
                   <trial>_StaticOptimization_activation.sto 활성도 (0~1)
                   <trial>_reserve_actuators.xml
                   <trial>_synergy_r.png / _l.png            시너지 그림
                   <trial>_synergy_r.npz / _l.npz            W, C, VAF 원자료
  <subject>/_concat_<조건>/             --concat 으로 이어붙인 분해
```

### `summary.csv`

| 열 | 뜻 |
| --- | --- |
| `subject` `trial` `model` | 어느 피험자·걸음·모델인지 |
| `mass_kg` | **선택된 모델의 체중.** 예상과 다르면 모델을 잘못 고른 것 |
| `t0` `t1` | 자동으로 정해진 계산 구간 (초) |
| `ankle_r/l` `knee_r/l` `hip_r/l` | 관절 모멘트 피크 (Nm/kg) |
| `grf_peak` | 수직 지면반력 피크 (체중 배수) |
| `residual` | 골반 잔차력, 피크 GRF 대비 % — **동역학 일관성 지표** |
| `cop_dist` | COP와 배정된 발의 평균 거리 (m) |
| `so_recon_max` | SO 근육 힘이 ID 모멘트를 재현한 오차, 최대 % |
| `so_reserve_max` | 보조 액추에이터가 떠맡은 비율, 최대 % |
| `syn_k_r` `syn_k_l` | 좌·우 시너지 개수 (`--criterion` 기준) |
| `syn_vaf_r` `syn_vaf_l` | 그때의 전체 VAF |
| `syn_kglobal_*` `syn_kdual_*` | 두 기준 각각의 개수. 못 찾았으면 비어 있다 |
| `id_verdict` `ma_verdict` `so_verdict` `syn_verdict` | 단계별 판정 |
| `error` | 실패했을 때 이유. trial 하나가 죽어도 나머지는 계속 돈다 |

### 첫 실행 때 로그에서 볼 세 줄

조용히 틀릴 수 있는 지점은 셋뿐이다. 나머지는 검산기가 본다.

```
모델 : dryrun/Rajagopal/Scale/subject_scaled_walk.osim  (체중 85.07 kg)   ← 그 사람 몸이 맞나
           r -> calcn_r   calcn_r=0.137  calcn_l=0.483   +0.346           ← 좌우가 안 바뀌었나
  trial: 1개  ['walk']                                                    ← 엉뚱한 파일끼리 안 묶였나
```

`[WARN] 짝을 확정할 수 없어 건너뜀`이 보이면 그 GRF는 **결과에서 빠진다.** 엉뚱한
IK와 묶어 계산하는 것보다 빠뜨리는 쪽이 낫기 때문이다.

---

## 검산 기준

판정은 세 단계다. 목적은 **자릿수**가 맞는지 보는 것이다.

- `PASS` 정상범위 안
- `WARN` 밴드를 벗어났지만 자릿수는 맞음 → 피험자·속도 차이일 수 있음
- `FAIL` 상한의 2배 초과 또는 하한의 절반 미만 → 파이프라인 오류를 의심

### ID — `validate.py`의 `BANDS`

| 항목 | 정상범위 | 벗어나면 의심할 것 |
| --- | --- | --- |
| 발목 저굴 모멘트 피크 | 1.2 ~ 1.6 Nm/kg | GRF-발 매핑 / 좌표계 / 단위 |
| 무릎 모멘트 피크 | 0.4 ~ 0.8 Nm/kg | 〃 |
| 고관절 모멘트 피크 | 0.7 ~ 1.2 Nm/kg | 〃 |
| 수직 GRF 피크 | 0.9 ~ 1.3 BW | GRF 단위(N vs BW) / 체중 |
| COP-발 거리 | < 0.30 m | COP 단위(mm vs m) / 좌우 반전 |
| 골반 잔차력 | < 5 % of peak GRF | 모델-모션 부정합 / 동역학 불일치 |

**한쪽 관절 모멘트만 0에 가깝게 나오면 매핑보다 먼저 `t0`~`t1` 구간 길이를 본다.**
구간이 한 걸음(약 0.8초)보다 짧으면 그쪽 다리의 push-off가 통째로 잘렸을 수 있다.
검산기가 이 경우를 알아서 구분해 안내한다 — NOTES.md의 "midstance 골이 지지 구간을
반으로 쪼갰다" 참고.

좌우를 따로 본다. 한쪽만 틀리면 매핑 문제다.

### 모멘트암 — `run_muscle_analysis.py`의 `ANTAGONISTS`

부호 규약을 하드코딩하지 않는다. 모델마다 축 방향이 다를 수 있어서다. 대신
**길항근 쌍의 부호가 반대인지**를 무리 전체로 보고, 크기가 0.02~0.08 m인지는
**대표 근육**(soleus/gasmed/gaslat, tibant, vasti, 햄스트링, psoas/iliacus, glmax)
으로만 본다. `fdl`·`perbrev` 같은 작은 원위근은 모멘트암이 원래 1 cm 남짓이라
함께 재면 밴드를 벗어난다. 주동근끼리 부호가 갈리거나 길항근과 부호가 같으면 FAIL.

### SO — `validate.validate_so`

SO가 배분한 근육 힘이 진짜로 ID 관절 모멘트를 만들어내는지 좌표마다 확인한다.

```
Σ(근육힘 × 모멘트암) + 보조액추에이터  =?  ID 모멘트
```

| 항목 | 정상범위 | 벗어나면 의심할 것 |
| --- | --- | --- |
| `recon/<좌표>` 재현 오차 (RMS / 피크모멘트) | < 5 % | SO 미수렴 / 시간축 불일치 / 근육 누락 |
| `reserve/<좌표>` 보조가 떠맡은 비율 | < 10 % | 근육만으로 모멘트 부족 — 모델·구간·GRF 확인 |

**두 번째가 진짜 품질 지표다.** 보조 액추에이터는 어떤 모멘트든 만들어낼 수 있어서
재현 오차는 언제나 0에 가깝게 나올 수 있다. 보조가 많이 개입했다면 숫자가 맞아도
"근육이 낸 힘"으로 해석하면 안 된다.

### 시너지 개수 — `synergy.py`

VAF(설명된 분산)로 고른다. 두 기준을 나란히 보고한다.

- **global** — 전체 VAF ≥ 90 %
- **dual** (기본) — 전체 VAF ≥ 90 % **그리고** 모든 근육의 개별 VAF ≥ 75 %

전체 VAF만 보면 활성도가 큰 근육 몇 개에 가려 조용한 근육이 전혀 설명되지 않아도
통과한다. `--kmax`까지 기준을 못 채우면 PASS를 주지 않고 경고를 찍는다. **표 위의
판정을 반드시 볼 것.**

**실제 데이터에서는 `global`을 본 결과로 쓴다.** `global`은 문헌값(4~6)과 일치하지만
`dual`은 흔들리고 일부는 상한에서도 못 채운다. 조건별로 trial을 이어붙여 표본을
3~4배로 늘려도(`--concat`) `dual`은 내려가지 않았다 — 표본 부족이 아니라는 뜻이다.
`dual`은 "이 근육은 개별 해석 금지"를 알려주는 진단 지표로 남긴다. 전말은 NOTES.md에.

---

## 파일

| 파일 | 역할 |
| --- | --- |
| `run.sh` | 어느 폴더에서든 스크립트 하나를 `osim` 환경 안에서 실행 |
| `check_env.py` | 환경 스모크 테스트 (모델 로드, 체중, 좌표·바디 존재, 그림 저장) |
| `b3d_to_opensim.py` | `.b3d` → `.osim`/`.mot` 변환 (`b3d` 환경에서) |
| `addb_io.py` | 파일 탐색·분류, `.mot` 파싱, GRF-발 자동 매핑, 지지 구간 탐지 |
| `run_trial.py` | trial 1개: 정합 확인 → 매핑 → 구간 결정 → ID → 검산 → 그래프 |
| `run_muscle_analysis.py` | 모멘트암 추출 + 길항근 부호 검산 |
| `run_static_optimization.py` | 근육별 힘 추정 + 보조 액추에이터 생성 |
| `synergy.py` | NMF 시너지 분해, VAF로 개수 선택, 표·그림 |
| `validate.py` | ID·SO 검산 밴드와 리포트, PNG |
| `batch.py` | subject × trial 전수 처리, `summary.csv` 집계 |
| `check_foot_height.py` | 발·골반 높이로 바닥 관통·단위 확인 |
| `models.py` | 모델별 설정 — 좌표·ROM·대표 모멘트암 (하지 배치 과제) |
| `actuator_model.py` | 인공근육 모듈 규격과 개념 모델 (하지 배치 과제) |
| `screen_modules.py` | Phase 0 — 모듈 개수 하한 스크리닝 (하지 배치 과제) |
| `animate.py` | Simbody 뷰어로 IK 모션 재생 — 눈으로 볼 때 (지면이 정확히 y=0) |
| `gait_video.py` | 보행 영상 파일 (하지 스틱피겨 + GRF 벡터). 화면 없이 돌아간다 |
| `dryrun/Rajagopal/` | 번들 예제 복사본. 원본 리소스는 건드리지 않는다 |
| `dryrun_addbfmt/` | 같은 예제를 **변환기 출력 형식**으로 다시 쓴 것. 회귀 시험용 |
| `report/make_ppt.py` | 9/8 요약 + subject10 병합 → 7명 확정본, 그림, 현황 PPT |
| `report/spanning_check.py` | 시너지 토크 방향벡터의 양의 스패닝 검사 (새 근육 k개 가능성 첫 관문) |
| `report/lumped_actuator_test.py` | 새 근육 k개(토크 비율 고정)로 측정 관절 모멘트를 재현하나 — 준정적 NNLS, 5 vs 8 |
| `report/moco_synergy_inverse.py` | 근육 80개에 SynergyController(신호 k개) 물리고 MocoInverse로 걸음 재현 — 근육 동역학 포함 |
| `report/moco_new_muscles.py` | 사람 근육 전부 떼고 **새 근육 k개/다리**(토크 비율 고정)로 걸음 재현 — 8개면 된다, 5개는 무릎이 안 된다. `--v-from SUBJ TRIAL|ALL`(방향 고정 교차 검증), `--reserve-weight 100`(비용 착시 제거) |
| `report/plot_new_muscles.py` | 위 결과 그림 (신호·토크·보조) |
| `report/plot_newmuscles_cross.py` | 방향 V 출처별 교차 검증 막대그림 (`fig_newmuscles_cross.png`) |
| `report/newmus_video.py` | 새 근육으로 걷는 모습 mp4 — 스틱피겨 + GRF + 신호 막대 + 토크 추적 (GUI 없이) |
| `report/emg_check.py` | 실측 EMG(LabValidation) 대 SO 활성도·시너지 검산 — soleus 0.78, tibant 0.10 |
| `raw_b3d/` | 내려받은 `.b3d`를 여기에 둔다 |
| `data/` | AddBiomechanics 다운로드를 여기에 푼다 |
| `out/` | 결과 |

---

## 하지 인공근육 배치 과제

`PLAN_lowerlimb_muscle_placement.md` 의 작업. **기존 보행 파이프라인은 건드리지 않고**
새 스크립트만 얹는다. Phase 마다 드라이런 회귀를 확인한다.

### Phase 0 — 모듈 개수 하한 스크리닝

```bash
python screen_modules.py --id-dir out/real --data data \
    --exclude-subject subject6 --exclude-subject subject9 \
    --out out/lowerlimb/phase0_no69 --bound both   # 7명 봉투 (subject10 포함)
```

| 옵션 | 뜻 |
| --- | --- |
| `--id-dir <폴더>` | 배치 결과. `*_id.sto` 를 훑어 **6 DOF 전부** 요구 토크를 뽑는다 |
| `--summary <csv>` | `summary.csv` 로도 되지만 발목·무릎·고관절 3개뿐이고 방향 구분이 없다 |
| `--data <폴더>` | subject 모델 위치. 체중을 알아야 `Nm/kg` 이 된다 |
| `--bound low\|high\|both` | 모듈 규격의 하한(보수적)·상한. 둘 차이가 곧 기계연에 물어볼 것 |
| `--exclude-subject <이름>` | 봉투에서 뺀다. 잔차 높은 subject6·9 의 영향을 보는 데 쓴다 |

결과는 `out/lowerlimb/phase0/summary.md` 와 `screening.csv`.

---

## 걷는 모습 보기

둘 다 하지 골격을 지면 y=0 기준으로 그린다.

### 창을 띄워 눈으로

```bash
python animate.py --model dryrun/Rajagopal/Scale/subject_scaled_walk.osim \
                  --ik dryrun/Rajagopal/IK/results_walk/ik_output_walk.mot
```

Simbody 뷰어가 뜬다. 재생 속도·일시정지는 창 안에서 조절하고, 창을 닫으면 끝난다.
**OpenSim GUI의 체크무늬 지면은 믿지 말 것** — 같은 모션인데 GUI에서는 몸이 지면
아래로 보인 적이 있다. 바닥 기준이 필요하면 이쪽을 쓴다.

### 파일로 (발표자료용)

```bash
python gait_video.py --model <osim> --ik <ik.mot> --grf <grf.mot> --out out/walk.mp4
```

시상면 스틱피겨에 **지면반력 벡터**(초록)와 그 순간의 vGRF(BW)를 함께 그린다.
오른다리 빨강, 왼다리 파랑. 카메라는 골반을 따라간다.

| 옵션 | 기본값 | 뜻 |
| --- | --- | --- |
| `--out <파일>` | — | `.mp4` 또는 `.gif` (ffmpeg 없으면 GIF로 넘어간다) |
| `--grf <파일>` | 없음 | 있으면 지면반력 화살표를 함께 그린다 |
| `--fps <N>` | `30` | |
| `--view sagittal\|frontal` | `sagittal` | 시상면 / 관상면 |
| `--no-follow` | 끔 | 카메라 고정 — 전체 궤적을 한 화면에 |
| `--start` / `--end` | 전체 | 구간 |

---

## 데이터 넣기 — `.b3d` 변환

AddBiomechanics 공개 데이터셋은 subject 하나당 **`.b3d` 바이너리 하나**로 배포된다.
스케일된 모델과 모든 trial의 관절 각도·지면반력이 다 들어있다.

### 1. 받기 (수동)

공유 드라이브에 정리본이 있다 — `공유 문서함 > addb_dataset_publication > {train,test} > With_Arm`.
subject 폴더마다 `.b3d` 하나다. 필요한 subject 폴더를 받아 `raw_b3d/`에 둔다.

**`With_Arm`을 쓴다.** 팔을 포함한 Rajagopal 모델(근육 80개)이라 SO·시너지까지 된다.

`app.addbiomechanics.org`에서 직접 받아도 되지만 **스크립트로는 못 받는다** —
익명 S3 리스팅이 막혀 있는 것을 확인했다.

### 2. 변환

```bash
conda activate b3d
python b3d_to_opensim.py raw_b3d --out data --walking-only
```

| 옵션 | 뜻 |
| --- | --- |
| `--out <폴더>` | 출력 위치 (기본 `data`) |
| `--inspect` | 변환하지 않고 안을 훑기만 — trial 목록, 근육 수, GRF 양호 프레임 수 |
| `--walking-only` | 보행/달리기로 표시된 trial만 |
| `--pass <N>` | processing pass 지정 (기본: `dynamics` 자동 선택) |
| `--keep-all-frames` | GRF 없는 구간도 그대로 (기본은 잘라냄) |
| `--min-frames <N>` | 이보다 짧은 GRF 구간은 건너뜀 (기본 50) |
| `--allow-manual-review` | `manualReview` 플래그를 무시 — 물리 판정은 우리 검산기가 한다 (아래 subject 선별 참고) |

인자로 폴더를 주면 `.b3d`를 재귀로 찾는다. 나오는 구조는 이렇다.

```
data/<subject>/
  <subject>_scaled.osim     스케일·질량보정된 모델
  <trial>_ik.mot            관절 각도 (rad, inDegrees=no)
  <trial>_grf.mot           발별 힘 / COP / 자유모멘트
```

**변환 로그에서 볼 두 줄.**

```
모델 근육 80개 — SO/시너지 가능        ← 없다고 나오면 ID까지만 된다
trial walk1  1204 프레임 (0.30~4.31s, 1204/1500 사용)   ← GRF 없는 구간이 잘린 양
```

기본값은 `missingGRFReason`이 정상인 가장 긴 연속 구간만 남기는 것이다. GRF가 없는
프레임을 ID에 넣으면 체중 전체가 골반 잔차로 잡혀 결과가 통째로 무의미해진다.

### 환경 만들기 (최초 1회)

```bash
conda create -n b3d python=3.11 "numpy<2" -y
conda activate b3d
pip install "nimblephysics==0.10.52.1"
```

nimblephysics는 PyPI에 **소스 배포본(sdist)이 아예 없다.** 미리 빌드된 휠뿐이라
인터프리터에 맞는 휠이 없으면 컴파일 폴백도 없이 그냥 실패한다.

|  | 되는 파이썬 (0.10.52.1 기준) |
| --- | --- |
| **macOS 애플실리콘** | **3.9, 3.11 (`universal2` 휠)** |
| macOS 인텔 | 3.8 ~ 3.11 |
| 리눅스 x86_64 | 3.8 ~ 3.12 (glibc 2.27 이상) |
| 리눅스 aarch64 | 3.8 ~ 3.11 |
| Windows | **없음** |

**3.13 이상은 어느 플랫폼에도 휠이 없다.** 버전 고정(`==0.10.52.1`)은 권한다 —
최신 0.10.52.2에는 macOS cp39 휠 하나만 올라와 있다.

---

## 결과 — Uhlrich2023 (7명 42 trial)

AddBiomechanics 공개 데이터셋 **`Uhlrich2023_Formatted_With_Arm`** (공유 드라이브
`test/With_Arm`). 팔을 포함한 Rajagopal 모델(37 DOF, 근육 80개)로 처리돼 있어
SO·시너지까지 전부 가능하다. 접지 바디는 `calcn_r` / `calcn_l`. 무릎 부하를 줄이는
보행 수정 중재 연구라 같은 사람이 평상(`walking1~4`)과 수정(`walkingTS1~4`)을 둘 다
걸었다. 평지 보행이라 trial 하나가 한 걸음(1.2~1.8 s)이다.

```bash
conda activate b3d
python b3d_to_opensim.py raw_b3d --inspect > out/inspect_uhlrich2023.txt
python b3d_to_opensim.py raw_b3d --out data --walking-only --allow-manual-review
conda activate osim
python batch.py --data data --out out/real --so --synergy --kmax 12
python report/make_ppt.py        # 7명으로 병합 + 그림 + PPT
```

### 데이터 출처 — OpenCap 논문

Uhlrich, Falisse, Kidziński 외, *OpenCap: Human movement dynamics from smartphone videos*,
PLOS Comput Biol 2023 (`인공근육휴머노이드/…/journal.pcbi.1011462.pdf`). 스마트폰 영상
도구를 검증하려고 찍은 실험실 데이터가 AddBiomechanics로 재처리된 것이다.

| | 논문 (p16, p18) | 우리 |
| --- | --- | --- |
| 피험자 | 건강한 성인 10명, 여6·남4, 27.7±3.8세, 69.2±11.6 kg (59.0~92.9), 1.74 m | subject2~11, 질량 범위 일치 |
| `walkingTS` | **trunk sway — 몸통을 디딤발 쪽으로 옆으로 기울여 걷기**, 평균 15° 더 | `lumbar_bending` 진폭 24° → 47° (2.0×), 골반·다리는 불변 |
| 목적 변수 | 초기 입각기 무릎 내전 모멘트(KAM)·내측 접촉력(MCF). trunk sway로 21~46 % 감소(10명 중 9명 감소) | 우리 모델 무릎은 굴곡 1 DOF라 **KAM 미계산** |
| 마커·힘판 | 100 Hz 마커 31+20개, Bertec 힘판 3장 2000 Hz | — |
| 필터 | 걷기 6 Hz (4차 영위상 Butterworth) | `LOWPASS = 6.0` 동일 |
| 근육 힘 | OpenSim 4.3 **Static Optimization** → Joint Reaction Analysis | 우리도 SO |
| EMG | 16개 하지 근육, **SimTK 원본 데이터셋에 공유** (simtk.org/projects/opencap) | `.b3d`에는 없음 (`emgSignals` 비어 있음) |

다른 수정 조건도 같은 틀이다: `squatsAsym`(왼발 힘 줄이기), `STSweakLegs`(몸통 앞으로
숙이며 일어서기), `DJAsym`(착지 때 왼발 힘 줄이기). 훈련이 아니라 **그 자리에서 지시받고
한 번 한 것**이다.

### subject 선별 — 이름표가 아니라 물리로

`--inspect`가 찍는 **AddBiomechanics 자체 골반 잔차**(dynamics pass, 프레임 평균)를
체중으로 정규화하면 subject가 셋으로 갈린다.

| subject | AddB 잔차 (%BW) | 우리 잔차 (% peak GRF) | 판정 |
| --- | --- | --- | --- |
| 5 · 4 · 11 · 3 | 0.3 ~ 1.2 | 2.7 ~ 3.5 | 사용 |
| **10** | **0.5 ~ 0.6** | 4.2 | **사용** — `manualReview`였지만 물리 품질은 subject3보다 좋다 |
| 7 · 8 | 0.9 ~ 4.3 / 1.0 ~ 28 (trial 한둘만 나쁨) | 3.4 / 4.3 | 사용 |
| 2 | 5 ~ 13 | — | 제외 — `manualReview`인 데다 잔차도 나쁘다 |
| **6** | **13 ~ 18** | 8.1 | **제외** |
| **9** | **16 ~ 30** | 7.5 | **제외** |

- **subject6·9의 높은 잔차는 우리 계산 문제가 아니다.** AddBiomechanics 스스로도
  이 둘의 동역학을 못 맞췄다 — 다른 사람의 20~40배, 순간 최대 2100~3300 N(체중의
  3~5배). 우리 검산기가 이 둘만 골라낸 것이 정상 작동의 증거다. SO 근육 힘도 신뢰할
  수 없으므로 뺀다.
- **`manualReview`는 물리 판정이 아니다.** 프레임별 `torqueDiscrepancy` 같은 판정과
  달리 trial 전체에 균일하게 찍히는 "사람이 검토하지 않음" 표시다. subject10은 이
  플래그가 붙었는데 AddB 잔차가 최상급이고, subject2는 같은 플래그인데 잔차도 실제로
  나쁘다. 그래서 `--allow-manual-review`로 플래그를 열고 우리 검산기로 거른다.
  두 사람을 같은 이유로 묶어 제외했던 예전 판단이 틀렸다.

### 전체 결과

| 단계 | 결과 |
| --- | --- |
| 모멘트암 | 36 PASS / 6 WARN |
| SO | 39 PASS / 3 WARN — `so_reserve_max` 최대 **2.06 %** |
| 시너지 | 40 PASS / 2 WARN |
| ID | 11 PASS / 31 WARN / **0 FAIL** — 잔차 평균 3.49 %, 최대 5.82 % |

**`so_reserve_max`가 전부 2.1 % 미만이다.** 근육 힘을 "근육이 낸 힘"으로 해석해도
된다. ID WARN은 대부분 `walkingTS`의 무릎 모멘트가 밴드 하단 아래로 내려간 것 —
중재 효과가 잡힌 것이지 오류가 아니다(아래 "WARN을 어떻게 읽었나").

### 시너지 개수: 5개, 조건 차이 없음

`global` 기준(전체 VAF ≥ 90 %), 84개 분해(42 trial × 좌우):

| k | 4 | **5** | 6 | 7 |
| --- | --- | --- | --- | --- |
| 분해 수 | 5 | **56** | 22 | 1 |

평균 5.23. 문헌의 보행 시너지 4~6개와 일치한다. walking 5.26 vs walkingTS 5.19,
subject별 차이 −0.5 ~ +0.3으로 방향이 뒤섞인다. **보행 수정은 제어 차원을 바꾸지
않는다.** n=7에 k가 정수값이라 개수로는 아무것도 주장할 수 없다 — 비교는 W·C에서.

k=5에서 subject를 넘나들며 같은 다섯 묶음이 반복된다.

| 시너지 | 묶이는 근육 |
| --- | --- |
| 고관절 외전·골반 안정 | glmed1-3, glmin1-3, tfl |
| push-off + 고관절 굴곡 | soleus, psoas, iliacus, fhl, perlong, tibpost |
| 고관절 신전 | glmax1-3, bflh, semimem, semiten, addmagIsch |
| 배굴 | tibant, edl, ehl (± vasti) |
| 내전 | addlong, addbrev, addmagProx/Mid, grac, piri |

### "예전엔 8개가 안정적이라 했는데?"

같은 데이터의 두 기준이다. `global`(전체 VAF ≥ 90 %)이 5, `dual`(+ 모든 근육 개별
VAF ≥ 75 %)이 평균 8.0(5~12)이다. 추가 3개는 `gaslat`/`gasmed`, `glmed1/2/3`처럼
모멘트암이 겹치는 근육 사이에 SO가 임의로 나눈 부하를 재현하는 데 쓰인다 — 근육
하나짜리 성분이 나오는 것이 증거고, 표본을 4배로 늘려도 안 내려간다. **본 결과는 5.**
`dual`은 "이 근육은 개별 해석 금지" 진단 지표다. 하드웨어 설계에 8을 쓰면 SO의 동전
던지기를 기계에 새기게 된다.

### WARN을 어떻게 읽었나

`walkingTS` trial에서만 무릎 모멘트가 밴드 하단(0.4 Nm/kg) 아래로 내려갔다.

| trial | knee_r | knee_l |
| --- | --- | --- |
| walking1/2/3 (평상) | 0.409 / 0.441 / 0.415 | 0.477 / 0.472 / 0.597 |
| walkingTS2/3/4 (수정) | 0.347 / 0.296 / 0.287 | 0.384 / 0.304 / 0.341 |

Uhlrich2023은 무릎 부하를 줄이는 보행 수정을 다루는 연구다. **평상 보행은 밴드
안이고 수정 조건에서만 낮다** — 파이프라인 오류가 아니라 중재 효과가 그대로 잡힌
것으로 본다. 검산 밴드는 자릿수를 보는 도구지 이 연구의 가설을 판정하는 기준이 아니다.

### 결과 파일

| 파일 | 뜻 |
| --- | --- |
| `out/inspect_uhlrich2023.txt` | `--inspect` 전체 로그 — AddB 자체 잔차·마커 RMS·GRF 사유 |
| `out/summary_2026-09-08_48trial.csv` | 9/8 실행 원본 (8명 48 trial, subject6·9 포함) |
| `out/summary_7subjects.csv` | 위에 subject10을 더하고 6·9를 뺀 **확정본** |
| `out/real/` | trial별 산출물 (9/8 실행 + 오늘 subject10) |
| `report/` | 현황 PPT·그림·생성 스크립트, 시너지 토크 서명 |

---

## 다음 단계

1. **중복 협응근 묶기**를 `synergy.py`에 옵션으로 — `gasmed`+`gaslat`, `glmed1-3`,
   `glmax1-3`, `addmag*`. 묶은 뒤 `dual`이 `global` 쪽으로 내려오면 SO 임의성이
   원인이라는 증명이자 W 비교의 전제가 된다. **이 실험 하나가 W 비교가 가능한
   이야기인지 아닌지를 가른다.**
2. **walking vs walkingTS를 W(구성)·C(타이밍)로 비교.** `--min-peak`로 빠지는 근육이
   조건마다 다르므로 두 조건의 유지 근육 교집합을 먼저 맞춘다. 같은 조건 trial끼리
   W가 얼마나 일치하는지(재현성)를 먼저 재서 숫자로 둔다 — 평지 보행이라 trial당 한
   걸음, `--concat`해도 조건당 3~4 stride뿐이다.
3. **시너지당 새 근육 하나(다리당 5개)로 걷기 — 그대로는 안 된다.**
   - 이론: 당김 전용 액추에이터로 6 DOF 토크 공간을 양의 결합으로 덮으려면 다리당
     최소 7개다(Phase 0, `actuator_model.positive_spanning_floor`). 5개는
     subtalar·hip_rotation을 잠근 4 DOF에서만 n+1로 턱걸이다.
   - 실측: 5개 시너지의 토크 방향벡터(단위 활성당 평균)가 4 DOF를 양으로 덮는지
     12개 분해(6 subject × 좌우)로 검사했더니 **1/12만 성립**(`report/spanning_check.py`).
     걷기에 필요한 방향만 있고 균형·외란 대응 방향은 못 낸다.
   - 되는 구조: **하드웨어 근육 수 K는 스패닝·스트로크·힘으로(≥7~8, Phase 0~5
     계획), 제어는 시너지 5개로**(`SynergyController`에 W). "8"은 근육 수, "5"는
     제어 차원 — 둘 다 맞다.
   - **했다 (9/14).** 사람 근육을 전부 떼고 새 근육 k개/다리를 달아 MocoInverse로 같은
     걸음을 재현시켰다: **8개 → 보조 최악 8 %(재현), 5개 → 무릎 23~48 %(실패).** 근육
     80개를 시너지 신호로 모는 버전도 같은 순위(8개 22 %, 5개 85 %). NOTES 3절
     "새 근육 8개를 실제로 달아서 걸었다" 참고.
   - **결정: 다리당 8개.** 방향 V 를 고정하고 다른 걸음을 걷게 하는 교차 검증(9/14 오후)
     에서 걸음 하나의 NMF 방향은 다른 걸음에 u 피크 46 이 필요해 실패했고, subject3
     **6걸음을 합쳐 뽑은 8개 방향**은 여섯 걸음 전부 0.2 % 이내, subject4 도 3 % 로 낸다
     (`report/fig_newmuscles_cross.png`, 설계안 `report/newmus_V_subject3_pooled_k8.npz`).
     Moco 보조 가중치는 100 으로 둔다(1 이면 비용 착시). 걷는 영상은
     `report/newmus_video.py` → `report/video_newmus_*.mp4`.
     **EMG 검산(9/14)**: SO 는 발목·햄스트링은 맞고(r 0.6~0.8) recfem·tibant 는 안 켠다(0.1~0.3),
     EMG 시너지는 2~3개 — NOTES 3절. V 를 ID 토크에서 직접 뽑으면 이 SO 오류에서 벗어난다.
     다음: ① V 를 NMF 에 맡기지 말고 필요한 u 최소화로 직접 최적화(7명 전체 + 스쿼트·점프·STS), ② hip_rotation
     포함 5 DOF, ③ 접촉 모델로 스스로 걷기(균형), ④ 그 토크 서명을 내는 실제 부착 경로.
4. **하지 배치 과제 Phase 1~5**는 아래 "하지 인공근육 배치 과제" 절과 NOTES 4-A.
   Phase 0의 열린 질문(어느 봉투를 쓰나)은 오늘 닫혔다 — subject6·9는 데이터 탓이므로
   **`phase0_no69` 봉투를 쓰고**, subject10이 들어왔으니 7명으로 다시 뽑는다.
   `PLAN_lowerlimb_muscle_placement.md`는 참조만 남고 **파일이 없다** — Phase 목록은
   NOTES 4-A와 스크립트 docstring에서 복원한다.

### Claude Code로 이어갈 때

```
~/Desktop/scHooL/의료로봇연구실/addb_opensim/README.md 읽고 이어서 진행해줘
```
