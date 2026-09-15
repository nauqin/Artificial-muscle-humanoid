# 설계 노트

README.md가 "어떻게 돌리는가"라면, 이 문서는 **"왜 이렇게 만들었는가"와
"무엇이 잘못될 수 있는가"**다. 값이 이상할 때 여기부터 본다.

여기 적힌 수치는 전부 번들 Rajagopal 보행 예제(85.07 kg, 1.24 m/s)로 실측한 것이다.
파이프라인을 고친 뒤 값이 달라지면 회귀를 의심한다.

---

## 1. 설계 결정

### GRF-발 매핑은 이름이 아니라 거리로

`addb_io.assign_plates_to_feet()`

힘판 데이터를 어느 발에 연결할지, 보통은 컬럼 이름(`ground_force_r_v` 등)을 보고
정한다. 그런데 데이터 제공자마다 이름 짓는 방식이 다르고, **좌우가 바뀌어도 계산은
에러 없이 그럴듯한 그래프를 내놓는다.**

그래서 이름을 안 본다. 수직력이 체중의 10%를 넘는 프레임에서 모델을 IK 자세로 세운 뒤
(`realizePosition`), COP와 `calcn_r`/`calcn_l`의 실제 수평 거리를 재서 가까운 쪽에
배정한다. 택배 주소 라벨을 믿는 대신 GPS로 확인하는 셈이다.

번들 예제에서 배포된 정답(`dryrun/Rajagopal/ID/grf_walk.xml`)과 정확히 일치했고,
여유가 34 cm라 헷갈릴 여지가 없다.

```
r -> calcn_r   calcn_r=0.137  calcn_l=0.483   (여유 +0.346 m)
l -> calcn_l   calcn_r=0.474  calcn_l=0.137   (여유 +0.338 m)
```

두 힘판이 같은 발을 고르면 총 거리가 작은 조합으로 재배정하고 경고를 찍는다.

### ID 구간은 세 가지를 교집합해서 자동으로 정한다

`run_trial.resolve_window()`

| 자르는 것 | 안 자르면 |
| --- | --- |
| IK ∩ GRF 시간 교집합 | 한쪽에만 데이터가 있는 구간이 들어간다 |
| **지지 구간** — 총 수직 GRF가 체중의 70% 이상인 최장 연속 구간 | 발이 힘판 위에 없는 동안 체중 전체가 골반 잔차로 잡힌다 |
| **양 끝 마진 = 1/컷오프** (6 Hz면 0.167초) | 필터 후 미분한 가속도가 경계에서 망가진다 |

실측으로 확인한 효과다.

| 구간 처리 | 잔차 (peak GRF 대비) |
| --- | --- |
| 전체 구간 그대로 (0.000~2.370 s) | **93.61 %** |
| + 지지 구간 (0.372~2.059 s) | 41.87 % |
| + 양 끝 마진 (0.539~1.892 s) | **7.81 %** |

세 가지를 다 적용하면 **사람이 손으로 좋은 구간을 고른 것과 같은 결과가 자동으로**
나온다. 데이터가 수백 개여도 일일이 구간을 볼 필요가 없다는 뜻이다.

### `.mot` 분류에 문자열 검색을 쓰지 않는다

`addb_io.classify_mot()`

컬럼 이름에 `force`가 들어있는지로 GRF를 판단하면 안 된다. states 파일은
`/forceset/soleus_r/activation` 같은 라벨을 갖기 때문에 GRF로 오분류된다. 실제로
정적 자세 파일(`scale_output_walk.mot`)이 이 방식에 걸려 GRF 후보로 올라왔다.

지금은 **힘 x/y/z 와 작용점 x/y/z 가 모두 있는 파일**만 GRF로 인정한다. IK 파일의
`pelvis_tx/ty/tz`도 x/y/z 삼총사를 이루지만 작용점이 없어서 걸러진다.

### SO 보조 액추에이터

`run_static_optimization.build_reserve_actuators()`

근육만으로 관절 회전력을 못 맞추면 최적화 문제가 아예 안 풀린다. 그래서 모든 자유
좌표에 `CoordinateActuator`를 하나씩 단다.

- `optimal_force = 1` — 활성도가 곧 실제 일반화 힘이 된다. 비용이 활성도 제곱합이라
  근육(활성도 ≤ 1)에 비해 수백~수만 배 비싸져서 최적화가 알아서 피한다.
- 한도 `±1000` — 그래도 정말 필요하면 쓸 수 있게 넉넉히.
- **종속 좌표는 제외** — `knee_angle_r_beta` 같은 구속된 좌표에 액추에이터를 달면
  SO가 깨진다. `Coordinate.isConstrained(state)`로 거른다.

드라이런에서 29개가 붙었고(종속·잠긴 좌표 10개 제외), 실제로 떠맡은 비율은
0.19 % 미만이었다.

### SO 검증에 밴드를 둘 둔 이유

`validate.validate_so()`

좌표마다 이 등식을 확인한다.

```
Σ(근육힘 × 모멘트암) + 보조액추에이터  =?  ID 모멘트
```

모멘트암은 MuscleAnalysis 결과를 쓰고, 시간축이 다르므로 ID 격자에 보간한다.

`recon`(재현 오차) 하나만 보면 안 된다. **보조 액추에이터는 어떤 모멘트든 만들어낼 수
있으므로 재현 오차는 언제나 0에 가깝게 나온다.** 근육이 전혀 일을 안 해도 통과한다.
그래서 `reserve`(보조가 떠맡은 비율)를 함께 본다. 이쪽이 실제 품질 지표다.

좌표 이름은 **정확히 일치**시켜 찾는다. 부분일치를 쓰면 `knee_angle_r`이
`knee_angle_r_beta`의 액추에이터까지 끌어온다.

### 모멘트암 크기는 대표 근육으로만 본다

`run_muscle_analysis.ANTAGONISTS`

부호는 길항근 무리 **전체**로 본다 — 주동근끼리 부호가 갈리거나 길항근과 부호가
같으면 FAIL. 반면 크기(0.02~0.08 m)는 **대표 근육**으로만 본다.

처음에 무리 전체의 중앙값으로 쟀더니 발목이 0.017 m로 나와 WARN이 떴다.
`fdl`(0.013), `perbrev`(0.009), `perlong`(0.012), `tibpost`(0.014) 같은 작은 원위근이
중앙값을 끌어내린 것이고, 정작 아킬레스건(soleus −0.051, gasmed −0.051, gaslat −0.052)
과 `tibant`(+0.050)는 문헌값 그대로였다. 밴드는 주요 근육을 겨냥한 것이므로
대표군(아킬레스건·tibant·vasti·햄스트링·psoas/iliacus·glmax)에만 적용한다.

### 시너지 개수에 개별 근육 기준을 추가한 이유

`synergy.choose_k()`

NMF는 Lee & Seung 곱셈 업데이트(Frobenius norm)를 직접 구현했다. scikit-learn 기본
solver와 같은 알고리즘이고, 의존성을 늘리지 않으려고 직접 짰다. 초기값에 따라 국소해에
빠지므로 20회 재시작해 최선을 고른다.

k를 전체 VAF ≥ 90%로만 고르면 **k=5에서 통과한다.** 그런데 그때 최저 근육 VAF가
50% 수준이라 조용한 근육은 사실상 설명되지 않은 상태다. 활성도가 큰 근육 몇 개가
전체 VAF를 끌어올려 가려버린 것이다.

| k | 전체 VAF | 최저 근육 VAF | 미달 근육 |
| --- | --- | --- | --- |
| 4 | 87.4 % | 48.1 % | 9개 |
| 5 | 90.4 % | 50.4 % | 4개 ← global |
| 6 | 93.2 % | 52.4 % | 3개 |
| 7 | 94.8 % | 66.8 % | 2개 |
| **8** | **96.1 %** | **85.4 %** | **0개** ← dual |

EMG 기반 보행 시너지 연구는 보통 4~5개를 보고하지만, 그건 채널이 8~16개이고 실측
EMG가 더 매끄럽기 때문이다. 여기서는 근육 37~38개의 SO 활성도를 쓰므로 k가 더 필요하다.

---

## 2. 드라이런 기준값

번들 Rajagopal 보행, 85.07 kg, 1.24 m/s, 구간 0.539~1.892 s(자동), 6 Hz.

```bash
python batch.py --data dryrun --out out/batchtest --force --so --synergy
```

```
ankle_r 1.643 / ankle_l 1.681 Nm/kg   WARN   둘 다 저굴, push-off 시점
knee_r  0.670 / knee_l  0.701 Nm/kg   PASS
hip_r   0.726 / hip_l   0.702 Nm/kg   PASS
grf_peak 1.215 BW                     PASS
residual 7.81 %                       WARN
cop_dist 0.137 m                      PASS

모멘트암  아킬레스건(soleus/gasmed/gaslat) -0.051 m, tibant +0.050 m      PASS
          vasti -0.045 m, hamstrings +0.042~0.064 m                       PASS
          psoas +0.031 / iliacus +0.043, glmax -0.056~-0.088 m            PASS

SO        보조 액추에이터 29개
          재현 오차 최대 0.072 %, 보조 액추에이터 최대 0.188 %             PASS
          근육 피크힘 gasmed 2742 N, soleus 2378 N, tibant 1338 N

시너지    global 좌우 k=5 / dual 오른쪽 k=8 (VAF 96.1 %), 왼쪽 k=7 (96.6 %)  PASS
```

### 이 값은 재현된 것이다

작업물이 소실돼 macOS에서 다시 만들었고, 예전 윈도우 실측값과 대조했다.
Scale·IK를 새로 돌렸는데도 관절 모멘트는 소수 셋째 자리까지 같았다.

| | 예전 (Windows) | 지금 (macOS) |
| --- | --- | --- |
| ankle r/l | 1.644 / 1.680 | 1.643 / 1.681 |
| knee r/l | 0.670 / 0.700 | 0.670 / 0.701 |
| hip r/l | 0.726 / 0.703 | 0.726 / 0.702 |
| grf_peak | 1.215 | 1.215 |
| residual | 7.82 % | 7.81 % |
| 구간 | 0.539~1.892 | 0.539~1.892 |
| 보조 액추에이터 | 29개 | 29개 |
| recon / reserve 최대 | 0.07 / 0.19 % | 0.072 / 0.188 % |
| 전체 구간 잔차 | 93.6 % | 93.61 % |
| 15 Hz 잔차 | 22.9 % | 22.49 % |

**다른 것은 `dual` 기준의 k뿐이다** (예전 좌우 7 / 지금 8·7). 전체 VAF는 k=7에서
94.8 %로 예전(94.7 %)과 같고, 갈린 것은 최저 근육 VAF(66.8 % vs 79.7 %)다. 이는 아래
3절에서 다루는 `dual`의 불안정성 그 자체이며, `global`은 양쪽 다 k=5로 흔들리지 않는다.

### WARN 두 개는 코드 문제가 아니다

**발목 1.64/1.68** — 밴드 상단 1.6을 5% 넘는다. 그런데 좌우가 거의 같고, 힘이 최대가
되는 시점도 push-off로 정확하다. 이 피험자가 85 kg으로 무거운 편이라 개인차로 본다.
자릿수가 맞으므로 검산의 목적은 달성됐다.

**잔차 7.8%** — 공식 예제는 RRA로 보정한 뒤 ID를 돌리지만, 우리는 AddBiomechanics
워크플로(스케일 모델 + IK → 바로 ID)를 재현하려고 그 단계를 일부러 건너뛴다.
AddBiomechanics 데이터는 자체 동역학 최적화를 거쳐 배포되므로 잔차가 **더 낮게**
나와야 정상이고, 그래서 이 값이 "받은 데이터가 쓸만한가"를 판별하는 잣대가 된다.

### 시너지 구성

좌우가 독립 분해인데 같은 묶음이 나왔다.

| 시너지 | 묶인 근육 | 피크 시점 |
| --- | --- | --- |
| 심부 저굴근 | soleus, tibpost, fhl, perlong | push-off |
| 고관절 신전 | glmax, semimem, semiten, bflh, addmagIsch | 입각 초기 |
| 배굴근 | tibant, edl, ehl | 유각기 |
| 내전근 | addlong, addbrev, addmagProx, grac | 유각 후반 |
| 체중 수용 | glmax, vasti, glmed, glmin | 입각 초기 |

아킬레스건 모멘트암 5.1 cm는 문헌값(~5 cm)과 맞는다.

---

## 3. 밟은 함정

전부 **에러 없이 조용히 틀린 답을 주는** 종류였다.

### 저역통과 6 Hz가 옳다

| 컷오프 | 잔차 |
| --- | --- |
| 6 Hz | 7.81 % |
| 15 Hz | 22.49 % |
| 무필터 | 74.70 % |

무필터는 좌표 미분 노이즈로 완전히 무너진다.

### 모델 자동 선택이 템플릿을 골랐다

폴더에서 `.osim`을 자동으로 고를 때 그 사람 몸에 맞춘 모형 대신 배포판 기본 모형
(`Rajagopal2016.osim`)을 집었다. 체중이 85 kg이어야 하는데 75 kg으로 나왔고 잔차가
**84%**까지 뛰었다. 검산기가 잡아냈다.

지금은 템플릿 이름(`unscaled`, `generic`, `Rajagopal2016` 등)을 뒤로 밀고, 선택한
모델의 체중을 로그에 찍는다. 이상하면 `--model`로 지정한다.

재구축하면서 확인했다 — `generic_unscaled.osim`은 지금도 75.337 kg으로 로드된다.

### IK 후보가 하나면 이름을 안 본 탓에 엉뚱한 GRF와 묶였다

`grf_run.mot`이 이름 토큰이 하나도 안 겹치는데도 walk IK와 짝지어졌다. 원인은
"IK 후보가 하나뿐이면 검사를 건너뛴다"는 지름길이었다. 지금은 **IK도 GRF도 하나뿐일
때만** 이름을 안 보고 묶고, GRF가 여럿이면 반드시 토큰을 확인한다. 확정 못 하면
경고를 찍고 건너뛴다 — 엉뚱한 조합으로 계산해 조용히 틀리는 것보다 **빠뜨리는 쪽이
낫다.** 로그의 `[WARN] 짝을 확정할 수 없어 건너뜀`을 볼 것.

### 활성도 파일에 근육이 아닌 것이 섞여 있었다

Rajagopal 모델의 팔은 근육이 아니라 토크 액추에이터다. 그래서 SO 활성도 파일에
`shoulder_flex_r`, `elbow_flex_r`, `wrist_flex_r`, `pro_sup_r` 이 근육과 나란히
들어있다. 이름만 보고 거르면 놓친다 — `_r`로 끝나고 `reserve_`로 시작하지도 않는다.

그대로 시너지에 넣었더니 **`shoulder_flex_l` 이 가중치 1.00인 성분**이 만들어졌고,
`elbow_flex`가 거의 모든 시너지에 끼어들었다. 근육 수도 41개로 부풀었다.

지금은 `--model`로 모델을 읽어 `model.getMuscles()`에 있는 이름만 쓴다. 38개(오른쪽),
37개(왼쪽)로 줄고 성분도 생리적으로 읽히게 바뀌었다. **`synergy.py`에 `--model`은
선택 사항이 아니다.**

### 시너지 `--kmax`를 낮게 잡으면 조용히 틀린 답이 나온다

처음에 6으로 두었더니 기준을 만족하는 k가 하나도 없는데도 **아무 말 없이 6을 골라**
그럴듯한 표를 그렸다. 지금은 못 찾으면 PASS를 주지 않고 WARN과 경고를 찍는다.

### 개별 근육 VAF 기준이 SO의 임의 배분을 재현하라고 요구하고 있었다

실제 데이터 48 trial에서 시너지 개수가 5~12로 흔들렸고 4개는 `--kmax` 상한 12에
닿았다. 문헌은 보행에서 4~6개를 보고한다.

먼저 두 기준을 갈라 보았다. 96개 분해(48 trial × 좌우) 전부에서

| 기준 | 범위 | 평균 |
| --- | --- | --- |
| `global` (전체 VAF ≥ 90 %) | **4 ~ 6** | 5.2 |
| `dual` (+ 모든 근육 ≥ 75 %) | 5 ~ 12 | 7.9 |

흔들림은 전적으로 개별 근육 조건이 만들고 있었다. 그런데 `global` 시점의 최저 근육
VAF가 평균 55.9 %, 최소 24 %였다 — 드라이런에서 `dual`을 도입하게 만든 그 실패
모드가 실제 데이터에서도 재현된다. 그러니 "`dual`이 과하다"로 끝낼 수 없었다.

**표본 부족인지 확인하는 실험을 했다.** trial 하나는 보행 1.3초 조각이라 NMF 표본이
90~140개뿐이다. 조건별로 이어붙여(`--concat`) 표본을 3~4배로 늘렸다.

|  | 표본 | `global` | `dual` |
| --- | --- | --- | --- |
| trial 개별 | ~90 | 4~6 | 6~10 |
| walking 이어붙임 | 279 | 5 / 5 | 8 / 7 |
| walkingTS 이어붙임 | 396 | 6 / 6 | 12 / 12(미달) |

**`dual`은 내려가지 않았다.** `global`만 5~6으로 요지부동이었다. 표본 부족이 아니다.

결정적 증거는 시너지 구성에 있었다. 이어붙인 분해에서 **근육 하나짜리 성분**이 나온다.

```
시너지 3   기여 15%,  근육 1개      gaslat_r 1.00
시너지 5   기여 10%,  근육 1개      recfem_l 1.00
```

시너지는 함께 켜지는 근육 묶음인데, 혼자인 성분은 그 근육이 **다른 어떤 근육과도
함께 움직이지 않는다**는 뜻이다. `gaslat`은 `gasmed`와 모멘트암이 거의 같아 생리적
으로 같이 켜져야 한다. 따로 노는 것은 **SO가 둘 사이에 부하를 임의로 나눈 흔적**이고,
그 임의성은 매 사이클 다르게 나오므로 표본을 늘려도 평균되지 않는다.

`dual`을 끌어올린 근육 빈도도 같은 이야기를 한다 — `gaslat`(23회), `addmagMid`(17회),
`glmed3`, `glmax2`, `addmagProx`, `addlong`. 전부 모멘트암이 겹치는 중복 협응근이다.

**결론: `global`을 본 결과로 쓰고 `dual`은 진단 지표로 남긴다.** 개별 근육 VAF 기준은
SO 활성도가 개별 근육 수준에서 신뢰할 만할 때만 의미가 있는데, 중복 근육에서는 그렇지
않다. `dual`이 못 채워지는 것 자체가 "이 근육은 개별 해석 금지"라는 신호로 쓰인다.
`synergy.py`는 성분에 표시할 근육이 하나뿐이면 `← 근육 1개짜리 성분`을 찍는다.

**조건 비교 시 주의.** `--min-peak`(0.05)로 빠지는 근육이 walking은 2개,
walkingTS는 4~5개로 달랐다. 근육 집합이 다르면 두 조건의 시너지를 직접 비교할 수
없다. 조건 간 비교를 할 때는 집합을 먼저 맞춰야 한다.

### midstance 골이 지지 구간을 반으로 쪼갰다

증상은 **한쪽 발목 모멘트만 0에 가깝게 나오는 것**이었다. 실제 데이터 48 trial 중
4개가 FAIL이었고, 넷 다 왼쪽이었다.

| trial | 발목 r / l | 무릎 r / l | 분석 구간 |
| --- | --- | --- | --- |
| subject9 walkingTS2 | 1.396 / **0.145** | 0.579 / 0.765 | 0.46 s |
| subject9 walkingTS3 | 1.313 / **0.250** | 0.481 / 0.951 | 0.44 s |
| subject6 walkingTS2 | 1.255 / **0.462** | 0.446 / **1.488** | 0.48 s |
| subject6 walkingTS3 | 1.155 / **0.529** | 0.493 / **1.335** | 0.52 s |

검산기는 `의심: GRF-발 매핑 / 좌표계 / 단위`를 안내했지만 **전부 헛다리였다.**
매핑은 정상이었다 (`calcn_l -> calcn_l`, 여유 +0.331 m).

진짜 원인은 FAIL 4건이 **구간이 가장 짧은 4개와 정확히 일치한다**는 데 있었다.
FAIL 아닌 44개는 전부 0.629 s 이상이고 FAIL 넷은 0.52 s 이하 — 그 사이가 비어 있다.

보행의 총 수직 GRF는 쌍봉이고 midstance에 골이 있다. 정상 trial은 이 골이
0.80~0.82 BW에서 멈추는데, **문제의 trial은 0.67~0.70 BW까지 내려가 0.7 BW 문턱을
스쳤다.** 그 순간 `grf_support_window()`가 연속 구간을 둘로 쪼갰고, 최장 구간을
고르니 반 걸음만 남았다. 왼발 push-off가 구간 밖으로 밀려나 발목 모멘트가 0이 된
것이다. 지면반력 데이터 자체는 0~1.3초 내내 멀쩡했다.

지금은 끊긴 두 구간 사이가 **짧고(≤0.25 s) 얕으면(≥0.5 BW)** 이어 붙인다. 생리적
골은 0.1초 남짓이고 0.5 BW 밑으로 내려가지 않는다.

**히스테리시스(들어갈 때 0.7, 나올 때 0.3)를 먼저 시도했는데 그게 틀렸다.** 골
문제는 고쳤지만 양 끝까지 같이 늘려서, 드라이런에서 구간 끝이 1.892 → 1.933으로
밀렸다. 늘어난 부분을 열어보니 midstance 골이 아니라 **trial 끝에서 피험자가 힘판을
빠져나가는 구간**(오른발 0.79 → 0.30 BW, 왼발 0)이었다. 원래 문턱이 일부러 잘라내던
곳이다. 부작용 없는 처방은 양 끝을 건드리지 않고 중간 골만 잇는 쪽이었다.

또 하나 — 검산기가 이 경우를 스스로 구분한다. 관절 모멘트가 0에 가깝고 **구간이 한
걸음(0.8 s)보다 짧으면** 매핑이 아니라 midstance 골을 의심하라고 안내한다.

### AddBiomechanics GRF의 `_mx` 컬럼이 힘으로 분류됐다

nimblephysics가 내보내는 GRF 컬럼은 `ground_force_calcn_r_vx`(힘) / `_px`(COP) /
`_mx`(자유모멘트) 꼴이다. `torque`·`moment`라는 **단어**만 찾으면 `_m`으로 끝나는
스템이 force로 분류되고, `plates[key]["force"]`에 나중에 덮어쓴다. 결과적으로
**자유모멘트를 지면반력으로 알고 ID를 돌린다.** 컬럼 순서가 v → p → m이라 항상
마지막에 덮인다. 에러는 나지 않는다.

지금은 컬럼 이름을 토큰으로 쪼개 뒤에서부터 `v`/`p`/`m`을 찾고, plate key를 만들 때
그 토큰을 지워 세 항목이 같은 plate로 묶이게 한다. 기존 규약(`ground_force_r_vx`,
`ground_torque_r_x`, `1_ground_force_vx`)은 그대로 통과한다.

**변환기 출력 형식을 실제로 검증했다.** 드라이런 데이터를 변환기와 같은 형식
(라디안 `inDegrees=no`, `ground_force_<body>_v/p/m`)으로 다시 써서 파이프라인에
넣었더니 `summary.csv`의 모든 항목이 원본 형식과 **완전히 같았다** — 관절 모멘트,
잔차 7.81 %, cop 0.137 m, SO recon/reserve, 시너지 개수까지. COP 거리 판정도 컬럼
이름과 무관하게 `calcn_r -> calcn_r`을 여유 0.346 m로 맞혔다.

### `manualReview`로 두 사람을 같이 버렸는데 한 명은 최상급이었다

subject2·subject10은 모든 trial의 모든 프레임이 `manualReview`라 변환기가 0 trial로
건너뛰었다. 이름표만 보고 "품질 보증 안 됨"으로 같이 제외했는데, `--inspect`로
AddBiomechanics 자체 잔차를 찍어보니 둘은 정반대였다.

| | AddB 잔차 (%BW) | |
| --- | --- | --- |
| subject10 | 0.46 ~ 0.6 | subject3(0.7~1.2)보다 좋다 |
| subject2 | 5.3 ~ 13.0 | 실제로 나쁘다 |

`manualReview`는 `torqueDiscrepancy`처럼 프레임별 물리 판정이 아니라 trial 전체에
균일하게 찍히는 표시다. 지금은 `--allow-manual-review`로 그 플래그만 눈감고
(`SOFT_REASONS`), 물리 판정은 우리 검산기가 한다. subject10은 6 trial 전부 잔차
3.8~5.0 %, reserve < 0.41 %로 통과했다.

### subject6·9의 잔차는 데이터 탓이었다 — `--inspect`가 갈랐다

예전 노트의 열어둔 숙제. 우리 잔차가 subject6 8.1 %, subject9 7.5 %로 나머지
(2.7~4.3 %)의 두 배였는데, 파이프라인 탓인지 데이터 탓인지 몰랐다. AddBiomechanics가
스스로 보고한 골반 잔차(`getTrialLinearResidualNorms(trial, pass)`)로 갈랐다.

| subject | AddB 평균 잔차 (%BW) | 순간 최대 |
| --- | --- | --- |
| 3 · 4 · 5 · 10 · 11 | 0.3 ~ 1.2 | 8 ~ 50 N |
| 6 | 13.4 ~ 18.2 | 2100 ~ 2600 N |
| 9 | 16.2 ~ 29.7 | 3000 ~ 3300 N |

20~40배 차이고, 순간 최대가 체중의 3~5배다. AddBiomechanics 자신도 이 둘의 동역학을
못 맞췄다는 뜻이다. 우리 파이프라인이 이 둘만 골라낸 것은 오히려 정상 작동의 증거고,
둘을 빼면 잔차 분포가 균질해진다(평균 3.49 %, 최대 5.82 %). 7·8은 trial 단위로
한두 개만 나쁘다(subject8 walkingTS1 28 %BW) — subject 단위 현상이 아니라 그대로 쓴다.

**`getTrial*` API는 pass 인자를 받는다.** `getTrialLinearResidualNorms(trial)`처럼
빠뜨리면 예외가 나는데, 처음 `--inspect`는 그걸 `try/except`로 삼켜서 잔차 열이
조용히 `-`로 나왔다. pass 0(kinematics)은 잔차 900 N, pass 2(dynamics)는 4 N이니
어느 pass를 보느냐가 전부다.

### 소실됐다던 `out/real/`이 남아 있었다 — 배치가 덮어쓰기 직전이었다

9/8에 돌린 48 trial 결과가 `out/real/`에 그대로 있었다. 오늘 subject3 하나만 있는
`data/`로 같은 `--out`에 배치를 띄웠는데, `--force`가 없어도 **`summary.csv`는 끝에
무조건 새로 쓴다.** 끝났으면 48 trial 요약이 1 trial짜리로 바뀌었을 것이다. 중단하고
`out/summary_2026-09-08_48trial.csv`로 백업했다. 같은 `--out`을 다시 쓸 때는 먼저
`summary.csv`를 복사해 둘 것.

### OpenSim 도구 자체의 함정

- 번들 `dryrun/Rajagopal/ID/id_setup_walk.xml`은 그대로 못 쓴다. 존재하지 않는 RRA
  출력과 CMC 조정 모델을 입력으로 하며, 우리 워크플로와도 다르다. `run_trial.py`가
  새로 구성한다.
- **IK 툴은 `results_directory`를 미리 만들어주지 않아** 첫 실행이 실패한다.
  스크립트에서 `makedirs`로 처리한다.
- MuscleAnalysis 중 뜨는 `Unable to achieve required assembly error tolerance` 경고는
  무시해도 된다 (달성 1.5e-8 vs 요구 1e-10, 구속을 완화해 재시도하고 성공한다).
- MuscleAnalysis는 trial 하나에 `.sto`를 90여 개 쏟아낸다. 쓰는 것은 모멘트암뿐이라
  기본으로 정리한다 (`--keep-all`로 끌 수 있다).

### `initSystem()` 전에 모델을 건드리면 트레이스백 없이 죽는다

`addb_io.Pose`

`model.getSimbodyEngine().convertDegreesToRadians(storage)` 를 `initSystem()` 전에
부르면 **파이썬이 세그폴트(exit 139)로 사라진다.** 예외도, 스택도 없다.

파이프라인 본체에서는 늘 `model.initSystem()` 을 먼저 부르고 있어서 안 보이다가,
`animate.py` 에서 모델을 새로 열고 바로 Pose 를 만들었더니 터졌다. 원인을 찾느라
`Storage.print` 와 `printResult` 를 의심했는데 둘 다 무죄였다 — 그 앞 줄이 이미
죽어 있어서 print 가 안 찍혔을 뿐이다. **출력이 하나도 없으면 의심 구간을 앞으로
당겨 잡을 것.**

지금은 `Pose.__init__` 이 맨 처음에 `initSystem()` 을 부른다.

### 뷰어가 근육 동역학 때문에 죽는다

`VisualizerUtilities.showMotion()` 은 운동학만 재생하는데도 Dynamics 단계까지
realize 한다. 그러면 Millard 근육의 섬유 속도 Newton 반복이 수렴하지 않는다.

```
Exception caught in Millard2012EquilibriumMuscle::calcFiberVelocityInfo from addlong_r
addlong_r Fiber velocity Newton method did not converge
```

재생에 근육 동역학은 필요 없으므로 `animate.py` 가 모든 근육의
`ignore_tendon_compliance` / `ignore_activation_dynamics` 를 켜서 강체건으로 바꾼다.
**이건 재생 전용 처리다** — SO나 ID 쪽 모델에는 손대지 않는다.

### `Storage.exportToTable()` 은 파이썬에서 못 쓴다

SWIG가 감싸주지 않아 `SwigPyObject` 가 돌아오고 `getNumRows()` 조차 없다.
`TimeSeriesTable` 이 필요하면 값을 numpy로 뽑아 `.mot` 으로 쓴 뒤 다시 읽는다
(`addb_io.write_mot`). `animate.py` 가 그렇게 한다.

### 화면을 믿지 말 것

- 이 피험자는 발 Y 스케일이 **1.66배**라, 발 메시가 calcn 원점 아래로 내려간다.
  GUI에서 발이 살짝 잠겨 보이는 건 이 때문이고 IK 결과는 정상이다.
  `check_foot_height.py`로 확인할 수 있다 (발 마커 최저 +7.6 mm, calcn_l 최저
  −4.8 mm, 골반 0.953~1.001 m, 힘판 표면 y=0).
- **OpenSim GUI의 체크무늬 지면은 신뢰하지 말 것.** 같은 모델·모션인데 GUI에서는 몸
  전체가 지면 아래로 보였고, Simbody 뷰어(`animate.py`)에서는 정상으로 보였다.
  바닥 기준이 필요하면 Simbody 뷰어를 쓴다 — 지면이 정확히 y=0에 그려진다.

---

## 4. 환경 (macOS)

| 항목 | 값 |
| --- | --- |
| 기계 | Apple Silicon (arm64), macOS |
| 분석 환경 `osim` | Python 3.11.16 + OpenSim 4.6 + numpy 2.4.6 + pandas + matplotlib |
| 변환 환경 `b3d` | Python 3.11.16 + nimblephysics 0.10.52.1 + numpy 1.26 (+ torch) |
| GUI / CLI | `/Applications/OpenSim 4.6` (`opensim-cmd`도 arm64 네이티브) |

```bash
conda create -n osim -c opensim-org -c conda-forge python=3.11 opensim=4.6 pandas matplotlib -y
conda create -n b3d python=3.11 "numpy<2" -y
conda activate b3d && python -m pip install "nimblephysics==0.10.52.1"
```

### 윈도우 시절 이야기는 대부분 폐기됐다

예전 노트에는 `shell.bat` / `run.bat` / `site-packages\opensim_dlls.pth` /
"환경을 활성화하지 않으면 matplotlib이 트레이스백도 없이 죽는다" 가 있었다.
전부 Windows에서 conda의 `Library\bin`이 PATH에 없어 DLL을 못 찾던 문제라
**맥에서는 해당 없다.** `conda activate osim` 또는 `./run.sh` 하나면 된다.

`.b3d`를 리눅스에서 변환해 결과 폴더만 옮기던 이원화도 사라졌다 — nimblephysics는
**Windows 휠만 없고 macOS 휠은 있다.** 애플실리콘은 cp39·cp311 `universal2` 휠을 쓴다.

### 왜 환경을 둘로 나눴나

OpenSim 4.6 conda 패키지는 `numpy 2.4.*`를 못 박고, nimblephysics 0.10.52.1은
numpy 1.x 시절 빌드다. 한 환경에 밀어 넣으면 언제 깨질지 모른다. 같은 기계 안이라
분리해도 비용이 없다 — `.b3d`를 `data/`로 변환해 두고 `osim` 환경으로 넘어가면 된다.

### 파이썬은 3.11이어야 한다

OpenSim 4.6의 osx-arm64 conda 빌드는 py3.11/3.12/3.13이 있지만, nimblephysics는
**3.11이 애플실리콘에서 쓸 수 있는 가장 높은 버전**이다(0.10.52.1의 universal2 휠은
cp39·cp311뿐). 두 환경을 같은 3.11로 맞춰 두면 헷갈릴 일이 없다.

nimblephysics는 PyPI에 **sdist가 하나도 없다.** 미리 빌드된 휠뿐이라 인터프리터에
맞는 휠이 없으면 소스 컴파일 폴백도 없이 그냥 실패한다. 3.13 이상은 어느 플랫폼에도
휠이 없다.

### `pip`이 환경 밖을 가리킨다

`conda activate b3d` 후에도 `pip install nimblephysics`가

```
ERROR: Could not find a version that satisfies the requirement (from versions: none)
```

로 실패했다. 원인은 패키지가 아니라 `pip` 자체였다 — `which pip`이
`/Library/Frameworks/Python.framework/Versions/3.14/bin/pip`을 가리키고 있었고,
`pip debug --verbose`가 **cp314 태그**를 찍고 있었다. 3.14용 휠은 존재하지 않으니
"없다"가 맞는 말이었다.

**`python -m pip`를 쓰면 된다.** 프레임워크 파이썬이 설치된 맥에서는 conda 환경을
활성화해도 `pip` 실행파일이 PATH 앞쪽의 다른 인터프리터 것일 수 있다.

### OpenSim 4.6에 `SynergyController`가 있다

`osim.SynergyController` / `osim.SynergyVector`. 액추에이터 집합에 시너지 벡터를
붙이면 k개 입력 신호(`synergy_excitation_<i>`)로 전체를 구동하고, Moco에 넘기면
자동 감지된다. 시너지로 걷기 시뮬레이션의 상한 검증(80근육·5신호)에 그대로 쓴다.
libcasadi가 설치돼 있어 Moco도 돈다.

### 드라이런 자산은 앱 번들 안에 있다

macOS OpenSim 설치본에는 `Models/` 폴더가 밖으로 나와 있지 않다. 예제는

```
/Applications/OpenSim 4.6/OpenSim 4.6.app/Contents/Resources/opensim/Resources.zip
```

안의 `Models/Rajagopal/` 이다. 이걸 `dryrun/Rajagopal/`로 풀어 쓴다.

**단, 스케일된 모델과 IK 결과는 들어있지 않다.** 설정 XML과 원자료(마커·GRF)만
있으므로 한 번은 직접 돌려야 한다.

```bash
cd dryrun/Rajagopal/Scale
opensim-cmd run-tool scale_setup_walk_scaleOnly.xml   # → subject_scaled_walk.osim
opensim-cmd run-tool scale_setup_walk.xml             # 마커 재배치
mkdir -p ../IK/results_walk
cd ../IK && opensim-cmd run-tool ik_setup_walk.xml    # → results_walk/ik_output_walk.mot
```

`scale_setup_walk_scaleOnly.xml`에 `<mass>85.0668</mass>`가 박혀 있다. 결과 모델이
85.07 kg이 아니면 뭔가 잘못된 것이다.

---

## 4-A. 하지 인공근육 배치 — Phase 0 에서 정한 것

### 좌표 부호 규약은 문서가 아니라 데이터에서 확인했다

일반화 모멘트가 양수면 좌표를 키우는 방향이다. 실제 ID 결과에서 좌표별로 +·− 중
어느 쪽이 큰지 세어 규약을 확정했다 (보행이므로 저굴·신전·외전이 커야 정상).

| 좌표 | 최대 + (Nm) | 최대 − (Nm) | 읽는 법 |
| --- | --- | --- | --- |
| ankle_angle | 11 | **189** | − 가 저굴 — push-off |
| knee_angle | 58 | **92** | − 가 신전 |
| hip_flexion | **103** | 74 | + 가 굴곡 |
| hip_adduction | 39 | **105** | − 가 외전 — 입각기 중둔근 |
| subtalar_angle | **29** | 6 | 방향 이름은 아직 추정 |
| mtp_angle | 0.2 | 0.3 | **다른 좌표의 1/300** |

`mtp_angle` 을 기본 잠금으로 두는 근거가 이 표다. 0.3 Nm 에 구동원을 배정할 이유가 없다.

### 수동 힘 모델은 그대로 두면 발산한다

계획 2.2 의 `f_passive = F_max · (l − L_M)/(ε_p · l₀)` 에 `ε_p = 0.10`,
스트로크 `ε_max = 0.40` 을 넣으면 **완전 신장 자세에서 `f_passive = 4·F_max`** 가 된다.
길항근이 자기 능동 최대의 4배로 버티는 셈이라 어떤 배치도 불가 판정이 난다.
물리적으로도 SMA 코일실은 그 전에 항복한다.

`actuator_model.PASSIVE_CAP_RATIO`(기본 1.0)로 상한을 두었다. 이 값은 **가정**이므로
결과에 항상 같이 찍는다. 기계연의 힘-변위 곡선이 오면 실측으로 바꾼다.

### 모듈 수의 이론적 바닥은 다리당 7개다

당김 전용 액추에이터의 모멘트암 열들이 n차원 토크 공간을 **양의 결합**으로 덮으려면
최소 n+1 개가 필요하다. 하지 6 DOF → 다리당 7개, 양다리 14개.
`actuator_model.positive_spanning_floor()`.

Phase 3 이 클러스터 수 K=8 을 후보로 두는데, 이건 이미 이론적 바닥 바로 위다.
K=8 이 실제로 통과하면 "거의 최소"라는 뜻이고, 통과 못 하면 바닥이 8~10 사이에 있다.

### 요구 토크 봉투는 잔차가 높은 두 사람이 지배한다

48 trial 봉투에서 **12개 방향 중 6개의 최대값이 subject6·subject9 에서 나왔다.**
이 둘은 잔차가 7.5~8.2 % 로 유일하게 높은 피험자다. 동역학이 안 맞는 데이터가
설계 요구값을 정하는 셈이라, 빼고도 계산해 민감도를 본다
(`screen_modules.py --exclude-subject`).

| 방향 | 8명 전체 | subject6·9 제외 | 차이 |
| --- | --- | --- | --- |
| 고관절 외전 | 1.587 | 1.091 | **−31 %** |
| 고관절 내회전 | 0.236 | 0.173 | −27 % |
| 무릎 신전 | 1.488 | 1.119 | **−25 %** |
| 고관절 굴곡 | 1.233 | 1.081 | −12 % |
| 고관절 신전 | 0.921 | 0.853 | −7 % |
| **발목 저굴** | **1.961** | **1.961** | **0 %** |
| 발목 배굴 · 무릎 굴곡 · 내전 · 외회전 · subtalar | 동일 | 동일 | 0 % |

무릎 신전 최대값 1.488 Nm/kg 은 subject6 walkingTS2 에서 나오는데, 이 trial 은
midstance 골 함정 표에 올라 있던 바로 그 trial 이다. 발목은 구간 수정으로
0.462 → 1.356 으로 회복됐지만 무릎 1.488 은 그대로다. 밴드(0.4~0.8)의 약 2배라
**이 값이 진짜인지는 아직 모른다.**

**어느 쪽을 쓸지는 아직 정하지 않는다.** subject6·9 의 잔차가 데이터 탓인지
우리 계산 탓인지 가르는 것이 먼저다 (`b3d_to_opensim.py --inspect`).
지금은 두 봉투를 나란히 보관한다 — `out/lowerlimb/phase0/`, `phase0_no69/`.

### Phase 0 의 결론 — 인체 부착점 그대로는 안 된다

ROM 전체를 덮으려면 이동량 `Δl = r·ΔΘ` 가 필요한데, 인체 모멘트암은 이미 알고 있다
(NOTES 2절 실측). 로봇 70 kg 기준:

| 방향 | 대표근육 | r (cm) | 필요 Δl (mm) | 필요 힘 (N) |
| --- | --- | --- | --- | --- |
| 발목 저굴 | soleus/gasmed/gaslat | 5.1 | 62 | 2692 |
| 발목 배굴 | tibant | 5.0 | 61 | 163 |
| 무릎 신전 | vasti | 4.5 | 106 | 1482 |
| 무릎 굴곡 | hamstrings | 5.3 | 125 | 791 |
| 고관절 굴곡 | psoas/iliacus | 3.7 | 90 | 2023 |
| 고관절 신전 | glmax | 7.2 | **176** | 749 |

모듈 **하한 규격**(최대 스트로크 40 mm)으로는 **여섯 방향 전부 불가**다.
**상한 규격**(최대 120 mm)이면 네 방향이 M1 하나로 되고, 고관절 신전(176 mm)과
무릎 굴곡(125 mm)만 남는다.

즉 이 과제의 답은 "모듈을 몇 개 붙이나"가 아니라 **부착점을 옮겨 모멘트암을 줄이는
것**이다 (`Δl = r·ΔΘ` 이므로 r 을 줄이면 변위가 준다). 대신 같은 토크를 내려면
힘이 `τ/r` 로 커진다 — 이 맞바꿈이 Phase 5 최적화의 본질이다.

그리고 하한·상한의 차이가 결론을 뒤집으므로, **기계연에 힘-변위 곡선을 요청할
정량적 근거**가 생겼다.

**이 결론은 요구 토크와 무관하다.** 막는 것은 힘이 아니라 변위이고,
`Δl = r · ΔΘ` 에는 τ 가 들어가지 않는다. 위의 subject6·9 민감도(최대 31 %)가
바뀌어도 Δl 열은 한 자리도 안 움직인다. 즉 **Phase 0 의 결론은 요구 토크 봉투가
확정되기 전에도 이미 확정이다.** τ 는 필요한 힘(F = τ/r)과 모듈 개수에만 영향을 준다.

### `walkingTS` 는 trunk sway 다 — 이름이 아니라 데이터로 확인했다

OpenCap 논문(p16)이 "trunk leaned laterally over stance leg" 라고 정의한다. 우리 IK 로
확인하면 요추 측굴(`lumbar_bending`) 진폭만 24° → 47° 로 두 배가 되고 골반 기울기·고관절
내전·무릎 굴곡은 0.93~1.10× 로 사실상 그대로다(7명 21+21 trial). 즉 **다리는 같은 걸음을
걷고 상체만 흔든다.**

이게 시너지 해석을 바꾼다. 다리 근육 시너지 개수가 조건 간에 안 변한 것은 "중재가
제어 차원을 못 바꿨다"가 아니라 **다리 운동학을 거의 안 건드린 수정이니 예상되는 결과**다.
조건 차이가 있다면 고관절 외전근(S1: glmed·glmin·tfl) 쪽이다 — 몸통이 디딤발 위로 오면
외전근이 골반을 받칠 필요가 준다. 그리고 이 수정이 겨냥한 무릎 내전 모멘트는 Rajagopal
무릎이 1 DOF 라 우리 파이프라인에 없다. 보려면 Joint Reaction Analysis 를 얹어야 한다
(논문이 Mocap 쪽에 쓴 방법과 같다: SO → JRA).

원본 데이터셋(SimTK)에는 하지 근육 16개의 EMG 가 있다. `.b3d` 에는 안 들어왔다
(`frame.emgSignals` 가 비어 있다). SO 활성도와 시너지를 실측 EMG 로 대조할 수 있는
자료이므로, 시너지 W 를 믿어야 하는 단계가 오면 받아서 붙일 가치가 있다.

### Phase 0 의 열린 질문은 닫혔다 — `phase0_no69` 봉투를 쓴다

"subject6·9 의 잔차가 데이터 탓인지 우리 계산 탓인지" — `--inspect` 로 AddBiomechanics
자체 잔차를 찍어 **데이터 탓으로 확정**했다(3절). 따라서 요구 토크 봉투는
`out/lowerlimb/phase0_no69/` 쪽이다. 무릎 신전 1.488 Nm/kg(subject6 walkingTS2)은
봉투에서 빠지고 1.119 가 된다. Phase 0 의 결론(막는 것은 변위)은 원래 τ 와 무관하므로
그대로다. subject10 이 들어왔으니 봉투는 7명으로 다시 뽑는다
(`screen_modules.py --id-dir out/real --data data --exclude-subject subject6
--exclude-subject subject9`).

### 시너지 5개를 새 근육 5개로 1:1 만들 수는 없다 — 양의 스패닝 검사

"시너지당 새 근육 하나씩 달아서 걷게 할 수 있나"에 대한 답.

이론부터 막힌다. 당김 전용 액추에이터 n+1 개가 R^n 을 양으로 덮어야 임의 방향 토크를
낼 수 있는데, 6 DOF 면 7개다. 5개는 subtalar·hip_rotation 을 잠근 4 DOF 에서만 n+1 이다.

그래서 4 DOF(고관절 굴곡·내전, 무릎, 발목)로 줄여서 실측했다. 시너지 i 의 토크
방향벡터 = 근육 j 의 SO 힘을 `C_i W_ij / Σ_k C_k W_kj` 로 귀속시켜 모멘트암을 곱한 뒤
`C_i` 로 가중 평균한 것. n+1 개 벡터의 양의 스패닝은 rank = n 이고 영공간 벡터의
부호가 전부 같은지로 판정한다 (`report/spanning_check.py`).

| | 결과 |
| --- | --- |
| 6 subject × 좌우 12개 분해 | **양의 스패닝 성립 1 / 12** |
| 성립한 하나 (subject3 walking1 r) | S1 [+5.5, −29.9, −5.7, −26.3] S2 [+35.7, −4.2, −0.5, −38.6] S3 [−29.7, +4.7, +15.5, +6.6] S4 [+20.8, −28.3, +6.9, −52.0] S5 [+9.6, +2.2, −5.9, +2.9] Nm |

S1·S2·S4 가 전부 저굴+외전 쪽을 향한다. 걷기에 필요한 방향은 촘촘하고, 걷기가
요구하지 않은 방향(무릎 신전+배굴, 내전+저굴 …)은 비어 있다. **시너지는 "걷는 데 쓰는
방향"이지 "관절을 임의로 움직이는 데 필요한 방향"이 아니다.** 5개 근육으로 만든 로봇은
그 걸음은 재현해도, 균형을 잃는 방향의 외란에는 낼 토크가 없다.

결론은 두 숫자를 나누는 것이다.

- **하드웨어 근육 수 K** — 양의 스패닝(6 DOF 면 ≥7)·스트로크·힘으로 정한다. Phase 0~5.
  Phase 3 의 K=8 후보는 바닥(7) 바로 위라 타당하다.
- **제어 차원 k=5** — 시너지. `SynergyController` 에 W 를 넣어 K 개 근육을 5개 신호로
  구동한다. "8개가 안정적"과 "5개"는 다른 것을 세고 있었고 둘 다 맞다.

걷기 재현만이 목적이면 5개짜리 시뮬레이션은 된다(SynergyController 상한 검증 → 5개
lumped 액추에이터 + Moco 추적). 단 걷기 부분공간 밖으로는 못 나간다는 걸 알고 해야 한다.
그리고 S2·S4 의 저굴 중복(soleus/gasmed 를 SO 가 갈라놓은 것)은 근육 묶기 뒤 다시 뽑는다.

### 새 근육 5개는 걸음을 3분의 1에서 못 낸다, 8개는 낸다 — 준정적 검사

양의 스패닝(위)은 "임의 방향"을 묻는 검사라 걷기만 놓고 보면 과할 수 있다. 그래서 걷기
자체를 물었다. 새 근육 i = 고정 토크 방향벡터 V_i(4 DOF) × 흥분 u_i ≥ 0 으로 두고,
프레임마다 `min ‖Σ u_i V_i − τ_ID‖, u ≥ 0` (NNLS)을 풀어 **남는 토크**를 쟀다.
남는 토크 = 보조 액추에이터가 떠맡아야 할 몫 = 하드웨어가 못 내는 부분
(`report/lumped_actuator_test.py`, 7명 42 trial × 좌우).

| 새 근육 | 남는 토크 RMS/피크 (hip_flex · hip_add · knee · ankle) | 최악 프레임 | 관절 하나라도 >10 % 인 분해 |
| --- | --- | --- | --- |
| **5개** | 3.2 · 5.1 · **7.7** · 2.3 % | **32 %** | **36 %** |
| **8개** | 0.2 · 0.6 · 0.7 · 0.2 % | 4 % | 0 % |

같은 사람의 첫 walking trial 로 만든 근육을 나머지 trial(walkingTS 포함)에 써도 거의 같다
(5개: 3.9 · 6.7 · 9.8 · 3.1 %, 최악 37 %, >10 % 39 % / 8개: 0.0 · 0.9 · 0.6 · 0.2 %, 4 %).
**사람 안에서는 설계가 옮겨진다** — 몸통을 흔들어도 다리 근육 설계는 같다.

무릎이 제일 못 낸다. 5개 시너지엔 무릎 신전이 "체중 수용" 성분에 묻혀 있어서 단독으로
못 뽑는다. 그래서 5 → 8 로 갈 때 얻는 3개가 대체로 무릎·고관절 내전 방향의 독립 성분이다.
Phase 3 의 K=8 후보와 정확히 맞아떨어진다.

한계: 준정적(근육 동역학 없음), 흥분 상한 없음, 4 DOF. 동역학까지 넣은 답은
`report/moco_synergy_inverse.py`(MocoInverse + SynergyController)가 준다.

### 새 근육 8개를 실제로 달아서 걸었다 — MocoInverse (2026-09-14)

준정적 검사(위)를 진짜 시뮬레이션으로 올렸다. subject3 walking1, 운동학·지면반력은 실측
고정, 다물체 동역학 전부. 두 종류를 돌렸다 (`report/moco/`).

**A. 사람 근육 80개 + 시너지 신호 k개** (`moco_synergy_inverse.py`, `SynergyController`, 근육
동역학 포함)

| | 수렴 | 보조 최악 | 좌표별 |
| --- | --- | --- | --- |
| 기준선: 근육 80개 자유 제어 | 43회 0.8분 | **0.3 %** | 전부 < 0.3 % |
| 신호 8개/다리 | 47회 1.7분 | 21.8 % (무릎 r) | 무릎 14~22, 고관절 굴곡 4~11, 나머지 3~7 % |
| 신호 5개/다리 | 70회 2.0분 | 85 % (무릎 r) | 무릎 59~85, 고관절 내전 10~33, 발목 14~19 % |

**B. 사람 근육 전부 떼고 새 근육 k개/다리** (`moco_new_muscles.py`, 토크 비율 고정 액추에이터)

| | 수렴 | 보조 최악 | 좌표별 | 신호 최대 (시너지 피크=1) |
| --- | --- | --- | --- | --- |
| **새 근육 8개** | 14회 0.1분 | **8.3 % (무릎 r)** | 나머지 전부 ≤ 2 % | 4.0× |
| 새 근육 5개 | 13회 0.1분 | 48 % (무릎 l) | 무릎 23~48, 고관절 5~10, 발목 4~6 % | 7.8× |

**새 근육 8개면 걸음이 재현된다. 5개는 무릎이 안 된다.** 준정적 NNLS 예측(5개: 무릎
부족, 8개: 됨)과 정확히 같은 그림이고, 신호 상한을 3 → 10으로 풀어도 안 바뀌므로 용량이
아니라 방향(스패닝)의 문제다. `report/fig_newmuscles_walk_w0_u10.png`.

A 가 B 보다 나쁜 이유: A 는 시너지 재구성 오차(VAF 90~96 %)에 근육 동역학(활성화 지연,
힘-길이-속도)까지 얹혀서 무릎처럼 작은 모멘트(27 Nm)가 상대적으로 크게 틀린다. B 는
토크 비율만 고정하고 동역학이 없어 순수하게 "방향이 충분한가"를 본다. 하드웨어 질문엔 B 가
맞는 검사고, "사람 근육을 5개 신호로 몰 수 있나"엔 A 가 맞는 검사다.

한계는 그대로다: 토크형(경로 없음), 4 DOF(hip_rotation·subtalar 는 보조가 담당, 계산에 안
넣음), 운동학·GRF 고정(균형 아님), V 는 같은 trial 로 만든 것(in-sample). 다음은 cross-trial
V 와 접촉 모델.

**Moco 가 걸린 지점 다섯** — 전부 조용히 '불가능'으로 끝난다.
1. AddBiomechanics 모델의 잠긴 좌표(subtalar·mtp·손목)는 Moco 가 거부한다 → 그 좌표의
   CoordinateActuator 를 먼저 떼고 관절을 WeldJoint 로 바꾼다(`strip_locked`).
2. nimble IK 에는 `knee_angle_*_beta` 가 없다 → 커플러 함수로 채운다(`complete_kinematics`).
3. 시너지 가중치의 비정규 값(4.9e-324)은 XML 왕복에서 속성을 통째로 비운다 → 1e-8 미만은 0.
4. **`Muscle.min_control` 기본값 0.01** — 시너지 가중치가 0인 근육은 흥분이 0 이라 이 하한을
   영원히 못 지켜 k 와 무관하게 같은 오차(0.833)로 실패한다 → `setMinControl(0)`, 그리고
   모든 근육을 시너지에 넣는다(`--min-peak 0`).
5. 새 근육(CoordinateActuator)의 출력 토크가 보조와 같은 제곱 비용으로 목적함수에 들어가면
   둘을 정확히 50/50 으로 나눈다 → 새 근육 출력의 가중치 0, 비용은 신호 u 와 보조에만.
6. 파이썬 임시객체(`TableProcessor | TabOp`)를 C++ 에 바로 넘기면 segfault → 변수에 담아 `append`.

### 새 근육 8개 — 방향을 고정하고 다른 걸음을 걷게 했다 (2026-09-14 오후)

위 B 의 V 는 walking1 자기 것(in-sample)이었다. 하드웨어는 걸음마다 바꿀 수 없으니 **V 를
고정하고 다른 걸음·다른 사람을 걷게** 했다 (`moco_new_muscles.py --v-from`, 결과
`report/moco/results_newmus.csv`, 그림 `report/fig_newmuscles_cross.png`). 두 가지가 나왔다.

**1. 보조 가중치 1 은 검사가 아니다.** 목적함수가 Σu² + Σ(보조 Nm)² 라서, 새 근육이 낼
수 있어도 u 가 커지면 보조에 떠넘기는 게 싸다. walking1 V 를 walking2/3 에 쓰면 22~26 %,
심지어 자기 V 로도 walkingTS2 는 34 % 가 나왔는데 — 준정적 NNLS 는 같은 V 로 잔여 **0.0 %**
였다. 8개 V 는 4 DOF 를 양의 생성하므로(rank 4, 부호 검사 통과) 어떤 토크든 u ≥ 0 으로
정확히 낸다. 남는 건 **필요한 u 크기** 뿐이다. → `--reserve-weight 100`: 보조는 정말 못
내는 몫만 남는다. walking1 자기 V 가 0.4 %, walking2/3 은 8~10 %, subject4 는 15 %.

**2. 걸음 하나의 NMF 방향은 조건이 나쁘다.** walking1 V 로 다른 걸음을 내려면 오른다리
u 피크가 46(시너지 피크의 46배)까지 필요했다 (왼다리는 5~10). 무릎 신전 방향이 큰
벡터 두 개의 상쇄로만 나오기 때문. subject3 **여섯 걸음을 이어 붙여**(`lat.load_pooled`,
`--v-from subject3 ALL`, 근육 집합 통일하려고 min_peak 0) V 를 뽑으니 u 피크가 10 안으로
들어왔다 (오른 8~10, 왼 3~6).

| V 출처 (k=8, u ≤ 10) | s3 walking1 | walking2 | walking3 | TS2 | TS3 | TS4 | s4 walking1 |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 자기 걸음, 보조 가중치 1 | 8.2 | 5.6 | 13.5 | 33.9 | – | – | 13.5 |
| walking1 고정, 가중치 1 | (8.2) | 21.9 | 25.9 | 17.8 | – | – | 21.9 |
| walking1 고정, 가중치 100 | 0.4 | 7.9 | 9.7 | 0.9 | – | – | 14.7 |
| **6걸음 합침, 가중치 100** | **0.2** | **0.1** | **0.2** | **0.2** | **0.2** | **0.2** | **3.0** |
| 합침 k=5, 가중치 100 | 97.3 | | | | | | |

(값 = 보조 액추에이터 최악, % of ID 피크.) **합친 방향 8개면 피험자3의 여섯 걸음(trunk sway
포함)을 전부 0.2 % 이내로, 다른 사람(subject4)도 3 % 로 낸다. 5개는 97 %.** 이 V 가 현재의
설계안이다 — `report/newmus_V_subject3_pooled_k8.npz` (V_r, V_l: 8×4, 단위 u 당 Nm,
열 = hip_flexion, hip_adduction, knee_angle, ankle_angle).

읽는 법: V 의 행 하나가 새 근육 하나의 "토크 서명"이다. 예컨대 V_r 2행 (28, −10, 3, −43) 은
고관절 굴곡 + 발목 저굴 — 사람 근육엔 없는 조합이지만 push-off 시너지가 그렇게 생겼다.
3행 (−1, 0.3, −2, −1) 처럼 노름 3 Nm 짜리 약한 행도 있다 — 근육 하나가 놀고 있는 셈이라,
방향을 NMF 에 맡기지 말고 **필요한 u 를 최소화하도록 직접 최적화**하면 더 좋아질 여지가 있다.

교훈: Moco 의 "보조 %" 는 비용 설계에 따라 10배 넘게 변한다. 실현 가능성 질문엔 보조
가중치를 크게 두고(100), 준정적 NNLS 의 잔여·u 피크와 같이 봐야 한다.

영상: OpenSim GUI 가 안 열려서 `report/newmus_video.py` 로 만들었다 — 스틱피겨 + GRF + 이
순간의 새 근육 신호 8개 막대 + 토크 추적, 1/3 배속 3회 반복, mp4 (ffmpeg 은
`imageio-ffmpeg` 가 들고 온 바이너리). `report/video_newmus_subject3_walking1_k8_vALL_rw100.mp4`,
`..._walkingTS2_...`. `conda run` 을 셸 for 루프 안에서 부르면 조용히 죽는다 — 스크립트
파일로 env 의 python 절대경로를 직접 부른다.

### 실측 EMG 로 SO 를 검산했다 — 발목·햄스트링은 맞고, 대퇴직근·전경골근은 틀리다 (2026-09-14)

OpenCap 논문 원자료 `LabValidation_withoutVideos`(`~/Desktop/scHooL/의료로봇연구실/`, 7.2 GB)
를 받았다. 피험자마다 `EMGData/<trial>_EMG.sto`: 16채널 100 Hz 포락선(왼쪽 9: soleus gasmed
tibant recfem vasmed vaslat semiten bflh glmed1, 오른쪽 7: tibant·recfem 없음), 시간축이 AddB
trial 과 같다(둘 다 0 부터, 132행). 그 밖에 `ForceData`(2000 Hz 원본 힘판), `MarkerData/Mocap`
(.trc), `OpenSimData/Mocap/{IK,ID,SO,JR,Model}` — 논문 팀이 돌린 결과. **걷기 외 동작
(DJ 드롭점프, squats, STS 앉았다서기)의 IK·ID 도 들어 있다.** 검산은 `report/emg_check.py`
(`report/emg_check.csv`, `emg_check_synergy.csv`, `fig_emg_subject3_walking1.png`,
`fig_emg_summary.png`). 7명 × 6걸음 × 16근육 = 672쌍, SO 창 안에서 둘 다 피크 정규화.

| 근육 | r (지연 0) | r (EMG 를 0~120 ms 앞당긴 최선) | 최선 지연 |
| --- | --- | --- | --- |
| soleus | 0.74 | 0.78 | 25 ms |
| semiten | 0.64 | 0.68 | 17 |
| bflh | 0.59 | 0.60 | 6 |
| glmed1 | 0.45 | 0.58 | 49 |
| gasmed | 0.50 | 0.57 | 41 |
| vaslat / vasmed | 0.44 / 0.42 | 0.50 / 0.48 | 38 |
| **recfem** | 0.04 | 0.28 | 74 |
| **tibant** | −0.01 | 0.10 | 59 |

읽는 법. 발목 저굴근과 햄스트링은 SO 가 실제 근육 활동을 따라간다(지연 6~25 ms 는 전기역학
지연으로 정상). 대퇴사두는 반쯤. **대퇴직근과 전경골근은 SO 가 사실상 안 켠다** — 유각기
발등굽힘 토크가 작아 SO 엔 tibant 가 필요 없고, 동시수축(recfem 은 stance 초 vasti 와 같이
켜진다)은 최소 비용 해에 없다. 그림에서 보이는 두 가지 SO 특유 모양: 스파이크(프레임별
독립 최적화라 이웃 프레임과 상관없이 튄다)와 **1.0 에서 평평한 포화**(gasmed_l 처럼 힘이
모자라 한계에 붙는다 = 임의 배분·용량 문제).

시너지 수. 같은 근육 부분집합(한쪽 7~9개)으로 EMG 는 VAF 90 % 에 **2~3개**, SO 는 3~5개
(84 분해). SO 가 실제보다 복잡하다 — 스파이크·포화가 시너지를 더 먹는다. 그래도 SO 전체
근육 k=5 의 시간곡선 C 와 EMG 시너지 C 를 최대상관으로 짝지으면 평균 0.73 — 시너지의
**타이밍**은 SO 에 살아 있다.

새 근육 설계에 대한 뜻. 방향 V 는 SO 힘 F 를 시너지에 귀속시켜 만든 것이라 recfem·tibant
의 오류가 들어가 있다. 다만 V·u 는 어차피 ID 토크(실측)를 맞춰야 하므로 SO 는 "토크를
어느 근육에 묶는가"에만 관여한다. **다음 단계인 "V 를 ID 토크에서 직접 최적화"로 가면 SO
의존이 사라지고**, EMG 는 그때 u(t) 의 타이밍이 사람 시너지와 맞는지 보는 잣대로 쓴다.

### 새 근육 8개의 방향을 ID 토크에서 직접 설계했다 (2026-09-14 저녁)

SO 를 거치지 않는다. `report/design_directions.py`: 단위 방향 8개와 용량 c_i 를 골라, 모든 프레임의
관절 토크를 0 ≤ u ≤ c 로 정확히 내면서 Σc_i 최소. 안쪽(V 고정)은 LP(HiGHS), 바깥(V 32개 값)은
Powell 다중 시작 5개(축 정렬, NMF, 무작위 3). 데이터 46,006 프레임 = 7명 걷기 42(우리 ID, 8,616)
+ 논문 Mocap ID 의 squats1(18,784)·STS1(18,606), 좌우 합침, τ ÷ (질량×키) × (63.5×1.69) 로
피험자3 크기. 결과 `report/newmus_V_design_k8.npz` (V 단위방향, c 용량, V_nm = V·c),
`design_directions.csv`, `fig_design_directions.png`.

| 방향 8개 | 총용량 Σc (Nm, 3,000 프레임 표본) |
| --- | --- |
| 축 정렬 (관절마다 굴근·신근) | 532 |
| NMF (SO 시너지, subject3 합침) | 513 |
| **직접 최적화 (시작 2)** | **394** (−26 %) |

과제별 필요 용량(설계 V 고정): 걷기 387, 스쿼트 281, STS 302 → **걷기가 용량을 정한다.**
설계 용량으로 42걸음 전 프레임을 프레임별 유계 최소제곱으로 검산하면 못 내는 토크 최대 5.8 Nm
(무릎 피크의 6 %). 근육별 용량 33~77 Nm, 모멘트암 5 cm 면 650~1,540 N (KIMM 초대형 1개 범위).
관절 방향별 "필요 vs 가능": 고관절 굴곡 +70/−83 vs +88/−105, 내전 +31/−66 vs +58/−80, 무릎
+42/−99 vs +72/−105, 발목 +22/−119 vs +38/−131 Nm — 전부 여유 있음.

설계된 방향(행 = 근육, 열 = 고관절 굴곡·내전·무릎·발목, Nm): M2 (51, −17, −32, 4) 는 고관절
굴곡 + 무릎 신전 = 대퇴직근 꼴, M3 (8, −11, 9, −55) 는 발목 저굴 = 비복근 꼴, M6 (−34, −12, −14,
−5) 는 고관절 신전 = 대둔근 꼴, M7 (14, 1, 42, −21) 은 무릎 굴곡 + 저굴 = 햄스트링 꼴. 토크만
보고 설계했는데 사람 근육과 닮은 게 나온다. M0 (−49, 45, −30, 24) 처럼 사람에 없는 것도 있다.

**Moco 검증에서 잡힌 것 — MocoInverse 내부 토크가 ID 와 어긋나는 노드가 있다.** 설계 V 로
7명 12걸음을 돌리자 3·4·7·8·11번은 재현되는데 5·10번이 왼다리에서 100 % 넘게 실패했다. 추적:
부호 규약 동일(모델 회전축 전 피험자 같음), 용량 충분(준정적 0 Nm 부족), 격자 0.01·미리 필터·
골반 보조 가중치 전부 무효. 결정적 실험: **용량 무제한 축 정렬 액추에이터**로 돌려도 10번 왼
고관절이 0.963~0.983 s 에서 +120/−32 Nm 로 튄다(ID 26 Nm). 운동학 고정이면 토크는 유일하므로
이건 Moco 의 내부 ID(스플라인·격자·외력 보간)가 그 노드에서 틀린 것이다. 오른발 접촉 시작
(힘 0 → 28 N, COP (0,0,0) → (1.5, 0, −0.2) 점프) 시각과 일치. 12걸음 중 5걸음에 이런 노드가
있다(최대 24/97). 그래서 최종 검증 지표는 **보조 % 가 아니라 "새 근육 토크 vs 우리 ID"**, 자유
액추에이터로도 5 Nm 넘게 어긋난 노드는 뺀다 (`report/design_validation.csv`).

| 걸음 (u ≤ 1, 몸 크기 비율) | 어긋난 노드 | 설계 오차 최악 (% of ID 피크) |
| --- | --- | --- |
| s3 walking1/2/3, TS3, TS4 | 0 | 0.2~0.9 |
| s3 walkingTS2 | 10/147 | 19.6 (무릎, 어긋난 구간 옆) |
| s4, s7, s8 walking1 | 0~1 | 6.7 / 7.1 / 9.7 |
| s10 walking1 | 11/85 | 6.6 |
| s5 walking1 (1.37배 몸) | 24/97 | 5.9 |
| s11 walking2 (1.59배 몸) | 7/97 | 4.0 |

**한 설계로 7명이 걷는다.** 보조 % 만 봤다면 5·10번을 설계 실패로 잘못 읽었을 것이다. 남은 것:
그 Moco 노드 튐의 정확한 원인(외력 보간 vs 무릎 커플러 제약 λ) — 접촉 모델로 넘어가기 전에 잡는다.

`conda run` 이 루프에서 조용히 죽던 원인도 찾았다: **zsh 는 `$var` 를 단어로 나누지 않는다.**
`set -- $run` 에서 `$1` 이 문장 전체가 되어 인자가 깨졌다. `${=run}` 을 쓴다.

### 케이블 경로 근육 8개로 걸었다 — 부착점 최적화 1차 (2026-09-14 밤)

`report/path_design.py`: 설계 방향 v_i 에 대해 케이블의 모멘트암 벡터 r(q) = −∂L/∂q(중앙차분,
δ = 0.01 rad)가 자세 표본 120개(걷기 6걸음 각 10 + 스쿼트 30 + STS 30)에서 v_i 와 나란하도록
부착점을 고른다. 뼈 변환(pelvis·femur·tibia·calcn, 자세 × 9 섭동)을 한 번 저장해 두고 길이는
numpy 로만 계산하므로 평가가 마이크로초 → differential evolution(popsize 25, 400세대) + Powell
다듬기, 근육당 18초. 건너는 관절은 |v 성분| > 0.12 로 정하고 바디마다 점 하나(시작·경유·끝),
허용 상자 `BOX` 안. 벌점: 스트로크 > 0.40, 필요 힘(1.1 c_i ÷ 최소 유효 모멘트암) > 5,000 N,
유효 모멘트암 < 1.5 cm. 결과 `report/newmus_paths_k8.json`, 모델
`report/moco/newmus_paths_subject3_model.osim`(PathActuator 8×2, 왼쪽은 z 뒤집음), 그림
`fig_path_design.png`(방향 vs 모멘트암), `fig_path_layout.png`(자세 3개에서의 경로).

| 근육 | 경로 | 방향 각도차 | 유효 모멘트암 | 필요 힘 | 스트로크 |
| --- | --- | --- | --- | --- | --- |
| M0 | pelvis→calcn (4관절 성분) | 0.3° | 2.0~3.9 cm | 4,219 N | 0.11 |
| M1 | pelvis→calcn | 0.8° | 3.3~3.8 | 1,153 | 0.03 |
| M2 | pelvis→tibia (대퇴직근 꼴) | 4.1° | 1.5~3.7 | 4,591 | 0.08 |
| M3 | pelvis→calcn (비복근 꼴) | 0.3° | 3.1~4.1 | 2,054 | 0.04 |
| M4 | pelvis→calcn | 0.7° | 3.7~5.0 | 1,244 | 0.13 |
| M5 | pelvis→calcn | 0.4° | 3.2~3.5 | 1,107 | 0.06 |
| M6 | pelvis→calcn (대둔근 꼴) | 1.6° | 2.1~3.8 | 2,043 | 0.11 |
| M7 | pelvis→calcn (햄스트링 꼴) | 0.8° | 4.4~7.3 | 1,238 | 0.18 |

전부 5,000 N 안, 스트로크 0.40 안. 힘이 큰 M0·M2 는 유효 모멘트암이 2 cm 아래로 떨어지는
자세가 있어서다.

**MocoInverse 검증** (`report/moco_path_muscles.py`, PathActuator 제어 0~1, 보조 가중치 100,
골반 보조 0.01, 지표 = 오른다리 Σ F_i u_i r_i(q) vs 우리 ID, 창 양 끝 0.05 s 제외):

| subject3 | walking1 | walking2 | walking3 | TS2 | TS3 | TS4 |
| --- | --- | --- | --- | --- | --- | --- |
| 새 근육 토크 vs ID 최악 | 6.7 % | 5.3 | 7.9 | 23.2 (Moco 튐 노드 10개 제외; 포함하면 104) | 7.9 | 10.1 |
| 보조 최악 | 4.8 % | 5.1 | 6.8 | 6.7 | 7.5 | 10.3 |

**진짜 케이블 근육(모멘트암이 자세 따라 변함)으로도 6걸음 중 5걸음이 10 % 안.** 방향 설계의
"고정 비율" 가정이 걷기 범위에서는 잘 맞는다는 뜻. TS2 는 Moco 튐 노드 옆이 남는다.

**1차의 한계 — 다음 반복에서 고칠 것.**
1. **뼈 간섭을 안 봤다.** 감싸는 면(wrap)이 없고 직선 구간뿐이라, `fig_path_layout.png` 의 스쿼트
   자세에서 케이블이 뼈를 가로지른다. 경유점이 무릎 중심 2 cm 안에 놓인 근육이 여럿(femur
   y ≈ −0.40~−0.42, tibia y ≈ −0.02) — 모멘트암을 작은 오프셋으로 만든 것이라 실제로는 관절
   축을 관통하는 배치다. 관절 중심에서 최소 거리(예: 3 cm) 제약과 wrap cylinder 를 넣어야 한다.
2. **7/8 이 4바디를 다 건넌다.** 문턱 0.12 가 낮아 작은 성분까지 건너게 했다. 건너는 관절 수에
   벌점을 주거나 문턱을 0.25 로 올리면 2~3바디 경로가 나올 것이고, 만들기 쉬워진다.
3. 좌우 대칭(z 뒤집기)만 했고 왼다리 모멘트암은 따로 검증하지 않았다(보조 % 로만).
4. 피험자3 뼈대 치수 기준. 로봇 뼈대가 오면 상자 `BOX` 를 그 치수로 바꾸고 다시 돌린다.

GUI 로 보기: `report/moco/newmus_paths_subject3_model.osim` + `report/moco/paths_subject3_walking1_ik_full.mot`.
PathActuator 라 색은 고정(근육별 색), 활성도 색은 Muscle 로 바꿔야 나온다.

### 케이블 근육으로 스쿼트도 했다 — 그리고 입력은 반드시 같은 모델 짝이어야 한다 (2026-09-14 밤)

설계에 스쿼트가 들어갔으니 케이블 근육으로 스쿼트를 재현해 봤다. **처음엔 460 % 로 실패**했다:
논문 팀의 IK·힘판 데이터(`LabValidation`)를 우리 AddB 스케일 모델에 그대로 넣었기 때문이다.
두 모델은 같은 사람이라도 스케일·발 기하가 달라 발 위치와 힘판 COP 가 몇 cm 어긋나고, 그게
발목 토크로 크게 튄다(발목 ID 피크 27 Nm 인데 오차 120 Nm). 해결: subject3 `.b3d` 에서
squats1 을 우리 모델 짝으로 변환(`b3d_to_opensim.py`, 1,184 프레임 전부 양호, 잔차 7 N) →
`run_trial.py` 로 ID → `moco_path_muscles.py --t0 2.6 --t1 4.4`(두 번째 반복). 결과 오른다리
케이블 토크 vs ID 최악 **5.6 %** (고관절 굴곡 2.7, 내전 5.6, 무릎 0.9, 발목 4.2), 왼쪽 고관절만
보조 22 %(반복 상한 도달, 미수렴). 영상 `report/video_paths_subject3_squats1.mp4`
(`report/paths_video.py`: 케이블을 신호 u 로 진하기·굵기 표시).

STS1(앉았다 서기)은 AddB 에서 절반이 `unmeasuredExternalForceDetected`(의자가 미는 힘이 힘판에
없음)라 잘렸지만, **남은 9.16~10.81 s 는 엉덩이가 의자에서 떨어진 뒤 완전히 일어서는 구간**(골반
높이 0.56 → 0.94 m, 잔차 1 %)이라 쓸 수 있다. 케이블 근육으로 재현: 오른다리 토크 오차 최악
**3.8 %**, 보조 4.2 %. 영상 `report/video_paths_subject3_STS1.mp4`. 설계 데이터에 넣은 논문 팀
STS ID 가 의자 힘을 포함한 값인지는 여전히 확인 필요(열린 질문).

GUI 로 보기 (2026-09-15): `report/gui_export.py` 가 PathActuator 를 같은 경로의 Millard 근육으로
바꾼 `report/moco/gui_muscles_model.osim`(GUI 는 Muscle 만 활성도 색을 칠한다)과, 관절 각도 +
`/forceset/M*_*/activation`(= Moco 신호 u) 을 합친 `gui_<name>.mot` 을 만든다. 걷기는 지지 구간만
이 아니라 시행 전체(0.03~1.27 s, `paths_subject3_walking1_full`)로 다시 돌렸다 — 보조 최악 8 %.
이때 "토크 vs ID" 지표는 ID 창(0.25~1.12) 밖에서 무의미하다(np.interp 가 끝값을 붙듦).

교훈: 운동학·힘판·ID 는 **한 모델에서 나온 짝**으로만 쓴다. 다른 출처를 섞으면 발목부터 깨진다.
`moco_path_muscles.py` 에 `--ik --ext --id --t0 --t1 --name` 을 넣어 임의 입력을 받게 했지만,
그 입력이 `--model` 과 짝인지는 사람이 확인해야 한다.

### 계획 파일이 없다

README 가 `PLAN_lowerlimb_muscle_placement.md` 를 참조하지만 프로젝트 어디에도 없다
(`~/.claude/plans/` 에도 9/8 자 파일이 없다). Phase 0~5 의 정의는 이 절과
`screen_modules.py`·`actuator_model.py`·`models.py` 의 docstring, 그리고
`인공근육휴머노이드/artificial_muscle_design_validation_flow.svg` 의 흐름
(SO → 시너지 → 모델 설계 변경 → CMC 추적 검증 → 성공/실패 → Moco 새 패턴 → 반복)에서
복원한다. 다시 쓸 때는 프로젝트 안에 둔다.

## 5. 기본값

바꾸려면 각 파일 상단의 상수를 고친다.

| 상수 | 값 | 파일 |
| --- | --- | --- |
| `LOWPASS` | 6.0 Hz | `run_trial.py`, `run_muscle_analysis.py`, `run_static_optimization.py` |
| `EXCLUDE` | `Muscles`, `Actuators` | `run_trial.py` (ID에서 제외할 힘) |
| `SUPPORT_THRESHOLD` | 0.70 BW | `addb_io.py` (지지 구간 문턱) |
| `SUPPORT_GAP_MAX` / `SUPPORT_GAP_FLOOR` | 0.25 s / 0.50 BW | `addb_io.py` (midstance 골 잇기) |
| `CONTACT_THRESHOLD` | 0.10 BW | `addb_io.py` (매핑에 쓸 프레임) |
| `DEFAULT_COORDS` | 발목·무릎·고관절 좌우 6개 | `run_muscle_analysis.py` |
| `ANTAGONISTS` | 길항근 쌍 정의 (전체 / 대표) | `run_muscle_analysis.py` |
| `RESERVE_OPTIMAL_FORCE` / `RESERVE_LIMIT` | 1.0 / 1000 | `run_static_optimization.py` |
| `VAF_TOTAL` / `VAF_MUSCLE` | 0.90 / 0.75 | `synergy.py` |
| `K_MAX` | 10 | `synergy.py` |
| `MIN_PEAK` | 0.05 | `synergy.py` (조용한 근육 제외 기준) |
| `WEIGHT_SHOW` | 0.30 | `synergy.py` (표에 표시할 가중치 하한) |
| `N_RESTARTS` / `MAX_ITER` | 20 / 3000 | `synergy.py` (NMF 재시작·반복) |
| `MIN_DURATION` | 0.5 s | `batch.py` (정적 자세 파일 걸러내기) |
| `WALKING_TOKENS` | `walk`, `run`, `gait`, `tread` | `b3d_to_opensim.py` |
| `SOFT_REASONS` | `manualReview` | `b3d_to_opensim.py` (`--allow-manual-review`가 눈감는 사유) |
| `BANDS` | 검산 정상범위 | `validate.py` |
| `EPS_MAX` / `EPS_PASSIVE` | 0.40 / 0.10 | `actuator_model.py` (스트로크·수동 강성 변형률) |
| `RHO_F_TARGET` | 15 N/g (민감도 5.4) | `actuator_model.py` (힘밀도) |
| `PASSIVE_CAP_RATIO` | 1.0 | `actuator_model.py` (수동 힘 상한, 발산 방지) |
| `MODULES` | 표준 모듈 6등급 | `actuator_model.py` (발표자료 p.35~36) |
| `ROM_TARGET_DEG` | AAOS 기준 임시값 | `models.py` (뼈대 팀과 확정 전) |
| `REPRESENTATIVE_ARM_M` | 대표 근육 실측 모멘트암 | `models.py` (NOTES 2절) |
| `ROBOT_MASSES_KG` | 40 / 55 / 70 | `screen_modules.py` (로봇 질량 미정) |
