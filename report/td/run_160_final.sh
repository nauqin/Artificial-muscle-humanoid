#!/usr/bin/env bash
# 2026-11-21 175 cm 시간 영역 최종 조합(report/td/best_set.md)을 160 cm 에 그대로 적용.
#   175 cm 최종: 하체 td_br20 · 몸통 td · 어깨 td_br20_rb · 팔꿈치 td_br50_rb.
#   160 cm 에 이미 있는 것: 어깨 td_br20_rb, 팔꿈치 td_br50_rb (report/td/run_160.sh, summary_160.md — Moco 통과). 다시 최적화하지 않는다.
#   160 cm 에 없는 것: 하체·몸통 → 여기서 175 cm 와 같은 설정(report/td/run_175.sh)으로 만든다. 하체는 비교용으로 td(브래킷 없음)도 같이.
#   bash report/td/run_160_final.sh design   # 설계 (하체 2 + 몸통 1, 병렬)
#   bash report/td/run_160_final.sh moco     # Moco 검증 (하체: 걷기·실험실 걷기·스쿼트 / 몸통: 걷기 2·들기 2)
cd "$(dirname "$0")/../.."
PY=${PY:-/c/Users/USER/miniforge3/envs/osim/python.exe}; export PYTHONUTF8=1 MKL_THREADING_LAYER=SEQUENTIAL
M=data/robot40_h160/robot40_h160.osim; T=report/upper/tasks_40kg_h160.csv; B=0.9877; L=report/td/logs; mkdir -p $L
C="--model $M --box-scale $B --force-mode time --refine-trust 0.01"
LOW="--space lower --design report/newmus_V_design_lower_k8_40kg_h160.npz --ref 40 1.60 $C --stage1-stride 4"
TRK="--space trunk --tasks $T --design report/newmus_V_design_trunk_k4_40kg_h160.npz $C --stage1-stride 4"
X="bash report/modulefirst/run_moco.sh"; R=report; MO=report/moco; F="\[warning\]|\[info\]"
case ${1:-design} in
  design)
    $PY report/path_design.py $LOW --out-tag _40kg_h160_td > $L/pd160_lower_td.log 2>&1 &
    $PY report/path_design.py $LOW --bracket 0.02 --out-tag _40kg_h160_td_br20 > $L/pd160_lower_td_br20.log 2>&1 &
    $PY report/path_design.py $TRK --out-tag _40kg_h160_td > $L/pd160_trunk_td.log 2>&1 &
    wait; grep -h "시간 영역:\|2단계(시간 영역) 끝" $L/pd160_lower_td.log $L/pd160_lower_td_br20.log $L/pd160_trunk_td.log ;;
  nom1)   # 2026-11-21: M1 을 빼고(등급 M2~M6) 하체 td_br20 · 어깨 td_br20_rb 를 다시 — 수축이 90 mm 를 넘어 M1 이 배정된 케이블을 없애려는 것. 꼬리표 _nom1
    G="--grades M2,M3,M4,M5,M6"
    SHO="--space shoulder_r --tasks $T --design report/newmus_V_design_shoulder_r_k8_40kg_h160.npz $C --robust-frac 0.25 --dt 0.04 --alloc-iters 150"
    $PY report/path_design.py $LOW $G --bracket 0.02 --out-tag _40kg_h160_td_br20_nom1 > $L/pd160_lower_td_br20_nom1.log 2>&1 &
    $PY report/path_design.py $SHO $G --bracket 0.02 --out-tag _40kg_h160_td_br20_rb_nom1 > $L/pd160_shoulder_td_br20_rb_nom1.log 2>&1 &
    wait; tail -10 $L/pd160_lower_td_br20_nom1.log | grep "  M[0-7]:"; tail -10 $L/pd160_shoulder_td_br20_rb_nom1.log | grep "  M[0-7]:"
    st=lower_k8_40kg_h160_td_br20_nom1; nm=td160_lower_td_br20_nom1
    $PY report/modulefirst/prepare_models.py lower $R/newmus_paths_$st.json 2>&1 | grep -Ev "$F"
    $X lower $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walk ${nm}_walk --prefilter &
    $X lower $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walklab ${nm}_walklab &
    $X lower $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim squat ${nm}_squat &
    st=shoulder_r_k8_40kg_h160_td_br20_rb_nom1; nm=td160_shoulder_td_br20_rb_nom1
    $PY report/modulefirst/prepare_models.py shoulder_r $R/newmus_paths_$st.json 2>&1 | grep -Ev "$F"
    $X shoulder_r $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walk ${nm}_walk --prefilter &
    $X shoulder_r $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walks9 ${nm}_walks9 --prefilter &
    $X shoulder_r $R/newmus_paths_$st.json $MO/newmus_paths_${st}_load3.79_model.osim lift ${nm}_lift &
    $X shoulder_r $R/newmus_paths_$st.json $MO/newmus_paths_${st}_load3.79_model.osim lifts6 ${nm}_lifts6 &
    wait; grep -h "_nom1_" $MO/results_paths.csv | tail -8 ;;
  moco)
    for v in td td_br20; do st=lower_k8_40kg_h160_$v; nm=td160_lower_$v
      $PY report/modulefirst/prepare_models.py lower $R/newmus_paths_$st.json 2>&1 | grep -Ev "$F"
      $X lower $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walk ${nm}_walk --prefilter &
      $X lower $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walklab ${nm}_walklab &
      $X lower $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim squat ${nm}_squat &
    done
    st=trunk_k4_40kg_h160_td; nm=td160_trunk_td
    $PY report/modulefirst/prepare_models.py trunk $R/newmus_paths_$st.json 2>&1 | grep -Ev "$F"
    $X trunk $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walk ${nm}_walk --prefilter &
    $X trunk $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walks9 ${nm}_walks9 --prefilter &
    $X trunk $R/newmus_paths_$st.json $MO/newmus_paths_${st}_load3.79_model.osim lift ${nm}_lift &
    $X trunk $R/newmus_paths_$st.json $MO/newmus_paths_${st}_load3.79_model.osim lifts6 ${nm}_lifts6 &
    wait; grep -h "td160_lower_td\|td160_trunk_td" $MO/results_paths.csv | tail -12 ;;
esac
