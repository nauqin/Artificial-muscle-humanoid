#!/usr/bin/env bash
# 2026-11-05 새 방식(모듈 먼저, _mf) + 현재 방식 신뢰 구간 민감도(_tr) 후처리와 Moco.
#   후처리: 어깨 힘 ×1.5 json·모델, 들기 손 하중 모델, stroke_check (report/modulefirst/stroke_check_*.csv)
#   Moco: 걷기 --prefilter (맥과 같은 조건), 스쿼트·들기는 그냥. 이름 mf_<부위>_new_* / mf_<부위>_curtr_*
cd "$(dirname "$0")/../.."; R=report; MO=report/moco; MF=report/modulefirst; X="bash $MF/run_moco.sh"
PY=${PY:-/c/Users/USER/miniforge3/envs/osim/python.exe}; export PYTHONUTF8=1 MKL_THREADING_LAYER=SEQUENTIAL
M=data/robot40_h160/robot40_h160.osim; T=report/upper/tasks_40kg_h160.csv; B=0.9877
F="\[warning\]|\[info\]|Warn|자세 제외"
for tag in mf s80_win_tr; do
  $PY $MF/prepare_models.py trunk   $R/newmus_paths_trunk_k4_40kg_h160_$tag.json 2>&1 | grep -Ev "$F"
  $PY $MF/prepare_models.py elbow_r $R/newmus_paths_elbow_r_k2_40kg_h160_$tag.json 2>&1 | grep -Ev "$F"
done
$PY $MF/prepare_models.py shoulder_r $R/newmus_paths_shoulder_r_k6_40kg_h160_v635_mf.json --force-scale 1.5 2>&1 | grep -Ev "$F"
$PY $MF/prepare_models.py shoulder_r $R/newmus_paths_shoulder_r_k6_40kg_h160_v635_s80_win_tr.json --force-scale 1.5 2>&1 | grep -Ev "$F"
for n in lower_k8_40kg_h160_mf lower_k8_40kg_h160_s120_ref_tr; do
  $PY report/stroke_check.py --space lower --paths $R/newmus_paths_$n.json --model $M --out $MF/stroke_check_$n.csv 2>&1 | grep -Ev "$F"; done
for x in "trunk trunk_k4_40kg_h160_mf" "elbow_r elbow_r_k2_40kg_h160_mf" "shoulder_r shoulder_r_k6_40kg_h160_v635_mf_f15" \
         "trunk trunk_k4_40kg_h160_s80_win_tr" "elbow_r elbow_r_k2_40kg_h160_s80_win_tr" "shoulder_r shoulder_r_k6_40kg_h160_v635_s80_win_tr_f15"; do set -- $x
  $PY report/stroke_check.py --space $1 --paths $R/newmus_paths_$2.json --tasks $T --model $M --box-scale $B --out $MF/stroke_check_$2.csv 2>&1 | grep -Ev "$F"; done
for pair in "new mf mf v635_mf_f15" "curtr s120_ref_tr s80_win_tr v635_s80_win_tr_f15"; do set -- $pair; w=$1; lt=$2; ut=$3; st=$4
  $X lower $R/newmus_paths_lower_k8_40kg_h160_$lt.json $MO/newmus_paths_lower_k8_40kg_h160_${lt}_model.osim walk mf_lower_${w}_walk --prefilter &
  $X lower $R/newmus_paths_lower_k8_40kg_h160_$lt.json $MO/newmus_paths_lower_k8_40kg_h160_${lt}_model.osim squat mf_lower_${w}_squat &
  $X trunk $R/newmus_paths_trunk_k4_40kg_h160_$ut.json $MO/newmus_paths_trunk_k4_40kg_h160_${ut}_model.osim walk mf_trunk_${w}_walk --prefilter &
  $X trunk $R/newmus_paths_trunk_k4_40kg_h160_$ut.json $MO/newmus_paths_trunk_k4_40kg_h160_${ut}_load3.79_model.osim lift mf_trunk_${w}_lift &
  $X elbow_r $R/newmus_paths_elbow_r_k2_40kg_h160_$ut.json $MO/newmus_paths_elbow_r_k2_40kg_h160_${ut}_model.osim walk mf_elbow_${w}_walk --prefilter &
  $X elbow_r $R/newmus_paths_elbow_r_k2_40kg_h160_$ut.json $MO/newmus_paths_elbow_r_k2_40kg_h160_${ut}_load3.79_model.osim lift mf_elbow_${w}_lift &
  $X shoulder_r $R/newmus_paths_shoulder_r_k6_40kg_h160_$st.json $MO/newmus_paths_shoulder_r_k6_40kg_h160_${st}_model.osim walk mf_shoulder_${w}_walk --prefilter &
  $X shoulder_r $R/newmus_paths_shoulder_r_k6_40kg_h160_$st.json $MO/newmus_paths_shoulder_r_k6_40kg_h160_${st}_load3.79_model.osim lift mf_shoulder_${w}_lift &
  wait
done
