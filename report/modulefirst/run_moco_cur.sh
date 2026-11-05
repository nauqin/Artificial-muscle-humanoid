#!/usr/bin/env bash
# 현재 방식(이 컴퓨터 재실행) Moco 8건. 걷기는 --prefilter (맥과 같은 조건, 하체 걷기로 확인)
cd "$(dirname "$0")/../.."; R=report; MO=report/moco; X="bash report/modulefirst/run_moco.sh"
$X lower $R/newmus_paths_lower_k8_40kg_h160_s120_ref.json $MO/newmus_paths_lower_k8_40kg_h160_s120_ref_model.osim walk mf_lower_cur_walk --prefilter &
$X lower $R/newmus_paths_lower_k8_40kg_h160_s120_ref.json $MO/newmus_paths_lower_k8_40kg_h160_s120_ref_model.osim squat mf_lower_cur_squat &
$X trunk $R/newmus_paths_trunk_k4_40kg_h160_s80_win.json $MO/newmus_paths_trunk_k4_40kg_h160_s80_win_model.osim walk mf_trunk_cur_walk --prefilter &
$X trunk $R/newmus_paths_trunk_k4_40kg_h160_s80_win.json $MO/newmus_paths_trunk_k4_40kg_h160_s80_win_load3.79_model.osim lift mf_trunk_cur_lift &
$X elbow_r $R/newmus_paths_elbow_r_k2_40kg_h160_s80_win.json $MO/newmus_paths_elbow_r_k2_40kg_h160_s80_win_model.osim walk mf_elbow_cur_walk --prefilter &
$X elbow_r $R/newmus_paths_elbow_r_k2_40kg_h160_s80_win.json $MO/newmus_paths_elbow_r_k2_40kg_h160_s80_win_load3.79_model.osim lift mf_elbow_cur_lift &
$X shoulder_r $R/newmus_paths_shoulder_r_k6_40kg_h160_v635_s80_win_f15.json $MO/newmus_paths_shoulder_r_k6_40kg_h160_v635_s80_win_f15_model.osim walk mf_shoulder_cur_walk --prefilter &
$X shoulder_r $R/newmus_paths_shoulder_r_k6_40kg_h160_v635_s80_win_f15.json $MO/newmus_paths_shoulder_r_k6_40kg_h160_v635_s80_win_f15_load3.79_model.osim lift mf_shoulder_cur_lift &
wait
