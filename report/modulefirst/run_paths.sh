#!/usr/bin/env bash
# 2026-11-05 모듈 먼저 실험: 현재 방식(부위 전체 상한) 경로 설계 4건을 이 컴퓨터에서 다시. 새 방식은 run_paths_final.sh. 로그는 report/modulefirst/logs/path_<이름>.log
cd "$(dirname "$0")/../.."
PY=${PY:-/c/Users/USER/miniforge3/envs/osim/python.exe}; export PYTHONUTF8=1 MKL_THREADING_LAYER=SEQUENTIAL   # 윈도 MKL 다중 스레드(OpenMP)가 행렬곱에서 죽는다 → 단일 스레드
M=data/robot40_h160/robot40_h160.osim; T=report/upper/tasks_40kg_h160.csv; B=0.9877; MF=report/modulefirst
run() { n=$1; shift; $PY -u report/path_design.py "$@" > $MF/logs/path_$n.log 2>&1; echo "[$n] $(grep '2단계 끝' $MF/logs/path_$n.log)"; }
# 현재 방식 (부위 전체 상한, TUTORIAL 9-6 과 같은 조건). 하체는 --ref 가 이제 토크에 반영된다
run lower_cur    --space lower --design report/newmus_V_design_lower_k8_40kg_h160.npz --model $M --box-scale $B --ref 40 1.60 --out-tag _40kg_h160_s120_ref --stroke-max 0.12 &
run trunk_cur    --space trunk --tasks $T --design report/newmus_V_design_trunk_k4_40kg_h160.npz --model $M --box-scale $B --out-tag _40kg_h160_s80_win --f-module 1800 --stroke-max 0.08 &
run elbow_cur    --space elbow_r --tasks $T --design report/newmus_V_design_elbow_r_k2_40kg_h160.npz --model $M --box-scale $B --out-tag _40kg_h160_s80_win --stroke-max 0.08 &
run shoulder_cur --space shoulder_r --tasks $T --design report/newmus_V_design_shoulder_r_k6_40kg_h160_v635.npz --model $M --box-scale $B --out-tag _40kg_h160_v635_s80_win --stroke-max 0.08 &
wait
