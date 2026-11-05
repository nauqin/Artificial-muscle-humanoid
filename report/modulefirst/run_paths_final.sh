#!/usr/bin/env bash
# 2026-11-05 최종: 새 방식 4건(--modules, --refine-trust 0.01) + 현재 방식에 같은 신뢰 구간을 켠 민감도 확인 4건(_tr)
cd "$(dirname "$0")/../.."
PY=${PY:-/c/Users/USER/miniforge3/envs/osim/python.exe}; export PYTHONUTF8=1 MKL_THREADING_LAYER=SEQUENTIAL
M=data/robot40_h160/robot40_h160.osim; T=report/upper/tasks_40kg_h160.csv; B=0.9877; MF=report/modulefirst; TR="--refine-trust 0.01"
run() { n=$1; shift; $PY -u report/path_design.py "$@" > $MF/logs/path_$n.log 2>&1; echo "[$n] $(grep '2단계 끝' $MF/logs/path_$n.log)"; grep "수선 M" $MF/logs/path_$n.log; }
L="--space lower --design report/newmus_V_design_lower_k8_40kg_h160.npz --model $M --box-scale $B --ref 40 1.60"
TK="--space trunk --tasks $T --design report/newmus_V_design_trunk_k4_40kg_h160.npz --model $M --box-scale $B"
E="--space elbow_r --tasks $T --design report/newmus_V_design_elbow_r_k2_40kg_h160.npz --model $M --box-scale $B"
S="--space shoulder_r --tasks $T --design report/newmus_V_design_shoulder_r_k6_40kg_h160_v635.npz --model $M --box-scale $B"
run lower_mf    $L  --out-tag _40kg_h160_mf      --modules $MF/modules_lower_k8_40kg_h160.json $TR &
run trunk_mf    $TK --out-tag _40kg_h160_mf      --modules $MF/modules_trunk_k4_40kg_h160.json $TR &
run elbow_mf    $E  --out-tag _40kg_h160_mf      --modules $MF/modules_elbow_r_k2_40kg_h160.json $TR &
run shoulder_mf $S  --out-tag _40kg_h160_v635_mf --modules $MF/modules_shoulder_r_k6_40kg_h160_v635_f15.json $TR &
run lower_cur_tr    $L  --out-tag _40kg_h160_s120_ref_tr --stroke-max 0.12 $TR &
run trunk_cur_tr    $TK --out-tag _40kg_h160_s80_win_tr --f-module 1800 --stroke-max 0.08 $TR &
run elbow_cur_tr    $E  --out-tag _40kg_h160_s80_win_tr --stroke-max 0.08 $TR &
run shoulder_cur_tr $S  --out-tag _40kg_h160_v635_s80_win_tr --stroke-max 0.08 $TR &
wait
