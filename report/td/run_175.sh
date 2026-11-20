#!/usr/bin/env bash
# 2026-11-19 175 cm · 40 kg, 시간 영역 설계.
#   모델: tools/rescale_body.py --mass 40 --height 1.75 --tag _h175sk --segment-scale (마디별 추가 배율) --per-source
#   배분비: 상체는 175 cm ID 로 design_directions 새로 (--starts 5), 하체는 160 cm V 그대로·c ×1.75/1.60 (선형, npz note 참고)
cd "$(dirname "$0")/../.."
PY=${PY:-/c/Users/USER/miniforge3/envs/osim/python.exe}; export PYTHONUTF8=1 MKL_THREADING_LAYER=SEQUENTIAL
M=data/robot40_h175sk/robot40_h175sk.osim; T=report/upper/tasks_40kg_h175sk.csv; B=1.0802; L=report/td/logs; mkdir -p $L
BS="--body-scale humerus_r=1,1.0974,1 --body-scale ulna_r=1,1.0495,1"
C="--model $M --box-scale $B $BS --force-mode time --refine-trust 0.01"
LOW="--space lower --design report/newmus_V_design_lower_k8_40kg_h175sk.npz --ref 40 1.75 $C --stage1-stride 4"
TRK="--space trunk --tasks $T --design report/newmus_V_design_trunk_k4_40kg_h175sk.npz $C --stage1-stride 4"
SHO="--space shoulder_r --tasks $T --design report/newmus_V_design_shoulder_r_k8_40kg_h175sk.npz $C --dt 0.04 --alloc-iters 150 --stage1-stride 6"
ELB="--space elbow_r --tasks $T --design report/newmus_V_design_elbow_r_k2_40kg_h175sk.npz $C --grades M3,M4,M5 --stage1-stride 2"
br() { awk "BEGIN{print $1/1000}"; }
$PY report/path_design.py $LOW --out-tag _40kg_h175sk_td > $L/pd175_lower_td.log 2>&1 &
for h in 20 40; do $PY report/path_design.py $LOW --bracket $(br $h) --out-tag _40kg_h175sk_td_br$h > $L/pd175_lower_td_br$h.log 2>&1 & done
$PY report/path_design.py $TRK --out-tag _40kg_h175sk_td > $L/pd175_trunk_td.log 2>&1 &
$PY report/path_design.py $SHO --out-tag _40kg_h175sk_td > $L/pd175_shoulder_td.log 2>&1 &
for h in 20 40; do $PY report/path_design.py $SHO --bracket $(br $h) --out-tag _40kg_h175sk_td_br$h > $L/pd175_shoulder_td_br$h.log 2>&1 & done
$PY report/path_design.py $ELB --out-tag _40kg_h175sk_td > $L/pd175_elbow_td.log 2>&1 &
for h in 20 30 40 50; do $PY report/path_design.py $ELB --bracket $(br $h) --out-tag _40kg_h175sk_td_br$h > $L/pd175_elbow_td_br$h.log 2>&1 & done
wait
