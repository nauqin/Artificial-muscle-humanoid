#!/usr/bin/env bash
# 2026-11-19 밤 (3판): 배정 불가 벌점을 바닥값으로 고친 뒤. 팔꿈치는 처음부터 다시(1·2단계, 방향 가중치 0, 꼬리표 _w0),
# 나머지는 run_160.sh·run_175.sh 의 설계(2단계 방향 가중치 3)가 끝나면, 그 json 에서 2단계만 방향 가중치 0 으로 다시 (꼬리표 _w0)
# → 부위가 다 끝나면 그 부위 Moco (SFX=_w0). 최종 보고는 _w0 (README 13절).
cd "$(dirname "$0")/../.."
PY=${PY:-/c/Users/USER/miniforge3/envs/osim/python.exe}; export PYTHONUTF8=1 MKL_THREADING_LAYER=SEQUENTIAL; L=report/td/logs
fin() { for f in "$@"; do grep -qE "그림 저장|Traceback" $f 2>/dev/null || return 1; done; }
waitfor() { until fin "$@"; do sleep 30; done; }
M6=data/robot40_h160/robot40_h160.osim; T6=report/upper/tasks_40kg_h160.csv
M7=data/robot40_h175sk/robot40_h175sk.osim; T7=report/upper/tasks_40kg_h175sk.csv; BS="--body-scale humerus_r=1,1.0974,1 --body-scale ulna_r=1,1.0495,1"
E6="--space elbow_r --tasks $T6 --design report/newmus_V_design_elbow_r_k2_40kg_h160.npz --model $M6 --box-scale 0.9877 --force-mode time --grades M3,M4,M5 --refine-trust 0.01"
S6="--space shoulder_r --tasks $T6 --design report/newmus_V_design_shoulder_r_k8_40kg_h160.npz --model $M6 --box-scale 0.9877 --force-mode time --refine-trust 0.01 --dt 0.04 --alloc-iters 150"
C7="--model $M7 --box-scale 1.0802 $BS --force-mode time --refine-trust 0.01"
LOW7="--space lower --design report/newmus_V_design_lower_k8_40kg_h175sk.npz --ref 40 1.75 $C7"
TRK7="--space trunk --tasks $T7 --design report/newmus_V_design_trunk_k4_40kg_h175sk.npz $C7"
SHO7="--space shoulder_r --tasks $T7 --design report/newmus_V_design_shoulder_r_k8_40kg_h175sk.npz $C7 --dt 0.04 --alloc-iters 150"
ELB7="--space elbow_r --tasks $T7 --design report/newmus_V_design_elbow_r_k2_40kg_h175sk.npz $C7 --grades M3,M4,M5"
brf() { case $1 in td) echo "";; td_br*) echo "--bracket $(awk "BEGIN{print ${1#td_br}/1000}")";; esac; }
re() {   # 공통플래그 json줄기 꼬리표(_40kg_…_v) 변형 로그이름
  local flags=$1 stem=$2 tag=$3 v=$4 log=$5
  $PY report/path_design.py $flags $(brf $v) --refine-from report/newmus_paths_$stem.json --out-tag ${tag}_w0 > $L/${log}_w0.log 2>&1
}
full() {   # 공통플래그 꼬리표 변형 로그이름 — 처음부터 (1단계 포함)
  local flags=$1 tag=$2 v=$3 log=$4
  $PY report/path_design.py $flags $(brf $v) --stage1-stride 2 --out-tag ${tag}_w0 > $L/${log}_w0.log 2>&1
}
egroup() {  # 로그접두 공통플래그 꼬리표접두 moco스크립트 변형…
  local lp=$1 flags=$2 tp=$3 ms=$4; shift 4
  for v in "$@"; do full "$flags" ${tp}_$v $v ${lp}_$v & done; wait
  SFX=_w0 bash report/td/$ms elbow
}
group() {   # 로그접두 공통플래그 줄기접두 꼬리표접두 moco스크립트 부위 변형…
  local lp=$1 flags=$2 sp=$3 tp=$4 ms=$5 reg=$6; shift 6
  waitfor $(for v in "$@"; do echo $L/${lp}_$v.log; done)
  for v in "$@"; do grep -q "그림 저장" $L/${lp}_$v.log && re "$flags" ${sp}_$v ${tp}_$v $v ${lp}_$v & done; wait
  SFX=_w0 bash report/td/$ms $reg
}
egroup pd_elbow "$E6" _40kg_h160 run_moco_160.sh td td_br0 td_br20 td_br30 td_br40 td_br50 &
group pd_shoulder "$S6" shoulder_r_k8_40kg_h160 _40kg_h160 run_moco_160.sh shoulder td td_br20 td_br40 &
egroup pd175_elbow "$ELB7" _40kg_h175sk run_moco_175.sh td td_br20 td_br30 td_br40 td_br50 &
group pd175_trunk "$TRK7" trunk_k4_40kg_h175sk _40kg_h175sk run_moco_175.sh trunk td &
group pd175_lower "$LOW7" lower_k8_40kg_h175sk _40kg_h175sk run_moco_175.sh lower td td_br20 td_br40 &
group pd175_shoulder "$SHO7" shoulder_r_k8_40kg_h175sk _40kg_h175sk run_moco_175.sh shoulder td td_br20 td_br40 &
wait
echo "끝"
