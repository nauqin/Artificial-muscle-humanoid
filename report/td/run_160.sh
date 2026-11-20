#!/usr/bin/env bash
# 2026-11-19 시간 영역 설계 — 160 cm 비교 (공식·브래킷 효과만 분리). 기존 160 cm 모델·ID·배분비를 그대로 쓴다.
#   팔꿈치: 허용 등급 M3·M4·M5 (교수님), 상자만 / 브래킷 0·20·30·40·50 mm
#   어깨 8: 모든 등급, ×1.5 없음, 상자만 / 브래킷 20·40 mm (위팔 쪽 점만 — 몸통은 상자)
cd "$(dirname "$0")/../.."
PY=${PY:-/c/Users/USER/miniforge3/envs/osim/python.exe}; export PYTHONUTF8=1 MKL_THREADING_LAYER=SEQUENTIAL
M=data/robot40_h160/robot40_h160.osim; T=report/upper/tasks_40kg_h160.csv; B=0.9877; L=report/td/logs; mkdir -p $L
E="--space elbow_r --tasks $T --design report/newmus_V_design_elbow_r_k2_40kg_h160.npz --model $M --box-scale $B --force-mode time --grades M3,M4,M5 --refine-trust 0.01 --stage1-stride 2"
S="--space shoulder_r --tasks $T --design report/newmus_V_design_shoulder_r_k8_40kg_h160.npz --model $M --box-scale $B --force-mode time --refine-trust 0.01 --dt 0.04 --alloc-iters 150 --stage1-stride 6"
$PY report/path_design.py $E --out-tag _40kg_h160_td > $L/pd_elbow_td.log 2>&1 &
for h in 0 20 30 40 50; do
  $PY report/path_design.py $E --bracket $(awk "BEGIN{print $h/1000}") --out-tag _40kg_h160_td_br$h > $L/pd_elbow_td_br$h.log 2>&1 &
done
$PY report/path_design.py $S --out-tag _40kg_h160_td > $L/pd_shoulder_td.log 2>&1 &
for h in 20 40; do
  $PY report/path_design.py $S --bracket $(awk "BEGIN{print $h/1000}") --out-tag _40kg_h160_td_br$h > $L/pd_shoulder_td_br$h.log 2>&1 &
done
wait
