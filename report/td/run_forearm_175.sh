#!/usr/bin/env bash
# 2026-11-19 낮: 팔꿈치 + 아래팔 회내외 (forearm_r, 케이블 4개) 175 cm 시간 영역 설계 → Moco (걷기 s3·s9, 들기 s1·s6).
#   배분비: python report/design_directions.py --space forearm_r --tasks report/upper/tasks_40kg_h175sk.csv --k 4 --starts 5 --ref-mass 40 --ref-height 1.75 --out-tag _40kg_h175sk
#   허용 등급 M3·M4·M5 (교수님: 팔꿈치는 M3·4·5 위주), 강건성 0.25, 상자만 / 브래킷 ≤ 20 mm
cd "$(dirname "$0")/../.."
PY=${PY:-/c/Users/USER/miniforge3/envs/osim/python.exe}; export PYTHONUTF8=1 MKL_THREADING_LAYER=SEQUENTIAL
M=data/robot40_h175sk/robot40_h175sk.osim; T=report/upper/tasks_40kg_h175sk.csv; L=report/td/logs; F="\[warning\]|\[info\]"
BS="--body-scale humerus_r=1,1.0974,1 --body-scale ulna_r=1,1.0495,1 --body-scale radius_r=1,1.0495,1"
FA="--space forearm_r --tasks $T --design report/newmus_V_design_forearm_r_k4_40kg_h175sk.npz --model $M --box-scale 1.0802 $BS --force-mode time --robust-frac 0.25 --grades M3,M4,M5 --refine-trust 0.01 --stage1-stride 2"
$PY report/path_design.py $FA --out-tag _40kg_h175sk_td_rb > $L/pd175_forearm_td_rb.log 2>&1 &
$PY report/path_design.py $FA --bracket 0.02 --out-tag _40kg_h175sk_td_br20_rb > $L/pd175_forearm_td_br20_rb.log 2>&1 &
wait
export B=1.0802 UP=out/upper40_h175sk LO=out/lower40_h175sk DR=data/robot40_h175sk; X="bash report/modulefirst/run_moco.sh"; R=report; MO=report/moco
for v in td_rb td_br20_rb; do
  st=forearm_r_k4_40kg_h175sk_$v; nm=td175_forearm_$v
  [ -f $R/newmus_paths_$st.json ] || continue
  $PY report/modulefirst/prepare_models.py forearm_r $R/newmus_paths_$st.json --box-scale $B --base $DR/robot40_h175sk.osim 2>&1 | grep -Ev "$F"
  $X forearm_r $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walk ${nm}_walk --prefilter &
  $X forearm_r $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walks9 ${nm}_walks9 --prefilter &
  $X forearm_r $R/newmus_paths_$st.json $MO/newmus_paths_${st}_load3.79_model.osim lift ${nm}_lift &
  $X forearm_r $R/newmus_paths_$st.json $MO/newmus_paths_${st}_load3.79_model.osim lifts6 ${nm}_lifts6 &
done
wait; echo "아래팔 끝"
