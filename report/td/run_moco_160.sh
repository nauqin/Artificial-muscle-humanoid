#!/usr/bin/env bash
# 2026-11-19 160 cm 비교의 Moco 검증: 걷기 subject3 T047_12 · subject9 T027_21 (--prefilter), 찬 물병 들기 subject1 PFWJ2 · subject6 PFWJ_2.
#   기준(지금 설계): 팔꿈치 _s80, 어깨 _s80d (×1.5 없음), 어깨 _s80d_f15 (×1.5) — 같은 창으로 윈도에서 다시
#   시간 영역: 팔꿈치 _td, _td_br{0,20,30,40,50}, 어깨 _td, _td_br{20,40}
cd "$(dirname "$0")/../.."
PY=${PY:-/c/Users/USER/miniforge3/envs/osim/python.exe}; export PYTHONUTF8=1 MKL_THREADING_LAYER=SEQUENTIAL
X="bash report/modulefirst/run_moco.sh"; R=report; MO=report/moco; F="\[warning\]|\[info\]"
moco4() {   # space json(경로 이름 줄기) 이름
  local sp=$1 st=$2 nm=$3
  $PY report/modulefirst/prepare_models.py $sp $R/newmus_paths_$st.json 2>&1 | grep -Ev "$F"
  $X $sp $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walk ${nm}_walk --prefilter &
  $X $sp $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walks9 ${nm}_walks9 --prefilter &
  $X $sp $R/newmus_paths_$st.json $MO/newmus_paths_${st}_load3.79_model.osim lift ${nm}_lift &
  $X $sp $R/newmus_paths_$st.json $MO/newmus_paths_${st}_load3.79_model.osim lifts6 ${nm}_lifts6 &
}
case ${1:-all} in
  elbow|all) moco4 elbow_r elbow_r_k2_40kg_h160_s80 td160_elbow_s80
             for v in td td_br0 td_br20 td_br30 td_br40 td_br50; do [ -f $R/newmus_paths_elbow_r_k2_40kg_h160_$v${SFX:-}.json ] && moco4 elbow_r elbow_r_k2_40kg_h160_$v${SFX:-} td160_elbow_$v${SFX:-}; done; wait ;;&
  shoulder|all) moco4 shoulder_r shoulder_r_k8_40kg_h160_s80d td160_shoulder_s80d
             moco4 shoulder_r shoulder_r_k8_40kg_h160_s80d_f15 td160_shoulder_s80d_f15
             for v in td td_br20 td_br40; do [ -f $R/newmus_paths_shoulder_r_k8_40kg_h160_$v${SFX:-}.json ] && moco4 shoulder_r shoulder_r_k8_40kg_h160_$v${SFX:-} td160_shoulder_$v${SFX:-}; done; wait ;;
esac
