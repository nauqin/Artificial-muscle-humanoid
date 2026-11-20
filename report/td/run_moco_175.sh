#!/usr/bin/env bash
# 2026-11-19 175 cm 시간 영역 설계의 Moco 검증.
#   상체: 걷기 subject3 T047_12 · subject9 T027_21 (--prefilter), 찬 물병 들기 subject1 PFWJ2 · subject6 PFWJ_2 (175 cm 뼈대로 다시 낸 ID)
#   하체: Moore 걷기 T047_12 (--prefilter), LabValidation 걷기, 스쿼트 (175 cm 뼈대 + 출처별 비율로 줄인 지면반력 ID)
cd "$(dirname "$0")/../.."
PY=${PY:-/c/Users/USER/miniforge3/envs/osim/python.exe}; export PYTHONUTF8=1 MKL_THREADING_LAYER=SEQUENTIAL
export B=1.0802 UP=out/upper40_h175sk LO=out/lower40_h175sk DR=data/robot40_h175sk
X="bash report/modulefirst/run_moco.sh"; R=report; MO=report/moco; F="\[warning\]|\[info\]"
up4() {   # space 경로이름줄기 이름
  local sp=$1 st=$2 nm=$3
  $PY report/modulefirst/prepare_models.py $sp $R/newmus_paths_$st.json --box-scale $B --base $DR/robot40_h175sk.osim 2>&1 | grep -Ev "$F"
  $X $sp $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walk ${nm}_walk --prefilter &
  $X $sp $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walks9 ${nm}_walks9 --prefilter &
  $X $sp $R/newmus_paths_$st.json $MO/newmus_paths_${st}_load3.79_model.osim lift ${nm}_lift &
  $X $sp $R/newmus_paths_$st.json $MO/newmus_paths_${st}_load3.79_model.osim lifts6 ${nm}_lifts6 &
}
lo3() {
  local st=$1 nm=$2
  $X lower $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walk ${nm}_walk --prefilter &
  $X lower $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim walklab ${nm}_walklab &
  $X lower $R/newmus_paths_$st.json $MO/newmus_paths_${st}_model.osim squat ${nm}_squat &
}
ok() { [ -f $R/newmus_paths_$1.json ]; }
case ${1:-all} in
  lower|all) for v in td td_br20 td_br40; do ok lower_k8_40kg_h175sk_$v${SFX:-} && lo3 lower_k8_40kg_h175sk_$v${SFX:-} td175_lower_$v${SFX:-}; done; wait ;;&
  trunk|all) ok trunk_k4_40kg_h175sk_td${SFX:-} && up4 trunk trunk_k4_40kg_h175sk_td${SFX:-} td175_trunk_td${SFX:-}; wait ;;&
  elbow|all) for v in td td_br20 td_br30 td_br40 td_br50; do ok elbow_r_k2_40kg_h175sk_$v${SFX:-} && up4 elbow_r elbow_r_k2_40kg_h175sk_$v${SFX:-} td175_elbow_$v${SFX:-}; done; wait ;;&
  shoulder|all) for v in td td_br20 td_br40; do ok shoulder_r_k8_40kg_h175sk_$v${SFX:-} && up4 shoulder_r shoulder_r_k8_40kg_h175sk_$v${SFX:-} td175_shoulder_$v${SFX:-}; done; wait ;;
esac
