#!/usr/bin/env bash
# Moco 검증 한 건: run_moco.sh <space> <json> <model.osim> <task: walk|squat|lift> <name> [추가 옵션...]
# 하체 walk = Moore T047_12 1~4 s (지면반력 비율 환산), squat = LabValidation squats1 2.6~4.4 s; 상체 walk = Moore T047_12 1~4 s 외력 없음, lift = CARRT subject1 PFWJ2 0.3~5.5 s
set -e
cd "$(dirname "$0")/../.."
PY=${PY:-/c/Users/USER/miniforge3/envs/osim/python.exe}; export PYTHONUTF8=1
# 윈도: pip 로 깐 opensim 은 OpenSim 4.6 GUI 의 DLL 을 쓴다. CasADi ipopt 플러그인은 의존 DLL(ipopt-3, libgfortran …)을 PATH 에서 찾는다
OSIM_BIN=${OSIM_BIN:-"/c/OpenSim 4.6/bin"}; export PATH="$OSIM_BIN:$PATH"
space=$1; json=$2; model=$3; task=$4; name=$5; shift 5
if [ "$space" = lower ]; then
  if [ "$task" = walk ]; then args="--ik data/robot40_h160/moore/T047_segment_12_ik.mot --ext out/lower40_h160/moore/T047_segment_12_external_loads_win.xml --id out/lower40_h160/moore/T047_segment_12_id.sto --t0 1 --t1 4"
  else args="--ik data/robot40_h160/lower/squats1_segment_0_ik.mot --ext out/lower40_h160/subject3/squats1_segment_0_external_loads_win.xml --id out/lower40_h160/subject3/squats1_segment_0_id.sto --t0 2.6 --t1 4.4"; fi
else
  if [ "$task" = walk ]; then args="--ik data/moore2015/subject3/T047_segment_12_ik.mot --id out/upper40_h160/moore2015/subject3/T047_segment_12/T047_segment_12_id.sto --ext none --t0 1 --t1 4"
  else args="--ik data/carrt/subject1/PFWJ2_ik.mot --id out/upper40_h160/carrt/subject1/PFWJ2/PFWJ2_id.sto --ext none --t0 0.3 --t1 5.5"; fi
fi
$PY -u report/moco_path_muscles.py --space $space --paths $json --model $model $args --box-scale 0.9877 --name $name "$@" > report/modulefirst/logs/moco_$name.log 2>&1
grep "\[$name\]" report/modulefirst/logs/moco_$name.log | tail -1
