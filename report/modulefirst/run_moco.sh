#!/usr/bin/env bash
# Moco 검증 한 건: run_moco.sh <space> <json> <model.osim> <task> <name> [추가 옵션...]
# 하체 walk = Moore T047_12 1~4 s (지면반력 비율 환산), squat = LabValidation squats1 2.6~4.4 s; 상체 walk = Moore T047_12 1~4 s 외력 없음, lift = CARRT subject1 PFWJ2 0.3~5.5 s
set -e
cd "$(dirname "$0")/../.."
PY=${PY:-/c/Users/USER/miniforge3/envs/osim/python.exe}; export PYTHONUTF8=1
# 윈도: pip 로 깐 opensim 은 OpenSim 4.6 GUI 의 DLL 을 쓴다. CasADi ipopt 플러그인은 의존 DLL(ipopt-3, libgfortran …)을 PATH 에서 찾는다
OSIM_BIN=${OSIM_BIN:-"/c/OpenSim 4.6/bin"}; export PATH="$OSIM_BIN:$PATH"
# 2026-11-19: 몸 크기를 환경 변수로 (기본 = 지금까지의 160 cm). 175 cm: B=1.0802 UP=out/upper40_h175sk LO=out/lower40_h175sk DR=data/robot40_h175sk
B=${B:-0.9877}; UP=${UP:-out/upper40_h160}; LO=${LO:-out/lower40_h160}; DR=${DR:-data/robot40_h160}
space=$1; json=$2; model=$3; task=$4; name=$5; shift 5
#  추가 (2026-11-05 2차): 하체 walklab = LabValidation walking1 (ID 창 전체), 상체 walk2 = Moore T048_5 (1.6 m/s) 1~4 s, lift2 = CARRT subject1 PB1 (상자 0.52 kg 양손) 0.3~4.8 s
#  추가 (2026-11-19): 상체 walks9 = Moore subject9 T027_21 (1.6 m/s) 1~4 s, lifts6 = CARRT subject6 PFWJ_2 (찬 물병 3.79 kg) 0.3~6.2 s  ※ 두 시행 모두 설계 과제 목록(tasks*.csv)에 들어 있다 — 따로 남긴 검증 시행이 아님 (2026-11-19 밤 정정)
if [ "$space" = lower ]; then
  if [ "$task" = walk ]; then args="--ik $DR/moore/T047_segment_12_ik.mot --ext $LO/moore/T047_segment_12_external_loads_win.xml --id $LO/moore/T047_segment_12_id.sto --t0 1 --t1 4"
  elif [ "$task" = walklab ]; then args="--ik $DR/lower/walking1_segment_0_ik.mot --ext $LO/subject3/walking1_segment_0_external_loads_win.xml --id $LO/subject3/walking1_segment_0_id.sto"
  else args="--ik $DR/lower/squats1_segment_0_ik.mot --ext $LO/subject3/squats1_segment_0_external_loads_win.xml --id $LO/subject3/squats1_segment_0_id.sto --t0 2.6 --t1 4.4"; fi
else
  if [ "$task" = walk ]; then args="--ik data/moore2015/subject3/T047_segment_12_ik.mot --id $UP/moore2015/subject3/T047_segment_12/T047_segment_12_id.sto --ext none --t0 1 --t1 4"
  elif [ "$task" = walk2 ]; then args="--ik data/moore2015/subject3/T048_segment_5_ik.mot --id $UP/moore2015/subject3/T048_segment_5/T048_segment_5_id.sto --ext none --t0 1 --t1 4"
  elif [ "$task" = walks9 ]; then args="--ik data/moore2015/subject9/T027_segment_21_ik.mot --id $UP/moore2015/subject9/T027_segment_21/T027_segment_21_id.sto --ext none --t0 1 --t1 4"
  elif [ "$task" = lift2 ]; then args="--ik data/carrt/subject1/PB1_ik.mot --id $UP/carrt/subject1/PB1/PB1_id.sto --ext none --t0 0.3 --t1 4.8"
  elif [ "$task" = lifts6 ]; then args="--ik data/carrt/subject6/PFWJ_2_ik.mot --id $UP/carrt/subject6/PFWJ_2/PFWJ_2_id.sto --ext none --t0 0.3 --t1 6.2"
  else args="--ik data/carrt/subject1/PFWJ2_ik.mot --id $UP/carrt/subject1/PFWJ2/PFWJ2_id.sto --ext none --t0 0.3 --t1 5.5"; fi
fi
LOGDIR=${LOGDIR:-report/modulefirst/logs}   # 2026-11-19: 최적화 반복 평가는 다른 곳에 (기본은 그대로)
$PY -u report/moco_path_muscles.py --space $space --paths $json --model $model $args --box-scale $B --name $name "$@" > $LOGDIR/moco_$name.log 2>&1
grep "\[$name\]" $LOGDIR/moco_$name.log | tail -1
