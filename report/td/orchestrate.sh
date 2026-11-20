#!/usr/bin/env bash
# 2026-11-19 밤: 설계(run_160.sh·run_175.sh 가 백그라운드로 돌림)가 부위별로 끝나는 대로 그 부위 Moco 를 돌린다.
cd "$(dirname "$0")/../.."; L=report/td/logs
fin() { for f in "$@"; do grep -qE "그림 저장|Traceback" $f 2>/dev/null || return 1; done; }
waitfor() { until fin "$@"; do sleep 30; done; }
( waitfor $L/pd_elbow_td*.log; bash report/td/run_moco_160.sh elbow ) &
( waitfor $L/pd175_elbow_td*.log; bash report/td/run_moco_175.sh elbow ) &
( waitfor $L/pd175_trunk_td.log; bash report/td/run_moco_175.sh trunk ) &
( waitfor $L/pd175_lower_td*.log; bash report/td/run_moco_175.sh lower ) &
( waitfor $L/pd_shoulder_td*.log; bash report/td/run_moco_160.sh shoulder ) &
( waitfor $L/pd175_shoulder_td*.log; bash report/td/run_moco_175.sh shoulder ) &
wait; echo "모든 설계·Moco 끝"
