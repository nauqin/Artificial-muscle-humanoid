#!/usr/bin/env bash
# 2026-11-19: 아래팔 회내외 전용 케이블 2개 (노뼈 드럼 = OpenSim WrapCylinder) + 팔꿈치 최종 세트 td_br50_rb → Moco (README 13-10).
#   1) 설계만 (변수 10개, 설계 쪽 numpy) + OpenSim 길이·모멘트암 비교 + 모델 → report/td/prosup_wrap_prosup_wrap/
#   2) 바깥 반복: Moco (lift s1 + walk s3) → 회내외 보조가 남은 순간을 설계 프레임으로 → 재설계 (통과하면 멈춤, 최대 3회)
#   3) 최종 Moco walk·walks9·lift·lifts6 → results_paths.csv 의 td175_forearm_prosup_wrap_* 줄
#   강건성 X = 0.25 N·m (사용자 결정, 2026-11-19), 허용 등급 M3·M4·M5, 자뼈 가이드 브래킷 ≤ 20 mm
cd "$(dirname "$0")/../.."
PY=${PY:-/c/Users/USER/miniforge3/envs/osim/python.exe}; export PYTHONUTF8=1 MKL_THREADING_LAYER=SEQUENTIAL
TAG=${TAG:-_prosup_wrap}
$PY -u report/td/prosup_wrap.py --robust-nm 0.25 --out-tag $TAG --outer 3 --final "$@" > report/td/logs/prosup_wrap$TAG.log 2>&1
grep -E "→ 반복|→ 최종|OpenSim M|저장" report/td/logs/prosup_wrap$TAG.log
