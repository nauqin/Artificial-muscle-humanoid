#!/bin/bash
# 어느 폴더에서 실행하든 `osim` 환경 안에서 스크립트 하나를 돌린다.
#   ./run.sh check_env.py
#   ./run.sh batch.py --data data --out out --so --synergy
set -e
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate osim
cd "$HERE"
exec python "$@"
