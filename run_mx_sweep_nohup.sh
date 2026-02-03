#!/bin/bash
set -euo pipefail

LOG_DIR="/home/abduraa/MX_project/logs/mx_ooc_data/sweep"
mkdir -p "$LOG_DIR"

STAMP=$(date +"%Y%m%d_%H%M%S")
OUT_LOG="$LOG_DIR/sweep_nohup_${STAMP}.log"

nohup env PLAIN_FP8_ONLY=1 QUIET=1 BUILD=${BUILD:-0} ./run_mx_sweep.sh > "$OUT_LOG" 2>&1 &

echo "Started sweep in background. Log: $OUT_LOG"