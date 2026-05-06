#!/usr/bin/env bash
# Sweep SFP8-E4M3 (mx_e4m3) in subtile_128 mode for weak-corr 20k matrix.
# One run per epsilon: 1e-5, 1e-6, 1e-7, 1e-8.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_SCRIPT="$SCRIPT_DIR/build_run.sh"

BIN_PATH="/home/abduraa/MX_project/logs/my_cov_weak_20k.bin"
EPS_LIST="1e-5 1e-6 1e-7 1e-8"
NB=2048
CORES=16

OUT_DIR="/home/abduraa/MX_project/logs/mx_ooc_data_sfp8_subtile_eps"
mkdir -p "$OUT_DIR"
SUMMARY_FILE="$OUT_DIR/summary_sfp8_subtile_weak20k.tsv"
MASTER_LOG="$OUT_DIR/sweep_master.log"

{
  echo "[START] $(date -Is)"
  echo "[CFG] format=e4m3  mode=block  subtile=128"
  echo "[CFG] bin=$BIN_PATH  eps=$EPS_LIST  nb=$NB  cores=$CORES"
} | tee "$MASTER_LOG"

for eps in $EPS_LIST; do
  stamp=$(date +"%Y%m%d_%H%M%S")
  run_log="$OUT_DIR/run_sfp8_subtile128_weak20k_eps_${eps//-/m}_${stamp}.log"

  echo "[RUN] eps=$eps" | tee -a "$MASTER_LOG"

  env \
    FORMAT="e4m3" \
    MX_MODE="block" \
    MX_MX_MODE="block" \
    MX_BLOCK_SUBTILE="128" \
    MX_SOURCE_EPSILON="$eps" \
    MX_SCALE_AWARE_EPSILON="1" \
    MX_SKIP_KL="1" \
    MX_UNDERFLOW_MODE="fz" \
    BUILD="0" \
    NB="$NB" \
    CORES="$CORES" \
    BIN_PATH="$BIN_PATH" \
    SUMMARY_FILE="$SUMMARY_FILE" \
    OUT_DIR="$OUT_DIR" \
    "$RUN_SCRIPT" > "$run_log" 2>&1

  echo "[DONE] eps=$eps  log=$run_log" | tee -a "$MASTER_LOG"
done

echo "[END] $(date -Is)" | tee -a "$MASTER_LOG"
echo "Summary: $SUMMARY_FILE"
