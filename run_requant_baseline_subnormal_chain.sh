#!/usr/bin/env bash
# Waits for the running orchestrator to exit, then launches the OCP-subnormal
# baseline FP8 sweep for 40k + 65k (32k already done as a standalone run).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WAIT_PID=${WAIT_PID:-3713997}
OUT_DIR=${OUT_DIR:-"/home/abduraa/MX_project/logs/mx_ooc_data"}
CHAIN_LOG="$OUT_DIR/requant_baseline_subnormal_chain.log"

{
  echo "[CHAIN START] $(date -Is)"
  echo "[CHAIN] waiting for PID $WAIT_PID (orchestrator) to exit"
} | tee -a "$CHAIN_LOG"

while [[ -d "/proc/$WAIT_PID" ]]; do
  sleep 60
done
echo "[CHAIN] PID $WAIT_PID exited at $(date -Is). Launching baseline subnormal sweep for 40k + 65k." | tee -a "$CHAIN_LOG"

# 32k is already complete; only run 40k + 65k. Skip-logic in build_run_relative
# isn't aware of "already done" at the sweep level, but the runner's bin loop
# explicitly lists only the bins we want.
BUILD_FIRST=0 OUT_DIR="$OUT_DIR" \
  BIN_LIST="/home/abduraa/MX_project/logs/my_cov_weak_40k.bin /home/abduraa/MX_project/logs/my_cov_weak_65k.bin" \
  "$SCRIPT_DIR/run_requant_baseline_fp8_subnormal_gt20k.sh" 2>&1 | tee -a "$CHAIN_LOG"

echo "[CHAIN END] $(date -Is)" | tee -a "$CHAIN_LOG"
