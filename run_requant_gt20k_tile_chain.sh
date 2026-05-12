#!/usr/bin/env bash
# Waits for the existing baseline_subnormal chain (which itself waits for the
# main orchestrator) to exit, then launches the tile-mode orchestrator.
# Effective chain order:
#   run_requant_gt20k_by_bin.sh                  (subnormal MX, 96 runs)
#     -> run_requant_baseline_subnormal_chain.sh (baseline_subn 40k+65k, 7 runs)
#     -> run_requant_gt20k_tile_by_bin.sh        (tile-mode MX, 48 runs)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WAIT_PID=${WAIT_PID:-3714005}
OUT_DIR=${OUT_DIR:-"/home/abduraa/MX_project/logs/mx_ooc_data"}
CHAIN_LOG="$OUT_DIR/requant_gt20k_tile_chain.log"

{
  echo "[TILE CHAIN START] $(date -Is)"
  echo "[TILE CHAIN] waiting for PID $WAIT_PID (baseline_subnormal chain) to exit"
} | tee -a "$CHAIN_LOG"

while [[ -d "/proc/$WAIT_PID" ]]; do
  sleep 60
done
echo "[TILE CHAIN] PID $WAIT_PID exited at $(date -Is). Launching tile orchestrator." | tee -a "$CHAIN_LOG"

OUT_DIR="$OUT_DIR" \
  "$SCRIPT_DIR/run_requant_gt20k_tile_by_bin.sh" 2>&1 | tee -a "$CHAIN_LOG"

echo "[TILE CHAIN END] $(date -Is)" | tee -a "$CHAIN_LOG"
