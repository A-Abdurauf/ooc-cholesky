#!/usr/bin/env bash
# Sequel chain: after the 32k SFP4 main orch finishes, runs the tile orch
# for 32k only (picks up SFP4 tile variants), then resumes the 65k pipeline.
# This time we DON'T do broad pkill — we wait on PIDs only.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR=${OUT_DIR:-"/home/abduraa/MX_project/logs/mx_ooc_data"}
LOG="$OUT_DIR/requant_sfp4_after.log"
WAIT_PID=${WAIT_PID:-3800872}

{
  echo "[SFP4-AFTER START] $(date -Is)"
  echo "[INFO] waiting for main_orch(32k) PID $WAIT_PID to exit"
} | tee -a "$LOG"

while [[ -d "/proc/$WAIT_PID" ]]; do
  sleep 30
done
echo "[INFO] PID $WAIT_PID exited at $(date -Is)" | tee -a "$LOG"

# Tile orch limited to 32k bin.
echo "[INFO] launching tile orch (32k only) for SFP4 tile sweeps" | tee -a "$LOG"
BIN_LIST="/home/abduraa/MX_project/logs/my_cov_weak_32k.bin" OUT_DIR="$OUT_DIR" \
  "$SCRIPT_DIR/run_requant_gt20k_tile_by_bin.sh" 2>&1 | tee -a "$LOG"

# Resume the 65k pipeline.
echo "[INFO] launching full main_orch to resume 65k" | tee -a "$LOG"
nohup "$SCRIPT_DIR/run_requant_gt20k_by_bin.sh" \
  >> "$OUT_DIR/requant_gt20k_by_bin.nohup.log" 2>&1 &
NEW_ORCH=$!
disown $NEW_ORCH 2>/dev/null
echo "[INFO] resumed main_orch PID=$NEW_ORCH" | tee -a "$LOG"

sed -i "s/^WAIT_PID=.*/WAIT_PID=\${WAIT_PID:-$NEW_ORCH}/" "$SCRIPT_DIR/run_requant_baseline_subnormal_chain.sh"
nohup "$SCRIPT_DIR/run_requant_baseline_subnormal_chain.sh" \
  >> "$OUT_DIR/requant_baseline_subnormal_chain.nohup.log" 2>&1 &
NEW_BCHAIN=$!
disown $NEW_BCHAIN 2>/dev/null
echo "[INFO] baseline_subnormal chain PID=$NEW_BCHAIN" | tee -a "$LOG"

sed -i "s/^WAIT_PID=.*/WAIT_PID=\${WAIT_PID:-$NEW_BCHAIN}/" "$SCRIPT_DIR/run_requant_gt20k_tile_chain.sh"
nohup "$SCRIPT_DIR/run_requant_gt20k_tile_chain.sh" \
  >> "$OUT_DIR/requant_gt20k_tile_chain.nohup.log" 2>&1 &
NEW_TCHAIN=$!
disown $NEW_TCHAIN 2>/dev/null
echo "[INFO] tile chain PID=$NEW_TCHAIN" | tee -a "$LOG"

echo "[SFP4-AFTER END] $(date -Is)" | tee -a "$LOG"
