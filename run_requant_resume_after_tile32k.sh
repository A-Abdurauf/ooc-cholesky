#!/usr/bin/env bash
# After the tile-32k orchestrator (PID supplied via WAIT_PID) finishes, resume
# the main orchestrator + baseline_subnormal chain + tile chain for 40k/65k.
# Skip-logic in each sub-script preserves already-done entries.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WAIT_PID=${WAIT_PID:-3445370}
OUT_DIR=${OUT_DIR:-"/home/abduraa/MX_project/logs/mx_ooc_data"}
RESUME_LOG="$OUT_DIR/requant_resume_after_tile32k.log"

{
  echo "[RESUME START] $(date -Is)"
  echo "[RESUME] waiting for PID $WAIT_PID (tile-32k orchestrator) to exit"
} | tee -a "$RESUME_LOG"

while [[ -d "/proc/$WAIT_PID" ]]; do
  sleep 60
done
echo "[RESUME] PID $WAIT_PID exited at $(date -Is). Relaunching main orchestrator + chains." | tee -a "$RESUME_LOG"

# 1. Main orchestrator (will skip already-done 32k entries, do 40k + 65k for 8 MX sweeps).
nohup "$SCRIPT_DIR/run_requant_gt20k_by_bin.sh" \
  >> "$OUT_DIR/requant_gt20k_by_bin.nohup.log" 2>&1 &
ORCH_PID=$!
disown $ORCH_PID 2>/dev/null
echo "[RESUME] orchestrator PID=$ORCH_PID" | tee -a "$RESUME_LOG"

# 2. Baseline subnormal chain (40k + 65k) — waits for orch.
sed -i "s/^WAIT_PID=.*/WAIT_PID=\${WAIT_PID:-$ORCH_PID}/" "$SCRIPT_DIR/run_requant_baseline_subnormal_chain.sh"
nohup "$SCRIPT_DIR/run_requant_baseline_subnormal_chain.sh" \
  >> "$OUT_DIR/requant_baseline_subnormal_chain.nohup.log" 2>&1 &
BPID=$!
disown $BPID 2>/dev/null
echo "[RESUME] baseline_subnormal chain PID=$BPID" | tee -a "$RESUME_LOG"

# 3. Tile chain (40k + 65k) — waits for baseline_subnormal chain.
sed -i "s/^WAIT_PID=.*/WAIT_PID=\${WAIT_PID:-$BPID}/" "$SCRIPT_DIR/run_requant_gt20k_tile_chain.sh"
nohup "$SCRIPT_DIR/run_requant_gt20k_tile_chain.sh" \
  >> "$OUT_DIR/requant_gt20k_tile_chain.nohup.log" 2>&1 &
TPID=$!
disown $TPID 2>/dev/null
echo "[RESUME] tile chain PID=$TPID" | tee -a "$RESUME_LOG"

echo "[RESUME END] $(date -Is)" | tee -a "$RESUME_LOG"
