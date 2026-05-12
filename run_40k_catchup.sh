#!/usr/bin/env bash
# Catch-up sequence: finish all 40k runs before resuming 65k.
# Order:
#   1. baseline_fp8_subnormal at 40k only (4 runs)
#   2. tile-mode 4 sweeps at 40k only (16 runs)
#   3. Restart main orchestrator (will skip 32k+40k, do 65k)
#   4. Restart baseline_subnormal chain (waits for orch, then runs 65k)
#   5. Restart tile chain (waits for baseline chain, then runs 65k for tiles)
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR=${OUT_DIR:-"/home/abduraa/MX_project/logs/mx_ooc_data"}
CATCHUP_LOG="$OUT_DIR/requant_40k_catchup.log"

echo "[CATCHUP START] $(date -Is)" | tee -a "$CATCHUP_LOG"

# Step 1: baseline_fp8_subnormal at 40k only.
echo "[CATCHUP] Step 1: baseline_fp8_subnormal at 40k only" | tee -a "$CATCHUP_LOG"
BUILD_FIRST=0 OUT_DIR="$OUT_DIR" \
  BIN_LIST="/home/abduraa/MX_project/logs/my_cov_weak_40k.bin" \
  "$SCRIPT_DIR/run_requant_baseline_fp8_subnormal_gt20k.sh" 2>&1 | tee -a "$CATCHUP_LOG"

# Step 2: tile-mode sweeps at 40k only.
echo "[CATCHUP] Step 2: tile-mode sweeps at 40k only" | tee -a "$CATCHUP_LOG"
BIN_LIST="/home/abduraa/MX_project/logs/my_cov_weak_40k.bin" OUT_DIR="$OUT_DIR" \
  "$SCRIPT_DIR/run_requant_gt20k_tile_by_bin.sh" 2>&1 | tee -a "$CATCHUP_LOG"

echo "[CATCHUP] 40k batch done. Resuming main orch + chains for 65k." | tee -a "$CATCHUP_LOG"

# Step 3: main orchestrator (will skip 32k+40k via skip-logic, run 65k).
nohup "$SCRIPT_DIR/run_requant_gt20k_by_bin.sh" \
  >> "$OUT_DIR/requant_gt20k_by_bin.nohup.log" 2>&1 &
ORCH_PID=$!
disown $ORCH_PID 2>/dev/null
echo "[CATCHUP] main orch PID=$ORCH_PID" | tee -a "$CATCHUP_LOG"

# Step 4: baseline_subnormal chain (40k done, 65k still needed) — waits for orch.
sed -i "s/^WAIT_PID=.*/WAIT_PID=\${WAIT_PID:-$ORCH_PID}/" "$SCRIPT_DIR/run_requant_baseline_subnormal_chain.sh"
nohup "$SCRIPT_DIR/run_requant_baseline_subnormal_chain.sh" \
  >> "$OUT_DIR/requant_baseline_subnormal_chain.nohup.log" 2>&1 &
BPID=$!
disown $BPID 2>/dev/null
echo "[CATCHUP] baseline_subnormal chain PID=$BPID" | tee -a "$CATCHUP_LOG"

# Step 5: tile chain (40k done, 65k still needed) — waits for baseline chain.
sed -i "s/^WAIT_PID=.*/WAIT_PID=\${WAIT_PID:-$BPID}/" "$SCRIPT_DIR/run_requant_gt20k_tile_chain.sh"
nohup "$SCRIPT_DIR/run_requant_gt20k_tile_chain.sh" \
  >> "$OUT_DIR/requant_gt20k_tile_chain.nohup.log" 2>&1 &
TPID=$!
disown $TPID 2>/dev/null
echo "[CATCHUP] tile chain PID=$TPID" | tee -a "$CATCHUP_LOG"

echo "[CATCHUP END] $(date -Is)" | tee -a "$CATCHUP_LOG"
