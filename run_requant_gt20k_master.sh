#!/usr/bin/env bash
# Wait for the currently in-flight 65k run (example_dpotrf_gpu --bin ...weak_65k.bin)
# to finish, then pause the 65k orchestrator, run the 4 missing 32k SFP4 sweeps
# standalone, then resume the orchestrator chain.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR=${OUT_DIR:-"/home/abduraa/MX_project/logs/mx_ooc_data"}
ORCH_LOG="$OUT_DIR/requant_sfp4_32k_after_65k.log"

# Track the PIDs we'll need to kill / wait on.
EXAMPLE_PID=${EXAMPLE_PID:-3778029}
ORCH_PID=${ORCH_PID:-3713997}
BASELINE_CHAIN_PID=${BASELINE_CHAIN_PID:-3714005}
TILE_CHAIN_PID=${TILE_CHAIN_PID:-3714019}

{
  echo "[SFP4-32K WATCHER START] $(date -Is)"
  echo "[INFO] waiting for example_dpotrf_gpu PID $EXAMPLE_PID to exit (current 65k run)"
} | tee -a "$ORCH_LOG"

# Wait for the currently in-flight binary to finish.
while [[ -d "/proc/$EXAMPLE_PID" ]]; do
  sleep 30
done
echo "[INFO] PID $EXAMPLE_PID exited at $(date -Is). Pausing orch + chains." | tee -a "$ORCH_LOG"

# Pause: kill the orch + chain wrappers, plus any new example_dpotrf_gpu that
# spawned in the seconds after we detected the exit.
kill -TERM "$ORCH_PID" "$BASELINE_CHAIN_PID" "$TILE_CHAIN_PID" 2>/dev/null || true
sleep 2
pkill -TERM -P "$ORCH_PID" 2>/dev/null || true
pkill -TERM -f "run_requant_gt20k_by_bin" 2>/dev/null || true
pkill -TERM -f "run_requant_baseline_subnormal_chain" 2>/dev/null || true
pkill -TERM -f "run_requant_gt20k_tile_chain" 2>/dev/null || true
pkill -TERM -f example_dpotrf_gpu 2>/dev/null || true
sleep 2
pkill -KILL -f "run_requant_" 2>/dev/null || true
pkill -KILL -f example_dpotrf_gpu 2>/dev/null || true
sleep 1
echo "[INFO] paused. GPU should be free now." | tee -a "$ORCH_LOG"
nvidia-smi --query-gpu=utilization.gpu,memory.used --format=csv 2>/dev/null | tee -a "$ORCH_LOG"

# Run the 4 missing SFP4 32k sweeps standalone via run_requant_gt20k_by_bin.sh
# limited to the 32k bin. Skip-logic preserves the existing 16 sweeps' 32k rows
# and only runs the new SFP4 ones (2 from main orch script + 2 from tile script).
echo "[INFO] launching main_orch (32k only) for SFP4 sweeps" | tee -a "$ORCH_LOG"
BIN_LIST="/home/abduraa/MX_project/logs/my_cov_weak_32k.bin" OUT_DIR="$OUT_DIR" \
  "$SCRIPT_DIR/run_requant_gt20k_by_bin.sh" 2>&1 | tee -a "$ORCH_LOG"

echo "[INFO] launching tile_orch (32k only) for SFP4 tile sweeps" | tee -a "$ORCH_LOG"
BIN_LIST="/home/abduraa/MX_project/logs/my_cov_weak_32k.bin" OUT_DIR="$OUT_DIR" \
  "$SCRIPT_DIR/run_requant_gt20k_tile_by_bin.sh" 2>&1 | tee -a "$ORCH_LOG"

echo "[INFO] 32k SFP4 sweeps done. Resuming main orch + chains for 65k." | tee -a "$ORCH_LOG"

# Resume — orchestrator skips done entries, picks up where 65k left off.
nohup "$SCRIPT_DIR/run_requant_gt20k_by_bin.sh" \
  >> "$OUT_DIR/requant_gt20k_by_bin.nohup.log" 2>&1 &
NEW_ORCH=$!
disown $NEW_ORCH 2>/dev/null
echo "[INFO] resumed main orch PID=$NEW_ORCH" | tee -a "$ORCH_LOG"

sed -i "s/^WAIT_PID=.*/WAIT_PID=\${WAIT_PID:-$NEW_ORCH}/" "$SCRIPT_DIR/run_requant_baseline_subnormal_chain.sh"
nohup "$SCRIPT_DIR/run_requant_baseline_subnormal_chain.sh" \
  >> "$OUT_DIR/requant_baseline_subnormal_chain.nohup.log" 2>&1 &
NEW_BCHAIN=$!
disown $NEW_BCHAIN 2>/dev/null
echo "[INFO] resumed baseline_subnormal chain PID=$NEW_BCHAIN" | tee -a "$ORCH_LOG"

sed -i "s/^WAIT_PID=.*/WAIT_PID=\${WAIT_PID:-$NEW_BCHAIN}/" "$SCRIPT_DIR/run_requant_gt20k_tile_chain.sh"
nohup "$SCRIPT_DIR/run_requant_gt20k_tile_chain.sh" \
  >> "$OUT_DIR/requant_gt20k_tile_chain.nohup.log" 2>&1 &
NEW_TCHAIN=$!
disown $NEW_TCHAIN 2>/dev/null
echo "[INFO] resumed tile chain PID=$NEW_TCHAIN" | tee -a "$ORCH_LOG"

echo "[SFP4-32K WATCHER END] $(date -Is)" | tee -a "$ORCH_LOG"
