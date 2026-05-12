#!/usr/bin/env bash
# Chain-runner: waits for the current gt20k master sweep (if any) to finish,
# then launches the ladder_scaled_block128 sweep.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR=${OUT_DIR:-"/home/abduraa/MX_project/logs/mx_ooc_data"}
WAIT_PID=${WAIT_PID:-2943389}
CHAIN_LOG="$OUT_DIR/requant_ladder_block128_chain.log"

{
  echo "[CHAIN START] $(date -Is)"
  echo "[CHAIN] waiting for PID $WAIT_PID (running master sweep) to exit before launching ladder_block128"
} | tee -a "$CHAIN_LOG"

# Block until WAIT_PID is gone. /proc/<pid> disappears when the process exits.
while [[ -d "/proc/$WAIT_PID" ]]; do
  sleep 60
done
echo "[CHAIN] PID $WAIT_PID has exited at $(date -Is). Launching ladder_block128." | tee -a "$CHAIN_LOG"

# Build is unnecessary -- the executable was rebuilt by the previous master's first run.
BUILD_FIRST=0 OUT_DIR="$OUT_DIR" \
  "$SCRIPT_DIR/run_requant_ladder_scaled_block128_gt20k.sh" 2>&1 | tee -a "$CHAIN_LOG"

echo "[CHAIN END] $(date -Is)" | tee -a "$CHAIN_LOG"
