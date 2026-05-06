#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_SCRIPT="$SCRIPT_DIR/build_run_relative.sh"

# Fixed target for now
BIN_PATH=${BIN_PATH:-"/home/abduraa/MX_project/logs/my_cov_weak_20k.bin"}

# Sweep dimensions (override as needed)
EPS_LIST=${EPS_LIST:-"1e-8 1e-7 1e-6 1e-5"}  # requested sweep range
MODE_LIST=${MODE_LIST:-"tile subtile128"}  # tile | subtile128

read -r -a EPS_ARR <<< "$EPS_LIST"
read -r -a MODE_ARR <<< "$MODE_LIST"

# Runtime knobs
CORES=${CORES:-32}
NB=${NB:-2048}
BUILD_FIRST=${BUILD_FIRST:-1}
MX_SKIP_KL=${MX_SKIP_KL:-1}
MX_SELECTION_CRITERIA=${MX_SELECTION_CRITERIA:-bound}
MX_BOUND_DEBUG=${MX_BOUND_DEBUG:-1}
MX_BOUND_LADDER=${MX_BOUND_LADDER:-full}
MX_MX_MODE=${MX_MX_MODE:-block}
MX_BLOCK_SUBTILE=${MX_BLOCK_SUBTILE:-128}
MX_UNDERFLOW_MODE=${MX_UNDERFLOW_MODE:-fz}
MX_SCALE_AWARE_EPSILON=${MX_SCALE_AWARE_EPSILON:-1}
MX_FP32_SCALE_BITS=${MX_FP32_SCALE_BITS:-11}
MX_MODE=${MX_MODE:-$MX_MX_MODE}

OUT_DIR=${OUT_DIR:-"/home/abduraa/MX_project/logs/mx_ooc_data"}
SWEEP_DIR="$OUT_DIR/bound_sweep_20k"
RUN_LOG_DIR="$SWEEP_DIR/run_logs"
SUMMARY_FILE=${SUMMARY_FILE:-"$SWEEP_DIR/summary_bound_20k.tsv"}
MASTER_LOG="$SWEEP_DIR/sweep_master.log"

mkdir -p "$RUN_LOG_DIR"

if [[ ! -x "$RUN_SCRIPT" ]]; then
  echo "[ERROR] Missing executable: $RUN_SCRIPT" >&2
  exit 1
fi
if [[ ! -f "$BIN_PATH" ]]; then
  echo "[ERROR] Missing input bin: $BIN_PATH" >&2
  exit 1
fi

{
  echo "[START] $(date -Is)"
  echo "[CFG] bin=$BIN_PATH"
  echo "[CFG] eps=$EPS_LIST"
  echo "[CFG] modes=$MODE_LIST"
  echo "[CFG] bound_ladder=$MX_BOUND_LADDER"
  echo "[CFG] summary=$SUMMARY_FILE"
} | tee "$MASTER_LOG"

run_idx=0
for eps in "${EPS_ARR[@]}"; do
  for mode in "${MODE_ARR[@]}"; do
      run_mode="$MX_MX_MODE"
      run_subtile="$MX_BLOCK_SUBTILE"
      mode_label="$mode"

      case "${mode,,}" in
        tile)
          run_mode="tile"
          run_subtile="0"
          mode_label="tile"
          ;;
        subtile128|block128|block_subtile128)
          run_mode="block"
          run_subtile="128"
          mode_label="subtile128"
          ;;
        block)
          run_mode="block"
          run_subtile="$MX_BLOCK_SUBTILE"
          mode_label="block"
          ;;
        *)
          echo "[WARN] Unknown mode '$mode' (expected tile/subtile128/block). Skipping." | tee -a "$MASTER_LOG"
          continue
          ;;
      esac

        run_idx=$((run_idx + 1))
        BUILD_VAL=0
        if [[ "$run_idx" -eq 1 && "$BUILD_FIRST" == "1" ]]; then
          BUILD_VAL=1
        fi

        stamp=$(date +"%Y%m%d_%H%M%S_%N")
        run_log="$RUN_LOG_DIR/run_20k_eps_${eps}_mode_${mode_label}_${stamp}.log"

        echo "[RUN] eps=$eps mode=$mode_label build=$BUILD_VAL" | tee -a "$MASTER_LOG"

        RUN_ENV=(
          MX_SKIP_KL="$MX_SKIP_KL"
          BUILD="$BUILD_VAL"
          MX_SELECTION_CRITERIA="$MX_SELECTION_CRITERIA"
          MX_BOUND_DEBUG="$MX_BOUND_DEBUG"
          MX_BOUND_LADDER="$MX_BOUND_LADDER"
          MX_MX_MODE="$run_mode"
          MX_BLOCK_SUBTILE="$run_subtile"
          MX_UNDERFLOW_MODE="$MX_UNDERFLOW_MODE"
          MX_SCALE_AWARE_EPSILON="$MX_SCALE_AWARE_EPSILON"
          MX_FP32_SCALE_BITS="$MX_FP32_SCALE_BITS"
          MX_MODE="$run_mode"
          NB="$NB"
          CORES="$CORES"
          BIN_PATH="$BIN_PATH"
          SUMMARY_FILE="$SUMMARY_FILE"
        )

        RUN_ENV+=(MX_SOURCE_EPSILON="$eps")

        if ! (cd "$SCRIPT_DIR" && env "${RUN_ENV[@]}" "$RUN_SCRIPT" > "$run_log" 2>&1); then
          echo "[FAIL] eps=$eps mode=$mode_label" | tee -a "$MASTER_LOG"
          continue
        fi

        # Validate that tile count strings are present in latest summary row
        python3 - <<PY >> "$MASTER_LOG"
from pathlib import Path
import csv

p = Path(r"""$SUMMARY_FILE""")
if not p.exists():
    print("[WARN] summary missing")
    raise SystemExit(0)
rows = list(csv.DictReader(p.open(), delimiter='\t'))
if not rows:
    print("[WARN] summary has no data rows")
    raise SystemExit(0)
last = rows[-1]
ok_lower = bool((last.get("tile_counts_lower") or "").strip())
ok_full = bool((last.get("tile_counts_full") or "").strip())
print(f"[OK] rel={last.get('rel_factor_error','')} abs={last.get('abs_factor_error','')} lower_counts={'yes' if ok_lower else 'no'} full_counts={'yes' if ok_full else 'no'}")
PY

        echo "[DONE] eps=$eps mode=$mode_label log=$run_log" | tee -a "$MASTER_LOG"
  done
done

echo "[END] $(date -Is)" | tee -a "$MASTER_LOG"
echo "Master log: $MASTER_LOG"
echo "Summary: $SUMMARY_FILE"
