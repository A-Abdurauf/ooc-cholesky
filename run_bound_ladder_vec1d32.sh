#!/usr/bin/env bash
# Bound-ladder sweep using true MX 1D 32-element groups (vec1d).
# Ladder (set in dpotrf_mixed_precision.cpp): e2m1 -> mx_e4m3 -> mx_fp16 -> fp32 -> FP64.
# Records relative factorization error in the summary TSV.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_SCRIPT="$SCRIPT_DIR/build_run_relative.sh"

# Ensure modern cmake is on PATH (system cmake is 3.3.2, project needs >=3.15)
export PATH="/tools/cmake-3.29.2/bin:$PATH"
export CC=/usr/bin/gcc-11
export CXX=/usr/bin/g++-11

BIN_PATH=${BIN_PATH:-"/home/abduraa/MX_project/logs/my_cov_weak_20k.bin"}
EPS_LIST=${EPS_LIST:-"1e-8 1e-7 1e-6 1e-5"}
read -r -a EPS_ARR <<< "$EPS_LIST"

CORES=${CORES:-32}
NB=${NB:-2048}
BUILD_FIRST=${BUILD_FIRST:-1}

OUT_DIR=${OUT_DIR:-"/home/abduraa/MX_project/logs/mx_ooc_data"}
SWEEP_DIR="$OUT_DIR/bound_ladder_vec1d32"
RUN_LOG_DIR="$SWEEP_DIR/run_logs"
SUMMARY_FILE=${SUMMARY_FILE:-"$SWEEP_DIR/summary_bound_ladder_vec1d32.tsv"}
MASTER_LOG="$SWEEP_DIR/sweep_master.log"
mkdir -p "$RUN_LOG_DIR"

if [[ ! -f "$BIN_PATH" ]]; then
  echo "[ERROR] Missing input bin: $BIN_PATH" >&2
  exit 1
fi

{
  echo "[START] $(date -Is)"
  echo "[CFG] bin=$BIN_PATH"
  echo "[CFG] eps=$EPS_LIST"
  echo "[CFG] mode=vec1d subtile=32 ladder=full(custom: e2m1->mx_e4m3->mx_fp16->fp32)"
  echo "[CFG] summary=$SUMMARY_FILE"
} | tee "$MASTER_LOG"

run_idx=0
for eps in "${EPS_ARR[@]}"; do
  run_idx=$((run_idx + 1))
  BUILD_VAL=0
  if [[ "$run_idx" -eq 1 && "$BUILD_FIRST" == "1" ]]; then
    BUILD_VAL=1
  fi

  stamp=$(date +"%Y%m%d_%H%M%S_%N")
  run_log="$RUN_LOG_DIR/run_eps_${eps}_${stamp}.log"
  echo "[RUN] eps=$eps build=$BUILD_VAL" | tee -a "$MASTER_LOG"

  RUN_ENV=(
    MX_SKIP_KL=1
    BUILD="$BUILD_VAL"
    MX_SELECTION_CRITERIA=bound
    MX_BOUND_DEBUG=1
    MX_BOUND_LADDER=full
    MX_MX_MODE=vec1d
    MX_BLOCK_SUBTILE=32
    MX_UNDERFLOW_MODE=fz
    MX_SCALE_AWARE_EPSILON=1
    MX_FP32_SCALE_BITS=11
    MX_MODE=vec1d
    NB="$NB"
    CORES="$CORES"
    BIN_PATH="$BIN_PATH"
    SUMMARY_FILE="$SUMMARY_FILE"
    MX_SOURCE_EPSILON="$eps"
  )

  if ! (cd "$SCRIPT_DIR" && env "${RUN_ENV[@]}" "$RUN_SCRIPT" > "$run_log" 2>&1); then
    echo "[FAIL] eps=$eps log=$run_log" | tee -a "$MASTER_LOG"
    continue
  fi

  REL=$(grep -m1 "^relative_error:" "$run_log" | awk '{print $NF}')
  ABS=$(grep -m1 "^error:" "$run_log" | awk '{print $NF}')
  RESID=$(grep -m1 "^relative_residual:" "$run_log" | awk '{print $NF}')
  echo "[DONE] eps=$eps rel=$REL abs=$ABS resid=$RESID log=$run_log" | tee -a "$MASTER_LOG"
done

echo "[END] $(date -Is)" | tee -a "$MASTER_LOG"
echo "Master log: $MASTER_LOG"
echo "Summary:    $SUMMARY_FILE"
