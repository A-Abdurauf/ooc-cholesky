#!/usr/bin/env bash
# Baseline plain-FP8 (E4M3) sweep with OCP-spec subnormal support (no FTZ).
# Identical selection logic to run_requant_baseline_fp8_gt20k.sh (legacy
# epsilon-cutoff, no bound, no ladder), but invokes the patched quantizeFp
# via MX_FP8_SUBNORMAL=1 so the low-tier FP8 path encodes subnormals per
# the OCP FP8 specification (matches real H100/MI300 hardware) instead of
# flushing them to zero.
# Runs only on bins with N > 20k (weak_32k, weak_40k, weak_65k).
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_SCRIPT="$SCRIPT_DIR/build_run_relative.sh"

export PATH="/tools/cmake-3.29.2/bin:$PATH"
export CC=/usr/bin/gcc-11
export CXX=/usr/bin/g++-11

BIN_LIST=${BIN_LIST:-"\
/home/abduraa/MX_project/logs/my_cov_weak_32k.bin \
/home/abduraa/MX_project/logs/my_cov_weak_40k.bin \
/home/abduraa/MX_project/logs/my_cov_weak_65k.bin"}
EPS_LIST=${EPS_LIST:-"1e-8 1e-7 1e-6 1e-5"}
read -r -a EPS_ARR <<< "$EPS_LIST"

CORES=${CORES:-32}
BUILD_FIRST=${BUILD_FIRST:-1}

# Per-bin NB lookup: 32k/40k -> 2048, 65k -> 4096.
nb_for_bin() {
  local bin="$1"
  case "$(basename "$bin")" in
    *65k*)  echo 4096 ;;
    *40k*)  echo 4096 ;;
    *32k*)  echo 2048 ;;
    *)      echo "${NB_DEFAULT:-2048}" ;;
  esac
}

OUT_DIR=${OUT_DIR:-"/home/abduraa/MX_project/logs/mx_ooc_data"}
SWEEP_DIR="$OUT_DIR/requant_baseline_fp8_subnormal_gt20k"
RUN_LOG_DIR="$SWEEP_DIR/run_logs"
SUMMARY_FILE=${SUMMARY_FILE:-"$SWEEP_DIR/summary_requant_baseline_fp8_subnormal_gt20k.tsv"}
MASTER_LOG="$SWEEP_DIR/sweep_master.log"
mkdir -p "$RUN_LOG_DIR"

{
  echo "[START] $(date -Is)"
  echo "[CFG] sweep=requant_baseline_fp8_subnormal_gt20k"
  echo "[CFG] bins=$BIN_LIST"
  echo "[CFG] eps=$EPS_LIST"
  echo "[CFG] selection=legacy_epsilon_cutoff (no bound, no ladder)"
  echo "[CFG] buckets: fp32->fp32 fp16->fp16 low->fp8_e4m3 (FP64 fallback)"
  echo "[CFG] mx_mode=tile (no block scaling for plain FP8)"
  echo "[CFG] underflow=OCP_SUBNORMAL (MX_FP8_SUBNORMAL=1)"
  echo "[CFG] NB: 32k/40k=2048  65k=4096"
  echo "[CFG] summary=$SUMMARY_FILE"
} | tee "$MASTER_LOG"

run_idx=0
for BIN_PATH in $BIN_LIST; do
  if [[ ! -f "$BIN_PATH" ]]; then
    echo "[WARN] missing bin, skipping: $BIN_PATH" | tee -a "$MASTER_LOG"
    continue
  fi
  base=$(basename "$BIN_PATH" .bin)
  NB=$(nb_for_bin "$BIN_PATH")

  for eps in "${EPS_ARR[@]}"; do
    run_idx=$((run_idx + 1))
    BUILD_VAL=0
    if [[ "$run_idx" -eq 1 && "$BUILD_FIRST" == "1" ]]; then
      BUILD_VAL=1
    fi

    stamp=$(date +"%Y%m%d_%H%M%S_%N")
    run_log="$RUN_LOG_DIR/run_${base}_eps_${eps}_${stamp}.log"
    echo "[RUN] bin=$base eps=$eps nb=$NB build=$BUILD_VAL" | tee -a "$MASTER_LOG"

    # Legacy epsilon-cutoff selector + OCP subnormal handling for FP8 low path.
    RUN_ENV=(
      MX_SKIP_KL=1
      BUILD="$BUILD_VAL"
      MX_MX_MODE=tile
      MX_UNDERFLOW_MODE=fz
      MX_FP8_SUBNORMAL=1
      MX_SCALE_AWARE_EPSILON=1
      MX_FP32_SCALE_BITS=11
      MX_MODE=tile
      MX_ERROR_LEGACY=0
      MX_BUCKET_FP32=fp32
      MX_BUCKET_FP16=fp16
      FORMAT=fp8_e4m3
      NB="$NB"
      CORES="$CORES"
      BIN_PATH="$BIN_PATH"
      SUMMARY_FILE="$SUMMARY_FILE"
      MX_SOURCE_EPSILON="$eps"
    )

    if ! (cd "$SCRIPT_DIR" && env "${RUN_ENV[@]}" "$RUN_SCRIPT" > "$run_log" 2>&1); then
      echo "[FAIL] bin=$base eps=$eps log=$run_log" | tee -a "$MASTER_LOG"
      continue
    fi

    REL=$(grep -m1 "^relative_error:" "$run_log" | awk '{print $NF}')
    ABS=$(grep -m1 "^error:" "$run_log" | awk '{print $NF}')
    RESID=$(grep -m1 "^relative_residual:" "$run_log" | awk '{print $NF}')
    echo "[DONE] bin=$base eps=$eps nb=$NB rel=$REL abs=$ABS resid=$RESID log=$run_log" | tee -a "$MASTER_LOG"
  done
done

echo "[END] $(date -Is)" | tee -a "$MASTER_LOG"
echo "Master log: $MASTER_LOG"
echo "Summary:    $SUMMARY_FILE"
