#!/usr/bin/env bash
# MX-staircase ladder sweep — vec1D-32 granularity.
# Ladder: e2m1 (MXFP4) -> mx_e4m3 (MXFP8) -> fp16 (plain IEEE) -> fp32 -> fp64.
# Mirrors the canonical Full-Ladder sweep but with the FP16 rung as plain IEEE FP16
# (no shared scale) instead of MXFP16. Runs both fz and gu underflow modes.
# Sizes: {32k, 40k, 65k} (gt20k); ε: {1e-8 .. 1e-5}.
#
# Output: standalone CSV at $OUT_DIR/results.csv (NOT appended to requant_gt20k_memory.csv).
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/examples/example_dpotrf_gpu"

BIN_LIST=${BIN_LIST:-"\
/home/abduraa/MX_project/logs/my_cov_weak_32k.bin \
/home/abduraa/MX_project/logs/my_cov_weak_40k.bin \
/home/abduraa/MX_project/logs/my_cov_weak_65k.bin"}
EPS_LIST=${EPS_LIST:-"1e-5 1e-6 1e-7 1e-8"}
CORES=${CORES:-32}

# Per-bin NB: 32k/40k -> 2048, 65k -> 4096.
nb_for_bin() {
  case "$(basename "$1")" in
    *65k*) echo 4096 ;;
    *40k*) echo 2048 ;;
    *32k*) echo 2048 ;;
    *)     echo "${NB_DEFAULT:-2048}" ;;
  esac
}

n_for_bin() {
  case "$(basename "$1")" in
    *32k*) echo 32768 ;;
    *40k*) echo 40960 ;;
    *65k*) echo 65536 ;;
    *)     echo "0" ;;
  esac
}

OUT_DIR=${OUT_DIR:-"$SCRIPT_DIR/mx_staircase_vec1d32_gt20k"}
RUN_LOG_DIR="$OUT_DIR/run_logs"
CSV="$OUT_DIR/results.csv"
MASTER_LOG="$OUT_DIR/sweep_master.log"
mkdir -p "$RUN_LOG_DIR"

if [[ ! -f "$CSV" ]]; then
  echo "sweep,bin,n,nb,mx_mode,underflow,source_epsilon,rel_factor_error,abs_factor_error,relative_residual,tile_breakdown" > "$CSV"
fi

{
  echo "[START] $(date -Is)"
  echo "[CFG] ladder=mx_staircase (e2m1 -> mx_e4m3 -> fp16(plain) -> fp32 -> fp64)"
  echo "[CFG] mx_mode=vec1d subtile=32"
  echo "[CFG] bins=$BIN_LIST"
  echo "[CFG] eps=$EPS_LIST"
  echo "[CFG] underflow_modes=fz gu"
  echo "[CFG] csv=$CSV"
} | tee "$MASTER_LOG"

run_one() {
  local bin="$1" underflow="$2" eps="$3"
  local base=$(basename "$bin" .bin)
  local nb=$(nb_for_bin "$bin")
  local n=$(n_for_bin "$bin")
  local sweep="mx_staircase_vec1d32_${underflow}"
  local stamp=$(date +"%Y%m%d_%H%M%S")
  local run_log="$RUN_LOG_DIR/${sweep}_${base}_eps${eps}_${stamp}.log"

  echo "[RUN] bin=$base underflow=$underflow eps=$eps nb=$nb log=$(basename "$run_log")" | tee -a "$MASTER_LOG"

  (
    cd "$SCRIPT_DIR/examples"
    exec env \
      MX_SKIP_KL=1 \
      MX_SELECTION_CRITERIA=bound \
      MX_BOUND_LADDER=mx_staircase \
      MX_BOUND_DEBUG=1 \
      MX_MX_MODE=vec1d \
      MX_BLOCK_SUBTILE=32 \
      MX_MODE=vec1d \
      MX_UNDERFLOW_MODE="$underflow" \
      MX_FP8_SUBNORMAL=0 \
      MX_SCALE_AWARE_EPSILON=1 \
      MX_FP32_SCALE_BITS=11 \
      MX_ERROR_LEGACY=0 \
      MX_SOURCE_EPSILON="$eps" \
      "$EXE" --nb "$nb" --cores "$CORES" --bin "$bin"
  ) >"$run_log" 2>&1
  local rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "[FAIL] rc=$rc see $run_log" | tee -a "$MASTER_LOG"
    return
  fi

  local rel abs res
  rel=$(grep -m1 "^relative_error:" "$run_log" | awk '{print $NF}')
  abs=$(grep -m1 "^error:" "$run_log" | awk '{print $NF}')
  res=$(grep -m1 "^relative_residual:" "$run_log" | awk '{print $NF}')

  local tb
  tb=$(grep -oE "^\[TILE_TARGET\] \([0-9]+, [0-9]+\) \S+" "$run_log" \
       | awk '{print $NF}' | sort | uniq -c | sort -rn \
       | awk '{printf "%s=%s;",$2,$1}')

  echo "[OK] bin=$base eps=$eps rel=$rel abs=$abs res=$res tiles=$tb" | tee -a "$MASTER_LOG"
  echo "$sweep,$bin,$n,$nb,vec1d,$underflow,$eps,${rel:-NA},${abs:-NA},${res:-NA},\"$tb\"" >> "$CSV"
}

for underflow in fz gu; do
  for bin in $BIN_LIST; do
    if [[ ! -f "$bin" ]]; then
      echo "[WARN] missing bin: $bin" | tee -a "$MASTER_LOG"; continue
    fi
    for eps in $EPS_LIST; do
      run_one "$bin" "$underflow" "$eps"
    done
  done
done

echo "[END] $(date -Is)" | tee -a "$MASTER_LOG"
echo "CSV: $CSV"
