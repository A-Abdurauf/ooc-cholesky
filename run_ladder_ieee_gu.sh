#!/usr/bin/env bash
# IEEE-only ladder (FP8_E4M3 -> FP16 -> FP32 -> FP64) under MX_UNDERFLOW_MODE=gu.
# Mirrors the original LADDER_IEEE_ENV from run_requant_gt20k_by_bin.sh, with
# the single change MX_UNDERFLOW_MODE: fz -> gu. Output to standalone CSV so the
# new rows do not collide with the existing FTZ master-CSV rows.
#
# Sizes: 32k, 40k, 65k.  Epsilons: 1e-5, 1e-6, 1e-7, 1e-8.
# 12 runs total.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/examples/example_dpotrf_gpu"

BIN_LIST=${BIN_LIST:-"\
/home/abduraa/MX_project/logs/my_cov_weak_32k.bin \
/home/abduraa/MX_project/logs/my_cov_weak_40k.bin \
/home/abduraa/MX_project/logs/my_cov_weak_65k.bin"}
EPS_LIST=${EPS_LIST:-"1e-5 1e-6 1e-7 1e-8"}
CORES=${CORES:-32}

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

OUT_DIR=${OUT_DIR:-"$SCRIPT_DIR/ladder_ieee_gu"}
RUN_LOG_DIR="$OUT_DIR/run_logs"
CSV="$OUT_DIR/results.csv"
MASTER_LOG="$OUT_DIR/sweep_master.log"
mkdir -p "$RUN_LOG_DIR"

if [[ ! -f "$CSV" ]]; then
  echo "sweep,bin,n,nb,mx_mode,underflow,source_epsilon,rel_factor_error,abs_factor_error,relative_residual,tile_breakdown" > "$CSV"
fi

{
  echo "[START] $(date -Is)"
  echo "[CFG] ladder=ieee_only (fp8_e4m3 -> fp16 -> fp32 -> fp64)"
  echo "[CFG] underflow=gu (denorm) only"
  echo "[CFG] mx_mode=tile  (matches original LADDER_IEEE_ENV)"
  echo "[CFG] bins=$BIN_LIST"
  echo "[CFG] eps=$EPS_LIST"
  echo "[CFG] csv=$CSV"
} | tee "$MASTER_LOG"

run_one() {
  local bin="$1" eps="$2"
  local base=$(basename "$bin" .bin)
  local nb=$(nb_for_bin "$bin")
  local n=$(n_for_bin "$bin")
  local sweep="ladder_ieee_gu"
  local stamp=$(date +"%Y%m%d_%H%M%S")
  local run_log="$RUN_LOG_DIR/${sweep}_${base}_eps${eps}_${stamp}.log"

  echo "[RUN] bin=$base eps=$eps nb=$nb log=$(basename "$run_log")" | tee -a "$MASTER_LOG"

  (
    cd "$SCRIPT_DIR/examples"
    exec env \
      MX_SKIP_KL=1 \
      MX_FP8_SUBNORMAL=1 \
      MX_SELECTION_CRITERIA=bound \
      MX_BOUND_LADDER=ieee_only \
      MX_BOUND_DEBUG=1 \
      MX_MX_MODE=tile \
      MX_MODE=tile \
      MX_UNDERFLOW_MODE=gu \
      MX_SCALE_AWARE_EPSILON=0 \
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
  echo "$sweep,$bin,$n,$nb,tile,gu,$eps,${rel:-NA},${abs:-NA},${res:-NA},\"$tb\"" >> "$CSV"
}

for bin in $BIN_LIST; do
  if [[ ! -f "$bin" ]]; then
    echo "[WARN] missing bin: $bin" | tee -a "$MASTER_LOG"; continue
  fi
  for eps in $EPS_LIST; do
    run_one "$bin" "$eps"
  done
done

echo "[END] $(date -Is)" | tee -a "$MASTER_LOG"
echo "CSV: $CSV"
