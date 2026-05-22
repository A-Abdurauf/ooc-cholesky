#!/usr/bin/env bash
# True MXFP4 drop-in at 32k/nb=2048: legacy bound evaluates on MXFP8
# (MX_BOUND_LOW_AS=mx_e4m3) but storage is MXFP4 (MX_FORCE_FORMAT=e2m1).
# Runs both UNDERFLOW_MODE=fz (FTZ) and =gu (denorm), 4 epsilons each.
#
# Writes a small CSV: sweep,bin,n,nb,mx_mode,underflow,source_epsilon,
#                     rel_factor_error,abs_factor_error,relative_residual,
#                     tile_breakdown
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/examples/example_dpotrf_gpu"

BIN=${BIN:-"/home/abduraa/MX_project/logs/my_cov_weak_32k.bin"}
NB=${NB:-2048}
CORES=${CORES:-32}
EPS_LIST=${EPS_LIST:-"1e-5 1e-6 1e-7 1e-8"}

OUT_DIR=${OUT_DIR:-"$SCRIPT_DIR/true_dropin_mxfp4_32k"}
RUN_LOG_DIR="$OUT_DIR/run_logs"
CSV="$OUT_DIR/results.csv"
mkdir -p "$RUN_LOG_DIR"

if [[ ! -f "$CSV" ]]; then
  echo "sweep,bin,n,nb,mx_mode,underflow,source_epsilon,rel_factor_error,abs_factor_error,relative_residual,tile_breakdown" > "$CSV"
fi

run_one() {
  local underflow="$1" eps="$2"
  local sweep="true_mxfp4_dropin_${underflow}"
  local stamp
  stamp=$(date +"%Y%m%d_%H%M%S")
  local run_log="$RUN_LOG_DIR/${sweep}_eps${eps}_${stamp}.log"

  echo "[RUN] underflow=$underflow eps=$eps  log=$(basename "$run_log")"
  (
    cd "$SCRIPT_DIR/examples"
    exec env \
      MX_SKIP_KL=1 \
      MX_SELECTION_CRITERIA=bound \
      MX_BOUND_LADDER=legacy \
      MX_FORCE_FORMAT=e2m1 \
      MX_BOUND_LOW_AS=mx_e4m3 \
      MX_BOUND_DEBUG=1 \
      MX_MX_MODE=vec1d \
      MX_BLOCK_SUBTILE=32 \
      MX_MODE=vec1d \
      MX_UNDERFLOW_MODE="$underflow" \
      MX_FP8_SUBNORMAL=1 \
      MX_SCALE_AWARE_EPSILON=0 \
      MX_FP32_SCALE_BITS=11 \
      MX_ERROR_LEGACY=0 \
      MX_BUCKET_FP32=fp32 \
      MX_BUCKET_FP16=mx_fp16 \
      MX_SOURCE_EPSILON="$eps" \
      "$EXE" --nb "$NB" --cores "$CORES" --bin "$BIN"
  ) >"$run_log" 2>&1
  local rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "[FAIL] rc=$rc see $run_log"
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

  echo "[OK]   eps=$eps rel=$rel abs=$abs res=$res  tiles=$tb"
  echo "$sweep,$BIN,32768,$NB,vec1d,$underflow,$eps,${rel:-NA},${abs:-NA},${res:-NA},\"$tb\"" >> "$CSV"
}

for underflow in fz gu; do
  for eps in $EPS_LIST; do
    run_one "$underflow" "$eps"
  done
done

echo
echo "Done. CSV: $CSV"
