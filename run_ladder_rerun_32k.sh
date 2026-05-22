#!/usr/bin/env bash
# Re-run Full-Ladder MX at 32k/nb=2048 on the patched binary, both FTZ and
# denorm modes, 4 epsilons each. Same envs as the original
# requant_ladder_scaled_vec1d32 sweep, but executed in one batch on a
# consistent binary so the comparison with the true MXFP4 drop-in sweep is
# apples-to-apples.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/examples/example_dpotrf_gpu"

BIN=${BIN:-"/home/abduraa/MX_project/logs/my_cov_weak_32k.bin"}
NB=${NB:-2048}
CORES=${CORES:-32}
EPS_LIST=${EPS_LIST:-"1e-5 1e-6 1e-7 1e-8"}

OUT_DIR=${OUT_DIR:-"$SCRIPT_DIR/ladder_rerun_32k"}
RUN_LOG_DIR="$OUT_DIR/run_logs"
CSV="$OUT_DIR/results.csv"
mkdir -p "$RUN_LOG_DIR"

if [[ ! -f "$CSV" ]]; then
  echo "sweep,bin,n,nb,mx_mode,underflow,source_epsilon,rel_factor_error,abs_factor_error,relative_residual,tile_breakdown" > "$CSV"
fi

run_one() {
  local underflow="$1" eps="$2"
  local sweep="ladder_full_${underflow}"
  local stamp
  stamp=$(date +"%Y%m%d_%H%M%S")
  local run_log="$RUN_LOG_DIR/${sweep}_eps${eps}_${stamp}.log"

  echo "[RUN] underflow=$underflow eps=$eps  log=$(basename "$run_log")"
  (
    cd "$SCRIPT_DIR/examples"
    exec env \
      MX_SKIP_KL=1 \
      MX_SELECTION_CRITERIA=bound \
      MX_BOUND_LADDER=full \
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
