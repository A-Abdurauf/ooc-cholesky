#!/usr/bin/env bash
# Fill ONLY the missing cell from run_ladder_full_gu_missing.sh:
#   ladder=full × vec1d-32 × gu × 40k × nb=4096   (4 eps values)
#
# The original sweep script (run_ladder_full_gu_missing.sh) hard-codes
# nb=2048 for 40k via nb_for_bin(), so the 40k @ nb=4096 cell needed to
# match the staircase sweep (mx_staircase_40k_nb4096_gu) was never run.
# This script reproduces the exact env of the original sweep for that
# one cell and appends rows to the same CSV.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/examples/example_dpotrf_gpu"

EPS_LIST=${EPS_LIST:-"1e-5 1e-6 1e-7 1e-8"}
CORES=${CORES:-32}

BIN="/home/abduraa/MX_project/logs/my_cov_weak_40k.bin"
N=40960
NB=4096
GRANULARITY="vec1d"
SWEEP="ladder_full_${GRANULARITY}_gu"

OUT_DIR=${OUT_DIR:-"$SCRIPT_DIR/ladder_full_gu_missing"}
RUN_LOG_DIR="$OUT_DIR/run_logs"
CSV="$OUT_DIR/results.csv"
MASTER_LOG="$OUT_DIR/sweep_master.log"
mkdir -p "$RUN_LOG_DIR"

if [[ ! -f "$CSV" ]]; then
  echo "sweep,bin,n,nb,mx_mode,underflow,source_epsilon,rel_factor_error,abs_factor_error,relative_residual,tile_breakdown" > "$CSV"
fi

{
  echo "[START-FILL] $(date -Is)"
  echo "[CFG] fill cell: vec1d × gu × 40k × nb=4096"
  echo "[CFG] eps=$EPS_LIST  csv=$CSV"
} | tee -a "$MASTER_LOG"

run_one() {
  local eps="$1"
  local base=$(basename "$BIN" .bin)
  local stamp=$(date +"%Y%m%d_%H%M%S")
  local run_log="$RUN_LOG_DIR/${SWEEP}_${base}_nb${NB}_eps${eps}_${stamp}.log"

  echo "[RUN] granularity=$GRANULARITY bin=$base eps=$eps nb=$NB log=$(basename "$run_log")" | tee -a "$MASTER_LOG"

  (
    cd "$SCRIPT_DIR/examples"
    exec env \
      MX_SKIP_KL=1 \
      MX_SELECTION_CRITERIA=bound \
      MX_BOUND_LADDER=full \
      MX_BOUND_DEBUG=1 \
      MX_MX_MODE=vec1d MX_BLOCK_SUBTILE=32 MX_MODE=vec1d \
      MX_UNDERFLOW_MODE=gu \
      MX_FP8_SUBNORMAL=0 \
      MX_SCALE_AWARE_EPSILON=1 \
      MX_FP32_SCALE_BITS=11 \
      MX_ERROR_LEGACY=0 \
      MX_SOURCE_EPSILON="$eps" \
      "$EXE" --nb "$NB" --cores "$CORES" --bin "$BIN"
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

  echo "[OK] eps=$eps rel=$rel abs=$abs res=$res tiles=$tb" | tee -a "$MASTER_LOG"
  echo "$SWEEP,$BIN,$N,$NB,$GRANULARITY,gu,$eps,${rel:-NA},${abs:-NA},${res:-NA},\"$tb\"" >> "$CSV"
}

if [[ ! -f "$BIN" ]]; then
  echo "[FATAL] missing bin: $BIN" | tee -a "$MASTER_LOG"; exit 1
fi

for eps in $EPS_LIST; do
  run_one "$eps"
done

echo "[END-FILL] $(date -Is)" | tee -a "$MASTER_LOG"
echo "CSV: $CSV"
