#!/usr/bin/env bash
# Small fill-in: MX-staircase at N=40k, nb=4096, gu (denorm) underflow only.
# Fixes the 40k apples-to-apples problem in the 3-ladder GU figure where
# baseline + IEEE-ladder + Full-Ladder-MX are at nb=4096 but the original
# staircase sweep ran at nb=2048.
#
# 4 runs total: 40k/nb=4096 x {1e-5, 1e-6, 1e-7, 1e-8}, underflow=gu.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/examples/example_dpotrf_gpu"

BIN="/home/abduraa/MX_project/logs/my_cov_weak_40k.bin"
NB=4096
CORES=${CORES:-32}
EPS_LIST=${EPS_LIST:-"1e-5 1e-6 1e-7 1e-8"}

OUT_DIR=${OUT_DIR:-"$SCRIPT_DIR/mx_staircase_40k_nb4096_gu"}
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
  echo "[CFG] N=40960  nb=$NB  underflow=gu"
  echo "[CFG] eps=$EPS_LIST"
  echo "[CFG] csv=$CSV"
} | tee "$MASTER_LOG"

for eps in $EPS_LIST; do
  sweep="mx_staircase_vec1d32_gu"
  stamp=$(date +"%Y%m%d_%H%M%S")
  run_log="$RUN_LOG_DIR/${sweep}_my_cov_weak_40k_nb${NB}_eps${eps}_${stamp}.log"

  echo "[RUN] eps=$eps nb=$NB log=$(basename "$run_log")" | tee -a "$MASTER_LOG"

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
      MX_UNDERFLOW_MODE=gu \
      MX_FP8_SUBNORMAL=0 \
      MX_SCALE_AWARE_EPSILON=1 \
      MX_FP32_SCALE_BITS=11 \
      MX_ERROR_LEGACY=0 \
      MX_SOURCE_EPSILON="$eps" \
      "$EXE" --nb "$NB" --cores "$CORES" --bin "$BIN"
  ) >"$run_log" 2>&1
  rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "[FAIL] rc=$rc see $run_log" | tee -a "$MASTER_LOG"
    continue
  fi

  rel=$(grep -m1 "^relative_error:" "$run_log" | awk '{print $NF}')
  abs=$(grep -m1 "^error:" "$run_log" | awk '{print $NF}')
  res=$(grep -m1 "^relative_residual:" "$run_log" | awk '{print $NF}')
  tb=$(grep -oE "^\[TILE_TARGET\] \([0-9]+, [0-9]+\) \S+" "$run_log" \
       | awk '{print $NF}' | sort | uniq -c | sort -rn \
       | awk '{printf "%s=%s;",$2,$1}')

  echo "[OK] eps=$eps rel=$rel abs=$abs res=$res tiles=$tb" | tee -a "$MASTER_LOG"
  echo "$sweep,$BIN,40960,$NB,vec1d,gu,$eps,${rel:-NA},${abs:-NA},${res:-NA},\"$tb\"" >> "$CSV"
done

echo "[END] $(date -Is)" | tee -a "$MASTER_LOG"
echo "CSV: $CSV"
