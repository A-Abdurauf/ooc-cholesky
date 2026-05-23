#!/usr/bin/env bash
# Fill-in sweep: close the 4-class x {FZ, GU} x {tile, vec1d32} error matrix
# at 40k (NB=4096) and 65k (NB=4096).  32k cells are already complete in the
# restored CSVs, so they are skipped.
#
# 28 runs total (each = one GPU job, ~minutes):
#   mx_staircase_tile_fz       @ 40k                : 4 eps
#   mx_staircase_tile_gu       @ 40k                : 4 eps
#   mx_staircase_vec1d32_fz    @ 40k                : 4 eps
#   ladder_full_tile_fz        @ 40k, 65k           : 8 eps
#   ladder_full_vec1d_fz       @ 40k, 65k           : 8 eps
#
# Skipped (already complete elsewhere):
#   - baseline / IEEE-FZ        : main CSV all N
#   - mx_staircase_vec1d32_gu   : stair_vec + stair_40n4k
#   - ladder_full_*_gu          : rerun_32k + gu_missing
#   - 32k cells generally       : restored from stash
#
# Resume-safe: skips (sweep, eps, N, nb) tuples whose CSV already contains them.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/examples/example_dpotrf_gpu"

CORES=${CORES:-32}
EPS_LIST=${EPS_LIST:-"1e-5 1e-6 1e-7 1e-8"}

BIN_40K="/home/abduraa/MX_project/logs/my_cov_weak_40k.bin"
BIN_65K="/home/abduraa/MX_project/logs/my_cov_weak_65k.bin"
NB=4096

OUT_DIR=${OUT_DIR:-"$SCRIPT_DIR/fill_4class_matrix"}
RUN_LOG_DIR="$OUT_DIR/run_logs"
CSV="$OUT_DIR/results.csv"
MASTER_LOG="$OUT_DIR/sweep_master.log"
mkdir -p "$RUN_LOG_DIR"

if [[ ! -f "$CSV" ]]; then
  echo "sweep,bin,n,nb,mx_mode,underflow,source_epsilon,rel_factor_error,abs_factor_error,relative_residual,tile_breakdown" > "$CSV"
fi

# Helper: skip if the row already exists in CSV.
have_row() {
  local sweep="$1" n="$2" nb="$3" eps="$4"
  awk -F, -v s="$sweep" -v n="$n" -v nb="$nb" -v e="$eps" '
    NR>1 && $1==s && $3==n && $4==nb && $7==e { found=1 }
    END { exit !found }
  ' "$CSV"
}

# Run one: (sweep, ladder, mx_mode, underflow, bin, n)
run_one() {
  local sweep="$1" ladder="$2" mx_mode="$3" underflow="$4" bin="$5" n="$6" eps="$7"

  if have_row "$sweep" "$n" "$NB" "$eps"; then
    echo "[SKIP] $sweep N=$n nb=$NB eps=$eps (already in CSV)" | tee -a "$MASTER_LOG"
    return
  fi

  local base=$(basename "$bin" .bin)
  local stamp=$(date +"%Y%m%d_%H%M%S")
  local run_log="$RUN_LOG_DIR/${sweep}_${base}_nb${NB}_eps${eps}_${stamp}.log"

  local subtile_env=""
  if [[ "$mx_mode" == "vec1d" ]]; then
    subtile_env="MX_BLOCK_SUBTILE=32"
  fi

  echo "[RUN] sweep=$sweep N=$n nb=$NB eps=$eps log=$(basename "$run_log")" | tee -a "$MASTER_LOG"

  (
    cd "$SCRIPT_DIR/examples"
    exec env \
      MX_SKIP_KL=1 \
      MX_SELECTION_CRITERIA=bound \
      MX_BOUND_LADDER="$ladder" \
      MX_BOUND_DEBUG=1 \
      MX_MX_MODE="$mx_mode" \
      $subtile_env \
      MX_MODE="$mx_mode" \
      MX_UNDERFLOW_MODE="$underflow" \
      MX_FP8_SUBNORMAL=0 \
      MX_SCALE_AWARE_EPSILON=1 \
      MX_FP32_SCALE_BITS=11 \
      MX_ERROR_LEGACY=0 \
      MX_SOURCE_EPSILON="$eps" \
      "$EXE" --nb "$NB" --cores "$CORES" --bin "$bin"
  ) >"$run_log" 2>&1
  local rc=$?
  if [[ $rc -ne 0 ]]; then
    echo "[FAIL] rc=$rc see $run_log" | tee -a "$MASTER_LOG"
    return
  fi

  local rel abs res tb
  rel=$(grep -m1 "^relative_error:" "$run_log" | awk '{print $NF}')
  abs=$(grep -m1 "^error:"          "$run_log" | awk '{print $NF}')
  res=$(grep -m1 "^relative_residual:" "$run_log" | awk '{print $NF}')
  tb=$(grep -oE "^\[TILE_TARGET\] \([0-9]+, [0-9]+\) \S+" "$run_log" \
       | awk '{print $NF}' | sort | uniq -c | sort -rn \
       | awk '{printf "%s=%s;",$2,$1}')

  echo "[OK]  sweep=$sweep N=$n nb=$NB eps=$eps rel=$rel res=$res" | tee -a "$MASTER_LOG"
  echo "$sweep,$bin,$n,$NB,$mx_mode,$underflow,$eps,${rel:-NA},${abs:-NA},${res:-NA},\"$tb\"" >> "$CSV"
}

{
  echo "[START] $(date -Is)"
  echo "[CFG] NB=$NB cores=$CORES eps=$EPS_LIST"
  echo "[CFG] CSV=$CSV"
} | tee "$MASTER_LOG"

# Configs: each row is (sweep, ladder, mx_mode, underflow, bin_list)
# We run all listed eps for each.
declare -a CFG=(
  # mx_staircase variants only at 40k (65k already done)
  "mx_staircase_tile_fz       mx_staircase tile  fz $BIN_40K"
  "mx_staircase_tile_gu       mx_staircase tile  gu $BIN_40K"
  "mx_staircase_vec1d32_fz    mx_staircase vec1d fz $BIN_40K"
  # ladder_full FZ at 40k + 65k for both granularities
  "ladder_full_tile_fz        full         tile  fz $BIN_40K"
  "ladder_full_tile_fz        full         tile  fz $BIN_65K"
  "ladder_full_vec1d_fz       full         vec1d fz $BIN_40K"
  "ladder_full_vec1d_fz       full         vec1d fz $BIN_65K"
)

for entry in "${CFG[@]}"; do
  read -r sweep ladder mx_mode underflow bin <<< "$entry"
  case "$bin" in
    *40k*) n=40960 ;;
    *65k*) n=65536 ;;
    *)     n=0    ;;
  esac
  for eps in $EPS_LIST; do
    run_one "$sweep" "$ladder" "$mx_mode" "$underflow" "$bin" "$n" "$eps"
  done
done

echo "[END] $(date -Is)" | tee -a "$MASTER_LOG"
echo "CSV: $CSV"
