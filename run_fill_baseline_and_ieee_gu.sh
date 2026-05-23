#!/usr/bin/env bash
# Fill-in: baseline (FP8 subnormal) + Ladder IEEE in GU underflow mode.
# Both modes already have FZ rows in the main CSV.  Adding GU rows so the
# bound selector's underflow handling is consistent across all configs in
# the merged CSV.
#
# 24 runs total: 2 sweeps x 3 bins x 4 eps.
#
# Resume-safe: skips (sweep, N, NB, eps) tuples already in the local CSV.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/examples/example_dpotrf_gpu"

CORES=${CORES:-32}
EPS_LIST=${EPS_LIST:-"1e-5 1e-6 1e-7 1e-8"}

BIN_32K="/home/abduraa/MX_project/logs/my_cov_weak_32k.bin"
BIN_40K="/home/abduraa/MX_project/logs/my_cov_weak_40k.bin"
BIN_65K="/home/abduraa/MX_project/logs/my_cov_weak_65k.bin"

OUT_DIR=${OUT_DIR:-"$SCRIPT_DIR/baseline_ieee_gu"}
RUN_LOG_DIR="$OUT_DIR/run_logs"
CSV="$OUT_DIR/results.csv"
MASTER_LOG="$OUT_DIR/sweep_master.log"
mkdir -p "$RUN_LOG_DIR"

if [[ ! -f "$CSV" ]]; then
  echo "sweep,bin,n,nb,mx_mode,underflow,source_epsilon,rel_factor_error,abs_factor_error,relative_residual,tile_breakdown" > "$CSV"
fi

have_row() {
  local sweep="$1" n="$2" nb="$3" eps="$4"
  awk -F, -v s="$sweep" -v n="$n" -v nb="$nb" -v e="$eps" '
    NR>1 && $1==s && $3==n && $4==nb && $7==e { found=1 }
    END { exit !found }
  ' "$CSV"
}

# Determine NB based on bin filename (matches the convention of other sweeps).
nb_for_bin() {
  case "$(basename "$1")" in
    *65k*) echo 4096 ;;
    *40k*) echo 4096 ;;
    *32k*) echo 2048 ;;
    *)     echo 2048 ;;
  esac
}
n_for_bin() {
  case "$(basename "$1")" in
    *65k*) echo 65536 ;;
    *40k*) echo 40960 ;;
    *32k*) echo 32768 ;;
    *)     echo 0 ;;
  esac
}

run_one() {
  local sweep="$1" ladder="$2" fp8_subnormal="$3" bin="$4" eps="$5"
  local nb n
  nb=$(nb_for_bin "$bin"); n=$(n_for_bin "$bin")

  if have_row "$sweep" "$n" "$nb" "$eps"; then
    echo "[SKIP] $sweep N=$n nb=$nb eps=$eps (already in CSV)" | tee -a "$MASTER_LOG"
    return
  fi

  local base=$(basename "$bin" .bin)
  local stamp=$(date +"%Y%m%d_%H%M%S")
  local run_log="$RUN_LOG_DIR/${sweep}_${base}_nb${nb}_eps${eps}_${stamp}.log"

  echo "[RUN] sweep=$sweep N=$n nb=$nb eps=$eps log=$(basename "$run_log")" | tee -a "$MASTER_LOG"

  (
    cd "$SCRIPT_DIR/examples"
    exec env \
      MX_SKIP_KL=1 \
      MX_SELECTION_CRITERIA=bound \
      MX_BOUND_LADDER="$ladder" \
      MX_BOUND_DEBUG=1 \
      MX_MX_MODE=tile \
      MX_MODE=tile \
      MX_UNDERFLOW_MODE=gu \
      MX_FP8_SUBNORMAL="$fp8_subnormal" \
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

  local rel abs res tb
  rel=$(grep -m1 "^relative_error:" "$run_log" | awk '{print $NF}')
  abs=$(grep -m1 "^error:"          "$run_log" | awk '{print $NF}')
  res=$(grep -m1 "^relative_residual:" "$run_log" | awk '{print $NF}')
  tb=$(grep -oE "^\[TILE_TARGET\] \([0-9]+, [0-9]+\) \S+" "$run_log" \
       | awk '{print $NF}' | sort | uniq -c | sort -rn \
       | awk '{printf "%s=%s;",$2,$1}')

  echo "[OK]  sweep=$sweep N=$n nb=$nb eps=$eps rel=$rel res=$res" | tee -a "$MASTER_LOG"
  echo "$sweep,$bin,$n,$nb,tile,gu,$eps,${rel:-NA},${abs:-NA},${res:-NA},\"$tb\"" >> "$CSV"
}

{
  echo "[START] $(date -Is)"
  echo "[CFG] underflow=gu cores=$CORES eps=$EPS_LIST"
  echo "[CFG] CSV=$CSV"
} | tee "$MASTER_LOG"

# Two sweeps:
#  - baseline_fp8_subnormal_gu : MX_BOUND_LADDER=legacy, MX_FP8_SUBNORMAL=1
#  - ladder_ieee_gu            : MX_BOUND_LADDER=ieee_only, MX_FP8_SUBNORMAL=0
declare -a CFG=(
  "baseline_fp8_subnormal_gu  legacy     1 $BIN_32K"
  "baseline_fp8_subnormal_gu  legacy     1 $BIN_40K"
  "baseline_fp8_subnormal_gu  legacy     1 $BIN_65K"
  "ladder_ieee_gu             ieee_only  0 $BIN_32K"
  "ladder_ieee_gu             ieee_only  0 $BIN_40K"
  "ladder_ieee_gu             ieee_only  0 $BIN_65K"
)

for entry in "${CFG[@]}"; do
  read -r sweep ladder fp8_sub bin <<< "$entry"
  for eps in $EPS_LIST; do
    run_one "$sweep" "$ladder" "$fp8_sub" "$bin" "$eps"
  done
done

echo "[END] $(date -Is)" | tee -a "$MASTER_LOG"
echo "CSV: $CSV"
