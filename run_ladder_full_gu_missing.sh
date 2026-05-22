#!/usr/bin/env bash
# Full-Ladder MX (mx_fp16 rung) — denorm/gu underflow only — fills the gaps
# in master-CSV coverage so apples-to-apples comparison with the new
# mx_staircase sweeps is possible across all (granularity, underflow, size)
# cells.
#
# Missing cells filled here:
#   vec1D-32 × gu × {40k, 65k}     ← 8 runs
#   tile     × gu × {32k, 40k, 65k} ← 12 runs
# Total: 20 runs
#
# Apples-to-apples note: the staircase sweeps were run on the post-MX_BOUND_LOW_AS
# patched binary; this sweep also uses that binary, so all gu rows are produced
# under identical solver state.
#
# Output: standalone CSV at $OUT_DIR/results.csv.
set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXE="$SCRIPT_DIR/examples/example_dpotrf_gpu"

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

BIN_32K="/home/abduraa/MX_project/logs/my_cov_weak_32k.bin"
BIN_40K="/home/abduraa/MX_project/logs/my_cov_weak_40k.bin"
BIN_65K="/home/abduraa/MX_project/logs/my_cov_weak_65k.bin"

# (granularity, bin) cells we want gu for:
#   vec1D-32 × {40k, 65k}   (32k already covered by ladder_rerun_32k)
#   tile     × {32k, 40k, 65k}
CELLS=(
  "vec1d:$BIN_40K"
  "vec1d:$BIN_65K"
  "tile:$BIN_32K"
  "tile:$BIN_40K"
  "tile:$BIN_65K"
)

OUT_DIR=${OUT_DIR:-"$SCRIPT_DIR/ladder_full_gu_missing"}
RUN_LOG_DIR="$OUT_DIR/run_logs"
CSV="$OUT_DIR/results.csv"
MASTER_LOG="$OUT_DIR/sweep_master.log"
mkdir -p "$RUN_LOG_DIR"

if [[ ! -f "$CSV" ]]; then
  echo "sweep,bin,n,nb,mx_mode,underflow,source_epsilon,rel_factor_error,abs_factor_error,relative_residual,tile_breakdown" > "$CSV"
fi

{
  echo "[START] $(date -Is)"
  echo "[CFG] ladder=full (e2m1 -> mx_e4m3 -> mx_fp16 -> fp32 -> fp64)"
  echo "[CFG] underflow=gu only (denorm)"
  echo "[CFG] cells: ${CELLS[*]}"
  echo "[CFG] eps=$EPS_LIST"
  echo "[CFG] csv=$CSV"
} | tee "$MASTER_LOG"

run_one() {
  local granularity="$1" bin="$2" eps="$3"
  local base=$(basename "$bin" .bin)
  local nb=$(nb_for_bin "$bin")
  local n=$(n_for_bin "$bin")
  local sweep="ladder_full_${granularity}_gu"
  local stamp=$(date +"%Y%m%d_%H%M%S")
  local run_log="$RUN_LOG_DIR/${sweep}_${base}_eps${eps}_${stamp}.log"

  echo "[RUN] granularity=$granularity bin=$base eps=$eps nb=$nb log=$(basename "$run_log")" | tee -a "$MASTER_LOG"

  local mode_envs=()
  if [[ "$granularity" == "vec1d" ]]; then
    mode_envs=(MX_MX_MODE=vec1d MX_BLOCK_SUBTILE=32 MX_MODE=vec1d)
  else
    mode_envs=(MX_MX_MODE=tile MX_MODE=tile)
  fi

  (
    cd "$SCRIPT_DIR/examples"
    exec env \
      MX_SKIP_KL=1 \
      MX_SELECTION_CRITERIA=bound \
      MX_BOUND_LADDER=full \
      MX_BOUND_DEBUG=1 \
      "${mode_envs[@]}" \
      MX_UNDERFLOW_MODE=gu \
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

  echo "[OK] granularity=$granularity bin=$base eps=$eps rel=$rel abs=$abs res=$res tiles=$tb" | tee -a "$MASTER_LOG"
  echo "$sweep,$bin,$n,$nb,$granularity,gu,$eps,${rel:-NA},${abs:-NA},${res:-NA},\"$tb\"" >> "$CSV"
}

for cell in "${CELLS[@]}"; do
  granularity="${cell%%:*}"
  bin="${cell##*:}"
  if [[ ! -f "$bin" ]]; then
    echo "[WARN] missing bin: $bin" | tee -a "$MASTER_LOG"; continue
  fi
  for eps in $EPS_LIST; do
    run_one "$granularity" "$bin" "$eps"
  done
done

echo "[END] $(date -Is)" | tee -a "$MASTER_LOG"
echo "CSV: $CSV"
