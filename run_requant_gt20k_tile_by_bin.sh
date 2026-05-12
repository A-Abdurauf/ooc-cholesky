#!/usr/bin/env bash
# Tile-mode (1 shared scale per NB×NB tile) MX sweep — coarsest granularity.
# Same per-tile selection as the canonical legacy / lowmidscale / lowscale / ladder
# sweeps; only difference is storage geometry (tile vs block 128 vs vec1D 32).
# Iterates (bin × sweep × eps) like run_requant_gt20k_by_bin.sh.
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
OUT_DIR=${OUT_DIR:-"/home/abduraa/MX_project/logs/mx_ooc_data"}
mkdir -p "$OUT_DIR"
ORCH_LOG="$OUT_DIR/requant_gt20k_tile_by_bin.log"

nb_for_bin() {
  local bin="$1"
  case "$(basename "$bin")" in
    *65k*)  echo 4096 ;;
    *40k*)  echo 4096 ;;
    *32k*)  echo 2048 ;;
    *)      echo "${NB_DEFAULT:-2048}" ;;
  esac
}

already_done() {
  local sum="$1" bin="$2" eps="$3"
  [[ -f "$sum" ]] || return 1
  awk -F'\t' -v b="$bin" -v e="$eps" 'NR>1 && $1==b && $7==e {found=1; exit} END{exit !found}' "$sum"
}

run_sweep_bin_eps() {
  local sweep_name="$1"; shift
  local bin="$1"; shift
  local eps="$1"; shift
  local nb="$1"; shift
  local sweep_dir="$OUT_DIR/$sweep_name"
  local run_log_dir="$sweep_dir/run_logs"
  local sum="$sweep_dir/summary_${sweep_name}.tsv"
  local master_log="$sweep_dir/sweep_master.log"
  mkdir -p "$run_log_dir"

  if already_done "$sum" "$bin" "$eps"; then
    echo "[SKIP] sweep=$sweep_name bin=$(basename $bin) eps=$eps (already in summary)" | tee -a "$ORCH_LOG"
    return 0
  fi

  local stamp run_log
  stamp=$(date +"%Y%m%d_%H%M%S_%N")
  run_log="$run_log_dir/run_$(basename "$bin" .bin)_eps_${eps}_${stamp}.log"
  echo "[RUN] sweep=$sweep_name bin=$(basename $bin) eps=$eps nb=$nb" | tee -a "$ORCH_LOG" "$master_log"

  local rc
  if (cd "$SCRIPT_DIR" && env "$@" \
        BUILD=0 \
        NB="$nb" \
        CORES="$CORES" \
        BIN_PATH="$bin" \
        SUMMARY_FILE="$sum" \
        MX_SOURCE_EPSILON="$eps" \
        "$RUN_SCRIPT" > "$run_log" 2>&1); then
    rc=0
  else
    rc=$?
  fi

  local rel abs resid
  rel=$(grep -m1 "^relative_error:" "$run_log" | awk '{print $NF}')
  abs=$(grep -m1 "^error:" "$run_log" | awk '{print $NF}')
  resid=$(grep -m1 "^relative_residual:" "$run_log" | awk '{print $NF}')
  if [[ "$rc" -eq 0 ]]; then
    echo "[DONE] sweep=$sweep_name bin=$(basename $bin) eps=$eps nb=$nb rel=$rel abs=$abs resid=$resid log=$run_log" | tee -a "$ORCH_LOG" "$master_log"
  else
    echo "[FAIL] sweep=$sweep_name bin=$(basename $bin) eps=$eps nb=$nb rc=$rc log=$run_log" | tee -a "$ORCH_LOG" "$master_log"
  fi
}

# Tile-mode env blocks. MX_MX_MODE=tile gives 1 shared scale per NB×NB tile.
# All four variants honour MX_FP8_SUBNORMAL=1 (OCP gradual underflow) and Path A
# (legacy ε-cutoff selector) so tile assignments match baseline cuts.
LEGACY_TILE_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_MX_MODE=tile
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=tile
  MX_ERROR_LEGACY=0
  MX_BUCKET_FP32=mx_fp32
  MX_BUCKET_FP16=mx_fp16
  FORMAT=mx_e4m3
)
LOWMIDSCALE_TILE_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_MX_MODE=tile
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=tile
  MX_ERROR_LEGACY=0
  MX_BUCKET_FP32=fp32
  MX_BUCKET_FP16=mx_fp16
  FORMAT=mx_e4m3
)
LOWSCALE_TILE_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_MX_MODE=tile
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=tile
  MX_ERROR_LEGACY=0
  MX_BUCKET_FP32=fp32
  MX_BUCKET_FP16=fp16
  FORMAT=mx_e4m3
)
# Same as lowscale_tile but with Scaled FP4 (E2M1) as the low tier.
LOWSCALE_E2M1_TILE_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_MX_MODE=tile
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=tile
  MX_ERROR_LEGACY=0
  MX_BUCKET_FP32=fp32
  MX_BUCKET_FP16=fp16
  FORMAT=e2m1
)
# Same as lowmidscale_tile but with Scaled FP4 (E2M1) as the low tier.
LOWMIDSCALE_E2M1_TILE_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_MX_MODE=tile
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=tile
  MX_ERROR_LEGACY=0
  MX_BUCKET_FP32=fp32
  MX_BUCKET_FP16=mx_fp16
  FORMAT=e2m1
)
# Ladder + tile: full ladder selection (bound formula), tile-mode storage.
LADDER_TILE_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_SELECTION_CRITERIA=bound
  MX_BOUND_DEBUG=1
  MX_BOUND_LADDER=full
  MX_MX_MODE=tile
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=1
  MX_FP32_SCALE_BITS=11
  MX_MODE=tile
  MX_ERROR_LEGACY=0
)

{
  echo "[ORCH START] $(date -Is)"
  echo "[ORCH] sweep set: tile-mode (legacy_scaled, lowmidscale, lowscale, ladder_scaled)"
  echo "[ORCH] BIN_LIST=$BIN_LIST"
  echo "[ORCH] EPS_LIST=$EPS_LIST"
  echo "[ORCH] order: bin -> sweep -> eps"
} | tee -a "$ORCH_LOG"

for bin in $BIN_LIST; do
  if [[ ! -f "$bin" ]]; then
    echo "[WARN] missing bin: $bin" | tee -a "$ORCH_LOG"
    continue
  fi
  nb=$(nb_for_bin "$bin")
  echo "[ORCH] >>> entering bin=$(basename $bin) nb=$nb" | tee -a "$ORCH_LOG"

  for eps in "${EPS_ARR[@]}"; do
    run_sweep_bin_eps "requant_legacy_scaled_tile_gt20k" "$bin" "$eps" "$nb" "${LEGACY_TILE_ENV[@]}"
  done
  for eps in "${EPS_ARR[@]}"; do
    run_sweep_bin_eps "requant_lowmidscale_tile_gt20k"   "$bin" "$eps" "$nb" "${LOWMIDSCALE_TILE_ENV[@]}"
  done
  for eps in "${EPS_ARR[@]}"; do
    run_sweep_bin_eps "requant_lowscale_tile_gt20k"      "$bin" "$eps" "$nb" "${LOWSCALE_TILE_ENV[@]}"
  done
  if [[ "$(basename "$bin")" == *32k* ]]; then
    for eps in "${EPS_ARR[@]}"; do
      run_sweep_bin_eps "requant_lowscale_e2m1_tile_gt20k"    "$bin" "$eps" "$nb" "${LOWSCALE_E2M1_TILE_ENV[@]}"
    done
    for eps in "${EPS_ARR[@]}"; do
      run_sweep_bin_eps "requant_lowmidscale_e2m1_tile_gt20k" "$bin" "$eps" "$nb" "${LOWMIDSCALE_E2M1_TILE_ENV[@]}"
    done
  fi
  for eps in "${EPS_ARR[@]}"; do
    run_sweep_bin_eps "requant_ladder_scaled_tile_gt20k" "$bin" "$eps" "$nb" "${LADDER_TILE_ENV[@]}"
  done

  echo "[ORCH] <<< finished bin=$(basename $bin)" | tee -a "$ORCH_LOG"
done

echo "[ORCH END] $(date -Is)" | tee -a "$ORCH_LOG"
