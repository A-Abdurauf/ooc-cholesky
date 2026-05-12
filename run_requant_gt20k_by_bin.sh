#!/usr/bin/env bash
# By-bin orchestrator: runs (bin × sweep × eps) so 32k completes for ALL sweeps
# before any 40k or 65k run starts. Already-completed (bin, eps) entries (detected
# in the per-sweep summary TSV) are skipped, so this resumes prior progress.
#
# Sweep order within a bin:
#   legacy_scaled_block128 -> legacy_scaled_vec1d32_FP32tile (32k only) ->
#   ladder_scaled_vec1d32  -> ladder_scaled_block128 ->
#   lowscale_block128      -> lowscale_vec1d32 ->
#   lowmidscale_block128   -> lowmidscale_vec1d32.
# "_FP32tile" suffix on the legacy vec1d sweep flags that MX_FP32 silently
# degrades to a single tile-wide shared scale in vec1d mode (apply_mx_quant_fp64
# lacks a vec1d branch). Once apply_mx_quant_fp64 is patched, re-run under the
# clean name "requant_legacy_scaled_vec1d32_gt20k".
# The two "lowscale" sweeps keep upper buckets as plain IEEE fp32/fp16 and only
# the low tier (MX_E4M3) is shared-scale, so they sidestep the MX_FP32 issue.
# Baseline FP8 is fully done already and is NOT re-run here.
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
ORCH_LOG="$OUT_DIR/requant_gt20k_by_bin.log"

nb_for_bin() {
  local bin="$1"
  case "$(basename "$bin")" in
    *65k*)  echo 4096 ;;
    *40k*)  echo 4096 ;;
    *32k*)  echo 2048 ;;
    *)      echo "${NB_DEFAULT:-2048}" ;;
  esac
}

# Already-done check: returns 0 if (bin, eps) row exists in $sum.
already_done() {
  local sum="$1" bin="$2" eps="$3"
  [[ -f "$sum" ]] || return 1
  awk -F'\t' -v b="$bin" -v e="$eps" 'NR>1 && $1==b && $7==e {found=1; exit} END{exit !found}' "$sum"
}

# Run one (bin, eps) for a given sweep config.
# Args: sweep_name eps env_string... (env vars passed inline as KEY=VALUE pairs).
run_sweep_bin_eps() {
  local sweep_name="$1"; shift
  local bin="$1"; shift
  local eps="$1"; shift
  local nb="$1"; shift
  # remaining args are KEY=VALUE pairs
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

# Sweep config blocks (env vars per sweep).
LEGACY_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_MX_MODE=block
  MX_BLOCK_SUBTILE=128
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=block
  MX_ERROR_LEGACY=0
  MX_BUCKET_FP32=mx_fp32
  MX_BUCKET_FP16=mx_fp16
  FORMAT=mx_e4m3
)
LEGACY_VEC1D32_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_MX_MODE=vec1d
  MX_BLOCK_SUBTILE=32
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=vec1d
  MX_ERROR_LEGACY=0
  MX_BUCKET_FP32=mx_fp32
  MX_BUCKET_FP16=mx_fp16
  FORMAT=mx_e4m3
)
# IEEE-only ladder: bound-based selection, no MX scaling anywhere.
# Hardcoded ladder steps in dpotrf_mixed_precision.cpp:1592:
#   fp8_e4m3 -> fp16 -> fp32 -> fp64 fallback
# Second reference baseline: same plain-IEEE buckets as baseline_fp8 but
# selection uses the bound formula instead of legacy epsilon-cutoff.
LADDER_IEEE_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_SELECTION_CRITERIA=bound
  MX_BOUND_DEBUG=1
  MX_BOUND_LADDER=ieee_only
  MX_MX_MODE=tile
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=tile
  MX_ERROR_LEGACY=0
)
LADDER_VEC1D32_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_SELECTION_CRITERIA=bound
  MX_BOUND_DEBUG=1
  MX_BOUND_LADDER=full
  MX_MX_MODE=vec1d
  MX_BLOCK_SUBTILE=32
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=1
  MX_FP32_SCALE_BITS=11
  MX_MODE=vec1d
  MX_ERROR_LEGACY=0
)
LADDER_BLOCK128_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_SELECTION_CRITERIA=bound
  MX_BOUND_DEBUG=1
  MX_BOUND_LADDER=full
  MX_MX_MODE=block
  MX_BLOCK_SUBTILE=128
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=1
  MX_FP32_SCALE_BITS=11
  MX_MODE=block
  MX_ERROR_LEGACY=0
)
# Lower-tier-only scaled: low (MX_E4M3) gets shared-scale; upper buckets
# are plain IEEE fp32/fp16 (no MX_FP32, no MX_FP16). Selection: bound + legacy.
LOWSCALE_BLOCK128_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_MX_MODE=block
  MX_BLOCK_SUBTILE=128
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=block
  MX_ERROR_LEGACY=0
  MX_BUCKET_FP32=fp32
  MX_BUCKET_FP16=fp16
  FORMAT=mx_e4m3
)
LOWSCALE_VEC1D32_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_MX_MODE=vec1d
  MX_BLOCK_SUBTILE=32
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=vec1d
  MX_ERROR_LEGACY=0
  MX_BUCKET_FP32=fp32
  MX_BUCKET_FP16=fp16
  FORMAT=mx_e4m3
)
# Same as lowscale_* but with Scaled FP4 (E2M1) as the low tier.
LOWSCALE_E2M1_BLOCK128_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_MX_MODE=block
  MX_BLOCK_SUBTILE=128
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=block
  MX_ERROR_LEGACY=0
  MX_BUCKET_FP32=fp32
  MX_BUCKET_FP16=fp16
  FORMAT=e2m1
)
LOWSCALE_E2M1_VEC1D32_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_MX_MODE=vec1d
  MX_BLOCK_SUBTILE=32
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=vec1d
  MX_ERROR_LEGACY=0
  MX_BUCKET_FP32=fp32
  MX_BUCKET_FP16=fp16
  FORMAT=e2m1
)
# Same as lowmidscale_* but with Scaled FP4 (E2M1) as the low tier.
# Low = SFP4, Mid = SFP16, FP32 plain.
LOWMIDSCALE_E2M1_BLOCK128_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_MX_MODE=block
  MX_BLOCK_SUBTILE=128
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=block
  MX_ERROR_LEGACY=0
  MX_BUCKET_FP32=fp32
  MX_BUCKET_FP16=mx_fp16
  FORMAT=e2m1
)
LOWMIDSCALE_E2M1_VEC1D32_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_MX_MODE=vec1d
  MX_BLOCK_SUBTILE=32
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=vec1d
  MX_ERROR_LEGACY=0
  MX_BUCKET_FP32=fp32
  MX_BUCKET_FP16=mx_fp16
  FORMAT=e2m1
)
# Low + FP16 tiers scaled: MX_E4M3 (low) + MX_FP16 (mid) get shared-scale;
# FP32 stays as plain IEEE (no MX_FP32, so no apply_mx_quant_fp64 vec1d issue).
LOWMIDSCALE_BLOCK128_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_MX_MODE=block
  MX_BLOCK_SUBTILE=128
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=block
  MX_ERROR_LEGACY=0
  MX_BUCKET_FP32=fp32
  MX_BUCKET_FP16=mx_fp16
  FORMAT=mx_e4m3
)
LOWMIDSCALE_VEC1D32_ENV=(
  MX_SKIP_KL=1
  MX_FP8_SUBNORMAL=1
  MX_MX_MODE=vec1d
  MX_BLOCK_SUBTILE=32
  MX_UNDERFLOW_MODE=fz
  MX_SCALE_AWARE_EPSILON=0
  MX_FP32_SCALE_BITS=11
  MX_MODE=vec1d
  MX_ERROR_LEGACY=0
  MX_BUCKET_FP32=fp32
  MX_BUCKET_FP16=mx_fp16
  FORMAT=mx_e4m3
)

{
  echo "[ORCH START] $(date -Is)"
  echo "[ORCH] BIN_LIST=$BIN_LIST"
  echo "[ORCH] EPS_LIST=$EPS_LIST"
  echo "[ORCH] order: bin -> sweep -> eps (32k full first, then 40k, then 65k)"
} | tee -a "$ORCH_LOG"

for bin in $BIN_LIST; do
  if [[ ! -f "$bin" ]]; then
    echo "[WARN] missing bin: $bin" | tee -a "$ORCH_LOG"
    continue
  fi
  nb=$(nb_for_bin "$bin")
  echo "[ORCH] >>> entering bin=$(basename $bin) nb=$nb" | tee -a "$ORCH_LOG"

  for eps in "${EPS_ARR[@]}"; do
    run_sweep_bin_eps "requant_legacy_scaled_block128_gt20k" "$bin" "$eps" "$nb" "${LEGACY_ENV[@]}"
  done
  # apply_mx_quant_fp64 now has a proper vec1d branch (block_mode || vec1d_mode),
  # so legacy_scaled_vec1d32 runs at full quality on all bins.
  for eps in "${EPS_ARR[@]}"; do
    run_sweep_bin_eps "requant_legacy_scaled_vec1d32_gt20k"  "$bin" "$eps" "$nb" "${LEGACY_VEC1D32_ENV[@]}"
  done
  for eps in "${EPS_ARR[@]}"; do
    run_sweep_bin_eps "requant_ladder_scaled_vec1d32_gt20k"  "$bin" "$eps" "$nb" "${LADDER_VEC1D32_ENV[@]}"
  done
  for eps in "${EPS_ARR[@]}"; do
    run_sweep_bin_eps "requant_ladder_scaled_block128_gt20k" "$bin" "$eps" "$nb" "${LADDER_BLOCK128_ENV[@]}"
  done
  for eps in "${EPS_ARR[@]}"; do
    run_sweep_bin_eps "requant_lowscale_block128_gt20k"      "$bin" "$eps" "$nb" "${LOWSCALE_BLOCK128_ENV[@]}"
  done
  for eps in "${EPS_ARR[@]}"; do
    run_sweep_bin_eps "requant_lowscale_vec1d32_gt20k"       "$bin" "$eps" "$nb" "${LOWSCALE_VEC1D32_ENV[@]}"
  done
  for eps in "${EPS_ARR[@]}"; do
    run_sweep_bin_eps "requant_lowmidscale_block128_gt20k"   "$bin" "$eps" "$nb" "${LOWMIDSCALE_BLOCK128_ENV[@]}"
  done
  for eps in "${EPS_ARR[@]}"; do
    run_sweep_bin_eps "requant_lowmidscale_vec1d32_gt20k"    "$bin" "$eps" "$nb" "${LOWMIDSCALE_VEC1D32_ENV[@]}"
  done
  # SFP4 (E2M1) low-tier drop-in variants — 32k probe only.
  if [[ "$(basename "$bin")" == *32k* ]]; then
    for eps in "${EPS_ARR[@]}"; do
      run_sweep_bin_eps "requant_lowscale_e2m1_block128_gt20k"    "$bin" "$eps" "$nb" "${LOWSCALE_E2M1_BLOCK128_ENV[@]}"
    done
    for eps in "${EPS_ARR[@]}"; do
      run_sweep_bin_eps "requant_lowscale_e2m1_vec1d32_gt20k"     "$bin" "$eps" "$nb" "${LOWSCALE_E2M1_VEC1D32_ENV[@]}"
    done
    for eps in "${EPS_ARR[@]}"; do
      run_sweep_bin_eps "requant_lowmidscale_e2m1_block128_gt20k" "$bin" "$eps" "$nb" "${LOWMIDSCALE_E2M1_BLOCK128_ENV[@]}"
    done
    for eps in "${EPS_ARR[@]}"; do
      run_sweep_bin_eps "requant_lowmidscale_e2m1_vec1d32_gt20k"  "$bin" "$eps" "$nb" "${LOWMIDSCALE_E2M1_VEC1D32_ENV[@]}"
    done
  fi
  # IEEE-only ladder (second baseline): bound-based selection, plain IEEE storage.
  for eps in "${EPS_ARR[@]}"; do
    run_sweep_bin_eps "requant_ladder_ieee_gt20k"            "$bin" "$eps" "$nb" "${LADDER_IEEE_ENV[@]}"
  done

  echo "[ORCH] <<< finished bin=$(basename $bin)" | tee -a "$ORCH_LOG"
done

echo "[ORCH END] $(date -Is)" | tee -a "$ORCH_LOG"
