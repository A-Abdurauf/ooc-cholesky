#!/bin/bash
set -u
set -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BUILD_RUN="$SCRIPT_DIR/build_run.sh"

# Sweep settings
MX_MODES=(tile block)
NBS=(32 64 128 256)
if [[ -n "${NBS_OVERRIDE:-}" ]]; then
  read -r -a NBS <<< "$NBS_OVERRIDE"
elif [[ -n "${NB_MAX:-}" ]]; then
  FILTERED=()
  for nb in "${NBS[@]}"; do
    if [[ "$nb" -le "$NB_MAX" ]]; then
      FILTERED+=("$nb")
    fi
  done
  if [[ ${#FILTERED[@]} -gt 0 ]]; then
    NBS=("${FILTERED[@]}")
  fi
fi
FORMATS=(e4m3 e5m2 e3m2 e2m3 e2m1)
if [[ "${PLAIN_FP8_ONLY:-0}" == "1" ]]; then
  FORMATS=(fp8_e4m3 fp8_e5m2)
elif [[ -n "${FORMATS_OVERRIDE:-}" ]]; then
  read -r -a FORMATS <<< "$FORMATS_OVERRIDE"
fi
FP32_BUCKETS=(fp32 mx_fp16)
FP16_BUCKETS=(fp16 mx_fp16)
EPS_LIST=(1e-8)
if [[ -n "${EPS_LIST_OVERRIDE:-}" ]]; then
  read -r -a EPS_LIST <<< "$EPS_LIST_OVERRIDE"
fi
if [[ -n "${FP32_BUCKETS_OVERRIDE:-}" ]]; then
  read -r -a FP32_BUCKETS <<< "$FP32_BUCKETS_OVERRIDE"
fi
if [[ -n "${FP16_BUCKETS_OVERRIDE:-}" ]]; then
  read -r -a FP16_BUCKETS <<< "$FP16_BUCKETS_OVERRIDE"
fi
SUBTILES=()
if [[ -n "${MX_BLOCK_SUBTILE_SIZES:-}" ]]; then
  read -r -a SUBTILES <<< "$MX_BLOCK_SUBTILE_SIZES"
fi
AUTO_SUBTILE_BY_NB=${AUTO_SUBTILE_BY_NB:-1}

calc_auto_subtile() {
  local nb="$1"
  local st=$((nb / 8))
  if [[ "$st" -lt 1 ]]; then
    st=1
  fi
  echo "$st"
}

LOGS_DIR="/home/abduraa/MX_project/logs"
BUILD=${BUILD:-0}
QUIET=${QUIET:-1}
SWEEP_LOG_DIR="/home/abduraa/MX_project/logs/mx_ooc_data/sweep"
mkdir -p "$SWEEP_LOG_DIR"
SUMMARY_FILE_PATH=${SUMMARY_FILE:-"/home/abduraa/MX_project/logs/mx_ooc_data/summary_rel_error.txt"}
SKIP_EXISTING=${SKIP_EXISTING:-0}
PARALLEL_JOBS=${PARALLEL_JOBS:-1}
if ! [[ "$PARALLEL_JOBS" =~ ^[0-9]+$ ]] || [[ "$PARALLEL_JOBS" -lt 1 ]]; then
  PARALLEL_JOBS=1
fi

if [[ ! -x "$BUILD_RUN" ]]; then
  echo "Missing executable build script: $BUILD_RUN"
  exit 1
fi

# Large weak sweep controls (requested mapping):
#   20k -> NB 2048
#   40k/65k -> NB 5120
#   >65k -> NB 10240
LARGE_WEAK_ONLY=${LARGE_WEAK_ONLY:-1}
AUTO_NB_BY_SIZE=${AUTO_NB_BY_SIZE:-1}
TARGET_WEAK_SIZES=${TARGET_WEAK_SIZES:-"20k 40k 65k 80k 100k"}

size_token_to_n() {
  local tok="${1,,}"
  if [[ "$tok" =~ ^([0-9]+)k$ ]]; then
    echo $((BASH_REMATCH[1] * 1000))
    return
  fi
  if [[ "$tok" =~ ^[0-9]+$ ]]; then
    echo "$tok"
    return
  fi
  echo ""
}

extract_n_from_bin() {
  local p="$1"
  local b
  b=$(basename "$p")
  local bl="${b,,}"
  if [[ "$bl" =~ weak_([0-9]+)k ]]; then
    echo $((BASH_REMATCH[1] * 1000))
    return
  fi
  if [[ "$bl" =~ weak_([0-9]+) ]]; then
    echo "${BASH_REMATCH[1]}"
    return
  fi
  echo ""
}

is_target_weak_size() {
  local n="$1"
  local t tn
  for t in $TARGET_WEAK_SIZES; do
    tn=$(size_token_to_n "$t")
    if [[ -n "$tn" && "$n" == "$tn" ]]; then
      return 0
    fi
  done
  return 1
}

select_nb_for_n() {
  local n="$1"
  if [[ "$n" -eq 20000 ]]; then
    echo 2048
  elif [[ "$n" -le 65000 ]]; then
    echo 5120
  else
    echo 10240
  fi
}

if [[ -n "${BIN_PATH:-}" ]]; then
  BIN_FILES=("$BIN_PATH")
elif [[ -n "${BIN_OVERRIDE:-}" ]]; then
  BIN_FILES=("$BIN_OVERRIDE")
elif [[ -n "${BIN_LIST:-}" ]]; then
  read -r -a BIN_FILES <<< "$BIN_LIST"
else
  mapfile -t BIN_FILES < <(
    find "$LOGS_DIR" -maxdepth 2 -type f -name "*weak*.bin" -printf "%f\t%p\n" \
      | sort -t$'\t' -k1,1V \
      | cut -f2-
  )
fi
if [[ -n "${BIN_FILTER_REGEX:-}" ]]; then
  mapfile -t BIN_FILES < <(printf '%s\n' "${BIN_FILES[@]}" | grep -E "$BIN_FILTER_REGEX")
fi
if [[ ${#BIN_FILES[@]} -eq 0 ]]; then
  echo "No weak-correlation bins found in $LOGS_DIR"
  exit 1
fi

if [[ "$LARGE_WEAK_ONLY" == "1" ]]; then
  FILTERED_BINS=()
  for _bin in "${BIN_FILES[@]}"; do
    _n=$(extract_n_from_bin "$_bin")
    if [[ -n "$_n" ]] && is_target_weak_size "$_n"; then
      FILTERED_BINS+=("$_bin")
    fi
  done
  BIN_FILES=("${FILTERED_BINS[@]}")
fi

if [[ ${#BIN_FILES[@]} -eq 0 ]]; then
  echo "No matching weak bins for TARGET_WEAK_SIZES='$TARGET_WEAK_SIZES'"
  exit 1
fi

echo "Found ${#BIN_FILES[@]} weak-correlation bin files."

ERROR_LOG="$SWEEP_LOG_DIR/sweep_errors.log"
: > "$ERROR_LOG"
SKIP_IF_NO_LOW=${SKIP_IF_NO_LOW:-0}

for BIN_PATH in "${BIN_FILES[@]}"; do
  BIN_N=$(extract_n_from_bin "$BIN_PATH")
  NB_RUNS=("${NBS[@]}")
  if [[ "$AUTO_NB_BY_SIZE" == "1" && -n "$BIN_N" && "$BIN_N" -ge 20000 ]]; then
    MAPPED_NB=$(select_nb_for_n "$BIN_N")
    NB_RUNS=("$MAPPED_NB")
    echo "[NB_MAP] $(basename "$BIN_PATH") -> n=$BIN_N, nb=$MAPPED_NB"
  fi
  for MX_MODE in "${MX_MODES[@]}"; do
    for NB in "${NB_RUNS[@]}"; do
      for EPS in "${EPS_LIST[@]}"; do
        for FORMAT in "${FORMATS[@]}"; do
          for FP32_BUCKET in "${FP32_BUCKETS[@]}"; do
            for FP16_BUCKET in "${FP16_BUCKETS[@]}"; do
              SKIP_REST_FORMATS=0
              SUBTILE_RUNS=(0)
              if [[ "$MX_MODE" == "block" ]]; then
                if [[ ${#SUBTILES[@]} -gt 0 ]]; then
                  SUBTILE_RUNS=(0 "${SUBTILES[@]}")
                elif [[ "$AUTO_SUBTILE_BY_NB" == "1" ]]; then
                  AUTO_SUBTILE=$(calc_auto_subtile "$NB")
                  SUBTILE_RUNS=(0 "$AUTO_SUBTILE")
                  echo "[SUBTILE_AUTO] nb=$NB -> subtile=$AUTO_SUBTILE"
                fi
              fi
              for SUBTILE in "${SUBTILE_RUNS[@]}"; do
                MODE_LABEL=$MX_MODE
                if [[ "$MX_MODE" == "block" && "$SUBTILE" -gt 0 ]]; then
                  MODE_LABEL="subtile_${SUBTILE}"
                fi
                if [[ "$SKIP_EXISTING" == "1" && -f "$SUMMARY_FILE_PATH" ]]; then
                  if awk -F'\t' -v bin="$BIN_PATH" -v fmt="$FORMAT" -v mode="$MODE_LABEL" -v nb="$NB" -v eps="$EPS" -v fp32="$FP32_BUCKET" -v fp16="$FP16_BUCKET" 'NR>1 && $1==bin && $2==fmt && $3==mode && $4==nb && $7==eps && $16==fp32 && $17==fp16 {found=1; exit} END{exit found?0:1}' "$SUMMARY_FILE_PATH"; then
                    echo "[SKIP_EXISTING] BIN=$BIN_PATH MX_MODE=$MODE_LABEL NB=$NB EPS=$EPS FORMAT=$FORMAT FP32_BUCKET=$FP32_BUCKET FP16_BUCKET=$FP16_BUCKET"
                    continue
                  fi
                fi
                echo "=== BIN=$BIN_PATH MX_MODE=$MODE_LABEL NB=$NB EPS=$EPS FORMAT=$FORMAT FP32_BUCKET=$FP32_BUCKET FP16_BUCKET=$FP16_BUCKET ==="
                RUN_ENV=(BUILD=$BUILD FP32_BUCKET=$FP32_BUCKET FP16_BUCKET=$FP16_BUCKET \
                         FORMAT=$FORMAT MX_MODE=$MX_MODE NB=$NB BIN_PATH="$BIN_PATH" MX_SOURCE_EPSILON=$EPS)
                if [[ -n "${CORES:-}" ]]; then
                  RUN_ENV+=(CORES=$CORES)
                  RUN_ENV+=(OMP_NUM_THREADS=$CORES MKL_NUM_THREADS=$CORES OPENBLAS_NUM_THREADS=$CORES BLIS_NUM_THREADS=$CORES)
                fi
                if [[ "$MX_MODE" == "block" && "$SUBTILE" -gt 0 ]]; then
                  RUN_ENV+=(MX_BLOCK_SUBTILE=$SUBTILE)
                else
                  RUN_ENV+=(MX_BLOCK_SUBTILE=0)
                fi
                STAMP=$(date +"%Y%m%d_%H%M%S_%N")
                RUN_LOG="$SWEEP_LOG_DIR/sweep_${STAMP}.log"
                if [[ "$PARALLEL_JOBS" -gt 1 ]]; then
                  while [[ $(jobs -rp | wc -l) -ge "$PARALLEL_JOBS" ]]; do
                    wait -n
                  done
                  (
                    if [[ "$QUIET" == "1" ]]; then
                      if ! (cd "$SCRIPT_DIR" && env "${RUN_ENV[@]}" "$BUILD_RUN") > "$RUN_LOG" 2>&1; then
                        echo "FAILED: BIN=$BIN_PATH MX_MODE=$MODE_LABEL NB=$NB EPS=$EPS FORMAT=$FORMAT FP32_BUCKET=$FP32_BUCKET FP16_BUCKET=$FP16_BUCKET" >> "$ERROR_LOG"
                      fi
                    else
                      if ! (cd "$SCRIPT_DIR" && env "${RUN_ENV[@]}" "$BUILD_RUN") | tee "$RUN_LOG"; then
                        echo "FAILED: BIN=$BIN_PATH MX_MODE=$MODE_LABEL NB=$NB EPS=$EPS FORMAT=$FORMAT FP32_BUCKET=$FP32_BUCKET FP16_BUCKET=$FP16_BUCKET" >> "$ERROR_LOG"
                      fi
                    fi
                  ) &
                else
                  if [[ "$QUIET" == "1" ]]; then
                    if ! (cd "$SCRIPT_DIR" && env "${RUN_ENV[@]}" "$BUILD_RUN") > "$RUN_LOG" 2>&1; then
                      echo "FAILED: BIN=$BIN_PATH MX_MODE=$MODE_LABEL NB=$NB EPS=$EPS FORMAT=$FORMAT FP32_BUCKET=$FP32_BUCKET FP16_BUCKET=$FP16_BUCKET" >> "$ERROR_LOG"
                    fi
                  else
                    if ! (cd "$SCRIPT_DIR" && env "${RUN_ENV[@]}" "$BUILD_RUN") | tee "$RUN_LOG"; then
                      echo "FAILED: BIN=$BIN_PATH MX_MODE=$MODE_LABEL NB=$NB EPS=$EPS FORMAT=$FORMAT FP32_BUCKET=$FP32_BUCKET FP16_BUCKET=$FP16_BUCKET" >> "$ERROR_LOG"
                    fi
                  fi
                fi
                if [[ "$SKIP_IF_NO_LOW" == "1" ]]; then
                  LOW_CNT=$(grep -m1 "^\[SUMMARY_COUNTS\] [0-9]" "$RUN_LOG" | awk -F'\t' '{print $NF}')
                  if [[ -n "$LOW_CNT" && "$LOW_CNT" -eq 0 ]]; then
                    SKIP_REST_FORMATS=1
                  fi
                fi
              done
              if [[ "$SKIP_REST_FORMATS" == "1" ]]; then
                break
              fi
            done
          done
        done
      done
    done
  done
done

if [[ "$PARALLEL_JOBS" -gt 1 ]]; then
  wait
fi

FAILS=0
if [[ -f "$ERROR_LOG" ]]; then
  FAILS=$(grep -c '^FAILED:' "$ERROR_LOG" || true)
fi

if [[ $FAILS -gt 0 ]]; then
  echo "Sweep completed with $FAILS failures. See $ERROR_LOG"
else
  echo "Sweep completed successfully."
fi
