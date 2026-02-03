#!/bin/bash
set -u
set -o pipefail

# Sweep settings
MX_MODES=(tile block)
NBS=(32 64 128 256)
FORMATS=(e4m3 e5m2 e3m2 e2m3 e2m1)
if [[ "${PLAIN_FP8_ONLY:-0}" == "1" ]]; then
  FORMATS=(fp8_e4m3 fp8_e5m2)
elif [[ -n "${FORMATS_OVERRIDE:-}" ]]; then
  read -r -a FORMATS <<< "$FORMATS_OVERRIDE"
fi
FP32_BUCKETS=(fp32 mx_fp16)
FP16_BUCKETS=(fp16 mx_fp16)
if [[ -n "${FP32_BUCKETS_OVERRIDE:-}" ]]; then
  read -r -a FP32_BUCKETS <<< "$FP32_BUCKETS_OVERRIDE"
fi
if [[ -n "${FP16_BUCKETS_OVERRIDE:-}" ]]; then
  read -r -a FP16_BUCKETS <<< "$FP16_BUCKETS_OVERRIDE"
fi

LOGS_DIR="/home/abduraa/MX_project/logs"
BUILD=${BUILD:-0}
QUIET=${QUIET:-1}
SWEEP_LOG_DIR="/home/abduraa/MX_project/logs/mx_ooc_data/sweep"
mkdir -p "$SWEEP_LOG_DIR"

if [[ -n "${BIN_LIST:-}" ]]; then
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

echo "Found ${#BIN_FILES[@]} weak-correlation bin files."

ERROR_LOG="$SWEEP_LOG_DIR/sweep_errors.log"
FAILS=0

for BIN_PATH in "${BIN_FILES[@]}"; do
  for MX_MODE in "${MX_MODES[@]}"; do
    for NB in "${NBS[@]}"; do
      for FORMAT in "${FORMATS[@]}"; do
        for FP32_BUCKET in "${FP32_BUCKETS[@]}"; do
          for FP16_BUCKET in "${FP16_BUCKETS[@]}"; do
            echo "=== BIN=$BIN_PATH MX_MODE=$MX_MODE NB=$NB FORMAT=$FORMAT FP32_BUCKET=$FP32_BUCKET FP16_BUCKET=$FP16_BUCKET ==="
            RUN_ENV=(BUILD=$BUILD FP32_BUCKET=$FP32_BUCKET FP16_BUCKET=$FP16_BUCKET \
                     FORMAT=$FORMAT MX_MODE=$MX_MODE NB=$NB BIN_PATH="$BIN_PATH")
            STAMP=$(date +"%Y%m%d_%H%M%S")
            RUN_LOG="$SWEEP_LOG_DIR/sweep_${STAMP}.log"
            if [[ "$QUIET" == "1" ]]; then
              if ! env "${RUN_ENV[@]}" ./build_run.sh > "$RUN_LOG" 2>&1; then
                echo "FAILED: BIN=$BIN_PATH MX_MODE=$MX_MODE NB=$NB FORMAT=$FORMAT FP32_BUCKET=$FP32_BUCKET FP16_BUCKET=$FP16_BUCKET" >> "$ERROR_LOG"
                FAILS=$((FAILS + 1))
              fi
            else
              if ! env "${RUN_ENV[@]}" ./build_run.sh | tee "$RUN_LOG"; then
                echo "FAILED: BIN=$BIN_PATH MX_MODE=$MX_MODE NB=$NB FORMAT=$FORMAT FP32_BUCKET=$FP32_BUCKET FP16_BUCKET=$FP16_BUCKET" >> "$ERROR_LOG"
                FAILS=$((FAILS + 1))
              fi
            fi
          done
        done
      done
    done
  done
done

if [[ $FAILS -gt 0 ]]; then
  echo "Sweep completed with $FAILS failures. See $ERROR_LOG"
else
  echo "Sweep completed successfully."
fi
