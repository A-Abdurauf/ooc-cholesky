#!/usr/bin/env bash
set -euo pipefail

# Sweep settings
EPS_LIST=${EPS_LIST:-"1e-8 1e-7 1e-6 1e-5"}
BIN_LIST=${BIN_LIST:-"/home/abduraa/MX_project/logs/my_cov_weak_8192.bin /home/abduraa/MX_project/logs/my_cov_medium_8192.bin"}
FP32_BUCKETS=${FP32_BUCKETS:-"fp32 mx_fp32 mx_fp16"}
FP16_BUCKETS=${FP16_BUCKETS:-"fp16 mx_fp16 e4m3"}
FORMATS=${FORMATS:-"e4m3 fp8_e4m3"}

# Runtime options
NB=${NB:-256}
CORES=${CORES:-16}
MX_MX_MODE=${MX_MX_MODE:-block}
MX_BLOCK_SUBTILE=${MX_BLOCK_SUBTILE:-16}
OUT_DIR=${OUT_DIR:-"/home/abduraa/MX_project/logs/mx_ooc_data"}

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
BUILD_RUN="$SCRIPT_DIR/build_run.sh"
mkdir -p "$OUT_DIR"

if [[ ! -x "$BUILD_RUN" ]]; then
  echo "[ERROR] Missing script: $BUILD_RUN" >&2
  exit 1
fi

SUMMARY_FILE=${SUMMARY_FILE:-"$OUT_DIR/summary_precision_sweep_8192.csv"}
if [[ ! -f "$SUMMARY_FILE" ]]; then
  echo "bin,format,mx_mode,nb,n,cores,source_epsilon,rel_factor_error,kl_divergence,total_tiles,fp64,fp32,fp16,bf16,mx_fp16,mx_fp32,mx_e4m3,mx_e5m2,fp8_e4m3,fp8_e5m2,e3m2,e2m3,e2m1,unknown,low_other,fp32_bucket,fp16_bucket" >> "$SUMMARY_FILE"
fi

for BIN_PATH in $BIN_LIST; do
  if [[ ! -f "$BIN_PATH" ]]; then
    echo "[WARN] Skipping missing bin: $BIN_PATH" >&2
    continue
  fi
  N_SUMMARY=""
  if [[ "$BIN_PATH" =~ _([0-9]+)\.bin$ ]]; then
    N_SUMMARY="${BASH_REMATCH[1]}"
  fi

  for EPS in $EPS_LIST; do
    for FP32_BUCKET in $FP32_BUCKETS; do
      for FP16_BUCKET in $FP16_BUCKETS; do
        for FORMAT in $FORMATS; do
          base_name=$(basename "$BIN_PATH")
          base_name=${base_name%.bin}
          eps_tag=${EPS//./p}
          eps_tag=${eps_tag//-/m}

          LOG_PATH="$OUT_DIR/run_${base_name}_eps_${eps_tag}_f32_${FP32_BUCKET}_f16_${FP16_BUCKET}_fmt_${FORMAT}.log"
          TILING_CSV="$OUT_DIR/tiling_${base_name}_eps_${eps_tag}_f32_${FP32_BUCKET}_f16_${FP16_BUCKET}_fmt_${FORMAT}.csv"

          MX_SOURCE_EPSILON="$EPS" \
          MX_SKIP_KL=1 \
          BUILD=0 \
          MX_ERROR_LEGACY=0 \
          CORES="$CORES" \
          CUDA_VISIBLE_DEVICES="${CUDA_VISIBLE_DEVICES:-0}" \
          MX_MX_MODE="$MX_MX_MODE" \
          MX_BLOCK_SUBTILE="$MX_BLOCK_SUBTILE" \
          MX_DEBUG_TILE_DUMP=0 \
          FP32_BUCKET="$FP32_BUCKET" \
          FP16_BUCKET="$FP16_BUCKET" \
          FORMAT="$FORMAT" \
          MX_FP32_SCALE_BITS="${MX_FP32_SCALE_BITS:-11}" \
          MX_MODE="$MX_MX_MODE" \
          NB="$NB" \
          BIN_PATH="$BIN_PATH" \
          "$BUILD_RUN" > "$LOG_PATH" 2>&1

          REL_ERR=$(grep -m1 "error:" "$LOG_PATH" | awk '{print $NF}')
          KL_VAL=$(grep -m1 "kl_divergence:" "$LOG_PATH" | awk '{print $NF}')

          COUNTS_LINE=$(python3 - <<PY
import re
import csv
from collections import Counter
from pathlib import Path

log_path = Path(r"""${LOG_PATH}""")
csv_path = Path(r"""${TILING_CSV}""")

nt = None
tiles = {}

re_nt = re.compile(r"\[DEBUG\] nb used: .*?, nt: (\d+)")
re_dtype = re.compile(r"\[TILE_DTYPE\] \((\d+), (\d+)\) (\S+)")
re_target = re.compile(r"\[TILE_TARGET\] \((\d+), (\d+)\) (\S+)")

for line in log_path.read_text().splitlines():
    if nt is None:
        m = re_nt.search(line)
        if m:
            nt = int(m.group(1))
            continue
    m = re_target.search(line)
    if m:
        r, c = int(m.group(1)), int(m.group(2))
        tiles[(r, c)] = m.group(3)
        continue
    m = re_dtype.search(line)
    if m:
        r, c = int(m.group(1)), int(m.group(2))
        tiles.setdefault((r, c), m.group(3))
        continue

if nt is None:
    nt = 0

csv_path.parent.mkdir(parents=True, exist_ok=True)
with csv_path.open("w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["row", "col", "dtype"])
    for r in range(nt):
        for c in range(r + 1):
            dtype = tiles.get((r, c), "unknown")
            writer.writerow([r, c, dtype])

rows = list(csv.DictReader(csv_path.open()))
types = [r["dtype"].lower() for r in rows]
cnt = Counter(types)

def get(*keys):
    return sum(cnt.get(k, 0) for k in keys)

total = len(types)
fp64 = get("fp64")
fp32 = get("fp32")
fp16 = get("fp16")
bf16 = get("bf16")
mx_fp16 = get("mx_fp16")
mx_fp32 = get("mx_fp32", "mx_f32")
mx_e4m3 = get("mx_e4m3", "e4m3")
mx_e5m2 = get("mx_e5m2", "e5m2")
fp8_e4m3 = get("fp8_e4m3", "fp8e4m3")
fp8_e5m2 = get("fp8_e5m2", "fp8e5m2")
e3m2 = get("e3m2")
e2m3 = get("e2m3")
e2m1 = get("e2m1")
unknown = get("unknown")
low_other = total - (fp64 + fp32 + fp16 + bf16 + mx_fp16 + mx_fp32 + mx_e4m3 + mx_e5m2 + fp8_e4m3 + fp8_e5m2 + e3m2 + e2m3 + e2m1 + unknown)

print(",".join(map(str, [
    total, fp64, fp32, fp16, bf16, mx_fp16, mx_fp32,
    mx_e4m3, mx_e5m2, fp8_e4m3, fp8_e5m2, e3m2, e2m3, e2m1, unknown, low_other
])))
PY
      )

          if [[ -n "$REL_ERR" ]]; then
            echo "${BIN_PATH},${FORMAT},${MX_MX_MODE},${NB},${N_SUMMARY},${CORES},${EPS},${REL_ERR},${KL_VAL},${COUNTS_LINE},${FP32_BUCKET},${FP16_BUCKET}" >> "$SUMMARY_FILE"
          fi
        done
      done
    done
  done
done
