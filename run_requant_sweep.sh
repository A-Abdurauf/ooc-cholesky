#!/usr/bin/env bash
# Runs the precision sweep with per-step MX re-quantization enabled.
# Results go to a separate directory so they can be compared against
# the baseline (no re-quantization) runs side-by-side.
set -euo pipefail

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

REQUANT_DIR="/home/abduraa/MX_project/logs/mx_ooc_data_requant"
mkdir -p "$REQUANT_DIR"

OUT_DIR="$REQUANT_DIR" \
SUMMARY_FILE="$REQUANT_DIR/summary_precision_sweep_8192_requant.csv" \
  "$SCRIPT_DIR/run_precision_sweep_8192.sh"

echo ""
echo "Re-quantized sweep complete."
echo "Logs  : $REQUANT_DIR"
echo "Summary: $REQUANT_DIR/summary_precision_sweep_8192_requant.csv"
