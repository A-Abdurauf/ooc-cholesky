#!/usr/bin/env python3
"""Derive a 'ladder-only' view of figures/all_error_runs.csv:

  - drop dropin_mxfp4 rows (only ladder + baseline modes kept)
  - rename granularity 'vec1d' -> 'MX'
  - normalise tile_breakdown format names:
      e2m1     -> MXFP4
      mx_e4m3  -> MXFP8
      mx_e5m2  -> MXFP8
      mx_fp16  -> MXFP16
      mx_fp32  -> FP32
      fp8_e4m3 -> FP8       (plain IEEE FP8)
      fp8_e5m2 -> FP8
      fp64/fp32/fp16  -> uppercase

Same delimiter is preserved (';' for lower-tri breakdowns, ',' for full-square
counts) and the tile_breakdown_kind column is left intact.
"""
import argparse
import csv
from pathlib import Path

INPUT  = "/home/abduraa/MX_project/ooc-cholesky/figures/all_error_runs.csv"
OUTPUT = "/home/abduraa/MX_project/ooc-cholesky/figures/all_error_runs_ladder_only.csv"

# Lower-case key -> public uppercase name.
FORMAT_RENAME = {
    "e2m1":     "MXFP4",
    "mx_e4m3":  "MXFP8",
    "mx_e5m2":  "MXFP8",
    "mx_fp16":  "MXFP16",
    "mx_fp32":  "FP32",
    "fp8_e4m3": "FP8",
    "fp8_e5m2": "FP8",
    "fp16":     "FP16",
    "fp64":     "FP64",
    "fp32":     "FP32",
}

GRANULARITY_RENAME = {"vec1d": "MX"}

DROP_MODES = {"dropin_mxfp4"}


def rename_breakdown(s, kind):
    """Rename format keys inside the tile_breakdown string.  Preserves the
    original separator (';' for lower-tri counts, ',' for full-square counts).
    """
    if not s:
        return s
    s = s.strip().strip('"')
    sep = ";" if kind == "lower" else ","
    parts = [p for p in s.split(sep) if "=" in p]
    out_parts = []
    # Accumulate by canonical name in case two source aliases collapse to the
    # same output (e.g. both fp8_e4m3 and fp8_e5m2 -> FP8).
    agg = {}
    order = []
    for p in parts:
        k, v = p.split("=", 1)
        k_in = k.strip().lower()
        new_k = FORMAT_RENAME.get(k_in, k.strip())
        try:
            count = int(v.strip())
        except ValueError:
            count = 0
        if new_k not in agg:
            agg[new_k] = 0
            order.append(new_k)
        agg[new_k] += count
    out_parts = [f"{k}={agg[k]}" for k in order]
    out = sep.join(out_parts)
    # Mirror the original trailing-separator habit of the lower-tri rows.
    if kind == "lower" and parts and s.endswith(";"):
        out += ";"
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in",  dest="inp", default=INPUT)
    ap.add_argument("--out", default=OUTPUT)
    args = ap.parse_args()

    rows_in = list(csv.DictReader(open(args.inp)))
    cols = list(rows_in[0].keys()) if rows_in else []

    kept = []
    dropped_modes = {}
    for r in rows_in:
        mode = r.get("mode", "")
        if mode in DROP_MODES:
            dropped_modes[mode] = dropped_modes.get(mode, 0) + 1
            continue
        # Rename granularity.
        g = r.get("granularity", "")
        r["granularity"] = GRANULARITY_RENAME.get(g, g)
        # Rename tile_breakdown keys.
        r["tile_breakdown"] = rename_breakdown(
            r.get("tile_breakdown", ""), r.get("tile_breakdown_kind", "lower"))
        kept.append(r)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in kept:
            w.writerow(r)

    print(f"Read     : {len(rows_in)} rows from {args.inp}")
    print(f"Dropped  : {sum(dropped_modes.values())} rows  ({dict(dropped_modes)})")
    print(f"Kept     : {len(kept)} rows -> {out}")


if __name__ == "__main__":
    main()
