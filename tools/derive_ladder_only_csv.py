#!/usr/bin/env python3
"""Derive a 'ladder-only' view of figures/all_error_runs.csv:

  - drop dropin_mxfp4 rows (only ladder + baseline modes kept)
  - rename mode column to human-readable labels:
      baseline             -> Baseline
      ladder_ieee          -> Ladder IEEE
      ladder_mx_staircase  -> Ladder MX (native)     (no MXFP16 rung)
      ladder_mx_full       -> Ladder MX+MXFP16       (with MXFP16 rung)
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
  - convert tile counts ALL to lower-triangular (halve off-diagonals, keep
    M FP64 diagonal tiles at full count).  After this every row reports
    lower-tri counts and the tile_breakdown_kind column is dropped.
  - order tile counts from highest precision DOWNWARD:
      FP64 -> FP32 -> FP16 -> MXFP16 -> FP8 -> MXFP8 -> MXFP4
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

MODE_RENAME = {
    "baseline":            "Baseline",
    "ladder_ieee":         "Ladder IEEE",
    "ladder_mx_staircase": "Ladder MX (native)",
    "ladder_mx_full":      "Ladder MX+MXFP16",
}

# Highest precision first, going downward.
FORMAT_ORDER = ["FP64", "FP32", "FP16", "MXFP16", "FP8", "MXFP8", "MXFP4"]

DROP_MODES = {"dropin_mxfp4"}


def to_lower_tri(counts, M):
    """Convert full-square tile counts to lower-tri counts.
    M FP64 diagonal tiles + (FP64_count - M)//2 off-diag FP64 tiles;
    every other format is fully off-diagonal so halve."""
    out = {}
    for fmt, v in counts.items():
        if fmt == "FP64":
            off = max(v - M, 0)
            out[fmt] = M + off // 2
        else:
            out[fmt] = v // 2
    return out


def format_breakdown(counts):
    """Render counts in canonical FORMAT_ORDER, semicolon-separated, lower-tri."""
    parts = []
    seen = set()
    for fmt in FORMAT_ORDER:
        if fmt in counts and counts[fmt] > 0:
            parts.append(f"{fmt}={counts[fmt]}")
            seen.add(fmt)
    # Append anything unknown last so we never silently drop counts.
    for fmt, n in counts.items():
        if fmt not in seen and n > 0:
            parts.append(f"{fmt}={n}")
    return ";".join(parts) + (";" if parts else "")


def parse_and_rename(s, kind):
    """Parse a tile_breakdown string, rename format keys, return dict."""
    if not s:
        return {}
    s = s.strip().strip('"')
    sep = ";" if kind == "lower" else ","
    parts = [p for p in s.split(sep) if "=" in p]
    out = {}
    for p in parts:
        k, v = p.split("=", 1)
        new_k = FORMAT_RENAME.get(k.strip().lower(), k.strip())
        try:
            count = int(v.strip())
        except ValueError:
            count = 0
        out[new_k] = out.get(new_k, 0) + count
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in",  dest="inp", default=INPUT)
    ap.add_argument("--out", default=OUTPUT)
    args = ap.parse_args()

    rows_in = list(csv.DictReader(open(args.inp)))
    # Output columns: drop tile_breakdown_kind (everything is now lower-tri).
    cols = [c for c in rows_in[0].keys() if c != "tile_breakdown_kind"] if rows_in else []

    kept = []
    dropped_modes = {}
    for r in rows_in:
        mode = r.get("mode", "")
        if mode in DROP_MODES:
            dropped_modes[mode] = dropped_modes.get(mode, 0) + 1
            continue
        # Rename mode + granularity.
        r["mode"]        = MODE_RENAME.get(mode, mode)
        g = r.get("granularity", "")
        r["granularity"] = GRANULARITY_RENAME.get(g, g)
        # Parse + rename tile-format keys.
        counts = parse_and_rename(r.get("tile_breakdown", ""),
                                  r.get("tile_breakdown_kind", "lower"))
        # Convert to lower-tri if needed.
        if r.get("tile_breakdown_kind", "lower") == "full":
            try:
                N  = int(r["N"]); nb = int(r["nb"])
                counts = to_lower_tri(counts, N // nb)
            except (TypeError, ValueError, KeyError):
                pass
        # Re-emit in canonical FP64-first order, semicolon-separated.
        r["tile_breakdown"] = format_breakdown(counts)
        r.pop("tile_breakdown_kind", None)
        kept.append(r)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in kept:
            w.writerow({c: r.get(c, "") for c in cols})

    print(f"Read     : {len(rows_in)} rows from {args.inp}")
    print(f"Dropped  : {sum(dropped_modes.values())} rows  ({dict(dropped_modes)})")
    print(f"Kept     : {len(kept)} rows -> {out}")
    print(f"Schema   : {cols}")


if __name__ == "__main__":
    main()
