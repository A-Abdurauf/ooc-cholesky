#!/usr/bin/env python3
"""Plot mixed-precision tile map with FP64 diagonal enforced.

Input: CSV with columns: row,col,dtype (dtype values: fp64, fp32, mx_e4m3, bf16, fp16, unknown)
Output: PNG heatmap.
"""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--csv", required=True, help="Path to tiling CSV")
    p.add_argument("--out", required=True, help="Output PNG path")
    p.add_argument("--nt", type=int, default=None, help="Tile grid size (optional)")
    p.add_argument("--bin", default="", help="Input bin path for metadata")
    p.add_argument("--n", default="", help="Matrix size")
    p.add_argument("--nb", default="", help="Tile size")
    p.add_argument("--rel-err", default="", help="Relative factorization error")
    p.add_argument("--format", default="", help="Format name")
    return p.parse_args()


def main():
    args = parse_args()
    csv_path = Path(args.csv)
    out_path = Path(args.out)

    tiles = {}
    max_idx = -1
    with csv_path.open() as f:
        reader = csv.DictReader(f)
        for row in reader:
            r = int(row["row"])
            c = int(row["col"])
            dtype = row.get("dtype", "unknown")
            tiles[(r, c)] = dtype
            max_idx = max(max_idx, r, c)

    nt = args.nt if args.nt is not None else (max_idx + 1 if max_idx >= 0 else 0)

    dtype_order = [
        "fp64",
        "fp32",
        "fp16",
        "bf16",
        "mx_fp16",
        "mx_e4m3",
        "mx_e5m2",
        "e3m2",
        "e2m3",
        "e2m1",
        "unknown",
    ]
    color_map = {
        "fp64": "#1f77b4",
        "fp32": "#ff7f0e",
        "fp16": "#2ca02c",
        "bf16": "#d62728",
        "mx_fp16": "#8c564b",
        "mx_e4m3": "#9467bd",
        "mx_e5m2": "#e377c2",
        "e3m2": "#bcbd22",
        "e2m3": "#17becf",
        "e2m1": "#aec7e8",
        "unknown": "#7f7f7f",
    }
    index_map = {dtype: idx for idx, dtype in enumerate(dtype_order)}
    cmap = ListedColormap([color_map[dtype] for dtype in dtype_order])

    grid = [[index_map["unknown"]] * nt for _ in range(nt)]
    for r in range(nt):
        grid[r][r] = index_map["fp64"]

    for (r, c), dtype in tiles.items():
        if r < nt and c < nt:
            grid[r][c] = index_map.get(dtype, index_map["unknown"])

    plt.figure(figsize=(7, 7))
    plt.imshow(grid, origin="lower", interpolation="nearest", cmap=cmap,
               vmin=0, vmax=len(dtype_order) - 1)
    plt.title("Mixed-Precision Tile Map (Diagonal FP64)")
    plt.xlabel("col")
    plt.ylabel("row")

    legend_handles = [
        Patch(facecolor=color_map[dtype], label=dtype)
        for dtype in dtype_order
    ]
    plt.legend(handles=legend_handles, loc="upper right", framealpha=0.9)

    bin_name = Path(args.bin).name if args.bin else ""
    corr = "unknown"
    lower_name = bin_name.lower()
    if "weak" in lower_name:
        corr = "weak"
    elif "strong" in lower_name:
        corr = "strong"

    info_lines = [
        f"bin: {bin_name or 'n/a'}",
        f"format: {args.format or 'auto'}",
        f"n: {args.n or 'n/a'}",
        f"tile: {args.nb or 'n/a'}",
        f"correlation: {corr}",
        f"rel error: {args.rel_err or 'n/a'}",
    ]
    info_text = " | ".join(info_lines)
    plt.gcf().text(0.5, 0.02, info_text, ha="center", va="bottom", fontsize=9,
                   bbox=dict(facecolor="white", edgecolor="black", alpha=0.8))
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path)


if __name__ == "__main__":
    main()
