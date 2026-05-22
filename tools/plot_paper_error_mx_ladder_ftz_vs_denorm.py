#!/usr/bin/env python3
"""Two-panel figure for MX-ladder FTZ vs denorm at N=32768.

Three bars per epsilon cluster:
  1. Plain FP8 baseline             (main CSV, requant_baseline_fp8_subnormal_gt20k)
  2. MX ladder vec1D-32, denorm     (ladder_rerun_32k, ladder_full_gu)
  3. MX ladder vec1D-32, FTZ        (ladder_rerun_32k, ladder_full_fz)

Panels:
  Left  - relative factorization error (log scale, per epsilon)
  Right - lower-triangular tile allocation, stacked by tile format

The right panel makes any shift in tile rungs between FTZ and denorm visible:
bar color = tile format, bar hatch = underflow handling (matches the left panel).

Only 32k is rendered because a clean FTZ run of the canonical full ladder
exists only at N=32768 (the rerun on the patched binary).  PDF only.
"""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_paper_tile_vs_mx_ftz_vs_denorm_32k_v2 import (  # noqa: E402
    parse_tile_breakdown, parse_tile_counts_full, tile_counts_lower,
)


mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9.5,
    "legend.title_fontsize": 10.5,
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
})

EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]
N         = 32768
NB        = 2048

MAIN_CSV  = "/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv"
RERUN_CSV = "/home/abduraa/MX_project/ooc-cholesky/ladder_rerun_32k/results.csv"

# (label, source-key, sweep-name, edge-color-for-error-bars, hatch)
BARS = [
    ("Plain FP8 baseline (no MX; subnormal-allowed)",
     "main",  "requant_baseline_fp8_subnormal_gt20k", "#4e79a7", ".."),
    ("MX ladder vec1D-32, denorm  (MXFP4 -> MXFP8 -> MXFP16 -> FP32 -> FP64)",
     "rerun", "ladder_full_gu",                       "#D62728", ""),
    ("MX ladder vec1D-32, FTZ      (MXFP4 -> MXFP8 -> MXFP16 -> FP32 -> FP64)",
     "rerun", "ladder_full_fz",                       "#D62728", "//"),
]

# Lower triangle tile formats (stack order = bottom -> top).  Colors match the
# 3-panel v2 script so figures share a visual language.
TILE_FORMATS = [
    ("FP64",            "fp64",     "#0072B2"),
    ("FP32",            "fp32",     "#E69F00"),
    ("MXFP16 / FP16",   "mx_fp16",  "#7A2C00"),
    ("FP8 plain (E4M3)","fp8_e4m3", "#56B4E9"),
    ("MXFP8 (E4M3)",    "mx_e4m3",  "#009E73"),
    ("MXFP4 (E2M1)",    "e2m1",     "#F0E442"),
]


def load_csv(csv_path, n_target, nb_target):
    rows = {}
    p = Path(csv_path)
    if not p.exists():
        return rows
    with p.open() as f:
        for r in csv.DictReader(f):
            try:
                if int(r["n"]) != n_target or int(r["nb"]) != nb_target:
                    continue
            except (ValueError, KeyError):
                continue
            rows[(r["sweep"], r["source_epsilon"])] = r
    return rows


def err(row):
    if row is None:
        return float("nan")
    try:
        return float(row["rel_factor_error"])
    except (TypeError, ValueError, KeyError):
        return float("nan")


def lower_tri_counts(row, source):
    """Return a {format: tile_count} dict in lower-triangular convention."""
    if row is None:
        return {}
    if source == "rerun":
        # rerun's tile_breakdown column is already lower-triangular (each
        # TILE_TARGET line counted once in the diagnostic phase).
        return parse_tile_breakdown(row.get("tile_breakdown", ""))
    # main-CSV rows store FULL-square counts in tile_counts_full -> halve.
    full = parse_tile_counts_full(row.get("tile_counts_full", ""))
    return tile_counts_lower(full, N, NB)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures")
    args = ap.parse_args()

    data = {
        "main":  load_csv(MAIN_CSV,  N, NB),
        "rerun": load_csv(RERUN_CSV, N, NB),
    }

    fig, (ax_err, ax_tile) = plt.subplots(
        1, 2, figsize=(15.5, 5.4),
        gridspec_kw={"width_ratios": [1.0, 1.05]},
    )

    n_eps  = len(EPS_ORDER)
    n_bars = len(BARS)
    group_w = 0.78
    bar_w   = group_w / n_bars
    x_centres = list(range(n_eps))

    # ---- LEFT: error.
    y_min, y_max = float("inf"), 0.0
    for bi, (label, src, sweep, color, hatch) in enumerate(BARS):
        xs, ys = [], []
        for ci, eps in zip(x_centres, EPS_ORDER):
            v = err(data[src].get((sweep, eps)))
            xs.append(ci - group_w/2 + (bi + 0.5) * bar_w)
            ys.append(v)
        ax_err.bar(xs, ys, bar_w * 0.92,
                   color=color, edgecolor="black", linewidth=0.4,
                   hatch=hatch, label=label)
        for x, y in zip(xs, ys):
            if np.isfinite(y) and y > 0:
                y_min = min(y_min, y); y_max = max(y_max, y)
                ax_err.text(x, y * 1.15, f"{y:.1e}",
                            ha="center", va="bottom",
                            fontsize=8, rotation=90, color="#222")

    for c, eps in zip(x_centres, EPS_ORDER):
        ev = float(eps)
        ax_err.hlines(ev, c - group_w/2 - 0.06, c + group_w/2 + 0.06,
                      colors="gray", linestyles="--", linewidth=0.9, alpha=0.7, zorder=0)
        y_min = min(y_min, ev)

    ax_err.set_xticks(x_centres)
    ax_err.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                            for e in EPS_ORDER])
    ax_err.set_yscale("log")
    ax_err.set_ylabel(r"Relative factorization error  $\|LL^\top - A\|_\infty / \|A\|_\infty$  (log)")
    ax_err.set_axisbelow(True)
    if y_min < float("inf"):
        ax_err.set_ylim(y_min / 4, max(y_max, max(float(e) for e in EPS_ORDER)) * 20)
    ax_err.text(0.985, 0.965, f"N = {N // 1024}k\nNB = {NB}\nError",
                transform=ax_err.transAxes, ha="right", va="top",
                fontsize=10.5,
                bbox=dict(boxstyle="round,pad=0.32",
                          facecolor="white", edgecolor="#888",
                          linewidth=0.8, alpha=0.92))

    # ---- RIGHT: stacked tile allocation (lower triangle).
    max_t = 0.0
    drawn_formats = set()
    for bi, (label, src, sweep, _color, hatch) in enumerate(BARS):
        xs = [c - group_w/2 + (bi + 0.5) * bar_w for c in x_centres]
        per_eps = []
        for eps in EPS_ORDER:
            r = data[src].get((sweep, eps))
            per_eps.append(lower_tri_counts(r, src))

        bottoms = [0.0] * n_eps
        for fmt_lbl, fmt_key, fmt_color in TILE_FORMATS:
            vals = [pe.get(fmt_key, 0) for pe in per_eps]
            if not any(v > 0 for v in vals):
                continue
            drawn_formats.add(fmt_lbl)
            ax_tile.bar(xs, vals, bar_w * 0.92, bottom=bottoms,
                        color=fmt_color, edgecolor="black", linewidth=0.3,
                        hatch=hatch)
            bottoms = [b + v for b, v in zip(bottoms, vals)]
        for x, t in zip(xs, bottoms):
            if t > 0:
                ax_tile.text(x, t * 1.01, f"{int(round(t))}",
                             ha="center", va="bottom",
                             fontsize=7.5, rotation=90, color="#222")
        max_t = max(max_t, max(bottoms) if bottoms else 0.0)

    ax_tile.set_xticks(x_centres)
    ax_tile.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                             for e in EPS_ORDER])
    ax_tile.set_xlim(-0.5, n_eps - 0.5)
    ax_tile.set_ylim(0, max_t * 1.18 if max_t > 0 else None)
    ax_tile.set_ylabel("Tile count (lower triangle, incl. diagonal)")
    ax_tile.set_axisbelow(True)
    ax_tile.xaxis.grid(False)
    ax_tile.yaxis.grid(True, alpha=0.3, linewidth=0.5)
    ax_tile.text(0.985, 0.965, f"N = {N // 1024}k\nNB = {NB}\nTile allocation",
                 transform=ax_tile.transAxes, ha="right", va="top",
                 fontsize=10.5,
                 bbox=dict(boxstyle="round,pad=0.32",
                           facecolor="white", edgecolor="#888",
                           linewidth=0.8, alpha=0.92))

    # ---- Compound legend at the bottom: bar-mode hatches + tile-format colors.
    bar_handles = [Patch(facecolor=color, edgecolor="black", hatch=hatch, label=label)
                   for label, _, _, color, hatch in BARS]
    fmt_handles = [Patch(facecolor=c, edgecolor="black", label=l)
                   for l, _, c in TILE_FORMATS if l in drawn_formats]

    leg1 = fig.legend(handles=bar_handles,
                      loc="lower center", bbox_to_anchor=(0.27, -0.13),
                      ncol=1, frameon=True, framealpha=0.95,
                      title="Bars (left -> right within each cluster)   .   grey dashed = source $\\varepsilon$",
                      handletextpad=0.6, labelspacing=0.35, borderpad=0.4)
    fig.add_artist(leg1)
    fig.legend(handles=fmt_handles,
               loc="lower center", bbox_to_anchor=(0.78, -0.13),
               ncol=2, frameon=True, framealpha=0.95,
               title="Tile format (right panel stack)",
               handletextpad=0.4, columnspacing=1.4, labelspacing=0.3,
               borderpad=0.4)

    fig.tight_layout(rect=[0, 0.06, 1, 1.0])

    out = Path(args.out_dir) / "paper_error_mx_ladder_ftz_vs_denorm_N32768.pdf"
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(out)
    plt.close(fig)


if __name__ == "__main__":
    main()
