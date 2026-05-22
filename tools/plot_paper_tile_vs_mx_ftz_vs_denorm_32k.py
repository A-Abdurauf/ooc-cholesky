#!/usr/bin/env python3
"""Three-panel figure: error · memory · tile allocation for the ladder at
32k/nb=2048 across {Tile, MX} × {denorm, FTZ}.

Panel 1 (left):   relative factor error per ε cluster, log y, 4 bars.
Panel 2 (middle): per-L-factor memory (lower-tri), stacked by datatype.
Panel 3 (right):  tile-allocation count (lower-tri), stacked by per-tile
                  format. Diagonal FP64 counted once, off-diagonals halved.

Reads:  /home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv
Writes: <out>.pdf and <out>.png
"""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

# Reuse the existing memory transform.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_paper_memory_progressive_vec1d32_32k import (  # noqa: E402
    to_lower_triangular,
    HALVE_COLS,
)


mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.titlesize": 11.5,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9.5,
    "legend.title_fontsize": 10,
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
})

EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]

# (label, sweep, color, hatch)
BARS = [
    ("Baseline IEEE", "requant_baseline_fp8_subnormal_gt20k",  "#555555", ".."),
    ("Tile · denorm", "requant_ladder_scaled_tile_gt20k",      "#1F77B4", ""),
    ("Tile · FTZ",    "requant_ladder_scaled_tile_FTZ_gt20k",  "#1F77B4", "//"),
    ("MX · denorm",   "requant_ladder_scaled_vec1d32_gt20k",   "#D62728", ""),
    ("MX · FTZ",      "requant_ladder_scaled_vec1d32_FTZ_gt20k", "#D62728", "//"),
]

DTYPE_BUCKETS = [
    ("FP64",         ["fp64_gb"],       "#0072B2"),
    ("FP32",         ["fp32_gb"],       "#E69F00"),
    ("FP16/MXFP16",  ["fp16_gb", "mx_fp16_gb"], "#7A2C00"),
    ("FP8 plain",    ["fp8_e4m3_gb", "fp8_e5m2_gb"], "#56B4E9"),
    ("MXFP8 (E4M3)", ["mx_e4m3_gb"],    "#009E73"),
    ("MXFP4 (E2M1)", ["e2m1_gb"],       "#F0E442"),
    ("Scale meta",   ["scale_meta_gb"], "#999999"),
]

# Tile-allocation stack: same colours as the datatype buckets where possible,
# but no "scale meta" (it's not a tile choice).
TILE_FORMATS = [
    ("FP64",         "fp64",     "#0072B2"),
    ("FP32",         "fp32",     "#E69F00"),
    ("FP16/MXFP16",  "mx_fp16",  "#7A2C00"),
    ("FP8 plain",    "fp8_e4m3", "#56B4E9"),
    ("MXFP8 (E4M3)", "mx_e4m3",  "#009E73"),
    ("MXFP4 (E2M1)", "e2m1",     "#F0E442"),
]

# Fold these aliases into the canonical key.
TILE_ALIASES = {
    "fp16":     "mx_fp16",   # plain FP16 stacks into the FP16/MXFP16 band
    "mx_fp32":  "fp32",      # MXFP32 (rarely appears) stacks into FP32
    "fp8_e5m2": "fp8_e4m3",  # both plain FP8 variants share the FP8-plain stack
}


def parse_tile_counts_full(s):
    out = {}
    for part in (s or "").split(","):
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        k = k.strip().lower()
        k = TILE_ALIASES.get(k, k)
        try:
            out[k] = out.get(k, 0) + int(v.strip())
        except ValueError:
            pass
    return out


def tile_counts_lower(counts_full, n, nb):
    """Convert full-square counts to lower-triangular counts.
    Diagonal tiles are FP64; off-diagonals are mirrored (halved)."""
    if n % nb:
        raise ValueError(f"n={n} not divisible by nb={nb}")
    M = n // nb
    out = {}
    for fmt, v in counts_full.items():
        if fmt == "fp64":
            off = max(v - M, 0)
            out[fmt] = M + off // 2
        else:
            out[fmt] = v // 2
    return out


def load(csv_path, n_target, nb_target):
    rows = {}
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            try:
                if int(r["n"]) != n_target or int(r["nb"]) != nb_target:
                    continue
            except (ValueError, KeyError):
                continue
            rows[(r["sweep"], r["source_epsilon"])] = r  # last wins
    return rows


def bucket_value(r, keys):
    if r is None:
        return 0.0
    return sum(float(r.get(k, 0) or 0) for k in keys)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv",
                    default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--n", type=int, default=32768)
    ap.add_argument("--nb", type=int, default=2048)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    data = load(args.csv, args.n, args.nb)
    if not data:
        sys.exit(f"No rows in {args.csv} for n={args.n} nb={args.nb}")

    # Lower-triangular memory transform.
    data_tri = {k: to_lower_triangular(r, args.n, args.nb) for k, r in data.items()}

    fig, (ax_err, ax_mem, ax_tile) = plt.subplots(
        1, 3, figsize=(20.5, 5.8),
        gridspec_kw={"width_ratios": [1, 1.05, 1.05]},
    )

    n_eps = len(EPS_ORDER)
    n_bars = len(BARS)
    group_width = 0.84
    bar_w = group_width / n_bars
    x_centres = list(range(n_eps))

    # ---- LEFT: relative-error bars (log y).
    max_err = 0.0
    for bi, (lbl, sweep, color, hatch) in enumerate(BARS):
        xs = [c - group_width / 2 + (bi + 0.5) * bar_w for c in x_centres]
        vals = []
        for eps in EPS_ORDER:
            r = data.get((sweep, eps))
            v = float(r["rel_factor_error"]) if r and r.get("rel_factor_error") else np.nan
            vals.append(v)
            if np.isfinite(v):
                max_err = max(max_err, v)
        bars = ax_err.bar(xs, vals, bar_w * 0.90, color=color, edgecolor="black",
                          linewidth=0.4, hatch=hatch, label=lbl)
        # Numeric label on top of each bar.
        for x, v in zip(xs, vals):
            if np.isfinite(v) and v > 0:
                ax_err.text(x, v * 1.15, f"{v:.1e}",
                            ha="center", va="bottom",
                            fontsize=7.2, rotation=90, color="#222")

    ax_err.set_yscale("log")
    ax_err.set_xticks(x_centres)
    ax_err.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                            for e in EPS_ORDER])
    ax_err.set_xlim(-0.5, n_eps - 0.5)
    ax_err.set_ylim(top=max_err * 30)
    ax_err.set_ylabel("Relative error")
    ax_err.set_title(f"Error  (ladder, N={args.n}, nb={args.nb})")
    ax_err.set_axisbelow(True)
    ax_err.xaxis.grid(False)
    ax_err.yaxis.grid(True, which="both", alpha=0.3, linewidth=0.5)
    ax_err.legend(loc="upper right", frameon=True, framealpha=0.95)

    # ---- RIGHT: stacked memory bars (lower-triangular).
    drawn = {lbl: False for lbl, _, _ in DTYPE_BUCKETS}
    max_y = 0.0

    for bi, (bar_lbl, sweep, color, hatch) in enumerate(BARS):
        xs = [c - group_width / 2 + (bi + 0.5) * bar_w for c in x_centres]
        bottoms = [0.0] * n_eps
        for dt_lbl, keys, dt_color in DTYPE_BUCKETS:
            vals = [bucket_value(data_tri.get((sweep, eps)), keys) for eps in EPS_ORDER]
            if not any(v > 0 for v in vals):
                continue
            drawn[dt_lbl] = True
            ax_mem.bar(xs, vals, bar_w * 0.90, bottom=bottoms,
                       color=dt_color, edgecolor="black", linewidth=0.3,
                       hatch=hatch)
            bottoms = [b + v for b, v in zip(bottoms, vals)]
        # Total above each bar.
        for x, t in zip(xs, bottoms):
            if t > 0:
                ax_mem.text(x, t * 1.01, f"{t:.2f}",
                            ha="center", va="bottom",
                            fontsize=7.2, rotation=90, color="#222")
        max_y = max(max_y, max(bottoms))

    ax_mem.set_xticks(x_centres)
    ax_mem.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                            for e in EPS_ORDER])
    ax_mem.set_xlim(-0.5, n_eps - 0.5)
    ax_mem.set_ylim(0, max_y * 1.20)
    ax_mem.set_ylabel("Memory per Cholesky L factor — lower triangle (GB)")
    ax_mem.set_title(f"Memory  (ladder, N={args.n}, nb={args.nb}, lower triangle)")
    ax_mem.set_axisbelow(True)
    ax_mem.xaxis.grid(False)
    ax_mem.yaxis.grid(True, alpha=0.3, linewidth=0.5)

    # ---- RIGHT: tile-allocation count (lower triangle, stacked by format).
    max_t = 0.0
    for bi, (bar_lbl, sweep, color, hatch) in enumerate(BARS):
        xs = [c - group_width / 2 + (bi + 0.5) * bar_w for c in x_centres]
        # Pull counts per ε for this sweep.
        per_eps = []
        for eps in EPS_ORDER:
            r = data.get((sweep, eps))
            full = parse_tile_counts_full(r["tile_counts_full"]) if r else {}
            per_eps.append(tile_counts_lower(full, args.n, args.nb))

        bottoms = [0.0] * n_eps
        for fmt_lbl, fmt_key, fmt_color in TILE_FORMATS:
            vals = [pe.get(fmt_key, 0) for pe in per_eps]
            if not any(v > 0 for v in vals):
                continue
            ax_tile.bar(xs, vals, bar_w * 0.90, bottom=bottoms,
                        color=fmt_color, edgecolor="black", linewidth=0.3,
                        hatch=hatch)
            bottoms = [b + v for b, v in zip(bottoms, vals)]
        for x, t in zip(xs, bottoms):
            if t > 0:
                ax_tile.text(x, t * 1.01, f"{int(round(t))}",
                             ha="center", va="bottom",
                             fontsize=7.2, rotation=90, color="#222")
        max_t = max(max_t, max(bottoms))

    ax_tile.set_xticks(x_centres)
    ax_tile.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                             for e in EPS_ORDER])
    ax_tile.set_xlim(-0.5, n_eps - 0.5)
    ax_tile.set_ylim(0, max_t * 1.18)
    ax_tile.set_ylabel("Tile count (lower triangle, incl. diagonal)")
    ax_tile.set_title(f"Tile allocation  (ladder, N={args.n}, nb={args.nb}, lower triangle)")
    ax_tile.set_axisbelow(True)
    ax_tile.xaxis.grid(False)
    ax_tile.yaxis.grid(True, alpha=0.3, linewidth=0.5)

    # Compound legend: datatype colours (left col) + bar-style key (right col).
    dtype_handles = [Patch(facecolor=c, edgecolor="black", label=l)
                     for l, _, c in DTYPE_BUCKETS if drawn[l]]
    style_handles = [
        Patch(facecolor="white", edgecolor="black", label="solid = denorm"),
        Patch(facecolor="white", edgecolor="black", hatch="//", label="hatched = FTZ"),
        Patch(facecolor="#555555", edgecolor="black", hatch="..",
              label="grey = Baseline IEEE"),
        Patch(facecolor="#1F77B4", edgecolor="black", label="blue = Tile"),
        Patch(facecolor="#D62728", edgecolor="black", label="red = MX (vec1D-32)"),
    ]
    leg1 = fig.legend(handles=style_handles,
                      loc="lower center", bbox_to_anchor=(0.22, -0.05),
                      ncol=2, frameon=True, framealpha=0.95, title="Bar style",
                      handletextpad=0.6, columnspacing=1.4, labelspacing=0.3,
                      borderpad=0.35)
    fig.add_artist(leg1)
    fig.legend(handles=dtype_handles,
               loc="lower center", bbox_to_anchor=(0.66, -0.05),
               ncol=3, frameon=True, framealpha=0.95,
               title="Datatype  /  Tile format",
               handletextpad=0.4, columnspacing=1.2, labelspacing=0.3,
               borderpad=0.35)

    fig.subplots_adjust(bottom=0.22)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pdf", ".png"):
        p = out.with_suffix(ext)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(p)
    plt.close(fig)


if __name__ == "__main__":
    main()
