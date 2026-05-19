#!/usr/bin/env python3
"""Paper-grade figure 3: Granularity comparison at N=32768.

For each of the 4 sweep families (low-only / low+mid / legacy 3-tier / ladder),
plot the 3 granularity variants (tile / block128 / vec1D32) on a single
log-error vs ε axis. 2×2 panel grid.
"""
import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt


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
EPS_FLOAT = [float(e) for e in EPS_ORDER]

# 4 sweep families; for each family, the 3 granularities.
FAMILIES = [
    ("Low-only scaled  (MX_E4M3 low tier, FP32/FP16 plain)", "lowscale", "lowonly"),
    ("Low + Mid scaled  (MX_E4M3 + MX_FP16, FP32 plain)",   "lowmid",  "lowmid"),
    ("All tiers scaled  (MX_FP32 + MX_FP16 + MX_E4M3)", "legacy",  "legacy"),
    ("Full ladder  (MXFP4 → MXFP8-E4M3 → MXFP16 → FP32)",   "ladder",  "ladder"),
]

# Family-key → (tile, block128, vec1d32) sweep names.
SWEEP_NAMES = {
    "lowscale": (
        "requant_lowscale_tile_gt20k",
        "requant_lowscale_block128_gt20k",
        "requant_lowscale_vec1d32_gt20k",
    ),
    "lowmid": (
        "requant_lowmidscale_tile_gt20k",
        "requant_lowmidscale_block128_gt20k",
        "requant_lowmidscale_vec1d32_gt20k",
    ),
    "legacy": (
        "requant_legacy_scaled_tile_gt20k",
        "requant_legacy_scaled_block128_gt20k",
        "requant_legacy_scaled_vec1d32_gt20k",
    ),
    "ladder": (
        "requant_ladder_scaled_tile_gt20k",
        "requant_ladder_scaled_block128_gt20k",
        "requant_ladder_scaled_vec1d32_gt20k",
    ),
}

GRANULARITY_LABELS = ["tile  (1 scale/tile)",
                      "block 128  (256 scales/tile)",
                      "vec1D 32  (131,072 scales/tile)"]
GRANULARITY_COLORS = ["#7B3294", "#0571B0", "#CA0020"]
GRANULARITY_MARKERS = ["s", "o", "^"]

BASELINES = [
    ("baseline IEEE (OCP subnormal)", "requant_baseline_fp8_subnormal_gt20k",
     "#000000", "x", "--"),
    ("ladder IEEE (no MX)",            "requant_ladder_ieee_gt20k",
     "#666666", "+", ":"),
]


def load(csv_path, n_target):
    by = {}
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            try:
                if int(r["n"]) != n_target: continue
            except: continue
            by[(r["sweep"], r["source_epsilon"])] = r
    return by


def series(data, sweep):
    xs, ys = [], []
    for eps in EPS_ORDER:
        r = data.get((sweep, eps))
        if not r: continue
        try:
            v = float(r["rel_factor_error"])
            if v > 0:
                xs.append(float(eps))
                ys.append(v)
        except: pass
    return xs, ys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--n", type=int, default=32768)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    data = load(args.csv, args.n)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True, sharey=True)
    axes = axes.flatten()
    for ax, (title, fam_key, _) in zip(axes, FAMILIES):
        tile_sw, block_sw, vec_sw = SWEEP_NAMES[fam_key]
        # Plot baselines first (in background).
        for lbl, sw, col, mk, ls in BASELINES:
            xs, ys = series(data, sw)
            if xs:
                ax.plot(xs, ys, marker=mk, linestyle=ls, color=col,
                        linewidth=1.4, markersize=7, alpha=0.65, label=lbl)
        # Plot family granularities.
        for sw, lbl, col, mk in zip((tile_sw, block_sw, vec_sw),
                                     GRANULARITY_LABELS,
                                     GRANULARITY_COLORS,
                                     GRANULARITY_MARKERS):
            xs, ys = series(data, sw)
            if xs:
                ax.plot(xs, ys, marker=mk, linestyle="-", color=col,
                        linewidth=1.8, markersize=8, markeredgecolor="black",
                        markeredgewidth=0.5, label=lbl)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(title, fontsize=10.5)
        ax.grid(True, which="both", alpha=0.3)

    # Shared axis labels
    for ax in axes[2:]:
        ax.set_xlabel(r"Source epsilon $\varepsilon$  (log)")
    for ax in axes[::2]:
        ax.set_ylabel(r"$\|LL^\top - A\|_\infty / \|A\|_\infty$  (log)")

    # Single legend at top
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center",
               ncol=3, bbox_to_anchor=(0.5, -0.04),
               frameon=True, framealpha=0.95,
               title="Storage geometry & references")

    fig.tight_layout(rect=[0, 0.04, 1, 0.99])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pdf", ".png"):
        p = out.with_suffix(ext)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(p)
    plt.close(fig)


if __name__ == "__main__":
    main()
