#!/usr/bin/env python3
"""Paper-grade companion figure: tile allocation vs eps at N=32768.

For each of the 4 sweep families (low-only / low+mid / legacy 3-tier / ladder)
plot a stacked bar per eps showing the per-format tile count. 2x2 panel grid,
matched styling to `plot_paper_granularity_32k.py`.

Uses the vec1d32 variant of each family (canonical OCP MX storage). For
non-bound families the allocation is essentially identical across
granularities, so the choice doesn't change the picture; for the ladder
family this is the post-vec1d-bound-fix selection.
"""
import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


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
    "grid.alpha": 0.25,
    "grid.linewidth": 0.5,
})

EPS_ORDER = ["1e-8", "1e-7", "1e-6", "1e-5"]

# Family-key -> (panel title, vec1d sweep name).
FAMILIES = [
    ("Low-only scaled  (MX_E4M3 low, FP32/FP16 plain)",
     "requant_lowscale_vec1d32_gt20k"),
    ("Low + Mid scaled  (MX_E4M3 + MX_FP16, FP32 plain)",
     "requant_lowmidscale_vec1d32_gt20k"),
    ("All tiers scaled  (MX_FP32 + MX_FP16 + MX_E4M3)",
     "requant_legacy_scaled_vec1d32_gt20k"),
    ("Full ladder  (MXFP4 → MXFP8-E4M3 → MXFP16 → FP32)",
     "requant_ladder_scaled_vec1d32_gt20k"),
]

# Format stack order (top of bar = highest precision; matches typical ladder reading).
# Display label, csv key (canonical), color.
FORMAT_STACK = [
    ("FP64",      "fp64",     "#08306B"),  # darkest
    ("FP32",      "fp32",     "#2171B5"),
    ("MXFP32",    "mx_fp32",  "#6BAED6"),
    ("FP16",      "fp16",     "#9ECAE1"),
    ("MXFP16",    "mx_fp16",  "#FDD49E"),
    ("BF16",      "bf16",     "#FDAE61"),
    ("MXFP8 E5M2","mx_e5m2",  "#FD8D3C"),
    ("MXFP8 E4M3","mx_e4m3",  "#E6550D"),
    ("MXFP4 E2M1","e2m1",     "#A63603"),  # darkest at bottom of ladder
]

CANON = {
    "mx_f16": "mx_fp16",
    "mx_f32": "mx_fp32",
    "mx_fp8_e4m3": "mx_e4m3",
    "mx_fp8_e5m2": "mx_e5m2",
    "fp8e4m3": "fp8_e4m3",
    "fp8e5m2": "fp8_e5m2",
    "fp6_e3m2": "e3m2",
    "fp6e3m2":  "e3m2",
    "fp6_e2m3": "e2m3",
    "fp6e2m3":  "e2m3",
    "fp4_e2m1": "e2m1",
    "fp4e2m1":  "e2m1",
}


def canon_fmt(n):
    n = (n or "").strip().lower()
    return CANON.get(n, n)


def parse_counts(s):
    out = {}
    for part in (s or "").split(","):
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        k = canon_fmt(k)
        try:
            out[k] = out.get(k, 0) + int(v.strip())
        except ValueError:
            pass
    return out


def load_counts(csv_path, sweep, n_target):
    """Return {eps: {fmt: count}} from `tile_counts_full` column."""
    by = {}
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            try:
                if int(r["n"]) != n_target:
                    continue
            except (ValueError, KeyError):
                continue
            if r.get("sweep") != sweep:
                continue
            eps = r.get("source_epsilon", "")
            counts = parse_counts(r.get("tile_counts_full", ""))
            if counts:
                by[eps] = counts
    return by


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv",
                    default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--n", type=int, default=32768)
    ap.add_argument("--out", required=True)
    ap.add_argument("--normalize", action="store_true",
                    help="Plot fractions [0,1] instead of raw tile counts.")
    args = ap.parse_args()

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True, sharey=True)
    axes = axes.flatten()

    bar_width = 0.62

    for ax, (title, sweep) in zip(axes, FAMILIES):
        per_eps = load_counts(args.csv, sweep, args.n)
        if not per_eps:
            ax.set_title(f"{title}\n[no data]", fontsize=10.5)
            continue

        x = np.arange(len(EPS_ORDER))

        # Build stack values for each format. y[k][i] = count at eps_i for format k.
        bottoms = np.zeros(len(EPS_ORDER))
        for label, key, color in FORMAT_STACK:
            vals = np.array([per_eps.get(e, {}).get(key, 0) for e in EPS_ORDER],
                            dtype=float)
            if args.normalize:
                totals = np.array([sum(per_eps.get(e, {}).values()) or 1
                                   for e in EPS_ORDER], dtype=float)
                vals = vals / totals
            if vals.sum() == 0:
                continue
            ax.bar(x, vals, bar_width, bottom=bottoms,
                   label=label, color=color,
                   edgecolor="black", linewidth=0.35)
            bottoms += vals

        ax.set_xticks(x)
        ax.set_xticklabels([r"$10^{-8}$", r"$10^{-7}$", r"$10^{-6}$", r"$10^{-5}$"])
        ax.set_title(title, fontsize=10.5)
        ax.set_axisbelow(True)
        ax.set_xlim(-0.5, len(EPS_ORDER) - 0.5)
        # Y-axis only grid (x-axis is categorical).
        ax.xaxis.grid(False)
        ax.yaxis.grid(True, alpha=0.25, linewidth=0.5)

    # Shared axis labels
    for ax in axes[2:]:
        ax.set_xlabel(r"Source epsilon $\varepsilon$")
    for ax in axes[::2]:
        ax.set_ylabel("Tile fraction" if args.normalize else "Tile count (full matrix)")

    # Single legend at bottom.
    handles, labels = axes[0].get_legend_handles_labels()
    # If panel 0 lacks some formats, union with others
    seen = set(labels)
    for ax in axes[1:]:
        for h, l in zip(*ax.get_legend_handles_labels()):
            if l not in seen:
                handles.append(h)
                labels.append(l)
                seen.add(l)

    # Reorder legend to match FORMAT_STACK order, top-down
    label_to_handle = dict(zip(labels, handles))
    ordered_labels = [lbl for lbl, _, _ in FORMAT_STACK if lbl in label_to_handle]
    ordered_handles = [label_to_handle[l] for l in ordered_labels]

    fig.legend(ordered_handles, ordered_labels, loc="lower center",
               ncol=min(len(ordered_labels), 5),
               bbox_to_anchor=(0.5, -0.04),
               frameon=True, framealpha=0.95,
               title=f"Storage format  (N={args.n}, vec1d-32)")

    fig.tight_layout(rect=[0, 0.06, 1, 0.99])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pdf", ".png"):
        p = out.with_suffix(ext)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(p)
    plt.close(fig)


if __name__ == "__main__":
    main()
