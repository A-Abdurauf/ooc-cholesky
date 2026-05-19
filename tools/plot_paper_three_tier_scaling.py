#!/usr/bin/env python3
"""Paper figure: 3-tier (SFP32 + SFP16 + SFP8) drop-in scaling vs N.

4 subplots, one per requested ε ∈ {1e-5, 1e-6, 1e-7, 1e-8}.
Each panel: x = matrix size N, y = total GB (--metric=mem) or rel
factorisation error (--metric=err). Three lines per panel — one per
granularity (tile / block 128 / vec1D 32) — plus the IEEE baseline as
a dashed reference.

Usage:
  python3 plot_paper_three_tier_scaling.py --out path/fig.pdf
  python3 plot_paper_three_tier_scaling.py --metric err --out path/err.pdf
"""
import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


mpl.rcParams.update({
    "font.family":         "serif",
    "font.size":           11,
    "axes.titlesize":      12,
    "axes.labelsize":      11,
    "xtick.labelsize":     10,
    "ytick.labelsize":     10,
    "legend.fontsize":     10,
    "legend.title_fontsize": 10,
    "figure.dpi":          120,
    "savefig.dpi":         300,
    "pdf.fonttype":        42,
    "ps.fonttype":         42,
    "axes.grid":           True,
    "grid.alpha":          0.3,
    "grid.linewidth":      0.5,
    "lines.linewidth":     2.0,
    "lines.markersize":    8,
})

EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]

# All-tiers-scaled family: SFP32 + SFP16 + SFP8 at three granularities.
GRAN_SWEEP = [
    ("tile",       "requant_legacy_scaled_tile_gt20k"),
    ("block 128",  "requant_legacy_scaled_block128_gt20k"),
    ("vec1D 32",   "requant_legacy_scaled_vec1d32_gt20k"),
]
GRAN_COLOR = {
    "tile":       "#f28e2b",
    "block 128":  "#59a14f",
    "vec1D 32":   "#b07aa1",
}
GRAN_MARKER = {"tile": "o", "block 128": "s", "vec1D 32": "^"}

BASELINE_SWEEP = "requant_baseline_fp8_subnormal_gt20k"


def load(csv_path):
    rows = {}
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            try:
                n = int(r["n"])
            except (ValueError, KeyError):
                continue
            rows[(r["sweep"], r["source_epsilon"], n)] = r
    return rows


def value(r, metric):
    if r is None:
        return None
    try:
        if metric == "mem":
            v = float(r.get("total_gb", 0) or 0)
        else:
            v = float(r.get("rel_factor_error", 0) or 0)
        return v if v > 0 else None
    except (ValueError, TypeError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv",
        default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--metric", choices=("mem", "err"), default="mem")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    data = load(args.csv)
    Ns = sorted({n for (_, _, n) in data})
    if not Ns:
        raise SystemExit("No N values found in CSV.")

    fig, axes = plt.subplots(1, 4, figsize=(14, 4.0), sharey=True,
                             constrained_layout=False)

    metric_label = {
        "mem": "GB per Cholesky factor",
        "err": r"$\|LL^\top - A\|_\infty\,/\,\|A\|_\infty$",
    }[args.metric]
    use_log = (args.metric == "err")

    y_min, y_max = float("inf"), 0.0

    for ax, eps in zip(axes, EPS_ORDER):
        # Three granularity curves
        for gran_lbl, sweep in GRAN_SWEEP:
            ys, xs = [], []
            for n in Ns:
                v = value(data.get((sweep, eps, n)), args.metric)
                if v is None:
                    continue
                xs.append(n); ys.append(v)
            if not xs:
                continue
            ax.plot(xs, ys,
                    color=GRAN_COLOR[gran_lbl],
                    marker=GRAN_MARKER[gran_lbl],
                    label=gran_lbl, zorder=3)
            y_min = min(y_min, min(ys))
            y_max = max(y_max, max(ys))

        # IEEE baseline as a dashed reference
        xs_b, ys_b = [], []
        for n in Ns:
            v = value(data.get((BASELINE_SWEEP, eps, n)), args.metric)
            if v is None:
                continue
            xs_b.append(n); ys_b.append(v)
        if xs_b:
            ax.plot(xs_b, ys_b, color="black", ls="--", marker="x",
                    lw=1.4, ms=7, label="IEEE baseline", zorder=2, alpha=0.9)
            y_min = min(y_min, min(ys_b))
            y_max = max(y_max, max(ys_b))

        ax.set_title(f"$\\varepsilon = {eps}$")
        ax.set_xlabel("Matrix size $N$")
        ax.set_xticks(Ns)
        ax.set_xticklabels([f"{n//1024}k" for n in Ns])
        ax.set_axisbelow(True)

    if use_log:
        for ax in axes:
            ax.set_yscale("log")
        # Headroom on log axis
        for ax in axes:
            ax.set_ylim(y_min / 3, y_max * 6)
    else:
        for ax in axes:
            ax.set_ylim(0, y_max * 1.10)

    axes[0].set_ylabel(metric_label)

    # Shared legend at the bottom
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels,
               loc="lower center",
               bbox_to_anchor=(0.5, -0.04),
               ncol=4, frameon=True, framealpha=0.95,
               handletextpad=0.5, columnspacing=1.4)

    fig.suptitle(
        "3-tier drop-in (SFP32 + SFP16 + SFP8) scaling with $N$",
        fontsize=13, fontweight="bold", y=1.02)
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
