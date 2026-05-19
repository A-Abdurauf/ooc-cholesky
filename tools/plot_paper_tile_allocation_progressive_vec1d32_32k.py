#!/usr/bin/env python3
"""Paper-grade figure: tile-allocation companion to the progressive-error plot.

Single panel, 4 eps clusters, 6 side-by-side stacked bars per cluster showing
per-format tile counts for the cumulative MX-scaling progression
(baseline IEEE -> ladder IEEE -> +SFP8 -> +SFP8+SFP16 -> +SFP8+SFP16+MXFP32 ->
Full ladder). All MX-scaled bars use vec1D-32 granularity.

Pairs 1:1 with `plot_paper_error_progressive_vec1d32_32k.py`: each error bar
in that figure has a corresponding stacked allocation bar here.
"""
import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch


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

EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]

# Progression of bars per eps cluster. (short label, sweep name)
BAR_ORDER = [
    ("Baseline IEEE",  "requant_baseline_fp8_subnormal_gt20k"),
    ("Ladder IEEE",    "requant_ladder_ieee_gt20k"),
    ("Full ladder",    "requant_ladder_scaled_vec1d32_gt20k"),
]

# Format stack order: highest precision on top (deep blue) → lowest on bottom.
FORMAT_STACK = [
    ("FP64",        "fp64",       "#08306B"),
    ("FP32",        "fp32",       "#2171B5"),
    ("MXFP32",      "mx_fp32",    "#6BAED6"),
    ("FP16",        "fp16",       "#9ECAE1"),
    ("MXFP16",      "mx_fp16",    "#FDD49E"),
    ("BF16",        "bf16",       "#FDAE61"),
    ("FP8 E4M3",    "fp8_e4m3",   "#F16913"),
    ("MXFP8 E5M2",  "mx_e5m2",    "#FD8D3C"),
    ("MXFP8 E4M3",  "mx_e4m3",    "#E6550D"),
    ("MXFP4 E2M1",  "e2m1",       "#A63603"),
]

CANON = {
    "mx_f16": "mx_fp16", "mx_f32": "mx_fp32",
    "mx_fp8_e4m3": "mx_e4m3", "mx_fp8_e5m2": "mx_e5m2",
    "fp8e4m3": "fp8_e4m3",   "fp8e5m2": "fp8_e5m2",
    "fp6_e3m2": "e3m2", "fp6e3m2": "e3m2",
    "fp6_e2m3": "e2m3", "fp6e2m3": "e2m3",
    "fp4_e2m1": "e2m1", "fp4e2m1": "e2m1",
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


def load_counts(csv_path, n_target):
    """Returns dict[(sweep, eps)] -> {fmt: count}."""
    by = {}
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            try:
                if int(r["n"]) != n_target: continue
            except (ValueError, KeyError):
                continue
            counts = parse_counts(r.get("tile_counts_full", ""))
            if counts:
                # Last row per (sweep, eps) wins (picks post-fix data).
                by[(r["sweep"], r["source_epsilon"])] = counts
    return by


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv",
                    default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--n", type=int, default=32768)
    ap.add_argument("--out", required=True)
    ap.add_argument("--normalize", action="store_true",
                    help="Plot fractions [0,1] instead of raw counts.")
    args = ap.parse_args()

    counts_by = load_counts(args.csv, args.n)

    fig, ax = plt.subplots(figsize=(13, 5.8))

    n_eps  = len(EPS_ORDER)
    n_bars = len(BAR_ORDER)
    group_width = 0.86
    bar_w = group_width / n_bars
    x_centres = list(range(n_eps))

    # Track which formats actually appear so the legend matches what's drawn.
    used = {key: False for _, key, _ in FORMAT_STACK}
    max_y = 0.0

    for bi, (bar_label, sweep) in enumerate(BAR_ORDER):
        xs = [c - group_width / 2 + (bi + 0.5) * bar_w for c in x_centres]
        bottoms = np.zeros(n_eps)
        for label, key, color in FORMAT_STACK:
            vals = np.zeros(n_eps)
            for i, eps in enumerate(EPS_ORDER):
                counts = counts_by.get((sweep, eps), {})
                vals[i] = counts.get(key, 0)
            if args.normalize:
                totals = np.zeros(n_eps)
                for i, eps in enumerate(EPS_ORDER):
                    totals[i] = sum(counts_by.get((sweep, eps), {}).values()) or 1
                vals = vals / totals
            if vals.sum() == 0:
                continue
            used[key] = True
            ax.bar(xs, vals, bar_w * 0.90, bottom=bottoms,
                   color=color, edgecolor="black", linewidth=0.25)
            bottoms += vals
        # Total tile count above each bar.
        for x, total in zip(xs, bottoms):
            if total > 0:
                ax.text(x, total * 1.01,
                        f"{int(round(total))}" if not args.normalize else f"{total:.2f}",
                        ha="center", va="bottom",
                        fontsize=7.0, rotation=90, color="#222")
        max_y = max(max_y, bottoms.max())

    # Group centre tick labels for eps, plus per-bar sub-labels rotated 45°.
    ax.set_xticks(x_centres)
    ax.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                        for e in EPS_ORDER])
    ax.tick_params(axis="x", which="major", pad=42)  # room for rotated sub-labels

    # Replace newlines with spaces so 45°-rotated labels stay on one line.
    for ci in x_centres:
        for bi, (lbl, _) in enumerate(BAR_ORDER):
            x = ci - group_width / 2 + (bi + 0.5) * bar_w
            ax.text(x, -0.018, lbl.replace("\n", " "),
                    ha="right", va="top",
                    fontsize=6.5, rotation=45,
                    transform=ax.get_xaxis_transform(), clip_on=False)

    ax.set_axisbelow(True)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, alpha=0.25, linewidth=0.5)
    ax.set_xlim(-0.5, n_eps - 0.5)
    # Y headroom for the rotated value labels on top of each bar.
    ax.set_ylim(0, max_y * 1.18)
    ax.set_ylabel("Tile fraction" if args.normalize else "Tile count (full matrix)")
    ax.set_xlabel("")
    ax.set_title(f"Per-tile format mix under progressive MX scaling  (N={args.n}, MX-native storage, 1×32 row groups)")

    # Build legend from drawn formats only.
    legend_entries = [(lbl, color) for lbl, key, color in FORMAT_STACK if used[key]]
    handles = [Patch(facecolor=color, edgecolor="black", label=lbl)
               for lbl, color in legend_entries]
    fig.legend(handles=handles,
               loc="lower center", bbox_to_anchor=(0.5, -0.07),
               ncol=min(len(handles), 5),
               frameon=True, framealpha=0.95,
               title=(f"Storage format  ·  bar order: "
                      f"Baseline IEEE → +MXFP8 → +MXFP16 → +MXFP32 → "
                      f"Ladder IEEE (FP8→FP16→FP32→FP64) → "
                      f"Full ladder (MXFP4→MXFP8→MXFP16→FP32→FP64)"),
               handletextpad=0.4, columnspacing=1.6, labelspacing=0.3,
               borderpad=0.35)

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
