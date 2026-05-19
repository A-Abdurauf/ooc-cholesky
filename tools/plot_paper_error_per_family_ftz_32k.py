#!/usr/bin/env python3
"""Paper-grade figure: per-family granularity comparison under FTZ at N=32768.

2×2 panel grid (one per sweep family). Each panel: 4 ε clusters, 3 bars per
cluster — baseline FTZ (no MX), block-128 FTZ, MX-native vec1d-32 FTZ. Y-axis
is log rel_factor_error.

This is the FTZ analog of `plot_paper_error_per_family_with_baselines_32k.py`:
all bars in the comparison use `MX_FP8_SUBNORMAL=0` (flush small FP8 values to
zero). No tile-mode FTZ sweep exists; this figure stops at the
block-128 → vec1d-32 drop-in.
"""
import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 10.5,
    "axes.titlesize": 11.5,
    "axes.labelsize": 10.5,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
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

BASELINE_FTZ_SWEEP = "requant_baseline_fp8_gt20k"

# (panel_title, [(granularity_label, sweep_name), ...]) per sweep family.
FAMILIES = [
    ("Low only scaled  (MXFP8 only; FP32 / FP16 plain)", [
        ("tile",        "requant_lowscale_tile_FTZ_gt20k"),
        ("block 128",   "requant_lowscale_block128_FTZ_gt20k"),
        ("MX",          "requant_lowscale_vec1d32_FTZ_gt20k"),
    ]),
    ("Low + Mid scaled  (MXFP8 + MXFP16; FP32 plain)", [
        ("tile",        "requant_lowmidscale_tile_FTZ_gt20k"),
        ("block 128",   "requant_lowmidscale_block128_FTZ_gt20k"),
        ("MX",          "requant_lowmidscale_vec1d32_FTZ_gt20k"),
    ]),
    ("All tiers scaled  (MXFP32 + MXFP16 + MXFP8)", [
        ("tile",        "requant_legacy_scaled_tile_FTZ_gt20k"),
        ("block 128",   "requant_legacy_scaled_block128_FTZ_gt20k"),
        ("MX",          "requant_legacy_scaled_vec1d32_FTZ_gt20k"),
    ]),
    ("Full ladder  (MXFP4 → MXFP8 → MXFP16 → FP32)", [
        ("tile",        "requant_ladder_scaled_tile_FTZ_gt20k"),
        ("block 128",   "requant_ladder_scaled_block128_FTZ_gt20k"),
        ("MX",          "requant_ladder_scaled_vec1d32_FTZ_gt20k"),
    ]),
]

# Bar order per ε cluster — baseline FP8 first, then tile→block-128→MX progression.
BAR_ORDER = [
    ("baseline IEEE FP8 FTZ  (no MX)", BASELINE_FTZ_SWEEP, "..",  None,        "#4e79a7"),
    ("tile FTZ",                       None,               "",    "tile",      "#f28e2b"),
    ("block 128 FTZ",                  None,               "///", "block 128", "#59a14f"),
    ("MX FTZ  (1×32 native)",          None,               "xx",  "MX",        "#b07aa1"),
]


def load(csv_path, n_target):
    rows = {}
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            try:
                if int(r["n"]) != n_target: continue
            except: continue
            rows[(r["sweep"], r["source_epsilon"])] = r
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--n", type=int, default=32768)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    data = load(args.csv, args.n)

    fig, axes = plt.subplots(2, 2, figsize=(12, 6.0), sharey=True)
    axes = axes.flatten()

    n_eps = len(EPS_ORDER)
    n_bars = len(BAR_ORDER)
    group_width = 0.82
    bar_w = group_width / n_bars
    x_centres = list(range(n_eps))

    y_min, y_max = float("inf"), 0.0

    for ax, (title, grans) in zip(axes, FAMILIES):
        gran_map = {lbl: sw for lbl, sw in grans}
        for bi, (lbl, fixed_sweep, hatch, gran_key, color) in enumerate(BAR_ORDER):
            sweep = fixed_sweep if fixed_sweep is not None else gran_map.get(gran_key)
            if sweep is None: continue
            xs = [c - group_width/2 + (bi + 0.5) * bar_w for c in x_centres]
            ys = []
            for eps in EPS_ORDER:
                r = data.get((sweep, eps))
                try:
                    v = float(r["rel_factor_error"]) if r else 0.0
                except (ValueError, TypeError):
                    v = 0.0
                ys.append(v)
            ax.bar(xs, ys, bar_w * 0.92,
                   color=color, edgecolor="black", linewidth=0.3, hatch=hatch)
            for x, y in zip(xs, ys):
                if y > 0:
                    y_min = min(y_min, y); y_max = max(y_max, y)
                    ax.text(x, y * 1.15, f"{y:.1e}",
                            ha="center", va="bottom", fontsize=6.5, rotation=90,
                            color="#222")

        # Dashed grey reference line per ε column at the requested ε target.
        for c, eps in zip(x_centres, EPS_ORDER):
            ev = float(eps)
            ax.hlines(ev, c - group_width/2 - 0.06, c + group_width/2 + 0.06,
                      colors="gray", linestyles="--", linewidth=0.8, alpha=0.6, zorder=0)

        ax.set_xticks(x_centres)
        ax.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                            for e in EPS_ORDER], fontsize=9)
        ax.set_title(title, fontsize=9.5)
        ax.set_yscale("log")
        ax.set_axisbelow(True)

    # Headroom for the rotated scientific labels.
    if y_min < float("inf"):
        for ax in axes:
            ax.set_ylim(y_min / 3, y_max * 14)
            ax.tick_params(axis="y", labelsize=8.5)

    # Y-axis label on left column only.
    for ax in axes[::2]:
        ax.set_ylabel(r"Relative factorization error  $\|LL^\top - A\|_\infty / \|A\|_\infty$  (log)",
                      fontsize=9)

    # Bar-order legend.
    bar_handles = [
        Patch(facecolor=color, edgecolor="black", hatch=hatch, label=label)
        for label, _, hatch, _, color in BAR_ORDER
    ]
    fig.legend(handles=bar_handles,
               loc="lower center", bbox_to_anchor=(0.5, -0.05),
               ncol=3, frameon=True, framealpha=0.95,
               title="Bar order within each ε cluster (left → right)   ·   grey dashed = source ε",
               handletextpad=0.4, columnspacing=1.5, labelspacing=0.25,
               borderpad=0.35)

    fig.tight_layout(rect=[0, 0.07, 1, 0.99])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pdf",):
        p = out.with_suffix(ext)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(p)
    plt.close(fig)


if __name__ == "__main__":
    main()
