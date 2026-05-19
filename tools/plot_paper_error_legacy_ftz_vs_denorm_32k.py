#!/usr/bin/env python3
"""Paper-grade figure: legacy_scaled (All tiers scaled) at N=32768, granularity
× FTZ-vs-denormal comparison.

Single panel, 4 ε clusters, 6 bars per cluster — paired adjacent (denorm, FTZ)
for each granularity (tile, block-128, MX/vec1d-32). Highlights the FTZ vs
subnormal swing within the All-tiers-scaled (MXFP32+MXFP16+MXFP8) family.
"""
import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
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
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
})

EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]

BASELINE_SWEEP = "requant_baseline_fp8_gt20k"
BASELINE_COLOR = "#4e79a7"  # blue

# Granularity → base color. Within each granularity: (denorm, FTZ) pair.
# denorm is solid fill, FTZ is hatched fill in the same color.
GRANS = [
    ("tile",        "#f28e2b",  # orange
     "requant_legacy_scaled_tile_gt20k",
     "requant_legacy_scaled_tile_FTZ_gt20k"),
    ("block-128",   "#59a14f",  # green
     "requant_legacy_scaled_block128_gt20k",
     "requant_legacy_scaled_block128_FTZ_gt20k"),
    ("MX",   "#b07aa1",  # purple
     "requant_legacy_scaled_vec1d32_gt20k",
     "requant_legacy_scaled_vec1d32_FTZ_gt20k"),
]

# Bars per ε cluster: baseline first, then (denorm, FTZ) pair per granularity.
def bar_specs():
    out = [("Baseline IEEE", BASELINE_COLOR, "..", BASELINE_SWEEP)]
    for gran_lbl, color, sweep_denorm, sweep_ftz in GRANS:
        out.append((f"{gran_lbl}  denorm", color, "",   sweep_denorm))
        out.append((f"{gran_lbl}  FTZ",    color, "xx", sweep_ftz))
    return out


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
    ap.add_argument("--csv",
                    default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--n", type=int, default=32768)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    data = load(args.csv, args.n)
    BARS = bar_specs()

    fig, ax = plt.subplots(figsize=(11, 5.6))
    n_eps  = len(EPS_ORDER)
    n_bars = len(BARS)
    group_width = 0.86
    bar_w = group_width / n_bars
    x_centres = list(range(n_eps))

    y_min, y_max = float("inf"), 0.0

    for bi, (label, color, hatch, sweep) in enumerate(BARS):
        xs, ys = [], []
        for ci, eps in zip(x_centres, EPS_ORDER):
            r = data.get((sweep, eps))
            try:
                v = float(r["rel_factor_error"]) if r else 0.0
            except (ValueError, TypeError):
                v = 0.0
            xs.append(ci - group_width/2 + (bi + 0.5) * bar_w)
            ys.append(v)
        ax.bar(xs, ys, bar_w * 0.92,
               color=color, edgecolor="black", linewidth=0.4, hatch=hatch)
        for x, y in zip(xs, ys):
            if y > 0:
                y_min = min(y_min, y); y_max = max(y_max, y)
                ax.text(x, y * 1.15, f"{y:.1e}",
                        ha="center", va="bottom", fontsize=6.5, rotation=90,
                        color="#222")

    # Dashed grey reference per ε column at source ε.
    for c, eps in zip(x_centres, EPS_ORDER):
        ev = float(eps)
        ax.hlines(ev, c - group_width/2 - 0.07, c + group_width/2 + 0.07,
                  colors="gray", linestyles="--", linewidth=0.9, alpha=0.7, zorder=0)

    ax.set_xticks(x_centres)
    ax.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                        for e in EPS_ORDER])
    ax.set_yscale("log")
    ax.set_ylabel(r"Relative factorization error  $\|LL^\top - A\|_\infty / \|A\|_\infty$  (log)")
    ax.set_axisbelow(True)
    ax.set_xlim(-0.5, n_eps - 0.5)
    ax.set_title(f"All tiers scaled (MXFP32+MXFP16+MXFP8)  ·  N={args.n}, NB=2048")
    if y_min < float("inf"):
        ax.set_ylim(y_min / 3, y_max * 15)

    # Custom legend: baseline + per-granularity colors + denorm/FTZ pattern hint.
    gran_handles = [Patch(facecolor=BASELINE_COLOR, edgecolor="black",
                          hatch="..", label="Baseline IEEE")] + [
        Patch(facecolor=color, edgecolor="black", label=lbl)
        for lbl, color, _, _ in GRANS
    ]
    pat_handles = [
        Patch(facecolor="white", edgecolor="black", hatch="",   label="denorm"),
        Patch(facecolor="white", edgecolor="black", hatch="xx", label="FTZ"),
    ]
    lg1 = fig.legend(handles=gran_handles, loc="lower center",
                     bbox_to_anchor=(0.30, -0.05), ncol=4,
                     frameon=True, framealpha=0.95, title="Variant")
    fig.add_artist(lg1)
    fig.legend(handles=pat_handles, loc="lower center",
               bbox_to_anchor=(0.80, -0.05), ncol=2,
               frameon=True, framealpha=0.95, title="Mode")

    fig.tight_layout(rect=[0, 0.05, 1, 0.99])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pdf",):
        p = out.with_suffix(ext)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(p)
    plt.close(fig)


if __name__ == "__main__":
    main()
