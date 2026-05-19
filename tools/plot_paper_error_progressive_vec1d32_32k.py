#!/usr/bin/env python3
"""Paper-grade figure: progressive MX-scaling story at N=32768, vec1d-32 only.

Single panel, 4 epsilon clusters, 6 bars per cluster showing accumulation of
scaled tiers (baseline IEEE -> ladder IEEE -> SFP8 only -> +MXFP16 -> +MXFP32
-> Full ladder with MXFP4). All MX-scaled bars use the vec1D-32 granularity.

This is the "tell the whole story in one panel" version of
fig_error_per_family_with_baselines_32k, with the granularity axis collapsed.
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

# Progression of bars per eps cluster: each step "MX-scales one more tier".
# (label, sweep_name, hatch, color)
BAR_ORDER = [
    ("Baseline IEEE  (no MX; static format per matrix)",
     "requant_baseline_fp8_subnormal_gt20k",  "..",  "#4e79a7"),  # blue
    ("+ MXFP8 scaled  (MX_E4M3 low; FP32/FP16 plain)",
     "requant_lowscale_vec1d32_gt20k",        "",    "#f1c40f"),  # yellow
    ("+ MXFP8 + MXFP16 scaled  (MX_FP16 mid added)",
     "requant_lowmidscale_vec1d32_gt20k",     "//",  "#f28e2b"),  # orange
    ("+ MXFP8 + MXFP16 + MXFP32 scaled  (all upper tiers MX)",
     "requant_legacy_scaled_vec1d32_gt20k",   "\\\\","#59a14f"),  # green
    ("Ladder IEEE  (FP8 E4M3 → FP16 → FP32 → FP64)",
     "requant_ladder_ieee_gt20k",             "||",  "#a0a0a0"),  # grey
    ("Full ladder  (MXFP4 → MXFP8 E4M3 → MXFP16 → FP32 → FP64)",
     "requant_ladder_scaled_vec1d32_gt20k",   "xx",  "#b07aa1"),  # purple
]


def load(csv_path, n_target):
    rows = {}
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            try:
                if int(r["n"]) != n_target: continue
            except: continue
            # Keep last row per (sweep, eps) — dedup that picks post-fix.
            rows[(r["sweep"], r["source_epsilon"])] = r
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--n", type=int, default=32768)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    data = load(args.csv, args.n)

    fig, ax = plt.subplots(figsize=(11, 5.5))

    n_eps  = len(EPS_ORDER)
    n_bars = len(BAR_ORDER)
    group_width = 0.84
    bar_w = group_width / n_bars
    x_centres = list(range(n_eps))

    y_min, y_max = float("inf"), 0.0

    for bi, (label, sweep, hatch, color) in enumerate(BAR_ORDER):
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
               color=color, edgecolor="black", linewidth=0.4,
               hatch=hatch, label=label)
        for x, y in zip(xs, ys):
            if y > 0:
                y_min = min(y_min, y); y_max = max(y_max, y)
                ax.text(x, y * 1.15, f"{y:.1e}",
                        ha="center", va="bottom", fontsize=7.0, rotation=90,
                        color="#222")

    # Dashed grey reference line per eps column at the requested target.
    for c, eps in zip(x_centres, EPS_ORDER):
        ev = float(eps)
        ax.hlines(ev, c - group_width/2 - 0.08, c + group_width/2 + 0.08,
                  colors="gray", linestyles="--", linewidth=0.9, alpha=0.7, zorder=0)

    ax.set_xticks(x_centres)
    ax.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                        for e in EPS_ORDER])
    ax.set_yscale("log")
    ax.set_ylabel(r"Relative factorization error  $\|LL^\top - A\|_\infty / \|A\|_\infty$  (log)")
    ax.set_xlabel("")
    ax.set_axisbelow(True)
    ax.set_title(f"Progressive MX scaling at N={args.n} (MX-native storage, 1×32 row groups)")

    # Y headroom for rotated scientific labels.
    if y_min < float("inf"):
        ax.set_ylim(y_min / 3, y_max * 15)

    # Bar-order legend reflecting the cumulative-tier story.
    handles = [Patch(facecolor=color, edgecolor="black", hatch=hatch, label=label)
               for label, _, hatch, color in BAR_ORDER]
    fig.legend(handles=handles,
               loc="lower center", bbox_to_anchor=(0.5, -0.18),
               ncol=2, frameon=True, framealpha=0.95,
               title="Bar order within each cluster (left → right)   ·   grey dashed = source ε",
               handletextpad=0.6, columnspacing=2.0, labelspacing=0.35,
               borderpad=0.4)

    fig.tight_layout(rect=[0, 0.0, 1, 0.99])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pdf", ".png"):
        p = out.with_suffix(ext)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(p)
    plt.close(fig)


if __name__ == "__main__":
    main()
