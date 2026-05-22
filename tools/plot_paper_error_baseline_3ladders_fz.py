#!/usr/bin/env python3
"""Baseline + 3-ladder error comparison (FTZ).

Three panels (N = 32k / 40k / 65k), four epsilon clusters per panel, four bars
per cluster in this order:

  1. Baseline IEEE     static per-matrix format, no per-tile decision
  2. IEEE ladder       FP8_E4M3 -> FP16 -> FP32 -> FP64
  3. MX native         MXFP4 -> MXFP8 -> FP16 (plain IEEE) -> FP32 -> FP64
                       == "MX staircase on GB200" (no shared scale at the FP16 rung)
  4. Full Ladder MX    MXFP4 -> MXFP8 -> MXFP16 -> FP32 -> FP64

All under MX_UNDERFLOW_MODE=fz (FTZ) for apples-to-apples comparison.
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
N_PANELS = [32768, 40960, 65536]
# (N, nb) — nb is the per-bin tile size we actually swept at.
PANEL_META = {
    32768: {"N_label": "N = 32{,}768", "nb": 2048},
    40960: {"N_label": "N = 40{,}960", "nb": 2048},
    65536: {"N_label": "N = 65{,}536", "nb": 4096},
}

# (label, source-key, hatch, color)
BARS = [
    ("Baseline IEEE  (no MX; static format per matrix)",
     "baseline",  "..", "#4e79a7"),
    ("IEEE ladder  (FP8 E4M3 $\\to$ FP16 $\\to$ FP32 $\\to$ FP64)",
     "ieee",      "||", "#a0a0a0"),
    ("MX native  (MXFP4 $\\to$ MXFP8 $\\to$ FP16 plain $\\to$ FP32 $\\to$ FP64)",
     "staircase", "++", "#F0E442"),
    ("Full ladder MX  (MXFP4 $\\to$ MXFP8 $\\to$ MXFP16 $\\to$ FP32 $\\to$ FP64)",
     "full_mx",   "xx", "#b07aa1"),
]


def load_master(csv_path):
    """Master CSV: baseline + IEEE ladder + Full Ladder MX rows.
    Keys are (n, nb, eps) so panels can pick a specific nb and never silently mix."""
    by_sweep = {"baseline": {}, "ieee": {}, "full_mx": {}}
    mapping = {
        "requant_baseline_fp8_subnormal_gt20k":  "baseline",
        "requant_ladder_ieee_gt20k":             "ieee",
        "requant_ladder_scaled_vec1d32_gt20k":   "full_mx",
    }
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            key = mapping.get(r.get("sweep", ""))
            if key is None:
                continue
            try:
                n = int(r["n"]); nb = int(r["nb"]); eps = r["source_epsilon"]
                rel = float(r["rel_factor_error"])
            except (KeyError, ValueError, TypeError):
                continue
            by_sweep[key][(n, nb, eps)] = rel
    return by_sweep


def load_staircase(csv_paths):
    """Standalone CSV(s) from staircase sweeps. Filter underflow=fz.
    Multiple paths are merged (so e.g. nb=2048 + nb=4096 sweeps can both be read)."""
    out = {}
    for csv_path in csv_paths:
        p = Path(csv_path)
        if not p.exists():
            continue
        with p.open() as f:
            for r in csv.DictReader(f):
                if r.get("underflow") != "fz":
                    continue
                try:
                    n = int(r["n"]); nb = int(r["nb"]); eps = r["source_epsilon"]
                    rel = float(r["rel_factor_error"])
                except (KeyError, ValueError, TypeError):
                    continue
                out[(n, nb, eps)] = rel
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master-csv",
                    default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--staircase-csv",
                    default="/home/abduraa/MX_project/ooc-cholesky/mx_staircase_vec1d32_gt20k/results.csv")
    ap.add_argument("--out", required=True,
                    help="Output base path; .pdf is appended.")
    args = ap.parse_args()

    master = load_master(args.master_csv)
    stair  = load_staircase(args.staircase_csv)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5), sharey=True)

    n_eps   = len(EPS_ORDER)
    n_bars  = len(BARS)
    group_w = 0.84
    bar_w   = group_w / n_bars
    x_centres = list(range(n_eps))

    y_min_global, y_max_global = float("inf"), 0.0

    for ax, N in zip(axes, N_PANELS):
        for bi, (label, key, hatch, color) in enumerate(BARS):
            xs, ys = [], []
            for ci, eps in zip(x_centres, EPS_ORDER):
                if key == "staircase":
                    v = stair.get((N, eps))
                else:
                    v = master[key].get((N, eps))
                v = v if v is not None else 0.0
                xs.append(ci - group_w/2 + (bi + 0.5) * bar_w)
                ys.append(v)
            ax.bar(xs, ys, bar_w * 0.92,
                   color=color, edgecolor="black", linewidth=0.4,
                   hatch=hatch, label=label)
            for x, y in zip(xs, ys):
                if y > 0:
                    y_min_global = min(y_min_global, y)
                    y_max_global = max(y_max_global, y)
                    ax.text(x, y * 1.15, f"{y:.1e}",
                            ha="center", va="bottom", fontsize=7.0,
                            rotation=90, color="#222")

        # Dashed grey reference line per eps column at the target.
        for c, eps in zip(x_centres, EPS_ORDER):
            ev = float(eps)
            ax.hlines(ev, c - group_w/2 - 0.08, c + group_w/2 + 0.08,
                      colors="gray", linestyles="--", linewidth=0.9,
                      alpha=0.7, zorder=0)

        ax.set_xticks(x_centres)
        ax.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                            for e in EPS_ORDER])
        ax.set_yscale("log")
        ax.set_axisbelow(True)

        # Panel annotation: N and nb in the top-right corner.
        meta = PANEL_META[N]
        ax.text(0.985, 0.965,
                rf"${meta['N_label']}$" "\n" rf"$n_b = {meta['nb']}$",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=10.5,
                bbox=dict(boxstyle="round,pad=0.28",
                          facecolor="white", edgecolor="#888",
                          linewidth=0.5, alpha=0.92))

    axes[0].set_ylabel(r"Relative factorization error  $\|LL^\top - A\|_\infty / \|A\|_\infty$  (log)")

    if y_min_global < float("inf"):
        for ax in axes:
            ax.set_ylim(y_min_global / 3, y_max_global * 15)

    handles = [Patch(facecolor=color, edgecolor="black", hatch=hatch, label=label)
               for label, _, hatch, color in BARS]
    fig.legend(handles=handles,
               loc="lower center", bbox_to_anchor=(0.5, -0.10),
               ncol=2, frameon=True, framealpha=0.95,
               handletextpad=0.6, columnspacing=2.0, labelspacing=0.35,
               borderpad=0.4)

    fig.tight_layout(rect=[0, 0.02, 1, 1.0])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    p = out.with_suffix(".pdf")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    print(p)
    plt.close(fig)


if __name__ == "__main__":
    main()
