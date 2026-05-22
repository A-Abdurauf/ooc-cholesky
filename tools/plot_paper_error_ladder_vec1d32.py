#!/usr/bin/env python3
"""Paper-grade error figure for ladder mode at vec1D-32 granularity.

Single figure, one row of subplots, one subplot per N (default 32k / 40k / 65k).
Each subplot has 4 epsilon clusters with 3 bars:

    1. Plain FP8 baseline   (subnormals allowed; no MX scaling)
    2. IEEE ladder          (FP8 E4M3 -> FP16 -> FP32 -> FP64, no MX)
    3. Full ladder vec1D-32 (MXFP4 -> MXFP8 E4M3 -> MXFP16 -> FP32 -> FP64)

Bars 1 and 2 come from the main requant_gt20k_memory.csv.  Bar 3 is read from
the same file for N != 32768; for N == 32768, the freshly re-run
ladder_rerun_32k/results.csv (GU = gradual-underflow / denorms preserved) is
used so the underflow handling matches the IEEE-ladder bar, which keeps
subnormals natively.
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
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
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

BAR_FP8_BASELINE = (
    "Plain FP8 baseline (no MX; subnormal-allowed)", "..", "#4e79a7",
)
BAR_IEEE_LADDER  = (
    "IEEE ladder (FP8 E4M3 -> FP16 -> FP32 -> FP64; no MX)", "//", "#a0a0a0",
)
BAR_FULL_LADDER  = (
    "Full ladder vec1D-32 (MXFP4 -> MXFP8 -> MXFP16 -> FP32 -> FP64)",
    "xx", "#b07aa1",
)

BARS = [BAR_FP8_BASELINE, BAR_IEEE_LADDER, BAR_FULL_LADDER]
SWEEPS_MAIN = [
    "requant_baseline_fp8_subnormal_gt20k",
    "requant_ladder_ieee_gt20k",
    "requant_ladder_scaled_vec1d32_gt20k",
]

NB_BY_N = {32768: 2048, 40960: 4096, 65536: 4096}


def load_main(csv_path):
    rows = {}
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            try:
                n  = int(r["n"])
                nb = int(r["nb"])
            except (TypeError, ValueError):
                continue
            if NB_BY_N.get(n) is not None and nb != NB_BY_N[n]:
                continue
            rows[(r["sweep"], n, r["source_epsilon"])] = r
    return rows


def load_rerun(csv_path):
    rows = {}
    p = Path(csv_path)
    if not p.exists():
        return rows
    with p.open() as f:
        for r in csv.DictReader(f):
            try:
                n = int(r["n"])
            except (TypeError, ValueError):
                continue
            rows[(r["sweep"], n, r["source_epsilon"])] = r
    return rows


def err(row):
    if row is None:
        return 0.0
    try:
        return float(row["rel_factor_error"])
    except (TypeError, ValueError, KeyError):
        return 0.0


def _fmt_n(n):
    return f"{n // 1024}k" if n % 1024 == 0 else str(n)


def plot_subplot(ax, n, main_rows, rerun_rows, use_rerun_for_ladder, show_ylabel):
    n_eps  = len(EPS_ORDER)
    n_bars = len(BARS)
    group_w = 0.78
    bar_w   = group_w / n_bars
    x_centres = list(range(n_eps))

    y_min, y_max = float("inf"), 0.0

    for bi, ((label, hatch, color), sweep) in enumerate(zip(BARS, SWEEPS_MAIN)):
        xs, ys = [], []
        for ci, eps in zip(x_centres, EPS_ORDER):
            row = main_rows.get((sweep, n, eps))
            if (use_rerun_for_ladder
                    and sweep == "requant_ladder_scaled_vec1d32_gt20k"):
                r_rerun = rerun_rows.get(("ladder_full_gu", n, eps))
                if r_rerun is not None:
                    row = r_rerun
            v = err(row)
            xs.append(ci - group_w/2 + (bi + 0.5) * bar_w)
            ys.append(v)
        ax.bar(xs, ys, bar_w * 0.92,
               color=color, edgecolor="black", linewidth=0.4,
               hatch=hatch, label=label)
        for x, y in zip(xs, ys):
            if y > 0:
                y_min = min(y_min, y); y_max = max(y_max, y)
                ax.text(x, y * 1.15, f"{y:.1e}",
                        ha="center", va="bottom", fontsize=7.2, rotation=90,
                        color="#222")

    for c, eps in zip(x_centres, EPS_ORDER):
        ev = float(eps)
        ax.hlines(ev, c - group_w/2 - 0.06, c + group_w/2 + 0.06,
                  colors="gray", linestyles="--", linewidth=0.9, alpha=0.7, zorder=0)
        y_min = min(y_min, ev)

    ax.set_xticks(x_centres)
    ax.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                        for e in EPS_ORDER])
    ax.set_yscale("log")
    if show_ylabel:
        ax.set_ylabel(r"Relative factorization error  $\|LL^\top - A\|_\infty / \|A\|_\infty$  (log)")
    ax.set_axisbelow(True)

    if y_min < float("inf"):
        ax.set_ylim(y_min / 4, max(y_max, max(float(e) for e in EPS_ORDER)) * 20)

    nb = NB_BY_N.get(n, "?")
    ax.text(0.985, 0.965, f"N = {_fmt_n(n)}\nNB = {nb}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=10.5,
            bbox=dict(boxstyle="round,pad=0.32",
                      facecolor="white", edgecolor="#888", linewidth=0.8,
                      alpha=0.92))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv",
                    default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--rerun-csv",
                    default="/home/abduraa/MX_project/ooc-cholesky/ladder_rerun_32k/results.csv")
    ap.add_argument("--ns", nargs="+", type=int, default=[32768, 40960, 65536],
                    help="Matrix sizes to emit, one subplot each (left to right).")
    ap.add_argument("--out",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures/paper_error_ladder_vec1d32")
    ap.add_argument("--rerun-for-n", type=int, default=32768,
                    help="N for which the rerun CSV is preferred for the ladder bar.")
    args = ap.parse_args()

    main_rows  = load_main(args.csv)
    rerun_rows = load_rerun(args.rerun_csv)

    n_panels = len(args.ns)
    fig, axes = plt.subplots(1, n_panels,
                             figsize=(4.4 * n_panels + 0.8, 5.0),
                             sharey=True)
    if n_panels == 1:
        axes = [axes]

    for ax, n in zip(axes, args.ns):
        plot_subplot(ax, n, main_rows, rerun_rows,
                     use_rerun_for_ladder=(n == args.rerun_for_n),
                     show_ylabel=(ax is axes[0]))

    handles = [Patch(facecolor=color, edgecolor="black", hatch=hatch, label=label)
               for label, hatch, color in BARS]
    fig.legend(handles=handles,
               loc="lower center", bbox_to_anchor=(0.5, -0.04),
               ncol=3, frameon=True, framealpha=0.95,
               title="Bars within each cluster (left -> right)   .   grey dashed = source $\\varepsilon$",
               handletextpad=0.6, columnspacing=2.0, labelspacing=0.35,
               borderpad=0.4)

    fig.tight_layout(rect=[0, 0.03, 1, 1.0])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    p = out.with_suffix(".pdf")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    print(p)
    plt.close(fig)


if __name__ == "__main__":
    main()
