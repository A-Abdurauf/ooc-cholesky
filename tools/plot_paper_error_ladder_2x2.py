#!/usr/bin/env python3
"""Paper-ready error figure for the ladder story, eps-as-subplot layout.

Layout (mirrors paper_memory_nb4096_sweep_by_eps_2x2):
  - 1 figure, 4 subplots in one row (one per epsilon)
  - x-axis of each subplot: 3 N groups (32k, 40k, 65k)
  - each N group: 4 bars in this order:
      1. Baseline                  -- requant_baseline_fp8_subnormal_gt20k (main CSV, NB-matched)
      2. Ladder IEEE               -- requant_ladder_ieee_gt20k (main CSV, NB-matched)
      3. Ladder MX  (no MXFP16)    -- mx_staircase_vec1d32_gu  (GU underflow,
                                     32k+65k from mx_staircase_vec1d32_gt20k,
                                     40k from mx_staircase_40k_nb4096_gu)
      4. Ladder MX+MXFP16          -- 32k: ladder_full_gu  (rerun, GU)
                                     65k: ladder_full_vec1d_gu  (gu_missing, GU)
                                     40k: requant_ladder_scaled_vec1d32_gt20k
                                          (main CSV, FZ -- no GU run yet at NB=4096)

Per-bar NB at each N (so the comparison is fair within each N group):
  32k -> NB=2048; 40k -> NB=4096; 65k -> NB=4096.

Y-axis: relative factorization error (log).
Paper-ready styling: 9pt body, panel labels (a)/(b)/(c)/(d), hatches for B&W,
combined legend below.  PDF only.
"""
import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch


PAPER_RC = {
    "font.family": "serif",
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "legend.title_fontsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.4,
}


EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]
EPS_DISPLAY = {"1e-5": r"10^{-5}", "1e-6": r"10^{-6}",
               "1e-7": r"10^{-7}", "1e-8": r"10^{-8}"}
N_ORDER = [32768, 40960, 65536]
NB_BY_N = {32768: 2048, 40960: 4096, 65536: 4096}


MAIN_CSV       = "/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv"
RERUN_CSV      = "/home/abduraa/MX_project/ooc-cholesky/ladder_rerun_32k/results.csv"
STAIRCASE_CSV  = "/home/abduraa/MX_project/ooc-cholesky/mx_staircase_vec1d32_gt20k/results.csv"
STAIRCASE40_CSV = "/home/abduraa/MX_project/ooc-cholesky/mx_staircase_40k_nb4096_gu/results.csv"
GU_MISSING_CSV = "/home/abduraa/MX_project/ooc-cholesky/ladder_full_gu_missing/results.csv"


# (label, color, hatch).  Hatch convention matches the memory figure.
BARS = [
    ("Baseline",          "#4e79a7", ""),
    ("Ladder IEEE",       "#a0a0a0", ""),
    ("Ladder MX",         "#e07b39", "xx"),   # no MXFP16 (staircase)
    ("Ladder MX+MXFP16",  "#b07aa1", "//"),
]


def load_csv(path):
    rows = []
    p = Path(path)
    if not p.exists():
        return rows
    with p.open() as f:
        for r in csv.DictReader(f):
            rows.append(r)
    return rows


def err_lookup(rows_main, rows_rerun, rows_stair, rows_stair40, rows_gu,
               n, eps):
    """Return the 4-bar values (in BARS order) for one (N, eps) panel."""
    out = [float("nan")] * 4
    nb = NB_BY_N[n]

    # 1. Baseline.
    for r in rows_main:
        if (r["sweep"] == "requant_baseline_fp8_subnormal_gt20k"
                and int(r["n"]) == n and int(r["nb"]) == nb
                and r["source_epsilon"] == eps):
            out[0] = float(r["rel_factor_error"]); break

    # 2. Ladder IEEE.
    for r in rows_main:
        if (r["sweep"] == "requant_ladder_ieee_gt20k"
                and int(r["n"]) == n and int(r["nb"]) == nb
                and r["source_epsilon"] == eps):
            out[1] = float(r["rel_factor_error"]); break

    # 3. Ladder MX (no MXFP16) -- staircase vec1d32 GU.
    #    40k at NB=4096 -> dedicated new sweep CSV; other N's -> the original
    #    staircase CSV (if still present).
    stair_pool = rows_stair40 if n == 40960 else rows_stair
    for r in stair_pool:
        if (r["sweep"] == "mx_staircase_vec1d32_gu"
                and int(r["n"]) == n and int(r["nb"]) == nb
                and r["source_epsilon"] == eps):
            out[2] = float(r["rel_factor_error"]); break

    # 4. Ladder MX+MXFP16.  Fallback chain (GU preferred):
    #    32k: rerun ladder_full_gu  ->  main CSV (FZ)
    #    65k: ladder_full_vec1d_gu  ->  main CSV (FZ)
    #    40k: main CSV (FZ -- no GU run at NB=4096)
    chain = []
    if n == 32768:
        chain.append((rows_rerun, "ladder_full_gu"))
    elif n == 65536:
        chain.append((rows_gu, "ladder_full_vec1d_gu"))
    chain.append((rows_main, "requant_ladder_scaled_vec1d32_gt20k"))
    for pool, sweep in chain:
        if np.isfinite(out[3]):
            break
        for r in pool:
            if (r["sweep"] == sweep
                    and int(r["n"]) == n and int(r["nb"]) == nb
                    and r["source_epsilon"] == eps):
                out[3] = float(r["rel_factor_error"]); break
    return out


def _fmt_n(n):
    return f"{n // 1024}k"


def draw_subplot(ax, eps, datasets):
    rows_main, rows_rerun, rows_stair, rows_stair40, rows_gu = datasets
    n_groups = len(N_ORDER)
    n_bars = len(BARS)
    group_w = 0.78
    bar_w   = group_w / n_bars
    x_centres = list(range(n_groups))

    y_min, y_max = float("inf"), 0.0
    for ni, N in enumerate(N_ORDER):
        vals = err_lookup(rows_main, rows_rerun, rows_stair, rows_stair40,
                          rows_gu, N, eps)
        for bi, ((label, color, hatch), v) in enumerate(zip(BARS, vals)):
            x = ni - group_w/2 + (bi + 0.5) * bar_w
            if not np.isfinite(v) or v <= 0:
                continue
            ax.bar([x], [v], bar_w * 0.92,
                   color=color, edgecolor="black", linewidth=0.3,
                   hatch=hatch)
            y_min = min(y_min, v); y_max = max(y_max, v)

    # eps reference line.
    ev = float(eps)
    ax.axhline(ev, color="gray", linestyle="--", linewidth=0.7, alpha=0.65,
               zorder=0)
    y_min = min(y_min, ev)

    # Light vertical separators between N groups.
    for i in range(1, n_groups):
        ax.axvline(i - 0.5, color="#bbbbbb", linewidth=0.5, alpha=0.6, zorder=0)

    ax.set_xticks(x_centres)
    ax.set_xticklabels([_fmt_n(N) for N in N_ORDER])
    ax.set_xlim(-0.6, n_groups - 0.4)
    ax.set_yscale("log")
    ax.set_axisbelow(True)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, alpha=0.25, linewidth=0.4, which="both")
    return y_min, y_max


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures/paper_error_ladder_2x2")
    args = ap.parse_args()

    rows_main    = load_csv(MAIN_CSV)
    rows_rerun   = load_csv(RERUN_CSV)
    rows_stair   = load_csv(STAIRCASE_CSV)
    rows_stair40 = load_csv(STAIRCASE40_CSV)
    rows_gu      = load_csv(GU_MISSING_CSV)
    datasets = (rows_main, rows_rerun, rows_stair, rows_stair40, rows_gu)

    with mpl.rc_context(PAPER_RC):
        fig, axes = plt.subplots(1, len(EPS_ORDER),
                                 figsize=(7.16, 2.7),
                                 sharey=True)

        y_min_all, y_max_all = float("inf"), 0.0
        for ai, (ax, eps) in enumerate(zip(axes, EPS_ORDER)):
            y_min, y_max = draw_subplot(ax, eps, datasets)
            y_min_all = min(y_min_all, y_min)
            y_max_all = max(y_max_all, y_max)
            ax.text(0.02, 0.965,
                    f"({chr(ord('a') + ai)})  $\\varepsilon = {EPS_DISPLAY[eps]}$",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=8.5, fontweight="bold")

        if y_min_all < float("inf"):
            for ax in axes:
                ax.set_ylim(y_min_all / 4, y_max_all * 20)

        for ax in axes:
            ax.set_xlabel("Matrix size  $N$")
        axes[0].set_ylabel(r"Relative factorization error  (log)")

        bar_handles = [Patch(facecolor=c, edgecolor="black", hatch=h,
                             label=f"{i+1}. {lbl}")
                       for i, (lbl, c, h) in enumerate(BARS)]
        eps_handle = plt.Line2D([0], [0], color="gray", linestyle="--",
                                linewidth=0.9, label="source $\\varepsilon$")
        fig.legend(handles=bar_handles + [eps_handle],
                   loc="lower center", bbox_to_anchor=(0.5, -0.06),
                   ncol=5, frameon=True, framealpha=0.95,
                   title="Bars (left $\\rightarrow$ right within each $N$ group)",
                   handletextpad=0.5, columnspacing=1.4,
                   labelspacing=0.3, borderpad=0.35)

        fig.tight_layout(rect=[0, 0.10, 1, 1.0])

        out = Path(args.out).with_suffix(".pdf")
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(out)
        plt.close(fig)


if __name__ == "__main__":
    main()
