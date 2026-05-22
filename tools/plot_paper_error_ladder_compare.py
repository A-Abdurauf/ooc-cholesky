#!/usr/bin/env python3
"""Side-by-side ladder-mode error comparison: tile vs vec1D-32.

Reuses the per-subplot drawing from plot_paper_error_ladder.py.  Produces:

  - For each N: one standalone figure with two subplots side-by-side, left =
    tile granularity, right = vec1D-32 granularity.  N / NB shown in each
    subplot's top-right corner; no figure title.
  - One combined figure: 2 rows x len(NS) cols (row 1 = tile, row 2 = vec1D-32).

The vec1D-32 bar at N=32768 prefers the ladder_rerun (GU underflow); the tile
bar always comes from the main CSV.
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt

from plot_paper_error_ladder import (
    GRAN, NB_BY_N, draw_subplot, legend_handles,
    load_main, load_rerun, _fmt_n,
)


def annotate(ax, n, gran_key):
    nb = NB_BY_N.get(n, "?")
    title_gran = GRAN[gran_key][5]
    ax.text(0.985, 0.965, f"N = {_fmt_n(n)}\nNB = {nb}\n{title_gran}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=10.0,
            bbox=dict(boxstyle="round,pad=0.32",
                      facecolor="white", edgecolor="#888", linewidth=0.8,
                      alpha=0.92))


def plot_standalone(n, gran_keys, main_rows, rerun_rows, out_path, rerun_for_n):
    fig, axes = plt.subplots(1, len(gran_keys),
                             figsize=(4.4 * len(gran_keys) + 0.6, 5.0),
                             sharey=True)
    if len(gran_keys) == 1:
        axes = [axes]

    last_bars = None
    for ax, gk in zip(axes, gran_keys):
        last_bars = draw_subplot(ax, n, gk, main_rows, rerun_rows,
                                 use_rerun=(n == rerun_for_n and gk == "vec1d32"),
                                 show_ylabel=(ax is axes[0]),
                                 value_fontsize=8)
        annotate(ax, n, gk)

    # Build a union legend across granularities (FP8/IEEE bars are shared;
    # MX bars differ per granularity).
    bars_union = []
    seen = set()
    for gk in gran_keys:
        for ax_used in [None]:  # dummy
            pass
        mx_label, mx_hatch, mx_color, *_ = GRAN[gk]
        # The first call to draw_subplot returned a `bars` list, but we want
        # to assemble a deterministic order: FP8, IEEE, then each MX bar.
    # Re-derive:
    from plot_paper_error_ladder import BAR_FP8_BASELINE, BAR_IEEE_LADDER
    bars_union = [BAR_FP8_BASELINE, BAR_IEEE_LADDER]
    for gk in gran_keys:
        mx_label, mx_hatch, mx_color, *_ = GRAN[gk]
        bars_union.append((mx_label, mx_hatch, mx_color))

    fig.legend(handles=legend_handles(bars_union),
               loc="lower center", bbox_to_anchor=(0.5, -0.06),
               ncol=2, frameon=True, framealpha=0.95,
               title="Bars within each cluster (left -> right)   .   grey dashed = source $\\varepsilon$",
               handletextpad=0.6, columnspacing=2.0, labelspacing=0.35,
               borderpad=0.4)
    fig.tight_layout(rect=[0, 0.04, 1, 1.0])

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    p = out.with_suffix(".pdf")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    print(p)
    plt.close(fig)


def plot_grid(ns, gran_keys, main_rows, rerun_rows, out_path, rerun_for_n):
    n_rows = len(gran_keys)
    n_cols = len(ns)
    fig, axes = plt.subplots(n_rows, n_cols,
                             figsize=(4.4 * n_cols + 0.6, 4.6 * n_rows),
                             sharey=True, sharex=True, squeeze=False)

    for ri, gk in enumerate(gran_keys):
        for ci, n in enumerate(ns):
            ax = axes[ri][ci]
            draw_subplot(ax, n, gk, main_rows, rerun_rows,
                         use_rerun=(n == rerun_for_n and gk == "vec1d32"),
                         show_ylabel=(ci == 0))
            annotate(ax, n, gk)

    from plot_paper_error_ladder import BAR_FP8_BASELINE, BAR_IEEE_LADDER
    bars_union = [BAR_FP8_BASELINE, BAR_IEEE_LADDER]
    for gk in gran_keys:
        mx_label, mx_hatch, mx_color, *_ = GRAN[gk]
        bars_union.append((mx_label, mx_hatch, mx_color))

    fig.legend(handles=legend_handles(bars_union),
               loc="lower center", bbox_to_anchor=(0.5, -0.03),
               ncol=2, frameon=True, framealpha=0.95,
               title="Bars within each cluster (left -> right)   .   grey dashed = source $\\varepsilon$",
               handletextpad=0.6, columnspacing=2.0, labelspacing=0.35,
               borderpad=0.4)
    fig.tight_layout(rect=[0, 0.03, 1, 1.0])

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    p = out.with_suffix(".pdf")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    print(p)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv",
                    default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--rerun-csv",
                    default="/home/abduraa/MX_project/ooc-cholesky/ladder_rerun_32k/results.csv")
    ap.add_argument("--ns", nargs="+", type=int, default=[32768, 40960, 65536])
    ap.add_argument("--granularities", nargs="+", default=["tile", "vec1d32"],
                    choices=list(GRAN))
    ap.add_argument("--rerun-for-n", type=int, default=32768)
    ap.add_argument("--out-dir",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures")
    ap.add_argument("--layout", choices=["separate", "combined", "both"],
                    default="both")
    args = ap.parse_args()

    main_rows  = load_main(args.csv)
    rerun_rows = load_rerun(args.rerun_csv)
    stem = "_vs_".join(args.granularities)

    if args.layout in ("separate", "both"):
        for n in args.ns:
            out = Path(args.out_dir) / f"paper_error_ladder_{stem}_N{n}"
            plot_standalone(n, args.granularities, main_rows, rerun_rows, out,
                            rerun_for_n=args.rerun_for_n)

    if args.layout in ("combined", "both"):
        out = Path(args.out_dir) / f"paper_error_ladder_{stem}_grid"
        plot_grid(args.ns, args.granularities, main_rows, rerun_rows, out,
                  rerun_for_n=args.rerun_for_n)


if __name__ == "__main__":
    main()
