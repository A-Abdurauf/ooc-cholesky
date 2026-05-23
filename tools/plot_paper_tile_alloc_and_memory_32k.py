#!/usr/bin/env python3
"""Two-row paper figure for N=32k: tile allocation (top) + memory (bottom).

Row 1 mirrors paper_tile_allocation_32k.pdf -- per-format tile counts.
Row 2 shows the same 4 bars but stacked by memory (GB) using lower-triangle
+ half-diagonal accounting (FP64 diagonals = nb*(nb+1)/2 elements; off-diag
= nb*nb; MX formats add 1-byte shared scale per 32 elements).

Same colors / hatches / fonts as the memory and tile-allocation figures.
PDF only.  N/NB/tile-total info goes in the figure caption.
"""
import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_paper_memory_nb4096_sweep import (  # noqa: E402
    EPS_ORDER, EPS_DISPLAY,
    STACK, STACK_LABEL, STACK_HATCH, PAPER_RC,
)
from plot_paper_tile_allocation_32k_2x2 import (  # noqa: E402
    N, NB, M, TILES_TOTAL, BARS, ALIASES,
    MAIN_CSV, RERUN_CSV, STAIR_CSV,
    load_csv, parse_breakdown, parse_counts_full,
    full_to_lower_tri, find, counts_for_bar,
)

# Bytes per element for each (memory-bearing) format + whether the format
# carries a 1-byte/32-elem shared scale.  Keys must match STACK.
FMT_BYTES = {
    "FP64":      (8.0, False),
    "FP32":      (4.0, False),
    "FP16":      (2.0, False),
    "MXFP16":    (2.0, True),
    "FP8_plain": (1.0, False),
    "MXFP8":     (1.0, True),
    "MXFP4":     (0.5, True),
}


def memory_breakdown_gb(counts):
    """Convert lower-tri tile counts to per-format GB + 'Scale' bucket.

    Diagonal accounting: the M FP64 tiles each store only nb*(nb+1)/2 elements;
    every other tile stores the full nb*nb.
    """
    tile_full = NB * NB
    tile_diag = NB * (NB + 1) // 2
    GB = 1.0 / (1024 ** 3)
    out = {fmt: 0.0 for fmt, _ in STACK}
    for fmt, (bpe, has_scale) in FMT_BYTES.items():
        cnt = counts.get(fmt, 0)
        if cnt <= 0:
            continue
        if fmt == "FP64":
            diag = min(M, cnt)
            off  = max(cnt - diag, 0)
            data_bytes = diag * tile_diag * bpe + off * tile_full * bpe
        else:
            data_bytes = cnt * tile_full * bpe
        out[fmt] = data_bytes * GB
        if has_scale:
            out["Scale"] += cnt * (tile_full // 32) * GB
    return out


def _draw_panel(ax, kind, rows_main, rows_rerun, rows_stair, drawn_fmts):
    """kind ∈ {'tile', 'mem'} -- chooses count or GB stack."""
    n_eps  = len(EPS_ORDER)
    n_bars = len(BARS)
    group_w = 0.82
    bar_w   = group_w / n_bars
    x_centres = list(range(n_eps))
    panel_max = 0.0
    for ci, eps in zip(x_centres, EPS_ORDER):
        for bi, (bar_lbl, src, sweep) in enumerate(BARS):
            x = ci - group_w/2 + (bi + 0.5) * bar_w
            counts = counts_for_bar(rows_main, rows_rerun, rows_stair,
                                    src, sweep, eps)
            if kind == "mem":
                stack_vals = memory_breakdown_gb(counts)
            else:
                stack_vals = {fmt: counts.get(fmt, 0) for fmt, _ in STACK}
            bottom = 0.0
            for fmt, color in STACK:
                if kind == "tile" and fmt == "Scale":
                    continue   # tile panel: no scale bucket
                v = stack_vals.get(fmt, 0)
                if v <= 0:
                    continue
                drawn_fmts.add(fmt)
                ax.bar([x], [v], bar_w * 0.92,
                       bottom=[bottom],
                       color=color, edgecolor="black",
                       linewidth=0.25, hatch=STACK_HATCH.get(fmt, ""))
                bottom += v
            if bottom > panel_max:
                panel_max = bottom

    for i in range(1, n_eps):
        ax.axvline(i - 0.5, color="#bbbbbb",
                   linewidth=0.5, alpha=0.6, zorder=0)
    ax.set_xticks(x_centres)
    ax.set_xticklabels([rf"$\varepsilon = {EPS_DISPLAY[e][1:-1]}$"
                        for e in EPS_ORDER])
    # Shrink the gap between the x-axis line and the eps tick labels.
    ax.tick_params(axis="x", pad=1.5)
    ax.set_xlim(-0.5, n_eps - 0.5)
    if panel_max > 0:
        ax.set_ylim(0, panel_max * 1.05)
    ax.set_axisbelow(True)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, alpha=0.25, linewidth=0.4)
    return panel_max


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures/paper_tile_alloc_and_memory_32k")
    args = ap.parse_args()

    rows_main  = load_csv(MAIN_CSV)
    rows_rerun = load_csv(RERUN_CSV)
    rows_stair = load_csv(STAIR_CSV)

    with mpl.rc_context(PAPER_RC):
        fig, (ax_tile, ax_mem) = plt.subplots(
            2, 1, figsize=(7.16, 3.2), sharex=True,
            gridspec_kw={"height_ratios": [1.0, 1.0]},
        )
        drawn = set()
        _draw_panel(ax_tile, "tile", rows_main, rows_rerun, rows_stair, drawn)
        _draw_panel(ax_mem,  "mem",  rows_main, rows_rerun, rows_stair, drawn)

        ax_tile.set_ylabel("Tile count")
        ax_mem.set_ylabel("Memory (GB)")

        # FP16 + MXFP16 forced adjacent in legend.
        legend_drawn = set(drawn)
        if "FP16" in legend_drawn or "MXFP16" in legend_drawn:
            legend_drawn.add("FP16"); legend_drawn.add("MXFP16")
        fmt_handles = [Patch(facecolor=c, edgecolor="black",
                             hatch=STACK_HATCH.get(f, ""),
                             label=STACK_LABEL[f])
                       for f, c in STACK if f in legend_drawn]
        bar_handles = [Patch(facecolor="white", edgecolor="black",
                             label=f"{i+1}. {lbl}")
                       for i, (lbl, _, _) in enumerate(BARS)]

        fig.subplots_adjust(left=0.095, right=0.985, top=0.97,
                            bottom=0.32, hspace=0.06)
        # Pull the legend up just below the eps tick band.
        leg1 = fig.legend(handles=fmt_handles,
                          loc="upper center", bbox_to_anchor=(0.27, 0.22),
                          ncol=4, frameon=True, framealpha=0.95,
                          title="Tile format (bottom $\\rightarrow$ top of stack)",
                          handletextpad=0.5, columnspacing=1.2,
                          labelspacing=0.25, borderpad=0.3)
        fig.add_artist(leg1)
        fig.legend(handles=bar_handles,
                   loc="upper center", bbox_to_anchor=(0.78, 0.22),
                   ncol=2, frameon=True, framealpha=0.95,
                   title="Bars (left $\\rightarrow$ right within each $\\varepsilon$)",
                   handletextpad=0.5, columnspacing=1.2,
                   labelspacing=0.25, borderpad=0.3)

        out = Path(args.out).with_suffix(".pdf")
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(out)
        plt.close(fig)


if __name__ == "__main__":
    main()
