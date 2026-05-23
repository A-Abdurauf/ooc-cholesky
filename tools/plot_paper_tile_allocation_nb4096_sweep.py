#!/usr/bin/env python3
"""Tile-allocation sweep figure across all matrix sizes (NB=4096).

Paper-ready 2x2 layout mirroring paper_memory_nb4096_sweep_by_eps_2x2.pdf:
  - 4 subplots, one per epsilon (default by=eps) -- can also do by=mode
  - x-axis: 7 N values (20k -> 120k); within each N, 4 stacked bars in the
    order Baseline / Ladder IEEE / Ladder MX / Ladder MX+MXFP16.
  - Stack: tile counts (not GB) by format, same color scheme as the memory
    figure.  The "Scale" bucket is dropped since it is memory-only.

The bars sum to M*(M+1)/2 lower-triangular tiles -- so the per-N totals are
identical across the 4 modes, which makes the format-shift story crystal
clear (same total tiles, different colour mix).

Same FP16 -> MXFP16 remap as the memory figure: for ladder_mx_staircase_mxfp16
the binary logs "FP16" for tiles that are actually MXFP16-with-scale-meta.

PDF only.
"""
import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# Reuse parsing + bookkeeping from the memory script -- guarantees the figures
# stay consistent (same N range, modes, colors, hatches, paper rcParams).
import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_paper_memory_nb4096_sweep import (  # noqa: E402
    parse_log, LOG_PATH,
    EPS_ORDER, EPS_DISPLAY, N_ORDER, N_LABEL, _xlabel,
    MODE_BARS, MODE_PANEL_LABEL_SHORT,
    STACK, STACK_LABEL, STACK_HATCH, PAPER_RC,
)


def tile_counts_by_format(counts):
    """Return {format_label: tile_count} matching the STACK keys (no Scale)."""
    out = {fmt: 0 for fmt, _ in STACK}
    for fmt in ("FP64", "FP32", "FP16", "MXFP16",
                "FP8_plain", "MXFP8", "MXFP4"):
        out[fmt] = counts.get(fmt, 0)
    return out


def _draw_paper_subplot(ax, by, key, data, ns):
    if by == "eps":
        slots = [(lbl, mk) for lbl, mk in MODE_BARS]
        slot_key = lambda slot: (slot[1], None, key)
    else:
        slots = [(EPS_DISPLAY[e], e) for e in EPS_ORDER]
        slot_key = lambda slot: (key, None, slot[1])

    n_slots = len(slots)
    group_w = 0.78
    bar_w   = group_w / n_slots
    x_centres = list(range(len(ns)))

    drawn_fmts = set()
    panel_max = 0
    for ni, N in enumerate(ns):
        for bi, slot in enumerate(slots):
            x = ni - group_w/2 + (bi + 0.5) * bar_w
            mkey = slot_key(slot)
            counts = data.get((mkey[0], N, mkey[2]))
            if counts is None:
                continue
            tc = tile_counts_by_format(counts)
            bottom = 0
            for fmt, color in STACK:
                if fmt == "Scale":
                    continue
                v = tc.get(fmt, 0)
                if v <= 0:
                    continue
                drawn_fmts.add(fmt)
                ax.bar([x], [v], bar_w * 0.92,
                       bottom=[bottom],
                       color=color, edgecolor="black",
                       linewidth=0.25, hatch=STACK_HATCH.get(fmt, ""))
                bottom += v
            panel_max = max(panel_max, bottom)

    # Light vertical separators between N groups.
    for i in range(1, len(ns)):
        ax.axvline(i - 0.5, color="#bbbbbb", linewidth=0.5, alpha=0.6, zorder=0)

    ax.set_xticks(x_centres)
    ax.set_xticklabels([_xlabel(N) for N in ns])
    ax.set_xlim(-0.6, len(ns) - 0.4)
    ax.set_axisbelow(True)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, alpha=0.25, linewidth=0.4)
    return drawn_fmts, panel_max


def render_paper_2x2(data, by, out_path):
    if by == "eps":
        panel_keys   = EPS_ORDER
        panel_titles = [f"$\\varepsilon = {EPS_DISPLAY[e][1:-1]}$" for e in EPS_ORDER]
        right_title  = "Bars (left $\\rightarrow$ right within each $N$)"
        right_items  = [f"{i+1}. {lbl}" for i, (lbl, _) in enumerate(MODE_BARS)]
    elif by == "mode":
        panel_keys   = [mk for _, mk in MODE_BARS]
        panel_titles = [MODE_PANEL_LABEL_SHORT[mk] for mk in panel_keys]
        right_title  = "Bars (left $\\rightarrow$ right within each $N$)"
        right_items  = [f"{i+1}. $\\varepsilon = {EPS_DISPLAY[e][1:-1]}$"
                        for i, e in enumerate(EPS_ORDER)]
    else:
        raise ValueError(by)

    with mpl.rc_context(PAPER_RC):
        fig, axes = plt.subplots(2, 2, figsize=(7.16, 3.7),
                                 sharex=True, sharey=True)
        axes_flat = axes.flatten()

        drawn_all = set()
        global_max = 0
        for ai, (ax, key, title) in enumerate(zip(axes_flat, panel_keys, panel_titles)):
            drawn, panel_max = _draw_paper_subplot(ax, by, key, data, N_ORDER)
            drawn_all |= drawn
            global_max = max(global_max, panel_max)
            ax.text(0.02, 0.965,
                    f"({chr(ord('a') + ai)})  {title}",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=8.5, fontweight="bold")

        if global_max > 0:
            for ax in axes_flat:
                ax.set_ylim(0, global_max * 1.10)

        # Lock figure to top, reserve bottom band for legends -- matches memory fig.
        fig.subplots_adjust(left=0.10, right=0.985, top=0.97,
                            bottom=0.28,
                            wspace=0.06, hspace=0.18)

        for ax in axes[-1, :]:
            ax.set_xlabel(r"Matrix size  $N$")
        for ax in axes[:, 0]:
            ax.set_ylabel("Tile count")

        # Keep FP16 + MXFP16 adjacent in the legend.
        legend_drawn = set(drawn_all)
        if "FP16" in legend_drawn or "MXFP16" in legend_drawn:
            legend_drawn.add("FP16")
            legend_drawn.add("MXFP16")
        fmt_handles = [Patch(facecolor=c, edgecolor="black",
                             hatch=STACK_HATCH.get(f, ""),
                             label=STACK_LABEL[f])
                       for f, c in STACK
                       if f in legend_drawn and f != "Scale"]
        bar_handles = [Patch(facecolor="white", edgecolor="black", label=lbl)
                       for lbl in right_items]

        leg1 = fig.legend(handles=fmt_handles,
                          loc="upper center", bbox_to_anchor=(0.27, 0.18),
                          ncol=4, frameon=True, framealpha=0.95,
                          title="Tile format (bottom $\\rightarrow$ top of stack)",
                          handletextpad=0.5, columnspacing=1.2,
                          labelspacing=0.25, borderpad=0.3)
        fig.add_artist(leg1)
        fig.legend(handles=bar_handles,
                   loc="upper center", bbox_to_anchor=(0.78, 0.18),
                   ncol=2, frameon=True, framealpha=0.95,
                   title=right_title,
                   handletextpad=0.5, columnspacing=1.2,
                   labelspacing=0.25, borderpad=0.3)

        out = Path(out_path).with_suffix(".pdf")
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(out)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", default=LOG_PATH)
    ap.add_argument("--by", choices=["eps", "mode", "both"], default="both")
    ap.add_argument("--out-dir",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures")
    args = ap.parse_args()

    data = parse_log(args.log)
    targets = ["eps", "mode"] if args.by == "both" else [args.by]
    for by in targets:
        out = Path(args.out_dir) / f"paper_tile_allocation_nb4096_by_{by}_2x2"
        render_paper_2x2(data, by, out)


if __name__ == "__main__":
    main()
