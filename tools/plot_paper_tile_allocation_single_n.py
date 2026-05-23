#!/usr/bin/env python3
"""Paper-ready per-format tile allocation figure for a single N.

Single panel, 4 epsilon clusters x 4 stacked bars per cluster:
    1. Baseline
    2. Ladder IEEE
    3. Ladder MX  (no MXFP16, staircase)
    4. Ladder MX+MXFP16

Reads from recreate-ooc-chol/build/combined_perf_data.csv so the same script
works for any N in {20k, 32k, 40k, 64k, 80k, 100k, 120k}.  The combiner uses
NB=2048 for N=20480 and 32768, NB=4096 elsewhere.

FP16 -> MXFP16 remap for the Full MX ladder (the binary always logs the
second-highest tier as "FP16" even when the run used MXFP16).
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
    parse_log, LOG_PATH,
    EPS_DISPLAY, STACK, STACK_LABEL, STACK_HATCH, PAPER_RC,
)


EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]

# (display label, source-mode key in runs_all.log).
MODES = [
    ("Baseline",         "fp8"),
    ("Ladder IEEE",      "ladder_ieee_only"),
    ("Ladder MX",        "ladder_mx_staircase"),
    ("Ladder MX+MXFP16", "ladder_mx_staircase_mxfp16"),
]

# Default labels: 32k=32k, 40k=40k, 65536->64k, 81920->80k, 98304->100k, 122880->120k.
N_LABEL_DEFAULT = {
    20480: "20k", 32768: "32k", 40960: "40k",
    65536: "64k", 81920: "80k", 98304: "100k", 122880: "120k",
}


def pick(data, mode_log, N, eps):
    """Look up authoritative tile counts from the parsed runs_all.log dict.

    parse_log() already applies the FP16 -> MXFP16 remap for the full MX
    ladder mode (i.e. ladder_mx_staircase_mxfp16).
    """
    r = data.get((mode_log, N, eps))
    if r is None:
        return None
    return {
        "FP64":      r.get("FP64", 0),
        "FP32":      r.get("FP32", 0),
        "FP16":      r.get("FP16", 0),
        "MXFP16":    r.get("MXFP16", 0),
        "FP8_plain": r.get("FP8_plain", 0),
        "MXFP8":     r.get("MXFP8", 0),
        "MXFP4":     r.get("MXFP4", 0),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=65536,
                    choices=list(N_LABEL_DEFAULT))
    ap.add_argument("--out-dir",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures")
    args = ap.parse_args()

    data  = parse_log(LOG_PATH)
    N     = args.n
    label = N_LABEL_DEFAULT[N]

    with mpl.rc_context(PAPER_RC):
        fig, ax = plt.subplots(figsize=(7.16, 3.0))

        n_eps  = len(EPS_ORDER)
        n_bars = len(MODES)
        group_w = 0.82
        bar_w   = group_w / n_bars
        x_centres = list(range(n_eps))

        drawn_fmts = set()
        panel_max = 0
        for ci, eps in zip(x_centres, EPS_ORDER):
            for bi, (lbl, mode_log) in enumerate(MODES):
                x = ci - group_w/2 + (bi + 0.5) * bar_w
                counts = pick(data, mode_log, N, eps)
                if counts is None:
                    continue
                bottom = 0
                for fmt, color in STACK:
                    if fmt == "Scale":
                        continue
                    v = counts.get(fmt, 0)
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
        ax.tick_params(axis="x", pad=1.5)
        ax.set_xlim(-0.5, n_eps - 0.5)
        ax.set_ylim(0, panel_max * 1.05 if panel_max else 1)
        ax.set_ylabel("Tile count")
        ax.set_axisbelow(True)
        ax.xaxis.grid(False)
        ax.yaxis.grid(True, alpha=0.25, linewidth=0.4)

        # Force FP16 + MXFP16 adjacent in the legend.
        legend_drawn = set(drawn_fmts)
        if "FP16" in legend_drawn or "MXFP16" in legend_drawn:
            legend_drawn.add("FP16"); legend_drawn.add("MXFP16")
        fmt_handles = [Patch(facecolor=c, edgecolor="black",
                             hatch=STACK_HATCH.get(f, ""),
                             label=STACK_LABEL[f])
                       for f, c in STACK
                       if f in legend_drawn and f != "Scale"]
        bar_handles = [Patch(facecolor="white", edgecolor="black",
                             label=f"{i+1}. {lbl}")
                       for i, (lbl, _) in enumerate(MODES)]

        fig.subplots_adjust(left=0.085, right=0.985, top=0.96, bottom=0.30)
        leg1 = fig.legend(handles=fmt_handles,
                          loc="upper center", bbox_to_anchor=(0.27, 0.15),
                          ncol=4, frameon=True, framealpha=0.95,
                          title="Tile format (bottom $\\rightarrow$ top of stack)",
                          handletextpad=0.5, columnspacing=1.2,
                          labelspacing=0.25, borderpad=0.3)
        fig.add_artist(leg1)
        fig.legend(handles=bar_handles,
                   loc="upper center", bbox_to_anchor=(0.78, 0.15),
                   ncol=2, frameon=True, framealpha=0.95,
                   title="Bars (left $\\rightarrow$ right within each $\\varepsilon$)",
                   handletextpad=0.5, columnspacing=1.2,
                   labelspacing=0.25, borderpad=0.3)

        out = Path(args.out_dir) / f"paper_tile_allocation_{label}"
        out_path = out.with_suffix(".pdf")
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        print(out_path)
        plt.close(fig)


if __name__ == "__main__":
    main()
