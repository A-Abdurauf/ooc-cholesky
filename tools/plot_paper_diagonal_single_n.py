#!/usr/bin/env python3
"""Single-N variant of paper_diagonal_error_memory.pdf.

Two panels, one matrix size: top = relative error (log), bottom = memory (GB)
stacked by tile format.  X-axis: 4 epsilons; 2 bars per epsilon
(tile+FZ vs MX+GU).  Default N = 65536 (64k), overridable via --n.
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
    EPS_DISPLAY, STACK, STACK_LABEL, STACK_HATCH, PAPER_RC,
)
from plot_paper_diagonal_error_memory import (  # noqa: E402
    MERGED_CSV, BARS, FMT_BYTES, FMT_ALIAS_IN,
    load, pick, parse_lower, memory_breakdown_gb,
)

EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]
NB_BY_N   = {32768: 2048, 40960: 4096, 65536: 4096}
N_LABEL   = {32768: "32k", 40960: "40k", 65536: "64k"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="Ladder MX+MXFP16")
    ap.add_argument("--n", type=int, default=65536, choices=[32768, 40960, 65536])
    ap.add_argument("--out",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures/paper_diagonal_single_n")
    args = ap.parse_args()

    rows = load()
    N  = args.n
    nb = NB_BY_N[N]

    with mpl.rc_context(PAPER_RC):
        fig, (ax_e, ax_m) = plt.subplots(2, 1, figsize=(4.4, 3.6),
                                         sharex=True,
                                         gridspec_kw={"height_ratios": [1.0, 1.0]})

        n_eps = len(EPS_ORDER); n_bars = len(BARS)
        group_w = 0.78
        bar_w   = group_w / n_bars
        x_centres = list(range(n_eps))

        drawn_fmts = set()
        err_min, err_max = float("inf"), 0.0
        mem_max = 0.0
        for ci, eps in enumerate(EPS_ORDER):
            for bi, (lbl, gran, uflow, color, hatch) in enumerate(BARS):
                x = ci - group_w/2 + (bi + 0.5) * bar_w
                r = pick(rows, args.mode, gran, uflow, N, nb, eps)
                if r is None:
                    continue
                try:    err = float(r["rel_factor_error"])
                except: err = None
                counts = parse_lower(r["tile_breakdown"])
                gb     = memory_breakdown_gb(counts, N, nb)

                if err and err > 0:
                    ax_e.bar([x], [err], bar_w * 0.92,
                             color=color, edgecolor="black",
                             linewidth=0.25, hatch=hatch)
                    err_min = min(err_min, err); err_max = max(err_max, err)

                bottom = 0.0
                for fmt, fcolor in STACK:
                    v = gb.get(fmt, 0.0)
                    if v <= 0: continue
                    drawn_fmts.add(fmt)
                    ax_m.bar([x], [v], bar_w * 0.92, bottom=[bottom],
                             color=fcolor, edgecolor="black",
                             linewidth=0.2, hatch=STACK_HATCH.get(fmt, ""))
                    bottom += v
                if bottom > mem_max: mem_max = bottom

            # ε reference line on the top panel
            ax_e.hlines(float(eps), ci - group_w/2 - 0.06, ci + group_w/2 + 0.06,
                        colors="gray", linestyles="--",
                        linewidth=0.7, alpha=0.65, zorder=0)

        for ax in (ax_e, ax_m):
            ax.set_xticks(x_centres)
            ax.set_xticklabels([rf"$\varepsilon = {EPS_DISPLAY[e][1:-1]}$"
                                for e in EPS_ORDER])
            ax.set_xlim(-0.5, n_eps - 0.5)
            ax.tick_params(axis="x", pad=1.5)
            ax.set_axisbelow(True)
            ax.xaxis.grid(False)
        ax_e.set_yscale("log")
        ax_e.yaxis.grid(True, alpha=0.25, linewidth=0.4, which="both")
        ax_m.yaxis.grid(True, alpha=0.25, linewidth=0.4)
        if err_min < float("inf"): ax_e.set_ylim(err_min/4, err_max*20)
        if mem_max > 0:            ax_m.set_ylim(0, mem_max * 1.07)
        ax_e.set_ylabel("Rel. error  (log)")
        ax_m.set_ylabel("Memory (GB)")

        legend_drawn = set(drawn_fmts)
        if "FP16" in legend_drawn or "MXFP16" in legend_drawn:
            legend_drawn.add("FP16"); legend_drawn.add("MXFP16")
        fmt_handles = [Patch(facecolor=c, edgecolor="black",
                             hatch=STACK_HATCH.get(f, ""),
                             label=STACK_LABEL[f])
                       for f, c in STACK if f in legend_drawn]
        bar_handles = [Patch(facecolor=color, edgecolor="black",
                             hatch=hatch, label=lbl)
                       for lbl, _, _, color, hatch in BARS]
        eps_handle = plt.Line2D([0], [0], color="gray", linestyle="--",
                                linewidth=0.9, label="source $\\varepsilon$")

        fig.subplots_adjust(left=0.155, right=0.985, top=0.97,
                            bottom=0.28, hspace=0.10)
        leg1 = fig.legend(handles=bar_handles + [eps_handle],
                          loc="upper center", bbox_to_anchor=(0.32, 0.18),
                          ncol=1, frameon=True, framealpha=0.95,
                          title=f"Bars within each $\\varepsilon$ ($N = {N_LABEL[N]}$)",
                          handletextpad=0.5, columnspacing=1.2,
                          labelspacing=0.25, borderpad=0.3)
        fig.add_artist(leg1)
        fig.legend(handles=fmt_handles,
                   loc="upper center", bbox_to_anchor=(0.78, 0.18),
                   ncol=2, frameon=True, framealpha=0.95,
                   title="Tile format (bottom row stack)",
                   handletextpad=0.5, columnspacing=1.2,
                   labelspacing=0.25, borderpad=0.3)

        out_path = Path(f"{args.out}_{N_LABEL[N]}").with_suffix(".pdf")
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=300, bbox_inches="tight")
        print(out_path)
        plt.close(fig)


if __name__ == "__main__":
    main()
