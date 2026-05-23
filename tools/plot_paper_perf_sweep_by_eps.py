#!/usr/bin/env python3
"""Paper-ready performance sweep figure, eps-as-subplot layout.

1x4 panels, one per source epsilon (10^-5..10^-8).  Within each panel:
  - x-axis: matrix size N on a linear axis (true N positions, round-number
    tick labels 20k/40k/60k/80k/100k/120k)
  - y-axis: factorisation perf (TFLOPS)
  - 4 colored lines: one per ladder mode
      Baseline / Ladder IEEE / Ladder MX / Ladder MX+MXFP16

Reads recreate-ooc-chol/build/combined_perf_data.csv.  Shared y-axis.
Paper-ready styling -- shorter than the by-mode variant since the eps
panel label is small.
"""
import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_paper_memory_nb4096_sweep import PAPER_RC, EPS_DISPLAY  # noqa: E402
from plot_paper_perf_sweep_by_mode import (  # noqa: E402
    CSV, MODES, XTICK_POS, XTICK_LABELS, _norm_eps, load,
)


EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]

# Color + marker per mode -- matches the memory figure's bar palette so the
# two figures share a visual language.
MODE_COLOR = {
    "Baseline":         "#4e79a7",
    "Ladder IEEE":      "#7f7f7f",
    "Ladder MX":        "#e07b39",
    "Ladder MX+MXFP16": "#b07aa1",
}
MODE_MARKER = {
    "Baseline":         "o",
    "Ladder IEEE":      "s",
    "Ladder MX":        "^",
    "Ladder MX+MXFP16": "D",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=CSV)
    ap.add_argument("--out",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures/paper_perf_sweep_by_eps")
    args = ap.parse_args()

    rows = load()

    with mpl.rc_context(PAPER_RC):
        fig, axes = plt.subplots(1, 4, figsize=(8.0, 2.7), sharey=True)

        perf_max = 0.0
        for ai, (ax, eps) in enumerate(zip(axes, EPS_ORDER)):
            for mode_key, mode_label in MODES:
                pts = [r for r in rows
                       if r.get("mode") == mode_key and r["eps"] == eps]
                pts.sort(key=lambda r: r["N"])
                if not pts: continue
                xs = [r["N"] for r in pts]
                ys = [r["perf"] for r in pts]
                ax.plot(xs, ys,
                        marker=MODE_MARKER[mode_key], markersize=4.0,
                        linewidth=1.3,
                        color=MODE_COLOR[mode_key],
                        markerfacecolor=MODE_COLOR[mode_key],
                        markeredgecolor="black", markeredgewidth=0.35,
                        label=mode_label)
                perf_max = max(perf_max, max(ys))

            ax.text(0.04, 0.965,
                    f"$\\varepsilon = {EPS_DISPLAY[eps][1:-1]}$",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=9, fontweight="bold")
            ax.set_axisbelow(True)
            ax.grid(True, alpha=0.22, linewidth=0.4)
            ax.set_xticks(XTICK_POS)
            ax.set_xticklabels(XTICK_LABELS, rotation=45, ha="right", fontsize=7)
            ax.set_xlim(18000, 126000)
            ax.tick_params(axis="x", pad=1.5)
            ax.set_xlabel(r"Matrix size  $N$")

        for ax in axes:
            ax.set_ylim(0, perf_max * 1.12)
        axes[0].set_ylabel("Perf (TFLOPS)")

        mode_handles = [Line2D([0], [0], marker=MODE_MARKER[mk], linestyle="-",
                               color=MODE_COLOR[mk],
                               markerfacecolor=MODE_COLOR[mk],
                               markeredgecolor="black", markeredgewidth=0.4,
                               markersize=5.5, label=lbl)
                        for mk, lbl in MODES]

        # Larger top portion for the y-axis; tight gap to the legend.
        fig.subplots_adjust(left=0.075, right=0.985, top=0.96,
                            bottom=0.32, wspace=0.08)
        fig.legend(handles=mode_handles,
                   loc="upper center", bbox_to_anchor=(0.5, 0.10),
                   ncol=4, frameon=True, framealpha=0.95,
                   handletextpad=0.5, columnspacing=1.6,
                   labelspacing=0.25, borderpad=0.3)

        out = Path(args.out).with_suffix(".pdf")
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(out)
        plt.close(fig)


if __name__ == "__main__":
    main()
