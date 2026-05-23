#!/usr/bin/env python3
"""Paper-ready performance sweep figure.

1x4 panels, one per ladder mode (Baseline / Ladder IEEE / Ladder MX /
Ladder MX+MXFP16).  Within each panel:
  - x-axis: matrix size N (20k -> 120k, log-friendly labels)
  - y-axis: factorisation perf in TFLOPS
  - 4 colored curves: one per epsilon (1e-5, 1e-6, 1e-7, 1e-8)

Reads from recreate-ooc-chol/build/combined_perf_data.csv (4 modes x 7 N x
4 eps, OK rows only, 28 rows per mode).

Same paper-ready styling as the memory + error figures: 9pt body / 7.5pt
ticks / 8.5pt bold panel labels, in-axes panel titles in the top-left,
decimal-k tick labels (100k for N=98304), shared y-axis.
"""
import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

import sys
sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_paper_memory_nb4096_sweep import (  # noqa: E402
    PAPER_RC, EPS_DISPLAY, N_LABEL, _xlabel,
)


CSV = "/home/abduraa/MX_project/recreate-ooc-chol/build/combined_perf_data.csv"

# Match the canonical mode names already in the CSV.
MODES = [
    ("Baseline",         "Baseline"),
    ("Ladder IEEE",      "Ladder IEEE"),
    ("Ladder MX",        "Ladder MX"),
    ("Ladder MX+MXFP16", "Ladder MX+MXFP16"),
]

# Tick positions are round multiples of 10k; data points sit at their TRUE N
# values (e.g. N=65536 plots between the "60k" tick and the "80k" tick).
XTICK_POS    = [20000, 40000, 60000, 80000, 100000, 120000]
XTICK_LABELS = ["20k", "40k", "60k", "80k", "100k", "120k"]
EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]
# Wong / Okabe-Ito (color-blind safe).
EPS_COLOR = {
    "1e-5": "#0072B2",
    "1e-6": "#E69F00",
    "1e-7": "#009E73",
    "1e-8": "#D55E00",
}
EPS_MARKER = {"1e-5": "o", "1e-6": "s", "1e-7": "^", "1e-8": "D"}


def _norm_eps(s):
    """combined_perf_data.csv writes '1e-08' etc; canonicalise."""
    s = (s or "").strip()
    return s.replace("1e-05", "1e-5").replace("1e-06", "1e-6") \
            .replace("1e-07", "1e-7").replace("1e-08", "1e-8")


def load():
    rows = []
    with open(CSV) as f:
        for r in csv.DictReader(f):
            if r.get("status") and r["status"] != "OK":
                continue
            try:
                r["N"]    = int(r["N"])
                r["perf"] = float(r["perf_TFLOPS"])
                r["eps"]  = _norm_eps(r["eps"])
            except (TypeError, ValueError):
                continue
            rows.append(r)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default=CSV)
    ap.add_argument("--out",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures/paper_perf_sweep_by_mode")
    args = ap.parse_args()

    rows = load()

    with mpl.rc_context(PAPER_RC):
        fig, axes = plt.subplots(1, 4, figsize=(8.0, 2.7), sharey=True)

        ns_seen = set()
        perf_max = 0.0
        for ai, (ax, (mode_key, mode_label)) in enumerate(zip(axes, MODES)):
            for eps in EPS_ORDER:
                pts = [r for r in rows
                       if r.get("mode") == mode_key and r["eps"] == eps]
                pts.sort(key=lambda r: r["N"])
                if not pts:
                    continue
                xs = [r["N"] for r in pts]
                ys = [r["perf"] for r in pts]
                ax.plot(xs, ys,
                        marker=EPS_MARKER[eps], markersize=4.5,
                        linewidth=1.4, color=EPS_COLOR[eps],
                        markerfacecolor=EPS_COLOR[eps],
                        markeredgecolor="black", markeredgewidth=0.4,
                        label=rf"$\varepsilon = {EPS_DISPLAY[eps][1:-1]}$")
                ns_seen.update(xs)
                perf_max = max(perf_max, max(ys))

            # In-axes panel label
            ax.text(0.025, 0.965,
                    f"({chr(ord('a') + ai)})  {mode_label}",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=8.5, fontweight="bold")
            ax.set_axisbelow(True)
            ax.grid(True, alpha=0.22, linewidth=0.4)

        # Shared y-limit + x-axis ticks.  Bumped headroom so the (d) panel
        # label doesn't overlap the topmost line at N=120k.
        for ax in axes:
            ax.set_ylim(0, perf_max * 1.18)
            ax.set_xlabel(r"Matrix size  $N$")
            ax.tick_params(axis="x", pad=1.5)
        axes[0].set_ylabel("Perf (TFLOPS)")

        for ax in axes:
            ax.set_xticks(XTICK_POS)
            ax.set_xticklabels(XTICK_LABELS,
                               rotation=45, ha="right", fontsize=7)
            # Reasonable x-limits with a little padding around the data range.
            ax.set_xlim(18000, 126000)

        eps_handles = [Line2D([0], [0], marker=EPS_MARKER[e], linestyle="-",
                              color=EPS_COLOR[e], markerfacecolor=EPS_COLOR[e],
                              markeredgecolor="black", markeredgewidth=0.4,
                              markersize=5.5,
                              label=rf"$\varepsilon = {EPS_DISPLAY[e][1:-1]}$")
                       for e in EPS_ORDER]

        # Reserve enough at the bottom for the x-label band, then place the
        # legend BELOW it so they don't overlap.
        fig.subplots_adjust(left=0.075, right=0.985, top=0.96,
                            bottom=0.32, wspace=0.08)
        fig.legend(handles=eps_handles,
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
