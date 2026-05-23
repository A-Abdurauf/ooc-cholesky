#!/usr/bin/env python3
"""Paper-ready tile-allocation figure for N=32k alone.

One figure, 4 subplots (one per epsilon).  Inside each subplot, 4 stacked
bars showing per-format tile counts for:

    1. Baseline
    2. Ladder IEEE
    3. Ladder MX  (no MXFP16, staircase)
    4. Ladder MX+MXFP16

NB-pinned to 2048 at 32k.  Storage convention: lower triangle including the
diagonal (M(M+1)/2 tiles per bar).

Colors / hatches / fonts / panel labels match
paper_memory_nb4096_sweep_by_eps_2x2.pdf for a consistent look.
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


N  = 32768
NB = 2048
M  = N // NB  # 16
TILES_TOTAL = M * (M + 1) // 2  # 136

MAIN_CSV  = "/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv"
RERUN_CSV = "/home/abduraa/MX_project/ooc-cholesky/ladder_rerun_32k/results.csv"
STAIR_CSV = "/home/abduraa/MX_project/ooc-cholesky/mx_staircase_vec1d32_gt20k/results.csv"

# (display label, source-csv key, sweep name).  Bars 3+4 use rerun/stair CSVs
# whose tile_breakdown column is already lower-tri.  Bars 1+2 use main CSV
# whose tile_counts_full is full-square -> halve off-diagonals.
BARS = [
    ("Baseline",         "main",  "requant_baseline_fp8_subnormal_gt20k"),
    ("Ladder IEEE",      "main",  "requant_ladder_ieee_gt20k"),
    ("Ladder MX",        "stair", "mx_staircase_vec1d32_gu"),
    ("Ladder MX+MXFP16", "rerun", "ladder_full_gu"),
]

# Tile-format aliases used in the tile_breakdown / tile_counts_full strings.
ALIASES = {
    "fp16":     "FP16",
    "mx_fp16":  "MXFP16",
    "mx_e4m3":  "MXFP8",
    "mx_e5m2":  "MXFP8",
    "e2m1":     "MXFP4",
    "fp8_e4m3": "FP8_plain",
    "fp8_e5m2": "FP8_plain",
    "fp64":     "FP64",
    "fp32":     "FP32",
    "mx_fp32":  "FP32",
}


def load_csv(path):
    return list(csv.DictReader(open(path))) if Path(path).exists() else []


def parse_breakdown(s):
    """tile_breakdown column format: 'fp32=50;mx_e4m3=33;e2m1=27;fp64=16;...'"""
    out = {}
    for part in (s or "").strip().strip('"').split(";"):
        if "=" not in part: continue
        k, v = part.split("=", 1)
        canon = ALIASES.get(k.strip().lower(), k.strip())
        try:
            out[canon] = out.get(canon, 0) + int(v.strip())
        except ValueError:
            pass
    return out


def parse_counts_full(s):
    """tile_counts_full column format: 'fp64=16,fp32=100,fp16=20,fp8_e4m3=120'"""
    out = {}
    for part in (s or "").split(","):
        if "=" not in part: continue
        k, v = part.split("=", 1)
        canon = ALIASES.get(k.strip().lower(), k.strip())
        try:
            out[canon] = out.get(canon, 0) + int(v.strip())
        except ValueError:
            pass
    return out


def full_to_lower_tri(full_counts):
    """Convert full-square counts to lower-tri counts.  M diagonals are FP64."""
    out = {}
    for fmt, v in full_counts.items():
        if fmt == "FP64":
            off = max(v - M, 0)
            out[fmt] = M + off // 2
        else:
            out[fmt] = v // 2
    return out


def find(rows, sweep, eps):
    for r in rows:
        try:
            n_ok = int(r["n"]) == N and int(r["nb"]) == NB
        except (TypeError, ValueError, KeyError):
            continue
        if (r.get("sweep") == sweep and n_ok
                and r.get("source_epsilon") == eps):
            return r
    return None


def counts_for_bar(rows_main, rows_rerun, rows_stair, source, sweep, eps):
    if source == "main":
        r = find(rows_main, sweep, eps)
        return full_to_lower_tri(parse_counts_full(r.get("tile_counts_full", ""))) if r else {}
    if source == "rerun":
        r = find(rows_rerun, sweep, eps)
        return parse_breakdown(r.get("tile_breakdown", "")) if r else {}
    if source == "stair":
        r = find(rows_stair, sweep, eps)
        return parse_breakdown(r.get("tile_breakdown", "")) if r else {}
    return {}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures/paper_tile_allocation_32k")
    args = ap.parse_args()

    rows_main  = load_csv(MAIN_CSV)
    rows_rerun = load_csv(RERUN_CSV)
    rows_stair = load_csv(STAIR_CSV)

    with mpl.rc_context(PAPER_RC):
        fig, ax = plt.subplots(figsize=(7.16, 3.0))

        n_eps  = len(EPS_ORDER)
        n_bars = len(BARS)
        group_w = 0.82
        bar_w   = group_w / n_bars
        x_centres = list(range(n_eps))

        drawn_fmts = set()
        panel_max = 0

        for ci, eps in zip(x_centres, EPS_ORDER):
            for bi, (bar_lbl, src, sweep) in enumerate(BARS):
                x = ci - group_w/2 + (bi + 0.5) * bar_w
                counts = counts_for_bar(rows_main, rows_rerun, rows_stair,
                                        src, sweep, eps)
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
                if bottom > 0:
                    panel_max = max(panel_max, bottom)

        # Light vertical separators between ε groups.
        for i in range(1, n_eps):
            ax.axvline(i - 0.5, color="#bbbbbb",
                       linewidth=0.5, alpha=0.6, zorder=0)

        # ε labels at cluster centres -- per-bar labels are conveyed by the
        # numbered indexed legend below the figure, not on the x-axis.
        ax.set_xticks(x_centres)
        ax.set_xticklabels([rf"$\varepsilon = {EPS_DISPLAY[e][1:-1]}$"
                            for e in EPS_ORDER])

        ax.set_xlim(-0.5, n_eps - 0.5)
        ax.set_ylim(0, panel_max * 1.05 if panel_max else 1)
        ax.set_ylabel("Tile count")
        ax.set_axisbelow(True)
        ax.xaxis.grid(False)
        ax.yaxis.grid(True, alpha=0.25, linewidth=0.4)

        # (N / NB / tile-total info goes in the figure caption, not in-axes.)

        # Single legend below (FP16 + MXFP16 forced adjacent).
        legend_drawn = set(drawn_fmts)
        if "FP16" in legend_drawn or "MXFP16" in legend_drawn:
            legend_drawn.add("FP16"); legend_drawn.add("MXFP16")
        fmt_handles = [Patch(facecolor=c, edgecolor="black",
                             hatch=STACK_HATCH.get(f, ""),
                             label=STACK_LABEL[f])
                       for f, c in STACK
                       if f in legend_drawn and f != "Scale"]

        # Indexed bar-order legend, paired with the format-stack legend.
        bar_handles = [Patch(facecolor="white", edgecolor="black",
                             label=f"{i+1}. {lbl}")
                       for i, (lbl, _, _) in enumerate(BARS)]

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

        out = Path(args.out).with_suffix(".pdf")
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(out)
        plt.close(fig)


if __name__ == "__main__":
    main()
