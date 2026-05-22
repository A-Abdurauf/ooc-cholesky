#!/usr/bin/env python3
"""Paper-grade error figures for ladder mode at a chosen MX granularity.

Emits a 3-bar / 4-eps cluster layout per N:

    1. Plain FP8 baseline   (subnormals allowed; no MX scaling)
    2. IEEE ladder          (FP8 E4M3 -> FP16 -> FP32 -> FP64, no MX)
    3. Full ladder, <gran>  (MXFP4 -> MXFP8 E4M3 -> MXFP16 -> FP32 -> FP64)

Outputs (with --layout=both, the default):
  - one standalone figure per N (with title)
  - one combined figure with one subplot per N, no title; N/NB shown in-axes.

Bar 3 (the MX bar) is read from the main CSV.  For granularity=vec1d32 at
N=32768, we instead use the freshly-re-run ladder_rerun_32k/results.csv
(`ladder_full_gu` = denorms preserved, matching the IEEE-ladder bar).  Other
granularities have no rerun and fall back to the main CSV everywhere.
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
NB_BY_N   = {32768: 2048, 40960: 4096, 65536: 4096}

BAR_FP8_BASELINE = (
    "Plain FP8 baseline (no MX; subnormal-allowed)", "..", "#4e79a7",
)
BAR_IEEE_LADDER  = (
    "IEEE ladder (FP8 E4M3 -> FP16 -> FP32 -> FP64; no MX)", "//", "#a0a0a0",
)

# Per-granularity: (mx-bar label, mx-bar hatch, mx-bar color, main-CSV sweep,
#                   filename-stem suffix, title-fragment)
GRAN = {
    "vec1d32": (
        "Full ladder vec1D-32 (MXFP4 -> MXFP8 -> MXFP16 -> FP32 -> FP64)",
        "xx", "#b07aa1",
        "requant_ladder_scaled_vec1d32_gt20k",
        "vec1d32",
        "vec1D-32",
    ),
    "tile": (
        "Full ladder tile (MXFP4 -> MXFP8 -> MXFP16 -> FP32 -> FP64)",
        "++", "#e07b39",
        "requant_ladder_scaled_tile_gt20k",
        "tile",
        "tile",
    ),
    "block128": (
        "Full ladder block-128 (MXFP4 -> MXFP8 -> MXFP16 -> FP32 -> FP64)",
        "OO", "#59a14f",
        "requant_ladder_scaled_block128_gt20k",
        "block128",
        "block-128",
    ),
}


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


def draw_subplot(ax, n, gran_keys, main_rows, rerun_rows,
                 use_rerun, show_ylabel, value_fontsize=7.2):
    if isinstance(gran_keys, str):
        gran_keys = [gran_keys]
    bars   = [BAR_FP8_BASELINE, BAR_IEEE_LADDER]
    sweeps = ["requant_baseline_fp8_subnormal_gt20k",
              "requant_ladder_ieee_gt20k"]
    mx_sweeps_by_index = {}
    for gk in gran_keys:
        mx_label, mx_hatch, mx_color, mx_sweep, _, _ = GRAN[gk]
        mx_sweeps_by_index[len(bars)] = (gk, mx_sweep)
        bars.append((mx_label, mx_hatch, mx_color))
        sweeps.append(mx_sweep)

    n_eps  = len(EPS_ORDER)
    n_bars = len(bars)
    group_w = 0.78
    bar_w   = group_w / n_bars
    x_centres = list(range(n_eps))

    y_min, y_max = float("inf"), 0.0

    for bi, ((label, hatch, color), sweep) in enumerate(zip(bars, sweeps)):
        xs, ys = [], []
        for ci, eps in zip(x_centres, EPS_ORDER):
            row = main_rows.get((sweep, n, eps))
            # Prefer rerun (GU) for the vec1d32 bar at the rerun N.
            gk_for_bar = mx_sweeps_by_index.get(bi, (None, None))[0]
            if use_rerun and gk_for_bar == "vec1d32":
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
                        ha="center", va="bottom",
                        fontsize=value_fontsize, rotation=90, color="#222")

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

    return bars


def legend_handles(bars):
    return [Patch(facecolor=color, edgecolor="black", hatch=hatch, label=label)
            for label, hatch, color in bars]


def plot_standalone(n, gran_keys, main_rows, rerun_rows, out_path, use_rerun):
    if isinstance(gran_keys, str):
        gran_keys = [gran_keys]
    width = 9.2 + 0.8 * max(0, len(gran_keys) - 1)
    fig, ax = plt.subplots(figsize=(width, 5.2))
    bars = draw_subplot(ax, n, gran_keys, main_rows, rerun_rows,
                        use_rerun=use_rerun, show_ylabel=True, value_fontsize=8)

    title_grans = " + ".join(GRAN[gk][5] for gk in gran_keys)
    nb = NB_BY_N.get(n, "?")
    tail = "  (32k vec1d32 uses ladder_rerun GU)" if (
        use_rerun and "vec1d32" in gran_keys) else ""
    ax.set_title(f"Ladder mode at {title_grans} granularity, N = {n}, NB = {nb}{tail}")

    legend_ncol = 2 if len(bars) >= 4 else 1
    fig.legend(handles=legend_handles(bars),
               loc="lower center", bbox_to_anchor=(0.5, -0.16),
               ncol=legend_ncol, frameon=True, framealpha=0.95,
               title="Bar order within each cluster (left -> right)   .   grey dashed = source $\\varepsilon$",
               handletextpad=0.6, columnspacing=2.0,
               labelspacing=0.35, borderpad=0.4)
    fig.tight_layout(rect=[0, 0.0, 1, 0.99])

    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    p = out.with_suffix(".pdf")
    fig.savefig(p, dpi=300, bbox_inches="tight")
    print(p)
    plt.close(fig)


def plot_combined(ns, gran_keys, main_rows, rerun_rows, out_path, rerun_for_n):
    if isinstance(gran_keys, str):
        gran_keys = [gran_keys]
    n_panels = len(ns)
    fig, axes = plt.subplots(1, n_panels,
                             figsize=(4.7 * n_panels + 0.8, 5.0),
                             sharey=True)
    if n_panels == 1:
        axes = [axes]

    bars = None
    for ax, n in zip(axes, ns):
        bars = draw_subplot(ax, n, gran_keys, main_rows, rerun_rows,
                            use_rerun=(n == rerun_for_n),
                            show_ylabel=(ax is axes[0]))
        nb = NB_BY_N.get(n, "?")
        ax.text(0.985, 0.965, f"N = {_fmt_n(n)}\nNB = {nb}",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=10.5,
                bbox=dict(boxstyle="round,pad=0.32",
                          facecolor="white", edgecolor="#888", linewidth=0.8,
                          alpha=0.92))

    # Two-row legend when there are 4+ bars (better structure than one long row).
    legend_ncol = 2 if len(bars) >= 4 else len(bars)
    fig.legend(handles=legend_handles(bars),
               loc="lower center", bbox_to_anchor=(0.5, -0.08),
               ncol=legend_ncol, frameon=True, framealpha=0.95,
               title="Bars within each cluster (left -> right)   .   grey dashed = source $\\varepsilon$",
               handletextpad=0.6, columnspacing=2.4, labelspacing=0.4,
               borderpad=0.4)
    fig.tight_layout(rect=[0, 0.06, 1, 1.0])

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
    ap.add_argument("--granularities", nargs="+", choices=list(GRAN),
                    required=True,
                    help="One or more MX granularities to add as bars "
                         "(in addition to FP8 baseline + IEEE ladder).")
    ap.add_argument("--ns", nargs="+", type=int, default=[32768, 40960, 65536])
    ap.add_argument("--rerun-for-n", type=int, default=32768)
    ap.add_argument("--out-dir",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures")
    ap.add_argument("--layout", choices=["separate", "combined", "both"],
                    default="both")
    args = ap.parse_args()

    main_rows  = load_main(args.csv)
    rerun_rows = load_rerun(args.rerun_csv)
    stem_suffix = "_".join(GRAN[gk][4] for gk in args.granularities)

    if args.layout in ("separate", "both"):
        for n in args.ns:
            out = Path(args.out_dir) / f"paper_error_ladder_{stem_suffix}_N{n}"
            plot_standalone(n, args.granularities, main_rows, rerun_rows, out,
                            use_rerun=(n == args.rerun_for_n))

    if args.layout in ("combined", "both"):
        out = Path(args.out_dir) / f"paper_error_ladder_{stem_suffix}"
        plot_combined(args.ns, args.granularities, main_rows, rerun_rows, out,
                      rerun_for_n=args.rerun_for_n)


if __name__ == "__main__":
    main()
