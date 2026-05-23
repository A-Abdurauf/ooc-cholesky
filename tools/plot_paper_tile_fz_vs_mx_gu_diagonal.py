#!/usr/bin/env python3
"""Diagonal comparison figure: Full MX ladder under tile+FZ vs vec1d32+GU.

Two rows, 4 epsilon subplots per row, 3 N groups per subplot, 2 bars per N
group (tile+FZ left, vec1d+GU right):

    Row 1 (top)    : relative factorisation error, log y
    Row 2 (bottom) : per-format tile-allocation stacked bars (lower triangle)

Reads everything from the merged figures/all_error_runs.csv.
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

MERGED_CSV = "/home/abduraa/MX_project/ooc-cholesky/figures/all_error_runs.csv"
EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]
N_ORDER   = [32768, 40960, 65536]
NB_BY_N   = {32768: 2048, 40960: 4096, 65536: 4096}
N_LABEL   = {32768: "32k", 40960: "40k", 65536: "64k"}

ALIASES = {"fp16":"FP16","mx_fp16":"MXFP16","mx_e4m3":"MXFP8","mx_e5m2":"MXFP8",
           "e2m1":"MXFP4","fp64":"FP64","fp32":"FP32",
           "fp8_e4m3":"FP8_plain","fp8_e5m2":"FP8_plain","mx_fp32":"FP32"}

# Two bars per cluster.
BARS = [
    ("tile + FZ",   "tile",  "fz", "#4e79a7", ""),
    ("MX + GU",     "vec1d", "gu", "#b07aa1", "//"),
]


def load():
    return list(csv.DictReader(open(MERGED_CSV)))


def pick(rows, mode, gran, uflow, n, nb, eps):
    for r in rows:
        if (r["mode"] == mode and r["granularity"] == gran
                and r["underflow"] == uflow
                and int(r["N"]) == n and int(r["nb"]) == nb
                and r["eps"] == eps):
            return r
    return None


def parse_bd(s, kind):
    out = {}
    sep = ";" if kind == "lower" else ","
    for part in (s or "").strip().strip('"').split(sep):
        if "=" not in part: continue
        k, v = part.split("=", 1)
        k = ALIASES.get(k.strip().lower(), k.strip())
        try: out[k] = out.get(k, 0) + int(v.strip())
        except: pass
    return out


def lower_tri(counts, M):
    out = {}
    for fmt, v in counts.items():
        if fmt == "FP64":
            off = max(v - M, 0); out[fmt] = M + off // 2
        else:
            out[fmt] = v // 2
    return out


def counts_for(r, n, nb):
    if r is None: return {}
    raw = parse_bd(r["tile_breakdown"], r["tile_breakdown_kind"])
    if r["tile_breakdown_kind"] == "full":
        return lower_tri(raw, n // nb)
    return raw


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures/paper_tile_fz_vs_mx_gu_diagonal")
    args = ap.parse_args()

    rows = load()

    with mpl.rc_context(PAPER_RC):
        fig, axes = plt.subplots(2, 4, figsize=(8.0, 4.6),
                                 sharex="col")

        n_groups = len(N_ORDER)
        n_bars   = len(BARS)
        group_w  = 0.78
        bar_w    = group_w / n_bars
        x_centres = list(range(n_groups))

        # Per-row min/max so we can share y-limits horizontally.
        err_min, err_max = float("inf"), 0.0
        tile_max = 0
        drawn_fmts = set()

        # Pass 1: compute error/tile values and store them
        cell_data = {}
        for ci, eps in enumerate(EPS_ORDER):
            for ni, N in enumerate(N_ORDER):
                nb = NB_BY_N[N]
                vals = []
                for (lbl, gran, uflow, color, hatch) in BARS:
                    r = pick(rows, "ladder_mx_full", gran, uflow, N, nb, eps)
                    if r is None:
                        vals.append((None, {}, color, hatch))
                        continue
                    e = float(r["rel_factor_error"])
                    c = counts_for(r, N, nb)
                    if e > 0:
                        err_min = min(err_min, e); err_max = max(err_max, e)
                    total = sum(c.values())
                    tile_max = max(tile_max, total)
                    vals.append((e, c, color, hatch))
                cell_data[(ci, ni)] = vals

        # Pass 2: render top (error) and bottom (tile stack)
        for ci, eps in enumerate(EPS_ORDER):
            ax_e = axes[0, ci]
            ax_t = axes[1, ci]
            for ni, N in enumerate(N_ORDER):
                for bi, (e_val, counts, color, hatch) in enumerate(cell_data[(ci, ni)]):
                    x = ni - group_w/2 + (bi + 0.5) * bar_w
                    # Top: single error bar
                    if e_val is not None and e_val > 0:
                        ax_e.bar([x], [e_val], bar_w * 0.92,
                                 color=color, edgecolor="black",
                                 linewidth=0.25, hatch=hatch)
                    # Bottom: stacked by format
                    bottom = 0
                    for fmt, fcolor in STACK:
                        if fmt == "Scale": continue
                        v = counts.get(fmt, 0)
                        if v <= 0: continue
                        drawn_fmts.add(fmt)
                        ax_t.bar([x], [v], bar_w * 0.92, bottom=[bottom],
                                 color=fcolor, edgecolor="black",
                                 linewidth=0.2, hatch=STACK_HATCH.get(fmt, ""))
                        bottom += v

            # Reference eps line + separators on top row
            ev = float(eps)
            ax_e.axhline(ev, color="gray", linestyle="--",
                         linewidth=0.7, alpha=0.6, zorder=0)
            for ax in (ax_e, ax_t):
                for i in range(1, n_groups):
                    ax.axvline(i - 0.5, color="#bbbbbb",
                               linewidth=0.5, alpha=0.5, zorder=0)
                ax.set_xticks(x_centres)
                ax.set_xticklabels([N_LABEL[N] for N in N_ORDER])
                ax.tick_params(axis="x", pad=1.5)
                ax.set_xlim(-0.6, n_groups - 0.4)
                ax.set_axisbelow(True)
                ax.xaxis.grid(False)

            ax_e.set_yscale("log")
            ax_e.yaxis.grid(True, alpha=0.25, linewidth=0.4, which="both")
            ax_t.yaxis.grid(True, alpha=0.25, linewidth=0.4)
            ax_e.text(0.03, 0.965,
                      f"$\\varepsilon = {EPS_DISPLAY[eps][1:-1]}$",
                      transform=ax_e.transAxes, ha="left", va="top",
                      fontsize=8.5, fontweight="bold")

        # Shared y-limits.
        if err_min < float("inf"):
            for ax in axes[0, :]:
                ax.set_ylim(err_min/4, err_max*20)
        if tile_max > 0:
            for ax in axes[1, :]:
                ax.set_ylim(0, tile_max*1.05)
        # Hide redundant tick labels on non-leftmost columns.
        for ax in axes[:, 1:].flatten():
            ax.tick_params(labelleft=False)

        axes[0, 0].set_ylabel("Rel. error  (log)")
        axes[1, 0].set_ylabel("Tile count")

        # Legends below.
        legend_drawn = set(drawn_fmts)
        if "FP16" in legend_drawn or "MXFP16" in legend_drawn:
            legend_drawn.add("FP16"); legend_drawn.add("MXFP16")
        fmt_handles = [Patch(facecolor=c, edgecolor="black",
                             hatch=STACK_HATCH.get(f, ""),
                             label=STACK_LABEL[f])
                       for f, c in STACK if f in legend_drawn and f != "Scale"]
        bar_handles = [Patch(facecolor=color, edgecolor="black",
                             hatch=hatch, label=lbl)
                       for lbl, _, _, color, hatch in BARS]
        eps_handle = plt.Line2D([0], [0], color="gray", linestyle="--",
                                linewidth=0.9, label="source $\\varepsilon$")

        fig.subplots_adjust(left=0.082, right=0.99, top=0.97,
                            bottom=0.24, wspace=0.08, hspace=0.13)
        leg1 = fig.legend(handles=bar_handles + [eps_handle],
                          loc="upper center", bbox_to_anchor=(0.25, 0.14),
                          ncol=3, frameon=True, framealpha=0.95,
                          title="Bars (left $\\rightarrow$ right within each $N$)",
                          handletextpad=0.5, columnspacing=1.4,
                          labelspacing=0.25, borderpad=0.3)
        fig.add_artist(leg1)
        fig.legend(handles=fmt_handles,
                   loc="upper center", bbox_to_anchor=(0.72, 0.14),
                   ncol=4, frameon=True, framealpha=0.95,
                   title="Tile format (bottom row stack)",
                   handletextpad=0.5, columnspacing=1.2,
                   labelspacing=0.25, borderpad=0.3)

        out = Path(args.out).with_suffix(".pdf")
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(out)
        plt.close(fig)


if __name__ == "__main__":
    main()
