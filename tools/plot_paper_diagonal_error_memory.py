#!/usr/bin/env python3
"""Diagonal comparison figure: Full MX ladder under tile+FZ vs MX+GU.

Same 2-row layout as paper_tile_fz_vs_mx_gu_diagonal.pdf, but the bottom
row is per-format MEMORY (GB) instead of tile counts.  Storage convention:
lower triangle + half-diagonal (FP64 diagonals = nb*(nb+1)/2 elements);
MX formats add 1-byte shared scale per 32 elements.

Reads from figures/all_error_runs_ladder_only.csv (already lower-tri).
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

MERGED_CSV = "/home/abduraa/MX_project/ooc-cholesky/figures/all_error_runs_ladder_only.csv"
EPS_ORDER_DEFAULT = ["1e-5", "1e-6", "1e-7", "1e-8"]
EPS_ORDER          = EPS_ORDER_DEFAULT  # overridden by --eps
N_ORDER    = [32768, 40960, 65536]
NB_BY_N    = {32768: 2048, 40960: 4096, 65536: 4096}
N_LABEL    = {32768: "32k", 40960: "40k", 65536: "64k"}

# Two bars per cluster.
BARS = [
    ("tile + FZ", "tile", "fz", "#4e79a7", ""),
    ("MX + GU",   "MX",   "gu", "#b07aa1", "//"),
]

# Bytes per element + scale-bearing flag.  Keys must match STACK from the
# memory script.
FMT_BYTES = {
    "FP64":      (8.0, False),
    "FP32":      (4.0, False),
    "FP16":      (2.0, False),
    "MXFP16":    (2.0, True),
    "FP8_plain": (1.0, False),
    "MXFP8":     (1.0, True),
    "MXFP4":     (0.5, True),
}
# The merged ladder-only CSV uses "FP8" not "FP8_plain"; alias for lookup.
FMT_ALIAS_IN = {"FP8": "FP8_plain"}


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


def parse_lower(s):
    """Parse a lower-tri tile_breakdown like 'FP64=16;FP32=50;...'."""
    out = {}
    for part in (s or "").rstrip(";").split(";"):
        if "=" not in part: continue
        k, v = part.split("=", 1)
        canon = FMT_ALIAS_IN.get(k.strip(), k.strip())
        try: out[canon] = out.get(canon, 0) + int(v.strip())
        except: pass
    return out


def memory_breakdown_gb(counts, n, nb):
    """Return {format: GB} including a 'Scale' bucket for shared scales."""
    M = n // nb
    tile_full = nb * nb
    tile_diag = nb * (nb + 1) // 2
    GB = 1 / 1024 ** 3
    out = {fmt: 0.0 for fmt, _ in STACK}
    for fmt, (bpe, has_scale) in FMT_BYTES.items():
        cnt = counts.get(fmt, 0)
        if cnt <= 0: continue
        if fmt == "FP64":
            d = min(M, cnt); o = max(cnt - d, 0)
            data_bytes = d * tile_diag * bpe + o * tile_full * bpe
        else:
            data_bytes = cnt * tile_full * bpe
        out[fmt] += data_bytes * GB
        if has_scale:
            out["Scale"] += cnt * (tile_full // 32) * GB
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="Ladder MX+MXFP16",
                    help="Mode column value (canonical full MX ladder by default).")
    ap.add_argument("--eps", nargs="+", default=EPS_ORDER_DEFAULT,
                    help="Which epsilons to plot, e.g. --eps 1e-6 1e-7")
    ap.add_argument("--out",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures/paper_diagonal_error_memory")
    args = ap.parse_args()

    global EPS_ORDER
    EPS_ORDER = list(args.eps)
    rows = load()

    with mpl.rc_context(PAPER_RC):
        n_cols = len(EPS_ORDER)
        # Width scales with number of eps panels; keep height same.
        fig_w = 2.1 * n_cols + 0.5
        fig, axes = plt.subplots(2, n_cols, figsize=(fig_w, 4.4),
                                 sharex="col")
        if n_cols == 1:
            axes = axes.reshape(2, 1)

        n_groups = len(N_ORDER)
        n_bars   = len(BARS)
        group_w  = 0.78
        bar_w    = group_w / n_bars
        x_centres = list(range(n_groups))

        err_min, err_max = float("inf"), 0.0
        mem_max = 0.0
        drawn_fmts = set()
        cell_data = {}

        # Pass 1: gather error + per-format GB for each (eps, N, bar).
        for ci, eps in enumerate(EPS_ORDER):
            for ni, N in enumerate(N_ORDER):
                nb = NB_BY_N[N]
                vals = []
                for (lbl, gran, uflow, color, hatch) in BARS:
                    r = pick(rows, args.mode, gran, uflow, N, nb, eps)
                    if r is None:
                        vals.append((None, {}, color, hatch))
                        continue
                    try:    err = float(r["rel_factor_error"])
                    except: err = None
                    counts = parse_lower(r["tile_breakdown"])
                    gb     = memory_breakdown_gb(counts, N, nb)
                    if err and err > 0:
                        err_min = min(err_min, err); err_max = max(err_max, err)
                    total_gb = sum(gb.values())
                    mem_max  = max(mem_max, total_gb)
                    vals.append((err, gb, color, hatch))
                cell_data[(ci, ni)] = vals

        # Pass 2: render top (error) + bottom (memory stacked).
        for ci, eps in enumerate(EPS_ORDER):
            ax_e = axes[0, ci]; ax_m = axes[1, ci]
            for ni, N in enumerate(N_ORDER):
                for bi, (e_val, gb_dict, color, hatch) in enumerate(cell_data[(ci, ni)]):
                    x = ni - group_w/2 + (bi + 0.5) * bar_w
                    # Top: solid error bar
                    if e_val and e_val > 0:
                        ax_e.bar([x], [e_val], bar_w * 0.92,
                                 color=color, edgecolor="black",
                                 linewidth=0.25, hatch=hatch)
                    # Bottom: stacked memory bar
                    bottom = 0.0
                    for fmt, fcolor in STACK:
                        v = gb_dict.get(fmt, 0.0)
                        if v <= 0: continue
                        drawn_fmts.add(fmt)
                        ax_m.bar([x], [v], bar_w * 0.92, bottom=[bottom],
                                 color=fcolor, edgecolor="black",
                                 linewidth=0.2, hatch=STACK_HATCH.get(fmt, ""))
                        bottom += v

            # eps reference line on top row
            ev = float(eps)
            ax_e.axhline(ev, color="gray", linestyle="--",
                         linewidth=0.7, alpha=0.6, zorder=0)
            # Vertical separators between N groups
            for ax in (ax_e, ax_m):
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
            ax_m.yaxis.grid(True, alpha=0.25, linewidth=0.4)
            ax_e.text(0.03, 0.965,
                      f"$\\varepsilon = {EPS_DISPLAY[eps][1:-1]}$",
                      transform=ax_e.transAxes, ha="left", va="top",
                      fontsize=8.5, fontweight="bold")

        if err_min < float("inf"):
            for ax in axes[0, :]: ax.set_ylim(err_min/4, err_max*20)
        if mem_max > 0:
            for ax in axes[1, :]: ax.set_ylim(0, mem_max*1.07)
        for ax in axes[:, 1:].flatten():
            ax.tick_params(labelleft=False)

        axes[0, 0].set_ylabel("Rel. error  (log)")
        axes[1, 0].set_ylabel("Memory (GB)")

        # Legends.
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
                   title="Tile format (bottom row memory stack)",
                   handletextpad=0.5, columnspacing=1.2,
                   labelspacing=0.25, borderpad=0.3)

        out = Path(args.out).with_suffix(".pdf")
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(out)
        plt.close(fig)


if __name__ == "__main__":
    main()
