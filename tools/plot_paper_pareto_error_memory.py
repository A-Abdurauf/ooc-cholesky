#!/usr/bin/env python3
"""Pareto scatter: error vs memory for the Full MX ladder, four
(granularity, underflow) combinations.

4 epsilon subplots in one row.  Within each panel:
  - x-axis: memory (GB) on log scale, lower-tri + half-diagonal accounting
  - y-axis: relative factor error on log scale
  - 4 colors: tile+FZ, tile+GU, MX+FZ, MX+GU
  - 3 marker shapes: 32k (circle), 40k (square), 64k (triangle)
  - faint lines connecting the 4 combos at the same N -- visualises the
    trade-off shape per N
Lower-left = better.

Reads the merged ladder-only CSV (already lower-tri, canonical format
order).  PDF only.
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

CSV = "/home/abduraa/MX_project/ooc-cholesky/figures/all_error_runs_ladder_only.csv"

EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]
NS  = [32768, 40960, 65536]
NB_BY_N = {32768: 2048, 40960: 4096, 65536: 4096}
N_LABEL = {32768: "32k", 40960: "40k", 65536: "64k"}

# (gran, uflow) -> (label, color, hatch)
COMBOS = [
    ("tile", "fz", "tile + FZ", "#4e79a7"),
    ("tile", "gu", "tile + GU", "#9ec5db"),
    ("MX",   "fz", "MX + FZ",   "#b07aa1"),
    ("MX",   "gu", "MX + GU",   "#7a3d72"),
]
N_MARKER = {32768: "o", 40960: "s", 65536: "^"}

BPE = {"FP64":8,"FP32":4,"FP16":2,"MXFP16":2,"FP8":1,"MXFP8":1,"MXFP4":0.5}
SCALE_BEARING = {"MXFP16","MXFP8","MXFP4"}


def parse_bd(s):
    out = {}
    for part in (s or "").rstrip(";").split(";"):
        if "=" not in part: continue
        k, v = part.split("=", 1)
        try: out[k.strip()] = out.get(k.strip(), 0) + int(v.strip())
        except: pass
    return out


def memory_gb(counts, n, nb):
    M = n // nb
    tile_full = nb * nb
    tile_diag = nb * (nb + 1) // 2
    GB = 1 / 1024 ** 3
    total = 0.0
    for f, cnt in counts.items():
        bpe = BPE.get(f, 0)
        if f == "FP64":
            d = min(M, cnt); o = max(cnt - d, 0)
            total += (d * tile_diag + o * tile_full) * bpe
        else:
            total += cnt * tile_full * bpe
        if f in SCALE_BEARING:
            total += cnt * (tile_full // 32)
    return total * GB


def load_rows():
    return list(csv.DictReader(open(CSV)))


def lookup(rows, mode, gran, uflow, n, nb, eps):
    for r in rows:
        if (r["mode"] == mode and r["granularity"] == gran
                and r["underflow"] == uflow
                and int(r["N"]) == n and int(r["nb"]) == nb
                and r["eps"] == eps):
            try:
                err = float(r["rel_factor_error"])
            except (TypeError, ValueError):
                return None
            mem = memory_gb(parse_bd(r["tile_breakdown"]), n, nb)
            return err, mem
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="Ladder MX+MXFP16",
                    help="Mode column value to plot (full MX ladder by default).")
    ap.add_argument("--out",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures/paper_pareto_error_memory")
    args = ap.parse_args()

    rows = load_rows()

    with mpl.rc_context(PAPER_RC):
        fig, axes = plt.subplots(1, 4, figsize=(8.5, 2.8),
                                 sharey=True)

        all_err = []
        all_mem = []
        for ai, (ax, eps) in enumerate(zip(axes, EPS_ORDER)):
            # First pass to find data
            for n in NS:
                nb = NB_BY_N[n]
                # Collect the 4 combo points for this N
                pts = []
                for gran, uflow, _, _ in COMBOS:
                    v = lookup(rows, args.mode, gran, uflow, n, nb, eps)
                    pts.append(v)
                # Draw the connecting line (only across present points)
                xs = [p[1] for p in pts if p]
                ys = [p[0] for p in pts if p]
                if len(xs) >= 2:
                    ax.plot(xs, ys, color="#bbbbbb", linewidth=0.6,
                            alpha=0.7, zorder=1)
                # Draw markers per combo
                for (gran, uflow, lbl, color), p in zip(COMBOS, pts):
                    if p is None:
                        continue
                    err, mem = p
                    ax.scatter([mem], [err], s=42,
                               marker=N_MARKER[n],
                               facecolor=color, edgecolor="black",
                               linewidth=0.5, zorder=3)
                    all_err.append(err); all_mem.append(mem)

            ax.set_xscale("log"); ax.set_yscale("log")
            ax.set_xlabel("Memory (GB)")
            ax.set_axisbelow(True)
            ax.grid(True, which="both", alpha=0.22, linewidth=0.4)
            ax.text(0.04, 0.965,
                    f"({chr(ord('a') + ai)})  $\\varepsilon = {EPS_DISPLAY[eps][1:-1]}$",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=8.5, fontweight="bold")

        axes[0].set_ylabel("Relative factor error")

        # Shared symmetric padding on x.
        if all_mem:
            xmin, xmax = min(all_mem), max(all_mem)
            for ax in axes:
                ax.set_xlim(xmin * 0.78, xmax * 1.28)
        if all_err:
            for ax in axes:
                ax.set_ylim(min(all_err) * 0.4, max(all_err) * 2.5)

        # Legends.
        combo_handles = [Line2D([0], [0], marker="o", linestyle="",
                                markerfacecolor=color, markeredgecolor="black",
                                markeredgewidth=0.5, markersize=7, label=lbl)
                         for (_, _, lbl, color) in COMBOS]
        n_handles = [Line2D([0], [0], marker=N_MARKER[n], linestyle="",
                            markerfacecolor="#aaaaaa",
                            markeredgecolor="black", markeredgewidth=0.5,
                            markersize=7, label=N_LABEL[n])
                     for n in NS]
        line_handle = Line2D([0], [0], color="#bbbbbb", linewidth=0.9,
                             label="connect 4 combos at same $N$")

        fig.subplots_adjust(left=0.075, right=0.985, top=0.95,
                            bottom=0.32, wspace=0.10)
        leg1 = fig.legend(handles=combo_handles,
                          loc="upper center", bbox_to_anchor=(0.27, 0.20),
                          ncol=4, frameon=True, framealpha=0.95,
                          title="(granularity, underflow)",
                          handletextpad=0.4, columnspacing=1.1,
                          labelspacing=0.25, borderpad=0.3)
        fig.add_artist(leg1)
        fig.legend(handles=n_handles + [line_handle],
                   loc="upper center", bbox_to_anchor=(0.78, 0.20),
                   ncol=4, frameon=True, framealpha=0.95,
                   title="marker shape = $N$",
                   handletextpad=0.4, columnspacing=1.1,
                   labelspacing=0.25, borderpad=0.3)

        out = Path(args.out).with_suffix(".pdf")
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(out)
        plt.close(fig)


if __name__ == "__main__":
    main()
