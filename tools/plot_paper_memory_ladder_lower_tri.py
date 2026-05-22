#!/usr/bin/env python3
"""Memory figure for the ladder mode at 3 matrix sizes -- one combined figure.

Independent of the NB=4096 Recreate sweep: this script reads the main
ladder/baseline CSV (NB=2048 at 32k, NB=4096 at 40k and 65k) and the
ladder_rerun_32k CSV (for the 32k vec1D-32 GU value).

Layout: one figure with 3 subplots side-by-side (32k / 40k / 65k).  Within
each subplot, 4 epsilon clusters (1e-5, 1e-6, 1e-7, 1e-8), each cluster
showing 4 stacked memory bars:

    1. Baseline                    (no MX; subnormal-allowed)
    2. Ladder IEEE                 (IEEE rungs only)
    3. Ladder MX+MXFP16 (tile)     (tile-level MX scaling)
    4. Ladder MX+MXFP16 (vec1D-32) (vec1D-32 MX scaling)

For 32k the vec1D-32 bar reads ladder_full_gu from the rerun (GU underflow,
patched binary).  For 40k / 65k the vec1D-32 bar comes from the main CSV.

Storage: lower triangle + half-diagonal (FP64 diagonals store nb*(nb+1)/2
elements each).  PDF only.
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
    "legend.fontsize": 9.5,
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

MAIN_CSV  = "/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv"
RERUN_CSV = "/home/abduraa/MX_project/ooc-cholesky/ladder_rerun_32k/results.csv"

BYTES_PER_ELEM = {
    "fp64":     8.0,
    "fp32":     4.0,
    "mx_fp16":  2.0,
    "fp16":     2.0,
    "mx_e4m3":  1.0,
    "mx_e5m2":  1.0,
    "fp8_e4m3": 1.0,
    "fp8_e5m2": 1.0,
    "e3m2":     1.0,
    "e2m3":     1.0,
    "e2m1":     0.5,
    "mx_fp32":  4.0,
}
HAS_SHARED_SCALE = {"mx_fp16", "mx_e4m3", "mx_e5m2", "mx_fp32", "e2m1", "e3m2", "e2m3"}

TILE_ALIASES = {
    "fp16":     "mx_fp16",
    "mx_fp32":  "fp32",
    "fp8_e5m2": "fp8_e4m3",
    "mx_e5m2":  "mx_e4m3",
}

DTYPE_STACK = [
    ("FP64",             "fp64",       "#0072B2"),
    ("FP32",             "fp32",       "#E69F00"),
    ("FP16 / MXFP16",    "mx_fp16",    "#7A2C00"),
    ("FP8 plain (E4M3)", "fp8_e4m3",   "#56B4E9"),
    ("MXFP8 (E4M3)",     "mx_e4m3",    "#009E73"),
    ("MXFP4 (E2M1)",     "e2m1",       "#F0E442"),
    ("Scale meta",       "scale_meta", "#999999"),
]


def bars_for_n(n):
    common = [
        ("Baseline",                     "main", "requant_baseline_fp8_subnormal_gt20k"),
        ("Ladder IEEE",                  "main", "requant_ladder_ieee_gt20k"),
        ("Ladder MX+MXFP16 (tile)",      "main", "requant_ladder_scaled_tile_gt20k"),
    ]
    if n == 32768:
        common.append(("Ladder MX+MXFP16 (vec1D-32)",
                       "rerun", "ladder_full_gu"))
    else:
        common.append(("Ladder MX+MXFP16 (vec1D-32)",
                       "main",  "requant_ladder_scaled_vec1d32_gt20k"))
    return common


def parse_tile_counts_full(s):
    out = {}
    for part in (s or "").split(","):
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        k = k.strip().lower()
        k = TILE_ALIASES.get(k, k)
        try:
            out[k] = out.get(k, 0) + int(v.strip())
        except ValueError:
            pass
    return out


def parse_tile_breakdown(s):
    out = {}
    for part in (s or "").strip().strip('"').split(";"):
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        k = k.strip().lower()
        k = TILE_ALIASES.get(k, k)
        try:
            out[k] = out.get(k, 0) + int(v.strip())
        except ValueError:
            pass
    return out


def load_main_rows(csv_path, n, nb):
    rows = {}
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            try:
                if int(r["n"]) != n or int(r["nb"]) != nb:
                    continue
            except (ValueError, KeyError):
                continue
            rows[(r["sweep"], r["source_epsilon"])] = r
    return rows


def load_rerun_rows(csv_path, n, nb):
    rows = {}
    p = Path(csv_path)
    if not p.exists():
        return rows
    with p.open() as f:
        for r in csv.DictReader(f):
            try:
                if int(r["n"]) != n or int(r["nb"]) != nb:
                    continue
            except (ValueError, KeyError):
                continue
            rows[(r["sweep"], r["source_epsilon"])] = r
    return rows


def lower_tri_counts(row, source, n, nb):
    if row is None:
        return {}
    M = n // nb
    if source == "rerun":
        return parse_tile_breakdown(row.get("tile_breakdown", ""))
    full = parse_tile_counts_full(row.get("tile_counts_full", ""))
    out = {}
    for fmt, v in full.items():
        if fmt == "fp64":
            off = max(v - M, 0)
            out[fmt] = M + off // 2
        else:
            out[fmt] = v // 2
    return out


def lower_tri_bytes_by_dtype(counts_lower, n, nb):
    M = n // nb
    tile_elems_full = nb * nb
    tile_elems_diag = nb * (nb + 1) // 2
    by_dtype = {fmt: 0.0 for fmt in BYTES_PER_ELEM}
    scale_bytes_total = 0.0
    for fmt, count in counts_lower.items():
        bpe = BYTES_PER_ELEM.get(fmt, 0.0)
        if fmt == "fp64":
            diag = min(M, count)
            off  = max(count - diag, 0)
            data_bytes = diag * tile_elems_diag * bpe + off * tile_elems_full * bpe
        else:
            data_bytes = count * tile_elems_full * bpe
        by_dtype[fmt] = by_dtype.get(fmt, 0.0) + data_bytes
        if fmt in HAS_SHARED_SCALE and bpe > 0:
            scale_bytes_total += count * (tile_elems_full // 32)
    GB = 1.0 / (1024 ** 3)
    out_gb = {fmt: by_dtype[fmt] * GB for fmt in by_dtype}
    out_gb["scale_meta"] = scale_bytes_total * GB
    return out_gb


def _fmt_n(n):
    return f"{n // 1024}k" if n % 1024 == 0 else str(n)


def draw_subplot(ax, n, main_rows, rerun_rows, drawn_dtypes, show_ylabel):
    nb = NB_BY_N[n]
    bars = bars_for_n(n)
    n_eps  = len(EPS_ORDER)
    n_bars = len(bars)
    group_w = 0.84
    bar_w   = group_w / n_bars
    x_centres = list(range(n_eps))

    max_total = 0.0
    for bi, (label, src, sweep) in enumerate(bars):
        for ci, eps in zip(x_centres, EPS_ORDER):
            x = ci - group_w/2 + (bi + 0.5) * bar_w
            row = (main_rows if src == "main" else rerun_rows).get((sweep, eps))
            counts = lower_tri_counts(row, src, n, nb)
            gb_by_dtype = lower_tri_bytes_by_dtype(counts, n, nb)
            bottom = 0.0
            for dt_lbl, dt_key, dt_color in DTYPE_STACK:
                v = gb_by_dtype.get(dt_key, 0.0)
                if v <= 0:
                    continue
                drawn_dtypes.add(dt_lbl)
                ax.bar([x], [v], bar_w * 0.92,
                       bottom=[bottom],
                       color=dt_color, edgecolor="black", linewidth=0.3)
                bottom += v
            if bottom > 0:
                ax.text(x, bottom * 1.012, f"{bottom:.2f}",
                        ha="center", va="bottom",
                        fontsize=6.8, rotation=90, color="#222")
            max_total = max(max_total, bottom)

    ax.set_xticks(x_centres)
    ax.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                        for e in EPS_ORDER])
    ax.set_xlim(-0.5, n_eps - 0.5)
    ax.set_ylim(0, max_total * 1.20 if max_total > 0 else 1.0)
    if show_ylabel:
        ax.set_ylabel("Memory  (GB, lower triangle + half-diagonal)")
    ax.set_axisbelow(True)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, alpha=0.3, linewidth=0.5)

    ax.text(0.985, 0.965, f"N = {n // 1024}k\nNB = {nb}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=10.5,
            bbox=dict(boxstyle="round,pad=0.32",
                      facecolor="white", edgecolor="#888",
                      linewidth=0.8, alpha=0.92))
    return bars


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ns", nargs="+", type=int, default=[32768, 40960, 65536])
    ap.add_argument("--out",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures/paper_memory_ladder_lower_tri")
    args = ap.parse_args()

    n_panels = len(args.ns)
    fig, axes = plt.subplots(1, n_panels,
                             figsize=(5.6 * n_panels + 0.6, 5.6),
                             sharey=False)
    if n_panels == 1:
        axes = [axes]

    drawn_dtypes = set()
    last_bars = None
    for ax, n in zip(axes, args.ns):
        nb = NB_BY_N[n]
        main_rows  = load_main_rows(MAIN_CSV, n, nb)
        rerun_rows = load_rerun_rows(RERUN_CSV, n, nb)
        last_bars = draw_subplot(ax, n, main_rows, rerun_rows,
                                 drawn_dtypes,
                                 show_ylabel=(ax is axes[0]))

    dtype_handles = [Patch(facecolor=c, edgecolor="black", label=l)
                     for l, _, c in DTYPE_STACK if l in drawn_dtypes]
    bar_labels = [lbl for (lbl, _, _) in last_bars]
    bar_handles = [Patch(facecolor="white", edgecolor="black",
                         label=f"{i+1}. {lbl}")
                   for i, lbl in enumerate(bar_labels)]

    leg1 = fig.legend(handles=bar_handles,
                      loc="lower center", bbox_to_anchor=(0.25, -0.09),
                      ncol=2, frameon=True, framealpha=0.95,
                      title="Bars within each cluster (left -> right)",
                      handletextpad=0.6, columnspacing=1.6,
                      labelspacing=0.3, borderpad=0.4)
    fig.add_artist(leg1)
    fig.legend(handles=dtype_handles,
               loc="lower center", bbox_to_anchor=(0.75, -0.09),
               ncol=4, frameon=True, framealpha=0.95,
               title="Tile-format stack (bottom -> top)",
               handletextpad=0.5, columnspacing=1.4,
               labelspacing=0.3, borderpad=0.4)

    fig.tight_layout(rect=[0, 0.06, 1, 1.0])

    out = Path(args.out).with_suffix(".pdf")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(out)
    plt.close(fig)


if __name__ == "__main__":
    main()
