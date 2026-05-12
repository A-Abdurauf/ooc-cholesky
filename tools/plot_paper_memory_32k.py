#!/usr/bin/env python3
"""Paper-grade figure 2: Memory breakdown by datatype at N=32768.

X-axis: ε (loosest left). For each ε, one stacked bar per sweep family showing
the GB of each data-type bucket. Lines comparing the 4 families ×  3 granularities
plus the 2 reference baselines.

Reads from the merged memory CSV produced by build_requant_gt20k_memory_csv.py.
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
    "xtick.labelsize": 9.5,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "legend.title_fontsize": 10,
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
})

DTYPE_BUCKETS = [
    ("FP64",                  ["fp64_gb"]),
    ("FP32",                  ["fp32_gb"]),
    ("Scaled FP32",           ["mx_fp32_gb"]),
    ("FP16 / BF16",           ["fp16_gb", "bf16_gb"]),
    ("Scaled FP16",           ["mx_fp16_gb"]),
    ("Scaled FP8 (E4M3)",     ["mx_e4m3_gb"]),
    ("FP8 plain (E4M3)",      ["fp8_e4m3_gb"]),
    ("Scaled FP4 (E2M1)",     ["e2m1_gb"]),
    ("Scale meta",            ["scale_meta_gb"]),
]
DTYPE_COLOR = {
    "FP64":                 "#0072B2",
    "FP32":                 "#E69F00",
    "Scaled FP32":          "#9B6814",
    "FP16 / BF16":          "#D55E00",
    "Scaled FP16":          "#7A2C00",
    "Scaled FP8 (E4M3)":    "#009E73",
    "FP8 plain (E4M3)":     "#56B4E9",
    "Scaled FP4 (E2M1)":    "#F0E442",
    "Scale meta":           "#999999",
}

EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]

# Sweeps in the figure (no FTZ, no probes, no archived).
SWEEPS = [
    ("baseline IEEE (subn)",          "requant_baseline_fp8_subnormal_gt20k"),
    ("ladder IEEE (no MX)",           "requant_ladder_ieee_gt20k"),
    ("low only · tile",               "requant_lowscale_tile_gt20k"),
    ("low only · block128",           "requant_lowscale_block128_gt20k"),
    ("low only · vec1D32",            "requant_lowscale_vec1d32_gt20k"),
    ("low+mid · tile",                "requant_lowmidscale_tile_gt20k"),
    ("low+mid · block128",            "requant_lowmidscale_block128_gt20k"),
    ("low+mid · vec1D32",             "requant_lowmidscale_vec1d32_gt20k"),
    ("legacy 3-tier · tile",          "requant_legacy_scaled_tile_gt20k"),
    ("legacy 3-tier · block128",      "requant_legacy_scaled_block128_gt20k"),
    ("legacy 3-tier · vec1D32",       "requant_legacy_scaled_vec1d32_gt20k"),
    ("ladder · tile",                 "requant_ladder_scaled_tile_gt20k"),
    ("ladder · block128",             "requant_ladder_scaled_block128_gt20k"),
    ("ladder · vec1D32",              "requant_ladder_scaled_vec1d32_gt20k"),
]


def load(csv_path, n_target):
    rows = {}
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            try:
                if int(r["n"]) != n_target: continue
            except: continue
            key = (r["sweep"], r["source_epsilon"])
            rows[key] = r
    return rows


def bucket_value(r, keys):
    if r is None: return 0.0
    return sum(float(r.get(k, 0) or 0) for k in keys)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--n", type=int, default=32768)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    data = load(args.csv, args.n)
    n_eps = len(EPS_ORDER)
    n_sw = len(SWEEPS)

    # Layout: x = ε (4 groups), within each group n_sw stacked bars.
    fig_w = max(15, 0.55 * n_sw * n_eps + 4)
    fig, ax = plt.subplots(figsize=(fig_w, 7.0))

    group_width = 0.88
    bar_w = group_width / n_sw
    x_centres = list(range(n_eps))

    seen_buckets = []
    max_total = 0.0
    for si, (label_sw, sweep) in enumerate(SWEEPS):
        xs = [c - group_width/2 + (si + 0.5) * bar_w for c in x_centres]
        bottoms = [0.0] * n_eps
        for dt_label, keys in DTYPE_BUCKETS:
            vals = []
            for eps in EPS_ORDER:
                r = data.get((sweep, eps))
                vals.append(bucket_value(r, keys))
            if not any(v > 0 for v in vals): continue
            if dt_label not in seen_buckets:
                seen_buckets.append(dt_label)
            ax.bar(xs, vals, bar_w * 0.94, bottom=bottoms,
                   color=DTYPE_COLOR[dt_label], edgecolor="black", linewidth=0.25)
            bottoms = [b + v for b, v in zip(bottoms, vals)]
        max_total = max(max_total, max(bottoms) if bottoms else 0)

    # Sweep label row underneath (rotated for legibility).
    tick_pos = []
    tick_lbl = []
    for c in x_centres:
        for si, (label_sw, _) in enumerate(SWEEPS):
            tick_pos.append(c - group_width/2 + (si + 0.5) * bar_w)
            tick_lbl.append(label_sw)
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(tick_lbl, rotation=90, ha="center", va="top", fontsize=7.5)

    # Major group labels (ε) at the top of each cluster.
    for c, eps in zip(x_centres, EPS_ORDER):
        ax.text(c, max_total * 1.10, rf"$\varepsilon = {eps}$",
                ha="center", va="bottom", fontsize=12, fontweight="bold")

    ax.set_ylim(0, max_total * 1.18)
    ax.set_ylabel(f"Memory footprint of Cholesky factor (GB), N={args.n}")
    ax.set_xlabel("")  # group labels already above; sweep labels below
    ax.set_axisbelow(True)
    ax.set_title(
        f"Memory breakdown by datatype  —  weak_32k covariance, NB=2048\n"
        f"4 source epsilons × {n_sw} sweep configurations"
    )

    # Datatype legend (colour swatches), placed below x-tick labels.
    handles = [Patch(facecolor=DTYPE_COLOR[b], edgecolor="black", label=b)
               for b in seen_buckets]
    fig.legend(handles=handles, loc="lower center", ncol=min(len(seen_buckets), 5),
               bbox_to_anchor=(0.5, -0.04), frameon=True, framealpha=0.95,
               title="Datatype")

    fig.tight_layout(rect=[0, 0.04, 1, 1])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pdf", ".png"):
        p = out.with_suffix(ext)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(p)
    plt.close(fig)


if __name__ == "__main__":
    main()
