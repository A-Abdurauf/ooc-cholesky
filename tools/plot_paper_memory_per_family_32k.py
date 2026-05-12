#!/usr/bin/env python3
"""Paper-grade figure 2 (refined): Memory breakdown per sweep family at N=32768.

2×2 panel grid. One panel per sweep family (low-only / low+mid / legacy 3-tier /
ladder). Each panel: x = ε (4 clusters), 3 stacked bars per cluster
(tile / block128 / vec1D32) decomposed by datatype. Baseline IEEE GB shown
as a horizontal reference line per ε.

Cleaner than the 14-sweep monolithic version; one family at a time keeps the
bar count <= 12 per panel.
"""
import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D


mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 10.5,
    "axes.titlesize": 11.5,
    "axes.labelsize": 10.5,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
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
    ("FP64",          ["fp64_gb"]),
    ("FP32",          ["fp32_gb"]),
    ("SFP32",         ["mx_fp32_gb"]),
    ("FP16",          ["fp16_gb"]),
    ("SFP16",         ["mx_fp16_gb"]),
    ("SFP8 (E4M3)",   ["mx_e4m3_gb"]),
    ("FP8 (E4M3)",    ["fp8_e4m3_gb"]),
    ("SFP4 (E2M1)",   ["e2m1_gb"]),
    ("Scale meta",    ["scale_meta_gb"]),
]
DTYPE_COLOR = {
    "FP64":          "#0072B2",
    "FP32":          "#E69F00",
    "SFP32":         "#9B6814",
    "FP16":          "#D55E00",
    "SFP16":         "#7A2C00",
    "SFP8 (E4M3)":   "#009E73",
    "FP8 (E4M3)":    "#56B4E9",
    "SFP4 (E2M1)":   "#F0E442",
    "Scale meta":    "#999999",
}

EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]

# (panel_title, [(granularity_label, sweep_name), ...]) per sweep family.
FAMILIES = [
    ("Low only scaled  (SFP8 only; FP32 / FP16 plain)", [
        ("tile",        "requant_lowscale_tile_gt20k"),
        ("block 128",   "requant_lowscale_block128_gt20k"),
        ("vec1D 32",    "requant_lowscale_vec1d32_gt20k"),
    ]),
    ("Low + Mid scaled  (SFP8 + SFP16; FP32 plain)", [
        ("tile",        "requant_lowmidscale_tile_gt20k"),
        ("block 128",   "requant_lowmidscale_block128_gt20k"),
        ("vec1D 32",    "requant_lowmidscale_vec1d32_gt20k"),
    ]),
    ("All tiers scaled  (SFP32 + SFP16 + SFP8)", [
        ("tile",        "requant_legacy_scaled_tile_gt20k"),
        ("block 128",   "requant_legacy_scaled_block128_gt20k"),
        ("vec1D 32",    "requant_legacy_scaled_vec1d32_gt20k"),
    ]),
    ("Full ladder  (SFP4 → SFP8 → SFP16 → FP32)", [
        ("tile",        "requant_ladder_scaled_tile_gt20k"),
        ("block 128",   "requant_ladder_scaled_block128_gt20k"),
        ("vec1D 32",    "requant_ladder_scaled_vec1d32_gt20k"),
    ]),
]

BASELINE_SWEEP = "requant_baseline_fp8_subnormal_gt20k"
LADDER_IEEE_SWEEP = "requant_ladder_ieee_gt20k"


def load(csv_path, n_target):
    rows = {}
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            try:
                if int(r["n"]) != n_target: continue
            except: continue
            rows[(r["sweep"], r["source_epsilon"])] = r
    return rows


def bucket_value(r, keys):
    if r is None: return 0.0
    return sum(float(r.get(k, 0) or 0) for k in keys)


def total_gb(r):
    return float(r.get("total_gb", 0) or 0) if r else 0.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--n", type=int, default=32768)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    data = load(args.csv, args.n)

    fig, axes = plt.subplots(2, 2, figsize=(10, 6.0), sharey=True)
    axes = axes.flatten()

    n_eps = len(EPS_ORDER)
    n_gran = 3
    group_width = 0.78
    bar_w = group_width / n_gran
    x_centres = list(range(n_eps))

    seen_buckets = []
    max_y = 0.0

    # Reference GB per ε for the two IEEE-only references.
    baseline_gb = []
    ladder_ieee_gb = []
    for eps in EPS_ORDER:
        baseline_gb.append(total_gb(data.get((BASELINE_SWEEP, eps))))
        ladder_ieee_gb.append(total_gb(data.get((LADDER_IEEE_SWEEP, eps))))

    GRAN_HATCH = {"tile": "", "block 128": "///", "vec1D 32": "xx"}
    for panel_idx, (ax, (title, grans)) in enumerate(zip(axes, FAMILIES)):
        # Last panel (Full ladder) pairs with ladder IEEE; panels 0-2 pair with baseline IEEE.
        is_ladder_panel = (panel_idx == 3)
        for gi, (gran_lbl, sweep) in enumerate(grans):
            xs = [c - group_width/2 + (gi + 0.5) * bar_w for c in x_centres]
            bottoms = [0.0] * n_eps
            for dt_lbl, keys in DTYPE_BUCKETS:
                vals = []
                for eps in EPS_ORDER:
                    r = data.get((sweep, eps))
                    vals.append(bucket_value(r, keys))
                if not any(v > 0 for v in vals): continue
                if dt_lbl not in seen_buckets:
                    seen_buckets.append(dt_lbl)
                ax.bar(xs, vals, bar_w * 0.93, bottom=bottoms,
                       color=DTYPE_COLOR[dt_lbl], edgecolor="black", linewidth=0.3,
                       hatch=GRAN_HATCH.get(gran_lbl, ""))
                bottoms = [b + v for b, v in zip(bottoms, vals)]
            max_y = max(max_y, max(bottoms))

        # IEEE reference per ε cluster — baseline IEEE for panels 0-2 (Path A),
        # ladder IEEE for panel 3 (Full ladder, bound-based).
        ref_gb    = ladder_ieee_gb if is_ladder_panel else baseline_gb
        ref_color = "#666666"      if is_ladder_panel else "black"
        ref_style = ":"            if is_ladder_panel else "--"
        ref_width = 1.4            if is_ladder_panel else 1.1
        for c, b in zip(x_centres, ref_gb):
            if b > 0:
                ax.hlines(b, c - group_width/2 - 0.05, c + group_width/2 + 0.05,
                          colors=ref_color, linestyles=ref_style, linewidth=ref_width,
                          alpha=0.92, zorder=4)

        ax.set_xticks(x_centres)
        ax.set_xticklabels([f"ε = {e}" for e in EPS_ORDER], fontsize=9)
        ax.set_title(title, fontsize=9.5)
        ax.set_axisbelow(True)

    # Y-axis labels on left column only.
    for ax in axes[::2]:
        ax.set_ylabel("GB per Cholesky factor", fontsize=9)
    # Y-axis headroom + tighter tick label size.
    for ax in axes:
        ax.set_ylim(0, max_y * 1.10)
        ax.tick_params(axis="y", labelsize=8.5)

    # Two-row legend layout: dtypes on top row, granularity + ladder IEEE on bottom row.
    # Always show every datatype label even when no bar segment exercises it (e.g.,
    # plain FP8 (E4M3) is in the ladder IEEE reference line but not in any MX bar).
    dtype_handles = [Patch(facecolor=DTYPE_COLOR[label], edgecolor="black", label=label)
                     for label, _ in DTYPE_BUCKETS]
    gran_handles = [
        Patch(facecolor="white", edgecolor="black", hatch="",    label="tile"),
        Patch(facecolor="white", edgecolor="black", hatch="///", label="block 128"),
        Patch(facecolor="white", edgecolor="black", hatch="xx",  label="vec1D 32"),
    ]
    ieee_handles = [
        Line2D([0], [0], color="black",   linestyle="--", linewidth=1.2, label="baseline IEEE"),
        Line2D([0], [0], color="#666666", linestyle=":",  linewidth=1.4, label="ladder IEEE"),
    ]

    # Row 1 — datatypes, full-width.
    lg_dtype = fig.legend(handles=dtype_handles,
                          loc="lower center",
                          bbox_to_anchor=(0.5, -0.02),
                          ncol=len(dtype_handles),
                          frameon=True, framealpha=0.95, title="Datatype",
                          handletextpad=0.3, columnspacing=0.9, labelspacing=0.25,
                          borderpad=0.35)
    fig.add_artist(lg_dtype)

    # Row 2 — granularity (3 horizontal) and ladder IEEE reference, with a visible gap.
    lg_gran = fig.legend(handles=gran_handles,
                         loc="lower center",
                         bbox_to_anchor=(0.36, -0.10),
                         ncol=3, frameon=True, framealpha=0.95,
                         title="Granularity (bars left → right per ε)",
                         handletextpad=0.3, columnspacing=1.0, labelspacing=0.25,
                         borderpad=0.35)
    fig.add_artist(lg_gran)
    fig.legend(handles=ieee_handles,
               loc="lower center",
               bbox_to_anchor=(0.78, -0.10),
               ncol=2, frameon=True, framealpha=0.95,
               title="IEEE references",
               handletextpad=0.3, columnspacing=1.0, borderpad=0.35, labelspacing=0.25)

    fig.tight_layout(rect=[0, 0.12, 1, 0.99])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pdf", ".png"):
        p = out.with_suffix(ext)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(p)
    plt.close(fig)


if __name__ == "__main__":
    main()
