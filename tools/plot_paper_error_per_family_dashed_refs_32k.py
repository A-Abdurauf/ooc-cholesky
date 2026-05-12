#!/usr/bin/env python3
"""Paper-grade figure: Relative factorisation error per sweep family at N=32768.

Mirrors plot_paper_memory_per_family_32k.py: IEEE references are drawn as
horizontal dashed lines per ε cluster rather than bars. Each ε cluster shows
3 granularity bars (tile / block 128 / vec1D 32) for that family.

2×2 panel grid (one per sweep family). Y-axis log-scaled rel_factor_error.
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

EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]

BASELINE_SWEEP    = "requant_baseline_fp8_subnormal_gt20k"
LADDER_IEEE_SWEEP = "requant_ladder_ieee_gt20k"

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

# Granularity bar colours (consistent with the error-with-baselines figure).
GRAN_COLOR = {
    "tile":      "#f28e2b",  # orange
    "block 128": "#59a14f",  # green
    "vec1D 32":  "#b07aa1",  # purple
}
GRAN_HATCH = {"tile": "", "block 128": "///", "vec1D 32": "xx"}


def load(csv_path, n_target):
    rows = {}
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            try:
                if int(r["n"]) != n_target: continue
            except: continue
            rows[(r["sweep"], r["source_epsilon"])] = r
    return rows


def err(r):
    if r is None: return None
    try:
        v = float(r["rel_factor_error"])
        return v if v > 0 else None
    except (ValueError, TypeError, KeyError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--n", type=int, default=32768)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    data = load(args.csv, args.n)

    fig, axes = plt.subplots(2, 2, figsize=(10, 4.5), sharey=True)
    axes = axes.flatten()

    n_eps = len(EPS_ORDER)
    n_gran = 3
    group_width = 0.78
    bar_w = group_width / n_gran
    x_centres = list(range(n_eps))

    # Reference error values per ε (used for dashed lines).
    baseline_err = [err(data.get((BASELINE_SWEEP,    e))) for e in EPS_ORDER]
    ladder_ieee_err = [err(data.get((LADDER_IEEE_SWEEP, e))) for e in EPS_ORDER]

    y_min, y_max = float("inf"), 0.0

    for panel_idx, (ax, (title, grans)) in enumerate(zip(axes, FAMILIES)):
        # Last panel (index 3) is the Full ladder: it pairs with ladder IEEE.
        # Panels 0-2 use Path A (ε-cutoff) so they pair with baseline IEEE.
        is_ladder_panel = (panel_idx == 3)
        # Granularity bars per ε cluster.
        for gi, (gran_lbl, sweep) in enumerate(grans):
            xs = [c - group_width/2 + (gi + 0.5) * bar_w for c in x_centres]
            ys = []
            for eps in EPS_ORDER:
                v = err(data.get((sweep, eps)))
                ys.append(v if v is not None else 0.0)
            ax.bar(xs, ys, bar_w * 0.92,
                   color=GRAN_COLOR[gran_lbl], edgecolor="black", linewidth=0.3,
                   hatch=GRAN_HATCH[gran_lbl])
            for x, y in zip(xs, ys):
                if y > 0:
                    y_min = min(y_min, y)
                    y_max = max(y_max, y)
                    ax.text(x, y * 1.18, f"{y:.1e}",
                            ha="center", va="bottom", fontsize=6.0, rotation=90,
                            color="#222")

        # IEEE reference dashed line — baseline IEEE for panels 0-2, ladder IEEE for panel 3.
        ref_err = ladder_ieee_err if is_ladder_panel else baseline_err
        ref_color = "#666666" if is_ladder_panel else "black"
        ref_style = ":"        if is_ladder_panel else "--"
        ref_width = 1.4        if is_ladder_panel else 1.1
        for c, v in zip(x_centres, ref_err):
            if v is not None and v > 0:
                ax.hlines(v, c - group_width/2 - 0.05, c + group_width/2 + 0.05,
                          colors=ref_color, linestyles=ref_style, linewidth=ref_width,
                          alpha=0.92, zorder=4)
                y_min = min(y_min, v)
                y_max = max(y_max, v)

        # Requested-ε target reference (very faint grey dashed).
        for c, eps in zip(x_centres, EPS_ORDER):
            ev = float(eps)
            ax.hlines(ev, c - group_width/2 - 0.05, c + group_width/2 + 0.05,
                      colors="lightgray", linestyles=":", linewidth=0.7, alpha=0.7,
                      zorder=0)

        ax.set_xticks(x_centres)
        ax.set_xticklabels([f"ε = {e}" for e in EPS_ORDER], fontsize=9)
        ax.set_title(title, fontsize=9.5)
        ax.set_yscale("log")
        ax.set_axisbelow(True)

    if y_min < float("inf"):
        for ax in axes:
            ax.set_ylim(y_min / 3, y_max * 12)
            ax.tick_params(axis="y", labelsize=8.5)

    # Y-axis labels on left column only.
    for ax in axes[::2]:
        ax.set_ylabel(r"$\|LL^\top - A\|_F\,/\,\|A\|_F$  (log)", fontsize=9)

    # Bottom legends — granularity hatches (3) and IEEE references (2 lines, stacked).
    gran_handles = [
        Patch(facecolor=GRAN_COLOR["tile"],      edgecolor="black", hatch="",    label="tile"),
        Patch(facecolor=GRAN_COLOR["block 128"], edgecolor="black", hatch="///", label="block 128"),
        Patch(facecolor=GRAN_COLOR["vec1D 32"],  edgecolor="black", hatch="xx",  label="vec1D 32"),
    ]
    ieee_handles = [
        Line2D([0], [0], color="black",   linestyle="--", linewidth=1.2, label="baseline IEEE"),
        Line2D([0], [0], color="#666666", linestyle=":",  linewidth=1.4, label="ladder IEEE"),
    ]

    lg_gran = fig.legend(handles=gran_handles,
                         loc="lower center",
                         bbox_to_anchor=(0.36, -0.04),
                         ncol=3, frameon=True, framealpha=0.95,
                         title="Granularity (bars left → right per ε)",
                         handletextpad=0.3, columnspacing=1.0, labelspacing=0.25,
                         borderpad=0.35)
    fig.add_artist(lg_gran)
    fig.legend(handles=ieee_handles,
               loc="lower center",
               bbox_to_anchor=(0.78, -0.04),
               ncol=2, frameon=True, framealpha=0.95,
               title="IEEE references",
               handletextpad=0.3, columnspacing=1.0, borderpad=0.35, labelspacing=0.25)

    fig.tight_layout(rect=[0, 0.05, 1, 0.99])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pdf", ".png"):
        p = out.with_suffix(ext)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(p)
    plt.close(fig)


if __name__ == "__main__":
    main()
