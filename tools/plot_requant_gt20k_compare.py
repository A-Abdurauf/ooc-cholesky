#!/usr/bin/env python3
"""Comparison plots for the requant_gt20k sweeps at a given bin.

Three views per metric:
  - LEGACY:   baseline_fp8 + (legacy + lowmidscale + lowscale)
              "Legacy is a drop-in replacement for IEEE bucketing."
  - LADDER:   baseline_fp8 + (ladder_vec1d32 + ladder_block128)
  - COMBINED: all 9 sweeps in one figure with the same colour palette.

Two metrics per view:
  - MEMORY (stacked GB, broken down by datatype)
  - ERROR  (single bars, rel_factor_error, log y)

X-axis = source_epsilon, ordered 1e-5 -> 1e-6 -> 1e-7 -> 1e-8 (loosest first).
Within each epsilon cluster, the reference sweep (baseline_fp8) is leftmost.
"""
import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# Datatype buckets for the stacked memory bars.
DTYPE_BUCKETS = [
    ("FP64",                 ["fp64"]),
    ("FP32",                 ["fp32"]),
    ("Scaled FP32",          ["mx_fp32"]),
    ("FP16 / BF16",          ["fp16", "bf16"]),
    ("Scaled FP16",          ["mx_fp16"]),
    ("Scaled FP8 (E4M3)",    ["mx_e4m3"]),
    ("Scaled FP8 (E5M2)",    ["mx_e5m2"]),
    ("FP8 plain (E4M3/E5M2)", ["fp8_e4m3", "fp8_e5m2"]),
    ("Scaled FP6 (E3M2/E2M3)", ["e3m2", "e2m3"]),
    ("Scaled FP4 (E2M1)",    ["e2m1"]),
    ("Scale meta",           ["scale_meta"]),
]
DTYPE_COLOR = {
    "FP64":                  "#4e79a7",
    "FP32":                  "#f28e2b",
    "Scaled FP32":           "#bd6f25",
    "FP16 / BF16":           "#e15759",
    "Scaled FP16":           "#962f30",
    "Scaled FP8 (E4M3)":     "#76b7b2",
    "Scaled FP8 (E5M2)":     "#3a7975",
    "FP8 plain (E4M3/E5M2)": "#17becf",
    "Scaled FP6 (E3M2/E2M3)":"#59a14f",
    "Scaled FP4 (E2M1)":     "#edc948",
    "Scale meta":            "#9c755f",
}

# Per-sweep stable colours (used for the error bars).
SWEEP_ORDER_ALL = [
    "requant_baseline_fp8_subnormal_gt20k",
    "requant_ladder_ieee_gt20k",
    "requant_legacy_scaled_tile_gt20k",
    "requant_legacy_scaled_block128_gt20k",
    "requant_legacy_scaled_vec1d32_gt20k",
    "requant_lowmidscale_tile_gt20k",
    "requant_lowmidscale_block128_gt20k",
    "requant_lowmidscale_vec1d32_gt20k",
    "requant_lowscale_tile_gt20k",
    "requant_lowscale_block128_gt20k",
    "requant_lowscale_vec1d32_gt20k",
    "requant_lowscale_e2m1_tile_gt20k",
    "requant_lowscale_e2m1_block128_gt20k",
    "requant_lowscale_e2m1_vec1d32_gt20k",
    "requant_lowmidscale_e2m1_tile_gt20k",
    "requant_lowmidscale_e2m1_block128_gt20k",
    "requant_lowmidscale_e2m1_vec1d32_gt20k",
    "requant_ladder_scaled_tile_gt20k",
    "requant_ladder_scaled_block128_gt20k",
    "requant_ladder_scaled_vec1d32_gt20k",
]
SWEEP_SHORT = {
    "requant_baseline_fp8_gt20k":                    "baseline IEEE (FTZ)",
    "requant_baseline_fp8_subnormal_gt20k":          "baseline IEEE (OCP subnormal)",
    "requant_ladder_ieee_gt20k":                     "ladder IEEE (no MX)",
    "requant_legacy_scaled_tile_gt20k":              "legacy 3-tier · tile",
    "requant_legacy_scaled_block128_gt20k":          "legacy 3-tier · block128",
    "requant_legacy_scaled_vec1d32_gt20k":           "legacy 3-tier · vec1D32",
    "requant_legacy_scaled_vec1d32_FP32tile_gt20k":  "legacy 3-tier · vec1D (FP32 tile*)",
    "requant_lowmidscale_tile_gt20k":                "low+mid scaled · tile",
    "requant_lowmidscale_block128_gt20k":            "low+mid scaled · block128",
    "requant_lowmidscale_vec1d32_gt20k":             "low+mid scaled · vec1D32",
    "requant_lowscale_tile_gt20k":                   "low only scaled · tile",
    "requant_lowscale_block128_gt20k":               "low only scaled · block128",
    "requant_lowscale_vec1d32_gt20k":                "low only scaled · vec1D32",
    "requant_lowscale_e2m1_tile_gt20k":              "low only scaled (SFP4) · tile",
    "requant_lowscale_e2m1_block128_gt20k":          "low only scaled (SFP4) · block128",
    "requant_lowscale_e2m1_vec1d32_gt20k":           "low only scaled (SFP4) · vec1D32",
    "requant_lowmidscale_e2m1_tile_gt20k":           "low+mid scaled (SFP4) · tile",
    "requant_lowmidscale_e2m1_block128_gt20k":       "low+mid scaled (SFP4) · block128",
    "requant_lowmidscale_e2m1_vec1d32_gt20k":        "low+mid scaled (SFP4) · vec1D32",
    "requant_ladder_scaled_tile_gt20k":              "ladder · tile",
    "requant_ladder_scaled_block128_gt20k":          "ladder · block128",
    "requant_ladder_scaled_vec1d32_gt20k":           "ladder · vec1D32",
}
SWEEP_COLOR = {
    "requant_baseline_fp8_gt20k":                    "#000000",  # reference (FTZ)
    "requant_baseline_fp8_subnormal_gt20k":          "#555555",  # reference (subnormal)
    "requant_ladder_ieee_gt20k":                     "#888888",  # reference (IEEE ladder)
    "requant_legacy_scaled_tile_gt20k":              "#08306b",
    "requant_legacy_scaled_block128_gt20k":          "#1f77b4",
    "requant_legacy_scaled_vec1d32_gt20k":           "#5ba1d2",
    "requant_legacy_scaled_vec1d32_FP32tile_gt20k":  "#9ec8e5",
    "requant_lowmidscale_tile_gt20k":                "#0c5c1c",
    "requant_lowmidscale_block128_gt20k":            "#2ca02c",
    "requant_lowmidscale_vec1d32_gt20k":             "#98df8a",
    "requant_lowscale_tile_gt20k":                   "#54278f",
    "requant_lowscale_block128_gt20k":               "#9467bd",
    "requant_lowscale_vec1d32_gt20k":                "#c5b0d5",
    "requant_lowscale_e2m1_tile_gt20k":              "#7c3a09",
    "requant_lowscale_e2m1_block128_gt20k":          "#c46406",
    "requant_lowscale_e2m1_vec1d32_gt20k":           "#fdc086",
    "requant_lowmidscale_e2m1_tile_gt20k":           "#4a1a40",
    "requant_lowmidscale_e2m1_block128_gt20k":       "#a83291",
    "requant_lowmidscale_e2m1_vec1d32_gt20k":        "#df7fc7",
    "requant_ladder_scaled_tile_gt20k":              "#7f0707",
    "requant_ladder_scaled_block128_gt20k":          "#d62728",
    "requant_ladder_scaled_vec1d32_gt20k":           "#ff9896",
}

# View definitions (which sweeps go into which figure). Reference first in each.
# Note: requant_legacy_scaled_vec1d32_FP32tile_gt20k (historical) is intentionally
# excluded — that data was collected before apply_mx_quant_fp64 had a vec1d branch
# and MX_FP32 silently degraded to a tile-wide single shared scale. The
# patched-binary sweep runs under the clean name requant_legacy_scaled_vec1d32_gt20k.
VIEWS = {
    # Focused 2- to 3-sweep subgroups (always baseline first as reference).
    "lowonly": [
        "requant_baseline_fp8_subnormal_gt20k",
        "requant_lowscale_block128_gt20k",
        "requant_lowscale_vec1d32_gt20k",
    ],
    "lowmid": [
        "requant_baseline_fp8_subnormal_gt20k",
        "requant_lowmidscale_block128_gt20k",
        "requant_lowmidscale_vec1d32_gt20k",
    ],
    "fullscaled": [
        "requant_baseline_fp8_subnormal_gt20k",
        "requant_legacy_scaled_block128_gt20k",
        "requant_legacy_scaled_vec1d32_gt20k",
    ],
    "ladder": [
        "requant_baseline_fp8_subnormal_gt20k",
        "requant_ladder_ieee_gt20k",
        "requant_ladder_scaled_block128_gt20k",
        "requant_ladder_scaled_vec1d32_gt20k",
    ],
    # Cross-cuts by scaling geometry.
    "block128_only": [
        "requant_baseline_fp8_subnormal_gt20k",
        "requant_lowscale_block128_gt20k",
        "requant_lowmidscale_block128_gt20k",
        "requant_legacy_scaled_block128_gt20k",
        "requant_ladder_scaled_block128_gt20k",
    ],
    "vec1d32_only": [
        "requant_baseline_fp8_subnormal_gt20k",
        "requant_lowscale_vec1d32_gt20k",
        "requant_lowmidscale_vec1d32_gt20k",
        "requant_legacy_scaled_vec1d32_gt20k",
        "requant_ladder_scaled_vec1d32_gt20k",
    ],
    # Big-picture views.
    "legacy": [
        "requant_baseline_fp8_subnormal_gt20k",
        "requant_legacy_scaled_block128_gt20k",
        "requant_legacy_scaled_vec1d32_gt20k",
        "requant_lowmidscale_block128_gt20k",
        "requant_lowmidscale_vec1d32_gt20k",
        "requant_lowscale_block128_gt20k",
        "requant_lowscale_vec1d32_gt20k",
    ],
    "all": [s for s in SWEEP_ORDER_ALL if s != "requant_legacy_scaled_vec1d32_FP32tile_gt20k"],
}

VIEW_TITLE = {
    "lowonly":      "Baseline vs Low-only Scaled  (Scaled FP8 (E4M3) low tier; FP32/FP16 plain)",
    "lowmid":       "Baseline vs Low+Mid Scaled  (Scaled FP8 (E4M3) + Scaled FP16; FP32 plain)",
    "fullscaled":   "Baseline vs Full Scaled  (Scaled FP32 + Scaled FP16 + Scaled FP8 (E4M3))",
    "ladder":       "Baseline vs Full ladder  (Scaled FP4 (E2M1) → Scaled FP8 (E4M3) → Scaled FP16 → FP32)",
    "block128_only":"All block-128 variants  (2D 128×128 shared scale)",
    "vec1d32_only": "All vec1D-32 variants  (32-elt row-vector shared scale)",
    "legacy":       "All legacy bucketing variants  (drop-in for IEEE FP32/FP16/FP8 cuts)",
    "all":          "All sweeps (legacy + ladder + baseline)",
}

# eps order on the x-axis: loosest first.
EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]


def eps_norm(e):
    """Normalise eps string for matching, e.g. '1e-08' -> '1e-8'."""
    try:
        return f"{float(e):g}".replace("0e-0", "0e-").replace("e-0", "e-")
    except Exception:
        return (e or "").strip()


def bucket_value(r, keys):
    v = 0.0
    for k in keys:
        if k == "scale_meta":
            v += float(r.get("scale_meta_gb", 0) or 0)
        else:
            v += float(r.get(f"{k}_gb", 0) or 0)
    return v


def select_rows(rows, n_target):
    out = []
    for r in rows:
        try:
            if int(r["n"]) == n_target:
                out.append(r)
        except Exception:
            continue
    return out


def find_row(rows, sweep, eps):
    eps_n = eps_norm(eps)
    for r in rows:
        if r["sweep"] == sweep and eps_norm(r["source_epsilon"]) == eps_n:
            return r
    return None


HATCHES = ["", "///", "xxx", "...", "\\\\", "|||", "+++", "OO", "**"]


def make_memory_view(rows, sweeps, out_path, title, n_target):
    eps_list = EPS_ORDER
    n_sw = len(sweeps)
    n_eps = len(eps_list)

    fig_w = max(12, 1.7 * n_sw + 5)
    fig, ax = plt.subplots(figsize=(fig_w, 7.4))
    group_width = 0.86
    bar_w = group_width / max(1, n_sw)
    x_centres = list(range(n_eps))

    seen_buckets = []
    sweep_totals = {}  # (si, ei) -> total

    for si, sw in enumerate(sweeps):
        xs = [c - group_width / 2 + (si + 0.5) * bar_w for c in x_centres]
        bottoms = [0.0] * n_eps
        hatch = HATCHES[si % len(HATCHES)]
        for label, keys in DTYPE_BUCKETS:
            vals = []
            for eps in eps_list:
                r = find_row(rows, sw, eps)
                vals.append(bucket_value(r, keys) if r else 0.0)
            if not any(v > 0 for v in vals):
                continue
            if label not in seen_buckets:
                seen_buckets.append(label)
            ax.bar(xs, vals, bar_w * 0.94, bottom=bottoms,
                   color=DTYPE_COLOR[label], edgecolor="black", linewidth=0.3,
                   hatch=hatch)
            bottoms = [b + v for b, v in zip(bottoms, vals)]

        for ei, (xi, tot) in enumerate(zip(xs, bottoms)):
            sweep_totals[(si, ei)] = tot
            if tot > 0:
                ax.text(xi, tot * 1.01, f"{tot:.2f}",
                        ha="center", va="bottom", fontsize=6.5, rotation=0)

    # Headroom so the GB labels and the legends don't overlap.
    if sweep_totals:
        y_max = max(sweep_totals.values())
        ax.set_ylim(0, y_max * 1.25)

    ax.set_xticks(x_centres)
    ax.set_xticklabels([f"ε={e}" for e in eps_list], fontsize=10)
    ax.set_xlabel("Source epsilon  (loosest → tightest)", fontsize=10)
    ax.set_ylabel(f"Memory footprint of N={n_target} Cholesky factor (GB)", fontsize=10)
    ax.grid(axis="y", alpha=0.25)

    # Datatype legend (colours).
    dtype_handles = [Patch(facecolor=DTYPE_COLOR[b], edgecolor="black", label=b)
                     for b in seen_buckets]
    lg1 = ax.legend(handles=dtype_handles, loc="upper right",
                    fontsize=8, title="Datatype", title_fontsize=9,
                    framealpha=0.95)
    ax.add_artist(lg1)

    # Sweep legend (hatch + colour swatch).
    sweep_handles = [
        Patch(facecolor="white", edgecolor="black",
              hatch=HATCHES[si % len(HATCHES)],
              label=SWEEP_SHORT.get(sw, sw))
        for si, sw in enumerate(sweeps)
    ]
    ax.legend(handles=sweep_handles, loc="upper left",
              fontsize=8, title="Sweep (left → right within each ε cluster)",
              title_fontsize=9, framealpha=0.95)

    fig.suptitle(f"Memory breakdown by datatype — {title}\nN={n_target}",
                 y=0.99, fontsize=11)
    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(out_path)


def make_error_view(rows, sweeps, out_path, title, n_target):
    eps_list = EPS_ORDER
    n_sw = len(sweeps)
    n_eps = len(eps_list)

    fig_w = max(12, 1.7 * n_sw + 5)
    fig, ax = plt.subplots(figsize=(fig_w, 7.4))
    group_width = 0.86
    bar_w = group_width / max(1, n_sw)
    x_centres = list(range(n_eps))

    for si, sw in enumerate(sweeps):
        xs = [c - group_width / 2 + (si + 0.5) * bar_w for c in x_centres]
        ys = []
        for eps in eps_list:
            r = find_row(rows, sw, eps)
            try:
                ys.append(float(r["rel_factor_error"]) if r else 0.0)
            except Exception:
                ys.append(0.0)
        ax.bar(xs, ys, bar_w * 0.94,
               color=SWEEP_COLOR.get(sw, "#888"),
               edgecolor="black", linewidth=0.3,
               hatch=HATCHES[si % len(HATCHES)],
               label=SWEEP_SHORT.get(sw, sw))
        for xi, yv in zip(xs, ys):
            if yv > 0:
                ax.text(xi, yv * 1.12, f"{yv:.1e}",
                        ha="center", va="bottom", fontsize=6.5, rotation=0)

    # Dashed grey reference line per epsilon column.
    eps_floats = [float(e) for e in eps_list]
    for xc, ev in zip(x_centres, eps_floats):
        ax.hlines(ev, xc - group_width / 2, xc + group_width / 2,
                  colors="gray", linestyles="--", linewidth=0.8, alpha=0.7)

    ax.set_yscale("log")
    ax.set_xticks(x_centres)
    ax.set_xticklabels([f"ε={e}" for e in eps_list], fontsize=10)
    ax.set_xlabel("Source epsilon  (loosest → tightest)", fontsize=10)
    ax.set_ylabel(r"$\|LL^\top - A\|_F \,/\, \|A\|_F$  (log scale)", fontsize=10)
    ax.grid(axis="y", which="both", alpha=0.25)
    ax.legend(fontsize=8, loc="upper right", title="Sweep",
              title_fontsize=9, framealpha=0.95)

    fig.suptitle(f"Relative factorization error — {title}\nN={n_target}  "
                 "(dashed grey = source-epsilon target per column)",
                 y=0.99, fontsize=11)
    fig.tight_layout(rect=[0, 0.02, 1, 0.92])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(out_path)


def make_tradeoff_view(rows, sweeps, out_path, title, n_target):
    """Tradeoff = memory savings vs achieved error.
       X = total GB (linear), Y = rel_factor_error (log). One marker per eps.
       Lines connect a sweep's 4 eps points so the user can read the curve.
    """
    eps_list = EPS_ORDER
    fig, ax = plt.subplots(figsize=(10, 7))
    markers = ["o", "s", "D", "^", "v", "P", "X", "*", "<", ">"]

    for si, sw in enumerate(sweeps):
        xs, ys, lbls = [], [], []
        for eps in eps_list:
            r = find_row(rows, sw, eps)
            if not r:
                continue
            try:
                gb = float(r["total_gb"])
                er = float(r["rel_factor_error"])
            except Exception:
                continue
            if gb <= 0 or er <= 0:
                continue
            xs.append(gb)
            ys.append(er)
            lbls.append(eps)
        if not xs:
            continue
        ax.plot(xs, ys, "-", marker=markers[si % len(markers)],
                color=SWEEP_COLOR.get(sw, "#888"), markersize=8,
                label=SWEEP_SHORT.get(sw, sw), alpha=0.9)
        for x, y, lb in zip(xs, ys, lbls):
            ax.annotate(lb, (x, y), fontsize=6.5, xytext=(4, 3),
                        textcoords="offset points", color=SWEEP_COLOR.get(sw, "#444"))

    ax.set_yscale("log")
    ax.set_xlabel("Total memory footprint (GB)", fontsize=10)
    ax.set_ylabel(r"$\|LL^\top - A\|_F\ /\ \|A\|_F$  (log scale)", fontsize=10)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(fontsize=8, loc="upper right")
    fig.suptitle(f"Accuracy ↔ Memory tradeoff — {title} — N={n_target}\n"
                 "(lower-left = better; labels are source-epsilon values per point)",
                 y=0.99, fontsize=11)
    fig.tight_layout(rect=[0, 0.02, 1, 0.93])
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
    print(out_path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--n", type=int, default=32768)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    with Path(args.csv).open() as f:
        rows = list(csv.DictReader(f))
    rows_n = select_rows(rows, args.n)
    if not rows_n:
        raise SystemExit(f"No rows with n={args.n} in {args.csv}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    for view, sweeps in VIEWS.items():
        # Only include sweeps that actually have rows for this N.
        sweeps_have = [s for s in sweeps if any(r["sweep"] == s for r in rows_n)]
        title = VIEW_TITLE[view]
        make_memory_view(rows_n, sweeps_have, out_dir / f"memory_{view}_{args.n}.png", title, args.n)
        make_error_view (rows_n, sweeps_have, out_dir / f"error_{view}_{args.n}.png",  title, args.n)
        make_tradeoff_view(rows_n, sweeps_have, out_dir / f"tradeoff_{view}_{args.n}.png", title, args.n)


if __name__ == "__main__":
    main()
