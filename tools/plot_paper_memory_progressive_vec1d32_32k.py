#!/usr/bin/env python3
"""Paper-grade figure: progressive MX-scaling memory story at N=32768, MX-native.

Single panel, 4 ε clusters, 6 stacked bars per cluster, decomposed by dtype
(GB per Cholesky factor). Same bar layout/order as the merged error figure
(`plot_paper_error_progressive_vec1d32_32k.py`): Baseline IEEE → +MXFP8 →
+MXFP8+MXFP16 → +MXFP8+MXFP16+MXFP32 → Ladder IEEE → Full ladder.

All MX-scaled bars use vec1D-32 (MX-native) storage.
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
    "axes.titlesize": 11.5,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9.5,
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

# Bars: cumulative MX-tier progression, ladder pairs at the right.
# (label, sweep, hatch)
BAR_ORDER = [
    ("Baseline IEEE  (no MX; static format)",
     "requant_baseline_fp8_subnormal_gt20k",      ".."),
    ("+ MXFP8 scaled  (MX_E4M3 low; FP32/FP16 plain)",
     "requant_lowscale_vec1d32_gt20k",            ""),
    ("+ MXFP8 + MXFP16 scaled  (MX_FP16 mid added)",
     "requant_lowmidscale_vec1d32_gt20k",         "//"),
    ("+ MXFP8 + MXFP16 + MXFP32 scaled  (all upper tiers MX)",
     "requant_legacy_scaled_vec1d32_gt20k",       "\\\\"),
    ("Ladder IEEE  (FP8 E4M3 → FP16 → FP32 → FP64)",
     "requant_ladder_ieee_gt20k",                 "||"),
    ("Full ladder  (MXFP4 → MXFP8 E4M3 → MXFP16 → FP32 → FP64)",
     "requant_ladder_scaled_vec1d32_gt20k",       "xx"),
    ("Full ladder FTZ  (subnormals → 0)",
     "requant_ladder_scaled_vec1d32_FTZ_gt20k",   "++"),
]

# Columns that scale linearly with off-diagonal tile count (no diagonal tiles
# for these formats — Cholesky's diagonal is always FP64). When converting
# full→lower-triangular, divide them by 2.
HALVE_COLS = [
    "fp32_gb", "fp16_gb", "bf16_gb",
    "mx_fp32_gb", "mx_fp16_gb", "mx_e4m3_gb", "mx_e5m2_gb",
    "fp8_e4m3_gb", "fp8_e5m2_gb",
    "e3m2_gb", "e2m3_gb", "e2m1_gb",
    "scale_meta_gb",
]

# Datatypes stacked bottom→top in the bar; legend shows top→bottom.
# Use OCP MX naming (MXFP*) in labels.
DTYPE_BUCKETS = [
    ("FP64",         ["fp64_gb"],        "#0072B2"),
    ("FP32",         ["fp32_gb"],        "#E69F00"),
    ("MXFP32",       ["mx_fp32_gb"],     "#9B6814"),
    ("FP16",         ["fp16_gb"],        "#D55E00"),
    ("MXFP16",       ["mx_fp16_gb"],     "#7A2C00"),
    ("MXFP8 (E4M3)", ["mx_e4m3_gb"],     "#009E73"),
    ("FP8 (E4M3)",   ["fp8_e4m3_gb"],    "#56B4E9"),
    ("MXFP4 (E2M1)", ["e2m1_gb"],        "#F0E442"),
    ("Scale meta",   ["scale_meta_gb"],  "#999999"),
]


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


def parse_tile_counts(s):
    out = {}
    for part in (s or "").split(","):
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        try:
            out[k.strip()] = int(v.strip())
        except ValueError:
            pass
    return out


def to_lower_triangular(r, n, nb):
    """Rewrite the per-format GB columns to lower-triangular storage.

    - Off-diagonal tiles are mirrored in a symmetric matrix → halve them.
    - Diagonal tiles are FP64 (POTRF), and the lower half of each diagonal
      tile (nb*(nb+1)/2 elements) is enough — drop the strictly-upper part.
    """
    if r is None:
        return None
    M = n // nb
    counts = parse_tile_counts(r.get("tile_counts_full", ""))
    diag_count = M  # one diagonal tile per block row

    new = dict(r)

    # FP64 split: diagonal (stored as lower-tri within the tile) + off-diag.
    full_fp64_tiles = counts.get("fp64", 0)
    diag_used = min(diag_count, full_fp64_tiles)
    off_fp64_tiles = (full_fp64_tiles - diag_used) // 2
    fp64_bytes = (off_fp64_tiles * nb * nb * 8
                  + diag_used * (nb * (nb + 1) // 2) * 8)
    new["fp64_gb"] = f"{fp64_bytes / 1e9:.9f}"

    # All other formats: off-diagonal-only, so they simply halve.
    for col in HALVE_COLS:
        try:
            v = float(r.get(col) or 0)
        except ValueError:
            v = 0.0
        new[col] = f"{v / 2:.9f}"

    return new


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv",
                    default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--n", type=int, default=32768)
    ap.add_argument("--nb", type=int, default=2048,
                    help="Tile size (only used for --triangular).")
    ap.add_argument("--out", required=True)
    ap.add_argument("--triangular", action="store_true",
                    help="Account for Cholesky's lower-triangular storage: "
                         "off-diagonal tiles counted once, diagonal FP64 tiles "
                         "stored as nb*(nb+1)/2 lower-triangular halves.")
    args = ap.parse_args()

    data = load(args.csv, args.n)

    if args.triangular:
        if args.n % args.nb != 0:
            ap.error(f"--triangular requires n % nb == 0 (n={args.n}, nb={args.nb})")
        data = {k: to_lower_triangular(r, args.n, args.nb) for k, r in data.items()}

    fig, ax = plt.subplots(figsize=(14.5, 6.0))

    n_eps = len(EPS_ORDER)
    n_bars = len(BAR_ORDER)
    group_width = 0.86
    bar_w = group_width / n_bars
    x_centres = list(range(n_eps))

    drawn = {lbl: False for lbl, _, _ in DTYPE_BUCKETS}
    max_y = 0.0

    for bi, bar in enumerate(BAR_ORDER):
        sweep, hatch = bar[1], bar[2]
        xs = [c - group_width / 2 + (bi + 0.5) * bar_w for c in x_centres]
        bottoms = [0.0] * n_eps
        for dt_lbl, keys, color in DTYPE_BUCKETS:
            vals = []
            for eps in EPS_ORDER:
                r = data.get((sweep, eps))
                vals.append(bucket_value(r, keys))
            if not any(v > 0 for v in vals):
                continue
            drawn[dt_lbl] = True
            ax.bar(xs, vals, bar_w * 0.90, bottom=bottoms,
                   color=color, edgecolor="black", linewidth=0.3,
                   hatch=hatch)
            bottoms = [b + v for b, v in zip(bottoms, vals)]
        # Annotate total GB on top of each bar (rotated for fit).
        for x, total in zip(xs, bottoms):
            if total > 0:
                ax.text(x, total * 1.01, f"{total:.2f}",
                        ha="center", va="bottom",
                        fontsize=7.0, rotation=90, color="#222")
        max_y = max(max_y, max(bottoms))

    ax.set_xticks(x_centres)
    ax.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                        for e in EPS_ORDER])
    ax.set_xlim(-0.5, n_eps - 0.5)
    # Extra headroom for the rotated value labels above each bar.
    ax.set_ylim(0, max_y * 1.22)
    if args.triangular:
        ax.set_ylabel("Memory per Cholesky L factor — lower triangle (GB)")
    else:
        ax.set_ylabel("Memory per Cholesky factor (GB)")
    ax.set_xlabel("")
    tri_note = "lower triangle" if args.triangular else "full matrix"
    ax.set_title(f"Progressive MX scaling: memory at N={args.n}, nb={args.nb} "
                 f"({tri_note}, MX-native storage, 1×32 row groups)")
    ax.set_axisbelow(True)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, alpha=0.3, linewidth=0.5)

    # Two-row legend (datatypes + bar order).
    dtype_handles = [Patch(facecolor=color, edgecolor="black", label=lbl)
                     for lbl, _, color in DTYPE_BUCKETS if drawn[lbl]]
    bar_handles = [Patch(facecolor="white", edgecolor="black", hatch=hatch, label=label)
                   for label, _, hatch in BAR_ORDER]

    lg_dtype = fig.legend(handles=dtype_handles,
                          loc="lower center",
                          bbox_to_anchor=(0.5, 0.02),
                          ncol=min(len(dtype_handles), 5),
                          frameon=True, framealpha=0.95, title="Datatype",
                          handletextpad=0.4, columnspacing=1.2, labelspacing=0.3,
                          borderpad=0.35)
    fig.add_artist(lg_dtype)

    fig.legend(handles=bar_handles,
               loc="lower center", bbox_to_anchor=(0.5, -0.14),
               ncol=2, frameon=True, framealpha=0.95,
               title="Bar order within each ε cluster (left → right)",
               handletextpad=0.6, columnspacing=2.0, labelspacing=0.35,
               borderpad=0.4)

    # Reserve room at the bottom for the two-row legend.
    fig.subplots_adjust(bottom=0.30)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pdf", ".png"):
        p = out.with_suffix(ext)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(p)
    plt.close(fig)


if __name__ == "__main__":
    main()
