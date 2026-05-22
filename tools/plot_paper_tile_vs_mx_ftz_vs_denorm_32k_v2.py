#!/usr/bin/env python3
"""Three-panel figure (v2): adds Drop-in MXFP4 + uses the post-archive
ladder_rerun_32k data for MX·FTZ / MX·denorm.

Bars:
  Baseline IEEE       (main CSV, requant_baseline_fp8_subnormal_gt20k)
  Tile · denorm       (main CSV, requant_ladder_scaled_tile_gt20k)
  Tile · FTZ          (main CSV, requant_ladder_scaled_tile_FTZ_gt20k)
  MX · denorm         (ladder_rerun_32k/results.csv, ladder_full_gu)
  MX · FTZ            (ladder_rerun_32k/results.csv, ladder_full_fz)
  Drop-in MXFP4       (true_dropin_mxfp4_32k/results.csv, true_mxfp4_dropin_fz)
                      (gu row is bit-identical, so we plot one bar)

Panels: error · memory · tile allocation.
"""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_paper_memory_progressive_vec1d32_32k import (  # noqa: E402
    to_lower_triangular,
)

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

# (label, source, sweep_name, color, hatch)
#   source ∈ {"main", "rerun", "dropin"}
BARS = [
    ("Baseline IEEE",  "main",   "requant_baseline_fp8_subnormal_gt20k", "#555555", ".."),
    ("Tile · denorm",  "main",   "requant_ladder_scaled_tile_gt20k",     "#1F77B4", ""),
    ("Tile · FTZ",     "main",   "requant_ladder_scaled_tile_FTZ_gt20k", "#1F77B4", "//"),
    ("MX · denorm",    "rerun",  "ladder_full_gu",                       "#D62728", ""),
    ("MX · FTZ",       "rerun",  "ladder_full_fz",                       "#D62728", "//"),
    ("Drop-in MXFP4",  "dropin", "true_mxfp4_dropin_fz",                 "#F0E442", "xx"),
]

DTYPE_BUCKETS = [
    ("FP64",         ["fp64_gb"],                        "#0072B2"),
    ("FP32",         ["fp32_gb"],                        "#E69F00"),
    ("FP16/MXFP16",  ["fp16_gb", "mx_fp16_gb"],          "#7A2C00"),
    ("FP8 plain",    ["fp8_e4m3_gb", "fp8_e5m2_gb"],     "#56B4E9"),
    ("MXFP8 (E4M3)", ["mx_e4m3_gb"],                     "#009E73"),
    ("MXFP4 (E2M1)", ["e2m1_gb"],                        "#F0E442"),
    ("Scale meta",   ["scale_meta_gb"],                  "#999999"),
]

TILE_FORMATS = [
    ("FP64",         "fp64",     "#0072B2"),
    ("FP32",         "fp32",     "#E69F00"),
    ("FP16/MXFP16",  "mx_fp16",  "#7A2C00"),
    ("FP8 plain",    "fp8_e4m3", "#56B4E9"),
    ("MXFP8 (E4M3)", "mx_e4m3",  "#009E73"),
    ("MXFP4 (E2M1)", "e2m1",     "#F0E442"),
]

TILE_ALIASES = {
    "fp16":     "mx_fp16",
    "mx_fp32":  "fp32",
    "fp8_e5m2": "fp8_e4m3",
}

# Bytes per element + per-32-element MX scale byte. Storage formula assumes
# vec1D-32 grouping for shared-scale formats (matches all sweeps below).
BYTES_PER_ELEM = {
    "fp64": 8.0,
    "fp32": 4.0,
    "mx_fp16": 2.0,
    "fp16": 2.0,
    "mx_e4m3": 1.0,
    "fp8_e4m3": 1.0,
    "e2m1": 0.5,
}
USES_SCALE = {
    "mx_fp16": True,
    "mx_e4m3": True,
    "e2m1": True,
}


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
    # "fp32=50;mx_fp16=10;e2m1=60;fp64=16;"
    out = {}
    for part in (s or "").strip().strip('"').split(";"):
        part = part.strip()
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


def synth_row(rerun_row, n, nb):
    """Turn a rerun/dropin CSV row into a main-CSV-like dict.

    Computes per-format _gb columns from the tile_breakdown using bytes-per-
    element + the vec1D-32 scale overhead. Returns a dict matching the
    columns used by `to_lower_triangular()` and `bucket_value()`.
    """
    counts = parse_tile_breakdown(rerun_row["tile_breakdown"])
    tile_elems = nb * nb
    scale_groups_per_tile = tile_elems // 32  # vec1D-32
    gb = {
        "fp64_gb":      0.0,
        "fp32_gb":      0.0,
        "fp16_gb":      0.0,
        "mx_fp16_gb":   0.0,
        "fp8_e4m3_gb":  0.0,
        "fp8_e5m2_gb":  0.0,
        "mx_e4m3_gb":   0.0,
        "e2m1_gb":      0.0,
        "scale_meta_gb": 0.0,
    }
    GB = 1.0 / (1024 ** 3)
    for fmt, cnt in counts.items():
        bpe = BYTES_PER_ELEM.get(fmt, 0.0)
        data_bytes = bpe * tile_elems * cnt
        scale_bytes = scale_groups_per_tile * cnt if USES_SCALE.get(fmt, False) else 0
        # Store into the right bucket.
        if fmt == "fp64":      gb["fp64_gb"]      += data_bytes * GB
        elif fmt == "fp32":    gb["fp32_gb"]      += data_bytes * GB
        elif fmt == "mx_fp16": gb["mx_fp16_gb"]   += data_bytes * GB
        elif fmt == "mx_e4m3": gb["mx_e4m3_gb"]   += data_bytes * GB
        elif fmt == "e2m1":    gb["e2m1_gb"]      += data_bytes * GB
        elif fmt == "fp8_e4m3":gb["fp8_e4m3_gb"]  += data_bytes * GB
        gb["scale_meta_gb"] += scale_bytes * GB

    # Build a row compatible with the rest of the pipeline.
    return {
        "sweep": rerun_row["sweep"],
        "n": rerun_row["n"],
        "nb": rerun_row["nb"],
        "source_epsilon": rerun_row["source_epsilon"],
        "rel_factor_error": rerun_row["rel_factor_error"],
        "tile_counts_full": ",".join(f"{k}={v*2}" for k, v in counts.items()
                                     if k != "fp64") +
                            ("," if counts else "") +
                            (f"fp64={counts.get('fp64', 0) * 2 - (n // nb) if counts.get('fp64', 0) else 0}"
                             if False else  # disabled; use lower-tri counts directly downstream
                             f"fp64={counts.get('fp64', 0)}"),
        **gb,
    }


def load_main(csv_path, n_target, nb_target):
    rows = {}
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            try:
                if int(r["n"]) != n_target or int(r["nb"]) != nb_target:
                    continue
            except (ValueError, KeyError):
                continue
            rows[(r["sweep"], r["source_epsilon"])] = r  # last wins
    return rows


def load_aux(csv_path, n_target, nb_target):
    """Load rerun/dropin CSV and synthesize main-CSV-like rows."""
    rows = {}
    if not Path(csv_path).exists():
        return rows
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            try:
                if int(r["n"]) != n_target or int(r["nb"]) != nb_target:
                    continue
            except (ValueError, KeyError):
                continue
            rows[(r["sweep"], r["source_epsilon"])] = synth_row(r, n_target, nb_target)
    return rows


def bucket_value(r, keys):
    if r is None:
        return 0.0
    return sum(float(r.get(k, 0) or 0) for k in keys)


def tile_counts_lower(counts_full, n, nb):
    if n % nb:
        raise ValueError(f"n={n} not divisible by nb={nb}")
    M = n // nb
    out = {}
    for fmt, v in counts_full.items():
        if fmt == "fp64":
            off = max(v - M, 0)
            out[fmt] = M + off // 2
        else:
            out[fmt] = v // 2
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv",
                    default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--rerun-csv",
                    default="/home/abduraa/MX_project/ooc-cholesky/ladder_rerun_32k/results.csv")
    ap.add_argument("--dropin-csv",
                    default="/home/abduraa/MX_project/ooc-cholesky/true_dropin_mxfp4_32k/results.csv")
    ap.add_argument("--n", type=int, default=32768)
    ap.add_argument("--nb", type=int, default=2048)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    data_main   = load_main(args.csv, args.n, args.nb)
    data_rerun  = load_aux(args.rerun_csv, args.n, args.nb)
    data_dropin = load_aux(args.dropin_csv, args.n, args.nb)

    def pick(source, sweep, eps):
        if source == "main":   return data_main.get((sweep, eps))
        if source == "rerun":  return data_rerun.get((sweep, eps))
        if source == "dropin": return data_dropin.get((sweep, eps))
        return None

    # For tile counts and lower-tri memory: main CSV uses tile_counts_full
    # ("a=A,b=B,...") with FULL-square counts; synthesized rows already store
    # lower-tri counts inside their tile_counts_full column. Differentiate via
    # a per-row 'tile_counts_already_lower' marker we add below.
    for r in list(data_rerun.values()) + list(data_dropin.values()):
        r["tile_counts_already_lower"] = "1"

    def lower_tri_counts_from_row(r):
        full = parse_tile_counts_full(r.get("tile_counts_full", ""))
        if r.get("tile_counts_already_lower") == "1":
            return full  # already lower-tri
        return tile_counts_lower(full, args.n, args.nb)

    # Memory: rows from rerun/dropin already store lower-tri-derived gb values,
    # so DO NOT pass them through to_lower_triangular. For main-CSV rows we
    # still want the lower-tri transform.
    def lower_tri_memory(r):
        if r is None: return None
        if r.get("tile_counts_already_lower") == "1":
            return r
        return to_lower_triangular(r, args.n, args.nb)

    fig, (ax_err, ax_mem, ax_tile) = plt.subplots(
        1, 3, figsize=(20.5, 5.8),
        gridspec_kw={"width_ratios": [1, 1.05, 1.05]},
    )

    n_eps = len(EPS_ORDER)
    n_bars = len(BARS)
    group_width = 0.84
    bar_w = group_width / n_bars
    x_centres = list(range(n_eps))

    # ---- LEFT: relative-error bars.
    max_err = 0.0
    for bi, (lbl, src, sweep, color, hatch) in enumerate(BARS):
        xs = [c - group_width / 2 + (bi + 0.5) * bar_w for c in x_centres]
        vals = []
        for eps in EPS_ORDER:
            r = pick(src, sweep, eps)
            v = float(r["rel_factor_error"]) if r and r.get("rel_factor_error") else np.nan
            vals.append(v)
            if np.isfinite(v):
                max_err = max(max_err, v)
        ax_err.bar(xs, vals, bar_w * 0.90, color=color, edgecolor="black",
                   linewidth=0.4, hatch=hatch, label=lbl)
        for x, v in zip(xs, vals):
            if np.isfinite(v) and v > 0:
                ax_err.text(x, v * 1.15, f"{v:.1e}",
                            ha="center", va="bottom",
                            fontsize=7.0, rotation=90, color="#222")

    ax_err.set_yscale("log")
    ax_err.set_xticks(x_centres)
    ax_err.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                            for e in EPS_ORDER])
    ax_err.set_xlim(-0.5, n_eps - 0.5)
    ax_err.set_ylim(top=max_err * 30)
    ax_err.set_ylabel("Relative error")
    ax_err.set_title(f"Error  (ladder, N={args.n}, nb={args.nb})")
    ax_err.set_axisbelow(True)
    ax_err.xaxis.grid(False)
    ax_err.yaxis.grid(True, which="both", alpha=0.3, linewidth=0.5)
    ax_err.legend(loc="upper right", frameon=True, framealpha=0.95)

    # ---- MIDDLE: stacked memory bars.
    drawn = {lbl: False for lbl, _, _ in DTYPE_BUCKETS}
    max_y = 0.0
    for bi, (bar_lbl, src, sweep, color, hatch) in enumerate(BARS):
        xs = [c - group_width / 2 + (bi + 0.5) * bar_w for c in x_centres]
        bottoms = [0.0] * n_eps
        for dt_lbl, keys, dt_color in DTYPE_BUCKETS:
            vals = [bucket_value(lower_tri_memory(pick(src, sweep, eps)), keys)
                    for eps in EPS_ORDER]
            if not any(v > 0 for v in vals):
                continue
            drawn[dt_lbl] = True
            ax_mem.bar(xs, vals, bar_w * 0.90, bottom=bottoms,
                       color=dt_color, edgecolor="black", linewidth=0.3,
                       hatch=hatch)
            bottoms = [b + v for b, v in zip(bottoms, vals)]
        for x, t in zip(xs, bottoms):
            if t > 0:
                ax_mem.text(x, t * 1.01, f"{t:.2f}",
                            ha="center", va="bottom",
                            fontsize=7.0, rotation=90, color="#222")
        max_y = max(max_y, max(bottoms))

    ax_mem.set_xticks(x_centres)
    ax_mem.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                            for e in EPS_ORDER])
    ax_mem.set_xlim(-0.5, n_eps - 0.5)
    ax_mem.set_ylim(0, max_y * 1.20)
    ax_mem.set_ylabel("Memory per Cholesky L factor — lower triangle (GB)")
    ax_mem.set_title(f"Memory  (ladder, N={args.n}, nb={args.nb}, lower triangle)")
    ax_mem.set_axisbelow(True)
    ax_mem.xaxis.grid(False)
    ax_mem.yaxis.grid(True, alpha=0.3, linewidth=0.5)

    # ---- RIGHT: tile allocation.
    max_t = 0.0
    for bi, (bar_lbl, src, sweep, color, hatch) in enumerate(BARS):
        xs = [c - group_width / 2 + (bi + 0.5) * bar_w for c in x_centres]
        per_eps = []
        for eps in EPS_ORDER:
            r = pick(src, sweep, eps)
            if r is None:
                per_eps.append({})
                continue
            per_eps.append(lower_tri_counts_from_row(r))

        bottoms = [0.0] * n_eps
        for fmt_lbl, fmt_key, fmt_color in TILE_FORMATS:
            vals = [pe.get(fmt_key, 0) for pe in per_eps]
            if not any(v > 0 for v in vals):
                continue
            ax_tile.bar(xs, vals, bar_w * 0.90, bottom=bottoms,
                        color=fmt_color, edgecolor="black", linewidth=0.3,
                        hatch=hatch)
            bottoms = [b + v for b, v in zip(bottoms, vals)]
        for x, t in zip(xs, bottoms):
            if t > 0:
                ax_tile.text(x, t * 1.01, f"{int(round(t))}",
                             ha="center", va="bottom",
                             fontsize=7.0, rotation=90, color="#222")
        max_t = max(max_t, max(bottoms))

    ax_tile.set_xticks(x_centres)
    ax_tile.set_xticklabels([rf"$\varepsilon = 10^{{{int(float(e).__format__('e').split('e')[1])}}}$"
                             for e in EPS_ORDER])
    ax_tile.set_xlim(-0.5, n_eps - 0.5)
    ax_tile.set_ylim(0, max_t * 1.18)
    ax_tile.set_ylabel("Tile count (lower triangle, incl. diagonal)")
    ax_tile.set_title(f"Tile allocation  (ladder, N={args.n}, nb={args.nb}, lower triangle)")
    ax_tile.set_axisbelow(True)
    ax_tile.xaxis.grid(False)
    ax_tile.yaxis.grid(True, alpha=0.3, linewidth=0.5)

    # Compound legend.
    dtype_handles = [Patch(facecolor=c, edgecolor="black", label=l)
                     for l, _, c in DTYPE_BUCKETS if drawn[l]]
    style_handles = [
        Patch(facecolor="white",   edgecolor="black", label="solid = denorm"),
        Patch(facecolor="white",   edgecolor="black", hatch="//", label="hatched = FTZ"),
        Patch(facecolor="#555555", edgecolor="black", hatch="..", label="grey = Baseline IEEE"),
        Patch(facecolor="#1F77B4", edgecolor="black", label="blue = Tile"),
        Patch(facecolor="#D62728", edgecolor="black", label="red = MX (vec1D-32)"),
        Patch(facecolor="#F0E442", edgecolor="black", hatch="xx",
              label="yellow xx = MXFP4 drop-in"),
    ]
    leg1 = fig.legend(handles=style_handles,
                      loc="lower center", bbox_to_anchor=(0.22, -0.08),
                      ncol=2, frameon=True, framealpha=0.95, title="Bar style",
                      handletextpad=0.6, columnspacing=1.4, labelspacing=0.3,
                      borderpad=0.35)
    fig.add_artist(leg1)
    fig.legend(handles=dtype_handles,
               loc="lower center", bbox_to_anchor=(0.66, -0.08),
               ncol=3, frameon=True, framealpha=0.95,
               title="Datatype  /  Tile format",
               handletextpad=0.4, columnspacing=1.2, labelspacing=0.3,
               borderpad=0.35)

    fig.subplots_adjust(bottom=0.22)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pdf", ".png"):
        p = out.with_suffix(ext)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(p)
    plt.close(fig)


if __name__ == "__main__":
    main()
