#!/usr/bin/env python3
"""Paper-grade figure for the IEEE-only FP8 baseline at N=32768.

Two panels:
  (top)    Memory-accuracy curve in (total_GB, rel_factor_error) space,
           parameterised by source_epsilon. Markers labelled with ε.
  (bottom) Bucket-distribution stacked bar chart vs ε (loosest left).

Datatype legend uses OCP-spec names; FP8 baseline uses E4M3 low tier.
"""
import argparse
import csv
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# --- Paper-grade style ---
mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9,
    "legend.title_fontsize": 10,
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,   # TrueType so editors can re-flow text
    "ps.fonttype": 42,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
})

# Colourblind-safe (Wong / Okabe-Ito).
COLORS = {
    "FP64":    "#0072B2",
    "FP32":    "#E69F00",
    "FP16":    "#D55E00",
    "FP8 E4M3":"#009E73",
}

EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]


def parse_counts(s):
    out = {}
    for part in (s or "").split(","):
        if "=" not in part: continue
        k, v = part.split("=", 1)
        k = k.strip().lower()
        try:
            out[k] = int(v)
        except ValueError:
            pass
    return out


N_FROM_NAME = {"32k": 32768, "40k": 40960, "65k": 65536, "80k": 81920,
               "100k": 98304, "120k": 122880, "16384": 16384, "8192": 8192,
               "2048": 2048, "1024": 1024}


def n_for_bin(bin_path):
    base = Path(bin_path).name
    n_str = base.rsplit("_", 1)[-1].rsplit(".", 1)[0]
    return N_FROM_NAME.get(n_str, 0)


def load_baseline(merged_csv, sweep_name, n_target):
    """Load baseline rows from the merged memory CSV (has total_gb pre-computed)."""
    rows = {}
    with Path(merged_csv).open() as f:
        for r in csv.DictReader(f):
            if r["sweep"] != sweep_name:
                continue
            try:
                if int(r["n"]) != n_target:
                    continue
            except Exception:
                continue
            rows[r["source_epsilon"]] = r
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv",
                    help="Merged memory CSV (built by build_requant_gt20k_memory_csv.py)")
    ap.add_argument("--sweep", default="requant_baseline_fp8_subnormal_gt20k",
                    help="IEEE baseline sweep name to read")
    ap.add_argument("--n", type=int, default=32768)
    ap.add_argument("--out", required=True, help="Output basename (will produce .pdf and .png)")
    args = ap.parse_args()

    rows = load_baseline(args.csv, args.sweep, args.n)
    if not rows:
        raise SystemExit(f"no baseline rows found at N={args.n} in {args.sweep}")

    # Build per-eps tuples
    series = []
    for eps in EPS_ORDER:
        r = rows.get(eps)
        if not r: continue
        gb = float(r["total_gb"])
        err = float(r["rel_factor_error"])
        counts = parse_counts(r["tile_counts_full"])
        series.append({"eps": eps, "gb": gb, "err": err, "counts": counts})

    # Single-panel Pareto: memory on x, error on y (log).
    fig, ax = plt.subplots(figsize=(7.0, 5.0))

    xs = [s["gb"] for s in series]
    ys = [s["err"] for s in series]
    ax.plot(xs, ys, marker="o", markersize=10, linewidth=2.0,
            color="#0072B2", markerfacecolor="#0072B2",
            markeredgecolor="black", markeredgewidth=0.9, zorder=3)
    # Label each point with its source-epsilon value.
    for s in series:
        ax.annotate(rf"$\varepsilon = {s['eps']}$", (s["gb"], s["err"]),
                    xytext=(10, 6), textcoords="offset points",
                    fontsize=10.5, color="#222222",
                    bbox=dict(boxstyle="round,pad=0.25", facecolor="white",
                              edgecolor="#888888", alpha=0.94, linewidth=0.5))
    ax.set_yscale("log")
    ax.set_xlabel("Memory footprint of Cholesky factor  (GB)")
    ax.set_ylabel(r"Relative factorisation error  $\|LL^\top - A\|_F\,/\,\|A\|_F$")
    # X-axis: a touch of headroom so labels don't clip
    if xs:
        xlo = min(xs) * 0.95
        xhi = max(xs) * 1.05
        ax.set_xlim(xlo, xhi)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pdf", ".png"):
        p = out.with_suffix(ext)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(p)
    plt.close(fig)


if __name__ == "__main__":
    main()
