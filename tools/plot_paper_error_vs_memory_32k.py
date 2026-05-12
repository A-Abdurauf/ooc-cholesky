#!/usr/bin/env python3
"""Paper-grade figure: Relative factorisation error vs memory footprint at N=32768.

Pareto-style accuracy-vs-memory scatter. One marker per (sweep, epsilon) for
14 sweeps x 4 epsilons. Marker shape encodes the sweep family; marker colour
encodes the granularity (or the two IEEE references). The 4 epsilon points
inside a sweep are connected by a line to show the per-sweep trajectory.
"""
import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
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


EPS_ORDER = ["1e-8", "1e-7", "1e-6", "1e-5"]

# Sweep -> (family_key, granularity_key, short_label).
# family_key:   baseline | ladder_ieee | legacy | lowmid | low | ladder
# granularity:  tile | block128 | vec1d32 | ieee_ref
SWEEPS = [
    # IEEE-only references
    ("requant_baseline_fp8_subnormal_gt20k", "baseline",    "ieee_ref", "baseline IEEE (subn)"),
    ("requant_ladder_ieee_gt20k",            "ladder_ieee", "ieee_ref", "ladder IEEE (no scaling)"),
    # legacy 3-tier
    ("requant_legacy_scaled_tile_gt20k",      "legacy", "tile",     "legacy 3-tier - tile"),
    ("requant_legacy_scaled_block128_gt20k",  "legacy", "block128", "legacy 3-tier - block 128"),
    ("requant_legacy_scaled_vec1d32_gt20k",   "legacy", "vec1d32",  "legacy 3-tier - vec1D 32"),
    # low+mid
    ("requant_lowmidscale_tile_gt20k",        "lowmid", "tile",     "low+mid - tile"),
    ("requant_lowmidscale_block128_gt20k",    "lowmid", "block128", "low+mid - block 128"),
    ("requant_lowmidscale_vec1d32_gt20k",     "lowmid", "vec1d32",  "low+mid - vec1D 32"),
    # low only
    ("requant_lowscale_tile_gt20k",           "low",    "tile",     "low only - tile"),
    ("requant_lowscale_block128_gt20k",       "low",    "block128", "low only - block 128"),
    ("requant_lowscale_vec1d32_gt20k",        "low",    "vec1d32",  "low only - vec1D 32"),
    # ladder scaled
    ("requant_ladder_scaled_tile_gt20k",      "ladder", "tile",     "ladder - tile"),
    ("requant_ladder_scaled_block128_gt20k",  "ladder", "block128", "ladder - block 128"),
    ("requant_ladder_scaled_vec1d32_gt20k",   "ladder", "vec1d32",  "ladder - vec1D 32"),
]

# Marker shape per family (6 distinct shapes; baseline + ladder_ieee are the 2 refs).
FAMILY_MARKER = {
    "baseline":    "X",   # reference
    "ladder_ieee": "P",   # reference
    "legacy":      "s",
    "lowmid":      "D",
    "low":         "^",
    "ladder":      "o",
}
FAMILY_LABEL = {
    "baseline":    "baseline IEEE (subn)",
    "ladder_ieee": "ladder IEEE (no scaling)",
    "legacy":      "legacy 3-tier (SFP32+SFP16+SFP8)",
    "lowmid":      "low+mid scaled (SFP8+SFP16)",
    "low":         "low only scaled (SFP8)",
    "ladder":      "ladder (SFP4-SFP8-SFP16-FP32)",
}

# Granularity palette (Wong / Okabe-Ito-inspired). Refs get neutral colours.
GRAN_COLOR = {
    "tile":     "#0072B2",   # blue
    "block128": "#D55E00",   # vermillion
    "vec1d32":  "#009E73",   # bluish green
    "ieee_ref": "#555555",   # dark grey (refs)
}
GRAN_LABEL = {
    "tile":     "tile",
    "block128": "block 128",
    "vec1d32":  "vec1D 32",
}


def load(csv_path, n_target):
    rows = {}
    with Path(csv_path).open() as f:
        for r in csv.DictReader(f):
            try:
                if int(r["n"]) != n_target:
                    continue
            except Exception:
                continue
            rows[(r["sweep"], r["source_epsilon"])] = r
    return rows


def num(r, key):
    if r is None:
        return None
    try:
        return float(r[key])
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", default="/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv")
    ap.add_argument("--n", type=int, default=32768)
    ap.add_argument("--out", default="/home/abduraa/MX_project/logs/mx_ooc_data/plots/paper/fig_error_vs_memory_32k")
    args = ap.parse_args()

    data = load(args.csv, args.n)

    fig, ax = plt.subplots(figsize=(10, 6.0))

    pareto_pts = []  # (sweep, eps, total_gb, rel_err)
    for sweep, fam, gran, label in SWEEPS:
        xs, ys, eps_used = [], [], []
        for eps in EPS_ORDER:
            r = data.get((sweep, eps))
            gb = num(r, "total_gb")
            er = num(r, "rel_factor_error")
            if gb is None or er is None or er <= 0:
                continue
            xs.append(gb)
            ys.append(er)
            eps_used.append(eps)
            pareto_pts.append((sweep, eps, gb, er, fam, gran))

        if not xs:
            continue

        color = GRAN_COLOR[gran]
        marker = FAMILY_MARKER[fam]
        is_ref = gran == "ieee_ref"

        # Connect the eps trajectory.
        ax.plot(xs, ys,
                linestyle=(":" if is_ref else "-"),
                linewidth=(1.4 if is_ref else 1.0),
                color=color, alpha=(0.85 if is_ref else 0.55),
                zorder=2)

        # Points.
        ax.scatter(xs, ys,
                   marker=marker,
                   s=(95 if is_ref else 70),
                   facecolor=("white" if is_ref else color),
                   edgecolor=color,
                   linewidth=(1.6 if is_ref else 1.1),
                   alpha=0.95,
                   zorder=5)

        # Annotate the loosest (1e-5) and tightest (1e-8) epsilon for orientation.
        # Only do this for one sweep per family to avoid clutter; pick tile when
        # available, else the first available.
        if gran in ("tile", "ieee_ref"):
            if "1e-8" in eps_used:
                i = eps_used.index("1e-8")
                ax.annotate(r"$\varepsilon{=}10^{-8}$",
                            (xs[i], ys[i]),
                            xytext=(6, -2),
                            textcoords="offset points",
                            fontsize=7.5, color=color, alpha=0.9)
            if "1e-5" in eps_used:
                i = eps_used.index("1e-5")
                ax.annotate(r"$\varepsilon{=}10^{-5}$",
                            (xs[i], ys[i]),
                            xytext=(6, 2),
                            textcoords="offset points",
                            fontsize=7.5, color=color, alpha=0.9)

    ax.set_yscale("log")
    ax.set_xlabel("Memory footprint per Cholesky factor [GB]")
    ax.set_ylabel("Relative factorisation error (log scale)")

    ax.set_axisbelow(True)

    # Dense FP64 reference (vertical line).
    dense_fp64_gb = (args.n * args.n * 8) / 1e9
    # Don't actually draw FP64 line (off scale at 8.59 GB might be near right
    # edge depending on data range), but include in suptitle.

    # Legends: two-column composite legend (family shapes, granularity colours).
    family_keys_in_order = ["baseline", "ladder_ieee", "legacy", "lowmid", "low", "ladder"]
    family_handles = [
        Line2D([0], [0],
               marker=FAMILY_MARKER[k],
               color="white",
               markerfacecolor=("white" if k in ("baseline", "ladder_ieee") else "#555555"),
               markeredgecolor="#222222",
               markeredgewidth=1.2,
               markersize=(10 if k in ("baseline", "ladder_ieee") else 8),
               label=FAMILY_LABEL[k])
        for k in family_keys_in_order
    ]
    gran_handles = [
        Line2D([0], [0], marker="o", color="white",
               markerfacecolor=GRAN_COLOR["tile"], markeredgecolor=GRAN_COLOR["tile"],
               markersize=9, label="tile"),
        Line2D([0], [0], marker="o", color="white",
               markerfacecolor=GRAN_COLOR["block128"], markeredgecolor=GRAN_COLOR["block128"],
               markersize=9, label="block 128"),
        Line2D([0], [0], marker="o", color="white",
               markerfacecolor=GRAN_COLOR["vec1d32"], markeredgecolor=GRAN_COLOR["vec1d32"],
               markersize=9, label="vec1D 32"),
        Line2D([0], [0], marker="o", color="white",
               markerfacecolor="white", markeredgecolor=GRAN_COLOR["ieee_ref"],
               markeredgewidth=1.5, markersize=9, label="IEEE reference"),
    ]

    leg1 = ax.legend(handles=family_handles,
                     title="Sweep family (marker shape)",
                     loc="upper left",
                     bbox_to_anchor=(1.01, 1.0),
                     frameon=True, framealpha=0.95,
                     handletextpad=0.5, labelspacing=0.35,
                     borderpad=0.4)
    ax.add_artist(leg1)
    ax.legend(handles=gran_handles,
              title="Granularity (marker colour)",
              loc="lower left",
              bbox_to_anchor=(1.01, 0.0),
              frameon=True, framealpha=0.95,
              handletextpad=0.5, labelspacing=0.35,
              borderpad=0.4)

    fig.tight_layout(rect=[0, 0, 0.78, 0.99])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".pdf", ".png"):
        p = out.with_suffix(ext)
        fig.savefig(p, dpi=300, bbox_inches="tight")
        print(p)
    plt.close(fig)

    # Print sanity summary: lowest-error point per family at N.
    print("\n# Lowest rel_factor_error per family at N =", args.n)
    by_fam = {}
    for sweep, eps, gb, er, fam, gran in pareto_pts:
        cur = by_fam.get(fam)
        if cur is None or er < cur[3]:
            by_fam[fam] = (sweep, eps, gb, er, gran)
    print(f"{'family':<14s} {'gran':<10s} {'eps':<6s} {'GB':>8s} {'rel_err':>12s}")
    for fam in ["baseline", "ladder_ieee", "legacy", "lowmid", "low", "ladder"]:
        if fam not in by_fam:
            continue
        sweep, eps, gb, er, gran = by_fam[fam]
        print(f"{fam:<14s} {gran:<10s} {eps:<6s} {gb:>8.3f} {er:>12.3e}")


if __name__ == "__main__":
    main()
