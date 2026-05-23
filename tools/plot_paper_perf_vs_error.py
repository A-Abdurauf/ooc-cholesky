#!/usr/bin/env python3
"""Per-N perf-vs-error scatter.

3 subplots, one per N (32k, 40k, 64k).  Within each subplot:
  - x-axis: relative factor error (log)
  - y-axis: performance (TFLOPS)
  - 4 colored shapes: one per mode (Baseline / Ladder IEEE / Ladder MX /
    Ladder MX+MXFP16); points within each mode are connected by a line
    across the 4 eps values.
  - vertical dashed grey lines at the 4 source-eps targets (10^-5..10^-8)
    divide the x-axis into eps regions; eps labels along the top edge.

Performance comes from recreate-ooc-chol/build/combined_perf_data.csv.
Error comes from figures/all_error_runs_ladder_only.csv.  Each (mode, N,
eps) cell uses a 'preferred' (granularity, underflow) row in the error CSV
so we have a single error value per perf point.
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

PERF_CSV  = "/home/abduraa/MX_project/recreate-ooc-chol/build/combined_perf_data.csv"
ERROR_CSV = "/home/abduraa/MX_project/ooc-cholesky/figures/all_error_runs_ladder_only.csv"

# Match the perf CSV's canonical mode column.
MODES = [
    ("Baseline",         "Baseline"),
    ("Ladder IEEE",      "Ladder IEEE"),
    ("Ladder MX",        "Ladder MX (native)"),       # error CSV name
    ("Ladder MX+MXFP16", "Ladder MX+MXFP16"),
]
EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]
NS_TO_PLOT = [32768, 40960, 65536]
NB_BY_N    = {32768: 2048, 40960: 4096, 65536: 4096}
N_LABEL    = {32768: "32k", 40960: "40k", 65536: "64k"}

MODE_COLOR = {
    "Baseline":         "#4e79a7",
    "Ladder IEEE":      "#7f7f7f",
    "Ladder MX":        "#e07b39",
    "Ladder MX+MXFP16": "#b07aa1",
}
MODE_MARKER = {
    "Baseline":         "o",
    "Ladder IEEE":      "s",
    "Ladder MX":        "^",
    "Ladder MX+MXFP16": "D",
}
# Distinct dash patterns per mode so crossing lines stay legible.
MODE_LINESTYLE = {
    "Baseline":         "-",        # solid
    "Ladder IEEE":      (0, (5, 2)),    # long dashes
    "Ladder MX":        (0, (1, 1.5)),  # dense dots
    "Ladder MX+MXFP16": (0, (3, 1, 1, 1)),  # dash-dot
}

# Preferred (granularity, underflow) per mode for picking ONE error row per
# (mode, N, eps).  Falls back to (tile, gu) -> (MX, fz) -> (tile, fz) if the
# preferred is missing.
PREFERRED_GU = {
    "Baseline":            ("n/a", "fz"),
    "Ladder IEEE":         ("n/a", "fz"),
    "Ladder MX (native)":  ("MX",  "gu"),
    "Ladder MX+MXFP16":    ("MX",  "gu"),
}
FALLBACKS = [("MX", "gu"), ("tile", "gu"), ("MX", "fz"), ("tile", "fz")]


def _norm_eps(s):
    s = (s or "").strip()
    return s.replace("1e-05","1e-5").replace("1e-06","1e-6") \
            .replace("1e-07","1e-7").replace("1e-08","1e-8")


def load_perf():
    out = []
    for r in csv.DictReader(open(PERF_CSV)):
        if r.get("status") and r["status"] != "OK": continue
        try:
            r["N"]    = int(r["N"]); r["nb"] = int(r["nb"])
            r["eps"]  = _norm_eps(r["eps"])
            r["perf"] = float(r["perf_TFLOPS"])
        except (TypeError, ValueError):
            continue
        out.append(r)
    return out


def load_error():
    return list(csv.DictReader(open(ERROR_CSV)))


def find_err(err_rows, mode_name_err, n, nb, eps):
    """Return rel_factor_error preferring (gran, uflow) per PREFERRED_GU."""
    pref = PREFERRED_GU.get(mode_name_err)
    chain = [pref] + [g for g in FALLBACKS if g != pref]
    for gran, uflow in chain:
        for r in err_rows:
            if (r["mode"] == mode_name_err and r["granularity"] == gran
                    and r["underflow"] == uflow
                    and int(r["N"]) == n and int(r["nb"]) == nb
                    and r["eps"] == eps):
                try: return float(r["rel_factor_error"])
                except: return None
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures/paper_perf_vs_error")
    args = ap.parse_args()

    perf  = load_perf()
    err   = load_error()

    with mpl.rc_context(PAPER_RC):
        # Each subplot gets its OWN y-axis so within-panel trends are clear.
        fig, axes = plt.subplots(1, len(NS_TO_PLOT), figsize=(8.0, 2.7),
                                 sharey=False)

        for ai, (ax, N) in enumerate(zip(axes, NS_TO_PLOT)):
            nb = NB_BY_N[N]
            panel_perf_max = 0.0
            panel_err = []
            for mode_perf, mode_err in MODES:
                pts = []
                for eps in EPS_ORDER:
                    p = next((r for r in perf
                              if r.get("mode") == mode_perf
                              and r["N"] == N and r["eps"] == eps), None)
                    if p is None: continue
                    e = find_err(err, mode_err, N, nb, eps)
                    if e is None or e <= 0: continue
                    pts.append((eps, e, p["perf"]))
                    panel_perf_max = max(panel_perf_max, p["perf"])
                    panel_err.append(e)
                if not pts: continue
                pts.sort(key=lambda x: float(x[0]), reverse=True)  # loose-eps first
                xs = [x[1] for x in pts]
                ys = [x[2] for x in pts]
                ax.plot(xs, ys,
                        linestyle=MODE_LINESTYLE[mode_perf],
                        marker=MODE_MARKER[mode_perf], markersize=7,
                        linewidth=1.7,
                        color=MODE_COLOR[mode_perf],
                        markerfacecolor=MODE_COLOR[mode_perf],
                        markeredgecolor="black", markeredgewidth=0.5,
                        label=mode_perf)

            # Log x-axis (errors span ~5 orders of magnitude, linear would
            # collapse the tight-eps points to the left edge).  Original
            # direction: small (tight) errors on the left, loose on right.
            ax.set_xscale("log")
            ax.set_axisbelow(True)
            ax.grid(True, alpha=0.22, linewidth=0.4, which="both")
            ax.text(0.04, 0.965,
                    f"$N = {N_LABEL[N]}$",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=9, fontweight="bold")
            ax.set_xlabel("Relative factor error  (log)")

            # Per-panel y-axis range, fitted tightly to its own perf data.
            # Pad ~15% above and below the panel min/max so the data band
            # fills the subplot (e.g. 32k: ~25-37 -> ylim ~ (20, 40)).
            panel_ys = [p[2] for mp in MODES for p in []]  # placeholder
            # Recover panel_ys from the axes' lines.
            panel_ys = []
            for ln in ax.get_lines():
                panel_ys.extend([y for y in ln.get_ydata()])
            if panel_ys:
                lo, hi = min(panel_ys), max(panel_ys)
                pad = max((hi - lo) * 0.15, 1.0)
                ax.set_ylim(lo - pad, hi + pad)

            # Vertical dashed lines at the source-eps targets.  Replace the
            # default log-axis decade ticks with explicit ε labels so each
            # dashed line is annotated and the reader can see which target
            # ε each cluster of points belongs to.
            eps_floats = [float(e) for e in EPS_ORDER]
            for ev in eps_floats:
                ax.axvline(ev, color="#999999", linestyle="--",
                           linewidth=0.7, alpha=0.6, zorder=0)
            ax.set_xticks(eps_floats)
            ax.set_xticklabels([rf"$\varepsilon = {EPS_DISPLAY[e][1:-1]}$"
                                for e in EPS_ORDER], fontsize=7)
            # Suppress matplotlib's minor decade tick labels.
            ax.minorticks_off()

        # Y-label only on the leftmost subplot; the other panels keep their
        # own tick numbers (per-panel y-range) but no repeated axis label.
        axes[0].set_ylabel("Perf (TFLOPS)")

        mode_handles = [Line2D([0], [0], marker=MODE_MARKER[mp],
                               linestyle=MODE_LINESTYLE[mp],
                               color=MODE_COLOR[mp],
                               markerfacecolor=MODE_COLOR[mp],
                               markeredgecolor="black", markeredgewidth=0.4,
                               markersize=6, label=mp)
                        for mp, _ in MODES]
        eps_handle = Line2D([0], [0], color="#999999", linestyle="--",
                            linewidth=0.9, label="source $\\varepsilon$")

        # Wider wspace so each panel's own y-tick numbers don't crowd the
        # plotting area of the panel to their left.
        fig.subplots_adjust(left=0.075, right=0.985, top=0.92,
                            bottom=0.36, wspace=0.22)
        fig.legend(handles=mode_handles + [eps_handle],
                   loc="upper center", bbox_to_anchor=(0.5, 0.13),
                   ncol=5, frameon=True, framealpha=0.95,
                   handletextpad=0.5, columnspacing=1.4,
                   labelspacing=0.25, borderpad=0.3)

        out = Path(args.out).with_suffix(".pdf")
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(out)
        plt.close(fig)


if __name__ == "__main__":
    main()
