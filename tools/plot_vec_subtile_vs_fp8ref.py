#!/usr/bin/env python3
"""Bar chart: relative error vs epsilon for SFP8/SFP6/SFP4 variants + plain FP8 reference.

SFP = Scaled Floating Point (MX format):
  SFP8  = MX-FP8 E4M3 (mx_e4m3)
  SFP6  = MX-FP6 E3M2 (e3m2)
  SFP4  = MX-FP4 E2M1 (e2m1)
  FP8   = plain FP8 E4M3 (no microscaling, fp8_e4m3)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT_DIR = Path("/home/abduraa/MX_project/logs/figures")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Data extracted from logs ─────────────────────────────────────────────────
# Keys: eps → relative_error.  "default" bucket = fp32/fp16.

# x-axis: loose → tight (1e-5 first)
EPS_LABELS = ["1e-5", "1e-6", "1e-7", "1e-8"]
EPS_VALS   = [1e-5,   1e-6,   1e-7,   1e-8  ]

DATA = {
    "weak": {
        "FP8-E4M3 (Baseline)":                  [3.61837e-4,  1.97798e-6, 1.56698e-6, 3.17131e-7],
        "SFP8-E4M3 (MX)":            [9.36875e-6,  1.26366e-6, 1.0874e-6,  2.14469e-7],
        "SFP8-E4M3 (sub_128)":        [1.06662e-5,  1.26394e-6, 1.08777e-6, 2.14469e-7],
        "SFP6-E3M2 (MX)":            [1.64174e-5,  1.31827e-6, 1.10092e-6, 2.1447e-7 ],
        "SFP6-E3M2 (sub_128)":       [3.43307e-5,  1.30752e-6, 1.15991e-6, 2.1447e-7 ],
        "SFP4-E2M1 (MX)":            [3.71574e-5,  1.39938e-6, 1.15916e-6, 2.1447e-7 ],
        "SFP4-E2M1 (sub_128)":       [8.10382e-5,  1.92331e-6, 1.3587e-6,  2.1447e-7 ],
    },
}

# ── Visual style ─────────────────────────────────────────────────────────────
SERIES_ORDER = [
    "FP8-E4M3 (Baseline)",
    "SFP8-E4M3 (MX)",
    "SFP8-E4M3 (sub_128)",
    "SFP6-E3M2 (MX)",
    "SFP6-E3M2 (sub_128)",
    "SFP4-E2M1 (MX)",
    "SFP4-E2M1 (sub_128)",
]

COLORS = {
    "FP8-E4M3 (Baseline)":                  "#666666",
    "SFP8-E4M3 (MX)":            "#1b7837",
    "SFP8-E4M3 (sub_128)":        "#74c476",
    "SFP6-E3M2 (MX)":            "#d73027",
    "SFP6-E3M2 (sub_128)":       "#fc8d59",
    "SFP4-E2M1 (MX)":            "#4393c3",
    "SFP4-E2M1 (sub_128)":       "#92c5de",
}

HATCHES = {
    "FP8-E4M3 (Baseline)":                  "xx",
    "SFP8-E4M3 (MX)":            "",
    "SFP8-E4M3 (sub_128)":        "//",
    "SFP6-E3M2 (MX)":            "",
    "SFP6-E3M2 (sub_128)":       "//",
    "SFP4-E2M1 (MX)":            "",
    "SFP4-E2M1 (sub_128)":       "//",
}

# ── Plot ─────────────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(10, 5.5))

n_series = len(SERIES_ORDER)
n_eps    = len(EPS_VALS)
group_w  = 0.82
bar_w    = group_w / n_series
x        = np.arange(n_eps)

for matrix in ["weak"]:
    mdata = DATA[matrix]

    for j, series in enumerate(SERIES_ORDER):
        vals = mdata[series]
        xpos = x - group_w / 2 + (j + 0.5) * bar_w
        bars = ax.bar(xpos, vals, width=bar_w * 0.92,
                      color=COLORS[series],
                      hatch=HATCHES[series],
                      label=series,
                      edgecolor="white", linewidth=0.4,
                      zorder=3)

ax.set_yscale("log")
ax.yaxis.set_major_formatter(matplotlib.ticker.LogFormatterSciNotation())
ax.set_xticks(x)
ax.set_xticklabels(EPS_LABELS, fontsize=11)
ax.set_xlabel("Source epsilon (threshold, loose → tight)", fontsize=11)
ax.set_ylabel("Relative error", fontsize=11)
ax.set_title("Weak correlation error", fontsize=13)
ax.grid(axis="y", which="both", alpha=0.25, zorder=0)
ax.set_axisbelow(True)

handles, labels = ax.get_legend_handles_labels()
fig.tight_layout(rect=[0, 0.12, 1, 1.0])

fig.legend(handles, labels,
           loc="lower center",
           bbox_to_anchor=(0.5, -0.04),
           ncol=3, fontsize=9.5, frameon=True)

out_png = OUT_DIR / "rel_error_sfp_bar_weak.png"
out_pdf = OUT_DIR / "rel_error_sfp_bar_weak.pdf"
fig.savefig(out_png, dpi=180, bbox_inches="tight")
fig.savefig(out_pdf, bbox_inches="tight", orientation="landscape")
print(f"Saved: {out_png}")
print(f"Saved: {out_pdf}")
