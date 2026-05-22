#!/usr/bin/env python3
"""Memory sweep figure across all matrix sizes (NB=4096).

Reads tile counts from the [TILE_HIST] lines in
  recreate-ooc-chol/build/runs_all.log
(authoritative; the runs_all.tsv aggregator has a parser bug that drops some
FP16 counts).

Subplots: one per epsilon (default by=eps) or one per mode (by=mode).
Bars within each N group: the 4 modes / 4 epsilons respectively in this
order:
  1. Baseline                  -> fp8
  2. Ladder IEEE               -> ladder_ieee_only
  3. Ladder MX  (no MXFP16)    -> ladder_mx_staircase
  4. Ladder MX+MXFP16          -> ladder_mx_staircase_mxfp16

Memory accounting: Cholesky's lower triangle of L is stored once.
  - Off-diagonal lower-tri tiles contribute full nb*nb element storage.
  - Diagonal tiles (FP64) contribute only nb*(nb+1)/2 elements.
  - MX formats add 1 byte per 32-element shared-scale group.

FP16 -> MXFP16 remap: the binary always logs the second-highest tier as
"FP16" even for the full MX ladder, where those tiles are actually MXFP16
(with shared-scale meta).  The remap fixes the per-format byte counts, the
legend, and the bar color for that mode.

PDF only.
"""
import argparse
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch


# ---------- Plot rcParams (legacy 1xN layout) ----------

mpl.rcParams.update({
    "font.family": "serif",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 9.5,
    "legend.title_fontsize": 10.5,
    "figure.dpi": 120,
    "savefig.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
})


LOG_PATH = "/home/abduraa/MX_project/recreate-ooc-chol/build/runs_all.log"
EPS_ORDER = ["1e-5", "1e-6", "1e-7", "1e-8"]
EPS_DISPLAY = {"1e-5": r"$10^{-5}$", "1e-6": r"$10^{-6}$",
               "1e-7": r"$10^{-7}$", "1e-8": r"$10^{-8}$"}
EPS_COLOR = {"1e-5": "#0072B2", "1e-6": "#E69F00",
             "1e-7": "#009E73", "1e-8": "#D55E00"}
N_ORDER = [20480, 32768, 40960, 65536, 81920, 98304, 122880]


MODE_BARS = [
    ("Baseline",         "fp8"),
    ("Ladder IEEE",      "ladder_ieee_only"),
    ("Ladder MX",        "ladder_mx_staircase"),
    ("Ladder MX+MXFP16", "ladder_mx_staircase_mxfp16"),
]
MODE_PANEL_LABEL = {
    "fp8":                        "Baseline\n(FP64 -> FP32 -> FP16 -> FP8 plain)",
    "ladder_ieee_only":           "Ladder IEEE\n(IEEE rungs only)",
    "ladder_mx_staircase":        "Ladder MX\n(MX rungs, no MXFP16)",
    "ladder_mx_staircase_mxfp16": "Ladder MX+MXFP16\n(... -> MXFP16 -> MXFP8 -> MXFP4)",
}
MODE_PANEL_LABEL_SHORT = {
    "fp8":                        "Baseline",
    "ladder_ieee_only":           "Ladder IEEE",
    "ladder_mx_staircase":        "Ladder MX",
    "ladder_mx_staircase_mxfp16": "Ladder MX+MXFP16",
}


# ---------- Format catalogue ----------

# Bytes per element + whether the format carries a 1-byte/32-elem shared scale.
FMT_BYTES = {
    "FP64":      (8.0, False),
    "FP32":      (4.0, False),
    "FP16":      (2.0, False),
    "MXFP16":    (2.0, True),
    "MXFP8":     (1.0, True),
    "MXFP4":     (0.5, True),
    "FP8_plain": (1.0, False),
}
# Stack order (bottom -> top) and colors.
STACK = [
    ("FP64",      "#0072B2"),
    ("FP32",      "#E69F00"),
    ("FP16",      "#7A2C00"),
    ("MXFP16",    "#9E5530"),
    ("FP8_plain", "#56B4E9"),
    ("MXFP8",     "#009E73"),
    ("MXFP4",     "#F0E442"),
    ("Scale",     "#999999"),   # synthesised at the end
]
STACK_LABEL = {
    "FP64":      "FP64",
    "FP32":      "FP32",
    "FP16":      "FP16",
    "MXFP16":    "MXFP16",
    "FP8_plain": "FP8 plain",
    "MXFP8":     "MXFP8",
    "MXFP4":     "MXFP4",
    "Scale":     "Scale",
}
# Hatch encoding for B&W printability. MX-scaled formats are hatched.
STACK_HATCH = {
    "FP64":      "",
    "FP32":      "",
    "FP16":      "",
    "MXFP16":    "//",
    "FP8_plain": "",
    "MXFP8":     "xx",
    "MXFP4":     "..",
    "Scale":     "++",
}


# ---------- Log parsing ----------

HEAD_RE = re.compile(r"^===\s+(\S+)\s+N=(\d+)\s+eps=(\S+)\s+")
HIST_RE = re.compile(
    r"^\[TILE_HIST\]\s+N=(\d+)\s+nb=(\d+)\s+eps=\S+\s+nt=\d+\s+total=\d+"
    r"\s+FP64=(\d+)\s+FP32=(\d+)\s+FP16=(\d+)\s+MXFP16=(\d+)"
    r"\s+MXFP8=(\d+)\s+MXFP4=(\d+)\s+FP8_plain=(\d+)"
)


def parse_log(path):
    """Return {(mode, N, eps): counts_dict} from runs_all.log."""
    rows = {}
    cur_key = None  # (mode, N, eps)
    with Path(path).open() as f:
        for line in f:
            m = HEAD_RE.match(line)
            if m:
                mode = m.group(1)
                N    = int(m.group(2))
                eps  = m.group(3).strip()
                eps = eps.replace("1e-05", "1e-5").replace("1e-06", "1e-6") \
                         .replace("1e-07", "1e-7").replace("1e-08", "1e-8")
                cur_key = (mode, N, eps)
                continue
            m = HIST_RE.match(line)
            if not m or cur_key is None:
                continue
            N_hist = int(m.group(1))
            nb     = int(m.group(2))
            if N_hist != cur_key[1]:
                continue
            counts = {
                "FP64":      int(m.group(3)),
                "FP32":      int(m.group(4)),
                "FP16":      int(m.group(5)),
                "MXFP16":    int(m.group(6)),
                "MXFP8":     int(m.group(7)),
                "MXFP4":     int(m.group(8)),
                "FP8_plain": int(m.group(9)),
                "nb":        nb,
            }
            # The TILE_HIST output always reports the second-highest tier as
            # "FP16" regardless of whether the binary used plain IEEE FP16 or
            # MXFP16.  For the full MX ladder, those tiles are actually MXFP16
            # (scale-meta -bearing); remap so memory + color + legend reflect
            # that.
            if cur_key[0] == "ladder_mx_staircase_mxfp16" and counts["FP16"] > 0:
                counts["MXFP16"] += counts["FP16"]
                counts["FP16"]   = 0
            rows[cur_key] = counts
    return rows


# ---------- Memory accounting ----------

def memory_breakdown_gb(counts, N):
    """Return {format_label: GB} including a 'Scale' bucket for shared scales.

    Diagonal accounting: the M FP64 tiles each store only nb*(nb+1)/2 elements
    (lower triangle + diagonal); off-diagonal lower-tri tiles store full nb*nb.
    """
    nb = counts["nb"]
    M = N // nb
    tile_full = nb * nb
    tile_diag = nb * (nb + 1) // 2
    GB = 1.0 / (1024 ** 3)
    out = {fmt: 0.0 for fmt, _ in STACK}

    for fmt, (bpe, has_scale) in FMT_BYTES.items():
        cnt = counts.get(fmt, 0)
        if cnt <= 0:
            continue
        if fmt == "FP64":
            diag = min(M, cnt)
            off  = max(cnt - diag, 0)
            data_bytes = diag * tile_diag * bpe + off * tile_full * bpe
        else:
            data_bytes = cnt * tile_full * bpe
        out[fmt] += data_bytes * GB
        if has_scale:
            out["Scale"] += cnt * (tile_full // 32) * GB
    return out


# ---------- Legacy 1xN renderer (kept for compatibility) ----------

def draw_subplot_by_eps(ax, eps, data, ns, show_ylabel):
    n_modes = len(MODE_BARS)
    group_w = 0.78
    bar_w   = group_w / n_modes
    x_centres = list(range(len(ns)))

    drawn_fmts = set()
    panel_max = 0.0
    for ni, N in enumerate(ns):
        for bi, (mode_lbl, mode_key) in enumerate(MODE_BARS):
            x = ni - group_w/2 + (bi + 0.5) * bar_w
            counts = data.get((mode_key, N, eps))
            if counts is None:
                continue
            gb = memory_breakdown_gb(counts, N)
            bottom = 0.0
            for fmt, color in STACK:
                v = gb.get(fmt, 0.0)
                if v <= 0:
                    continue
                drawn_fmts.add(fmt)
                ax.bar([x], [v], bar_w * 0.92,
                       bottom=[bottom],
                       color=color, edgecolor="black", linewidth=0.3)
                bottom += v
            panel_max = max(panel_max, bottom)

    ax.set_xticks(x_centres)
    ax.set_xticklabels([f"{N // 1024}k" for N in ns])
    ax.set_xlim(-0.6, len(ns) - 0.4)
    if show_ylabel:
        ax.set_ylabel("Memory  (GB, lower triangle + half-diagonal)")
    ax.set_xlabel(r"Matrix size  $N\ (\times 10^{3})$")
    ax.set_axisbelow(True)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, alpha=0.3, linewidth=0.5)
    return drawn_fmts, panel_max


def draw_subplot_by_mode(ax, mode_key, data, ns, show_ylabel):
    n_eps = len(EPS_ORDER)
    group_w = 0.78
    bar_w   = group_w / n_eps
    x_centres = list(range(len(ns)))

    drawn_fmts = set()
    panel_max = 0.0
    for ni, N in enumerate(ns):
        for bi, eps in enumerate(EPS_ORDER):
            x = ni - group_w/2 + (bi + 0.5) * bar_w
            counts = data.get((mode_key, N, eps))
            if counts is None:
                continue
            gb = memory_breakdown_gb(counts, N)
            bottom = 0.0
            for fmt, color in STACK:
                v = gb.get(fmt, 0.0)
                if v <= 0:
                    continue
                drawn_fmts.add(fmt)
                ax.bar([x], [v], bar_w * 0.92,
                       bottom=[bottom],
                       color=color, edgecolor="black", linewidth=0.3)
                bottom += v
            if bottom > 0 and ni == len(ns) - 1:
                ax.plot([x], [bottom], marker="v", markersize=4,
                        color=EPS_COLOR[eps], markeredgecolor="black",
                        markeredgewidth=0.4, zorder=10)
            panel_max = max(panel_max, bottom)

    ax.set_xticks(x_centres)
    ax.set_xticklabels([f"{N // 1024}k" for N in ns])
    ax.set_xlim(-0.6, len(ns) - 0.4)
    if show_ylabel:
        ax.set_ylabel("Memory  (GB, lower triangle + half-diagonal)")
    ax.set_xlabel(r"Matrix size  $N\ (\times 10^{3})$")
    ax.set_axisbelow(True)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, alpha=0.3, linewidth=0.5)
    return drawn_fmts, panel_max


def render_legacy(data, by, out_path):
    if by == "eps":
        panel_keys  = EPS_ORDER
        panel_label = lambda k: f"$\\varepsilon = {EPS_DISPLAY[k][1:-1]}$"
        right_title = "Bars within each N group (left -> right)"
        right_items = [(f"{i+1}. {lbl}", "white")
                       for i, (lbl, _) in enumerate(MODE_BARS)]
    elif by == "mode":
        panel_keys  = [mk for _, mk in MODE_BARS]
        panel_label = lambda k: MODE_PANEL_LABEL[k]
        right_title = "Bars within each N group (left -> right)"
        right_items = [(f"{i+1}. $\\varepsilon = {EPS_DISPLAY[e][1:-1]}$", "white")
                       for i, e in enumerate(EPS_ORDER)]
    else:
        raise ValueError(by)

    n_panels = len(panel_keys)
    fig, axes = plt.subplots(1, n_panels,
                             figsize=(4.7 * n_panels + 0.8, 5.0),
                             sharey=True)
    if n_panels == 1:
        axes = [axes]

    drawn_all = set()
    global_max = 0.0
    for ax, key in zip(axes, panel_keys):
        if by == "eps":
            drawn, panel_max = draw_subplot_by_eps(
                ax, key, data, N_ORDER, show_ylabel=(ax is axes[0]))
        else:
            drawn, panel_max = draw_subplot_by_mode(
                ax, key, data, N_ORDER, show_ylabel=(ax is axes[0]))
        drawn_all |= drawn
        global_max = max(global_max, panel_max)
        ax.text(0.025, 0.965, panel_label(key),
                transform=ax.transAxes, ha="left", va="top",
                fontsize=10.5 if by == "mode" else 11.0,
                bbox=dict(boxstyle="round,pad=0.34",
                          facecolor="white", edgecolor="#888",
                          linewidth=0.8, alpha=0.92))

    if global_max > 0:
        for ax in axes:
            ax.set_ylim(0, global_max * 1.12)

    fmt_handles = [Patch(facecolor=c, edgecolor="black", label=STACK_LABEL[f])
                   for f, c in STACK if f in drawn_all]
    right_handles = [Patch(facecolor=c, edgecolor="black", label=lbl)
                     for lbl, c in right_items]

    legend_ncol = 2 if len(fmt_handles) >= 4 else len(fmt_handles)
    leg1 = fig.legend(handles=fmt_handles,
                      loc="lower center", bbox_to_anchor=(0.5, -0.08),
                      ncol=legend_ncol, frameon=True, framealpha=0.95,
                      title="Tile-format stack (bottom -> top)",
                      handletextpad=0.6, columnspacing=2.4,
                      labelspacing=0.4, borderpad=0.4)
    fig.add_artist(leg1)
    fig.legend(handles=right_handles,
               loc="lower center", bbox_to_anchor=(0.78, -0.08),
               ncol=2, frameon=True, framealpha=0.95,
               title=right_title,
               handletextpad=0.6, columnspacing=1.4,
               labelspacing=0.3, borderpad=0.4)

    fig.tight_layout(rect=[0, 0.06, 1, 1.0])

    out = Path(out_path).with_suffix(".pdf")
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=300, bbox_inches="tight")
    print(out)
    plt.close(fig)


# ---------- Paper-ready 2x2 renderer ----------

PAPER_RC = {
    "font.family": "serif",
    "font.size": 9,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "xtick.labelsize": 7.5,
    "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5,
    "legend.title_fontsize": 8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linewidth": 0.4,
}


def _draw_paper_subplot(ax, by, key, data, ns):
    if by == "eps":
        slots = [(lbl, mk) for lbl, mk in MODE_BARS]
        slot_key = lambda slot: (slot[1], None, key)
    else:
        slots = [(EPS_DISPLAY[e], e) for e in EPS_ORDER]
        slot_key = lambda slot: (key, None, slot[1])

    n_slots = len(slots)
    group_w = 0.78
    bar_w   = group_w / n_slots
    x_centres = list(range(len(ns)))

    drawn_fmts = set()
    panel_max = 0.0
    for ni, N in enumerate(ns):
        for bi, slot in enumerate(slots):
            x = ni - group_w/2 + (bi + 0.5) * bar_w
            mkey = slot_key(slot)
            counts = data.get((mkey[0], N, mkey[2]))
            if counts is None:
                continue
            gb = memory_breakdown_gb(counts, N)
            bottom = 0.0
            for fmt, color in STACK:
                v = gb.get(fmt, 0.0)
                if v <= 0:
                    continue
                drawn_fmts.add(fmt)
                ax.bar([x], [v], bar_w * 0.92,
                       bottom=[bottom],
                       color=color, edgecolor="black",
                       linewidth=0.25, hatch=STACK_HATCH.get(fmt, ""))
                bottom += v
            panel_max = max(panel_max, bottom)

    # Light vertical separators between N groups.
    for i in range(1, len(ns)):
        ax.axvline(i - 0.5, color="#bbbbbb", linewidth=0.5, alpha=0.6, zorder=0)

    ax.set_xticks(x_centres)
    ax.set_xticklabels([f"{N // 1024}k" for N in ns])
    ax.set_xlim(-0.6, len(ns) - 0.4)
    ax.set_axisbelow(True)
    ax.xaxis.grid(False)
    ax.yaxis.grid(True, alpha=0.25, linewidth=0.4)
    return drawn_fmts, panel_max


def render_paper_2x2(data, by, out_path):
    if by == "eps":
        panel_keys   = EPS_ORDER
        panel_titles = [f"$\\varepsilon = {EPS_DISPLAY[e][1:-1]}$" for e in EPS_ORDER]
        right_title  = "Bars (left $\\rightarrow$ right within each $N$)"
        right_items  = [f"{i+1}. {lbl}" for i, (lbl, _) in enumerate(MODE_BARS)]
    elif by == "mode":
        panel_keys   = [mk for _, mk in MODE_BARS]
        panel_titles = [MODE_PANEL_LABEL_SHORT[mk] for mk in panel_keys]
        right_title  = "Bars (left $\\rightarrow$ right within each $N$)"
        right_items  = [f"{i+1}. $\\varepsilon = {EPS_DISPLAY[e][1:-1]}$"
                        for i, e in enumerate(EPS_ORDER)]
    else:
        raise ValueError(by)

    with mpl.rc_context(PAPER_RC):
        fig, axes = plt.subplots(2, 2, figsize=(7.16, 4.2),
                                 sharex=True, sharey=True)
        axes_flat = axes.flatten()

        drawn_all = set()
        global_max = 0.0
        for ai, (ax, key, title) in enumerate(zip(axes_flat, panel_keys, panel_titles)):
            drawn, panel_max = _draw_paper_subplot(ax, by, key, data, N_ORDER)
            drawn_all |= drawn
            global_max = max(global_max, panel_max)
            ax.text(0.02, 0.965,
                    f"({chr(ord('a') + ai)})  {title}",
                    transform=ax.transAxes, ha="left", va="top",
                    fontsize=8.5, fontweight="bold")

        if global_max > 0:
            for ax in axes_flat:
                ax.set_ylim(0, global_max * 1.10)

        for ax in axes[-1, :]:
            ax.set_xlabel(r"Matrix size  $N$")
        for ax in axes[:, 0]:
            ax.set_ylabel("Memory (GB)")

        # Force FP16+MXFP16 to be adjacent in the legend.
        legend_drawn = set(drawn_all)
        if "FP16" in legend_drawn or "MXFP16" in legend_drawn:
            legend_drawn.add("FP16")
            legend_drawn.add("MXFP16")
        fmt_handles = [Patch(facecolor=c, edgecolor="black",
                             hatch=STACK_HATCH.get(f, ""),
                             label=STACK_LABEL[f])
                       for f, c in STACK if f in legend_drawn]
        bar_handles = [Patch(facecolor="white", edgecolor="black", label=lbl)
                       for lbl in right_items]

        leg1 = fig.legend(handles=fmt_handles,
                          loc="lower center", bbox_to_anchor=(0.27, -0.10),
                          ncol=4, frameon=True, framealpha=0.95,
                          title="Tile format (bottom $\\rightarrow$ top of stack)",
                          handletextpad=0.5, columnspacing=1.2,
                          labelspacing=0.3, borderpad=0.35)
        fig.add_artist(leg1)
        fig.legend(handles=bar_handles,
                   loc="lower center", bbox_to_anchor=(0.78, -0.10),
                   ncol=2, frameon=True, framealpha=0.95,
                   title=right_title,
                   handletextpad=0.5, columnspacing=1.2,
                   labelspacing=0.3, borderpad=0.35)

        fig.tight_layout(rect=[0, 0.07, 1, 1.0])

        out = Path(out_path).with_suffix(".pdf")
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300, bbox_inches="tight")
        print(out)
        plt.close(fig)


# ---------- Main ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", default=LOG_PATH)
    ap.add_argument("--by", choices=["eps", "mode", "both"], default="both",
                    help="Subplot dimension for the legacy 1xN layout.")
    ap.add_argument("--paper-2x2", action="store_true", default=True,
                    help="Emit paper-ready 2x2 variants (*_2x2.pdf).  Default: on.")
    ap.add_argument("--no-paper-2x2", dest="paper_2x2", action="store_false")
    ap.add_argument("--legacy", action="store_true", default=True,
                    help="Emit the legacy 1xN layout PDFs.  Default: on.")
    ap.add_argument("--no-legacy", dest="legacy", action="store_false",
                    help="Skip the legacy 1xN PDFs (only emit *_2x2.pdf).")
    ap.add_argument("--out-dir",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures")
    args = ap.parse_args()

    data = parse_log(args.log)
    targets = ["eps", "mode"] if args.by == "both" else [args.by]

    if args.legacy:
        for by in targets:
            out = Path(args.out_dir) / f"paper_memory_nb4096_sweep_by_{by}"
            render_legacy(data, by, out)

    if args.paper_2x2:
        for by in targets:
            out = Path(args.out_dir) / f"paper_memory_nb4096_sweep_by_{by}_2x2"
            render_paper_2x2(data, by, out)


if __name__ == "__main__":
    main()
