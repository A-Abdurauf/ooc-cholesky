#!/usr/bin/env python3
import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch


def canon_fmt(name: str) -> str:
    n = (name or "").strip().lower()
    alias = {
        "mx_f16": "mx_fp16",
        "mx_f32": "mx_fp32",
        "mx_fp8_e4m3": "mx_e4m3",
        "mx_fp8_e5m2": "mx_e5m2",
        "fp8e4m3": "fp8_e4m3",
        "fp8e5m2": "fp8_e5m2",
        "fp6_e3m2": "e3m2",
        "fp6e3m2": "e3m2",
        "fp6_e2m3": "e2m3",
        "fp6e2m3": "e2m3",
        "fp4_e2m1": "e2m1",
        "fp4e2m1": "e2m1",
    }
    return alias.get(n, n)


def bits_from_fmt(fmt: str) -> int:
    f = canon_fmt(fmt)
    if f == "fp64":
        return 64
    if f in ("fp32", "mx_fp32"):
        return 32
    if f in ("fp16", "bf16", "mx_fp16"):
        return 16
    if f in ("mx_e4m3", "mx_e5m2", "fp8_e4m3", "fp8_e5m2"):
        return 8
    if f in ("e3m2", "e2m3"):
        return 6
    if f == "e2m1":
        return 4
    return 8


def is_mx(fmt: str) -> bool:
    return canon_fmt(fmt).startswith("mx_")


def is_mx_fp32(fmt: str) -> bool:
    return canon_fmt(fmt) == "mx_fp32"


def infer_subtile_from_mode(mode: str) -> int:
    m = (mode or "").strip().lower()
    if m.startswith("subtile_"):
        try:
            return int(m.split("_", 1)[1])
        except Exception:
            return 0
    return 0


def scale_groups_per_tile(nb: int, granularity: str, mode: str, subtile_override: int) -> int:
    if granularity == "tile":
        return 1
    subtile = subtile_override if subtile_override > 0 else infer_subtile_from_mode(mode)
    if granularity == "subtile" and subtile <= 0:
        raise ValueError("--scale-granularity=subtile requires mode=subtile_* or --subtile-size")
    if granularity == "auto" and subtile <= 0:
        return 1
    nblocks = int(math.ceil(float(nb) / float(subtile)))
    return nblocks * nblocks


def parse_counts(s: str) -> dict:
    out = {}
    for part in (s or "").split(","):
        if "=" not in part:
            continue
        k, v = part.split("=", 1)
        k = canon_fmt(k)
        try:
            out[k] = out.get(k, 0) + int(v)
        except ValueError:
            pass
    return out


def eps_sort_key(eps: str):
    try:
        return float(eps)
    except Exception:
        return float("inf")


def classify(fmt: str) -> str:
    f = canon_fmt(fmt)
    if f == "fp64":
        return "FP64"
    if f in ("fp32", "mx_fp32"):
        return "FP32"
    if f in ("fp16", "bf16", "mx_fp16"):
        return "FP16/BF16"
    if f in ("mx_e4m3", "mx_e5m2", "fp8_e4m3", "fp8_e5m2"):
        return "MXFP8/FP8"
    if f in ("e3m2", "e2m3"):
        return "MXFP6"
    if f == "e2m1":
        return "MXFP4"
    return "OTHER"


def build_group_rows(raw_lines, header, rng_expr, nb, mode):
    picked = []
    for piece in [x.strip() for x in rng_expr.split("+") if x.strip()]:
        a, b = piece.split("-", 1)
        lo, hi = min(int(a), int(b)), max(int(a), int(b))
        for ln in range(lo, hi + 1):
            idx = ln - 1
            if 0 <= idx < len(raw_lines):
                line = raw_lines[idx].strip()
                if line:
                    picked.append(line)
    rows = list(csv.DictReader([header] + picked, delimiter="\t"))
    out = []
    for r in rows:
        eps = (r.get("source_epsilon") or "").strip()
        if not eps or eps.lower() == "na":
            continue
        try:
            rnb = int((r.get("nb") or "0").strip())
        except ValueError:
            continue
        if nb and rnb != nb:
            continue
        if mode:
            if (r.get("mx_mode") or "").strip() != mode:
                continue
        out.append(r)
    return out


def main():
    ap = argparse.ArgumentParser(description="Stacked memory breakdown by datatype across epsilon")
    ap.add_argument("--summary", required=True)
    ap.add_argument("--compare-ranges", required=True, help='e.g. "44-51+60-67,52-59,36-43"')
    ap.add_argument("--compare-labels", default="baseline,ieee,auto")
    ap.add_argument("--nb", type=int, default=1024)
    ap.add_argument("--mode", default="")
    ap.add_argument("--use-full", action="store_true")
    ap.add_argument("--mx-scale-bits", type=int, default=8)
    ap.add_argument("--mx-fp32-scale-bits", type=int, default=11)
    ap.add_argument("--scale-granularity", choices=["auto", "tile", "subtile"], default="auto")
    ap.add_argument("--subtile-size", type=int, default=0)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    summary = Path(args.summary)
    raw = summary.read_text().splitlines()
    if not raw:
        raise SystemExit("Empty summary")
    header = raw[0]
    ranges = [s.strip() for s in args.compare_ranges.split(",") if s.strip()]
    labels = [s.strip() for s in args.compare_labels.split(",") if s.strip()]
    if len(labels) != len(ranges):
        labels = [f"group{i+1}" for i in range(len(ranges))]

    count_col = "tile_counts_full" if args.use_full else "tile_counts_lower"
    bucket_order = ["FP64", "FP32", "FP16/BF16", "MXFP8/FP8", "MXFP6", "MXFP4", "OTHER", "SCALE_META"]
    colors = {
        "FP64": "#4e79a7",
        "FP32": "#f28e2b",
        "FP16/BF16": "#e15759",
        "MXFP8/FP8": "#76b7b2",
        "MXFP6": "#59a14f",
        "MXFP4": "#edc948",
        "OTHER": "#b07aa1",
        "SCALE_META": "#9c755f",
    }

    all_eps = set()
    per_group = {}
    tile_elements = args.nb * args.nb

    for rng, label in zip(ranges, labels):
        rows = build_group_rows(raw, header, rng, args.nb, args.mode)
        eps_to_acc = defaultdict(lambda: defaultdict(float))
        eps_to_n = defaultdict(int)
        for r in rows:
            eps = (r.get("source_epsilon") or "").strip()
            counts = parse_counts(r.get(count_col, ""))
            gpt = scale_groups_per_tile(args.nb, args.scale_granularity,
                                        (r.get("mx_mode") or "").strip(),
                                        args.subtile_size)
            for fmt, cnt in counts.items():
                bits = bits_from_fmt(fmt)
                gb = (cnt * tile_elements * bits / 8.0) / 1e9
                eps_to_acc[eps][classify(fmt)] += gb
                if is_mx(fmt):
                    sb = args.mx_fp32_scale_bits if is_mx_fp32(fmt) else args.mx_scale_bits
                    eps_to_acc[eps]["SCALE_META"] += (cnt * sb * gpt / 8.0) / 1e9
            eps_to_n[eps] += 1
            all_eps.add(eps)

        # average if merged ranges contribute multiple rows per epsilon
        out = {}
        for eps, comp in eps_to_acc.items():
            n = float(max(1, eps_to_n[eps]))
            out[eps] = {k: v / n for k, v in comp.items()}
        per_group[label] = out

    eps_order = sorted(all_eps, key=eps_sort_key)
    n_groups = len(labels)
    fig, ax = plt.subplots(figsize=(12, 5.8))
    x = list(range(len(eps_order)))
    group_width = 0.82
    bar_w = group_width / max(1, n_groups)
    hatches = ["", "//", "xx", "..", "\\", "++"]

    for gi, label in enumerate(labels):
        xpos = [xi - group_width / 2 + (gi + 0.5) * bar_w for xi in x]
        bottoms = [0.0] * len(eps_order)
        for b in bucket_order:
            vals = [per_group.get(label, {}).get(eps, {}).get(b, 0.0) for eps in eps_order]
            if not any(v > 0 for v in vals):
                continue
            ax.bar(
                xpos,
                vals,
                bottom=bottoms,
                color=colors[b],
                width=bar_w * 0.95,
                hatch=hatches[gi % len(hatches)],
                edgecolor="black",
                linewidth=0.2,
            )
            bottoms = [bb + vv for bb, vv in zip(bottoms, vals)]
        ypad = 0.01 * max(bottoms) if bottoms and max(bottoms) > 0 else 0.001
        for xi, tot in enumerate(bottoms):
            if tot > 0:
                ax.text(xpos[xi], tot + ypad, f"{tot:.3f}", ha="center", va="bottom", fontsize=7, rotation=90)

    ax.set_xticks(list(range(len(eps_order))))
    ax.set_xticklabels(eps_order)
    ax.set_xlabel("source_epsilon")
    ax.set_ylabel("GB")
    ax.grid(axis="y", alpha=0.25)

    dtype_legend = [Patch(facecolor=colors[k], edgecolor="black", label=k) for k in bucket_order if any(
        per_group.get(lbl, {}).get(eps, {}).get(k, 0.0) > 0.0 for lbl in labels for eps in eps_order
    )]
    group_legend = [Patch(facecolor="white", edgecolor="black", hatch=hatches[i % len(hatches)], label=labels[i]) for i in range(n_groups)]
    lg1 = ax.legend(handles=dtype_legend, loc="upper right", title="Datatype")
    ax.add_artist(lg1)
    ax.legend(handles=group_legend, loc="upper left", title="Group")
    mode_title = f", mode={args.mode}" if args.mode else ""
    fig.suptitle(f"Memory breakdown by datatype (NB={args.nb}, {count_col}{mode_title})", y=0.99)
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    print(str(out))


if __name__ == "__main__":
    main()
