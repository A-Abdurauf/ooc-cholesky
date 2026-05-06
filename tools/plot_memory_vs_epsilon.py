#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
import math
from pathlib import Path

import matplotlib.pyplot as plt


def bits_from_token(token: str, default_bits: int = 8) -> int:
    t = (token or "").strip().lower()
    if not t:
        return default_bits
    if "fp64" in t:
        return 64
    if "fp32" in t or "mx_fp32" in t:
        return 32
    if "mxfp6" in t or "mx_fp6" in t or "fp6" in t or "e3m2" in t or "e2m3" in t:
        return 6
    if "mxfp4" in t or "mx_fp4" in t or "fp4" in t or "e2m1" in t:
        return 4
    if "fp16" in t or "bf16" in t or "mx_fp16" in t:
        return 16
    if "fp8" in t or "e4m3" in t or "e5m2" in t:
        return 8
    return default_bits


def token_is_mx(token: str) -> bool:
    t = (token or "").strip().lower()
    return "mx_" in t


def token_is_mx_fp32(token: str) -> bool:
    t = (token or "").strip().lower()
    return "mx_fp32" in t or "mxfp32" in t


def parse_counts(s: str):
    out = {}
    if not s:
        return out
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
    for part in s.split(','):
        part = part.strip()
        if not part or '=' not in part:
            continue
        k, v = part.split('=', 1)
        k = alias.get(k.strip().lower(), k.strip().lower())
        try:
            out[k] = out.get(k, 0) + int(v.strip())
        except ValueError:
            pass
    return out


def eps_sort_key(eps: str):
    try:
        return float(eps)
    except Exception:
        return float("inf")


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
        raise ValueError("scale granularity 'subtile' requested but subtile size is not available")
    if granularity == "auto" and subtile <= 0:
        return 1
    blocks = int(math.ceil(float(nb) / float(subtile)))
    return blocks * blocks


def calc_memory_bytes(counts: dict, nb: int, mx_scale_bits: int, mx_fp32_scale_bits: int,
                      groups_per_tile: int):
    tile_elements = int(nb) * int(nb)
    data_bits = 0
    scale_bits = 0
    for fmt, cnt in counts.items():
        bits = bits_from_token(fmt)
        data_bits += cnt * tile_elements * bits
        if token_is_mx(fmt):
            per_tile_scale_bits = mx_fp32_scale_bits if token_is_mx_fp32(fmt) else mx_scale_bits
            scale_bits += cnt * per_tile_scale_bits * groups_per_tile
    return (data_bits + scale_bits) / 8.0


def calc_dp64_bytes_for_same_tiles(counts: dict, nb: int):
    tile_elements = int(nb) * int(nb)
    total_tiles = sum(int(v) for v in counts.values())
    return total_tiles * tile_elements * 8


def main():
    ap = argparse.ArgumentParser(description="Compare mixed-memory vs epsilon using summary tile counts.")
    ap.add_argument("--summary", required=True)
    ap.add_argument("--compare-ranges", required=True, help='e.g. "36-43,44-51,52-59"')
    ap.add_argument("--compare-labels", default="auto,legacy,ieee")
    ap.add_argument("--nb", type=int, default=1024)
    ap.add_argument("--mode", default="", help="Optional mx_mode filter, e.g. tile or subtile_128")
    ap.add_argument("--use-full", action="store_true", help="Use tile_counts_full (default uses tile_counts_lower)")
    ap.add_argument("--mx-scale-bits", type=int, default=8)
    ap.add_argument("--mx-fp32-scale-bits", type=int, default=11)
    ap.add_argument("--scale-granularity", choices=["auto", "tile", "subtile"], default="auto",
                    help="Scale metadata accounting granularity: auto/tile/subtile (default: auto)")
    ap.add_argument("--subtile-size", type=int, default=0,
                    help="Optional subtile size for scale granularity (e.g., 128). If 0, infer from mode=subtile_*")
    ap.add_argument("--out", required=True, help="Output figure path")
    ap.add_argument("--report-out", default="", help="Optional TSV report output")
    args = ap.parse_args()

    summary = Path(args.summary)
    raw_lines = summary.read_text().splitlines()
    if not raw_lines:
        raise SystemExit("Summary is empty.")
    header = raw_lines[0]

    ranges = [s.strip() for s in args.compare_ranges.split(',') if s.strip()]
    labels = [s.strip() for s in args.compare_labels.split(',') if s.strip()]
    if len(ranges) < 2:
        raise SystemExit("Need at least two ranges in --compare-ranges")
    if len(labels) != len(ranges):
        labels = [f"group{i+1}" for i in range(len(ranges))]

    count_col = "tile_counts_full" if args.use_full else "tile_counts_lower"

    def pick_rows(rng: str):
        a, b = rng.split('-', 1)
        lo, hi = min(int(a), int(b)), max(int(a), int(b))
        selected = []
        for ln in range(lo, hi + 1):
            idx = ln - 1
            if 0 <= idx < len(raw_lines):
                line = raw_lines[idx].strip()
                if line:
                    selected.append(line)
        rows = list(csv.DictReader([header] + selected, delimiter='\t'))
        out = []
        for r in rows:
            eps = (r.get("source_epsilon") or "").strip()
            if not eps or eps.lower() == "na":
                continue
            try:
                nb = int((r.get("nb") or "0").strip())
            except ValueError:
                continue
            if args.nb and nb != args.nb:
                continue
            if args.mode:
                mode_val = (r.get("mx_mode") or "").strip()
                if mode_val != args.mode:
                    continue
            out.append(r)
        return out

    by_group = {}
    eps_vals = set()
    report_rows = []
    dp_ref_candidates = []
    for rng, label in zip(ranges, labels):
        rows = pick_rows(rng)
        g_multi = defaultdict(list)
        for r in rows:
            eps = (r.get("source_epsilon") or "").strip()
            counts = parse_counts(r.get(count_col, ""))
            gpt = scale_groups_per_tile(args.nb, args.scale_granularity,
                                        (r.get("mx_mode") or "").strip(),
                                        args.subtile_size)
            mem_b = calc_memory_bytes(counts, args.nb, args.mx_scale_bits,
                                      args.mx_fp32_scale_bits, gpt)
            dp_b = calc_dp64_bytes_for_same_tiles(counts, args.nb)
            g_multi[eps].append(mem_b)
            dp_ref_candidates.append(dp_b)
            eps_vals.add(eps)
            report_rows.append({
                "group": label,
                "epsilon": eps,
                "nb": args.nb,
                "count_col": count_col,
                "memory_bytes": mem_b,
                "memory_MB": mem_b / 1e6,
                "memory_GB": mem_b / 1e9,
                "dp64_ref_bytes": dp_b,
                "dp64_ref_MB": dp_b / 1e6,
                "dp64_ref_GB": dp_b / 1e9,
                "reduction_vs_dp64_pct": 100.0 * (1.0 - (mem_b / dp_b)) if dp_b > 0 else float("nan"),
                "scale_groups_per_tile": gpt,
            })
        by_group[label] = {
            eps: (sum(vals) / len(vals)) for eps, vals in g_multi.items()
        }

    eps_order = sorted(eps_vals, key=eps_sort_key)
    x = list(range(len(eps_order)))
    ngrp = len(labels)
    group_width = 0.82
    bar_w = group_width / max(1, ngrp)

    fig, ax = plt.subplots(figsize=(12, 5.2))
    cmap = plt.get_cmap("tab10")
    for i, label in enumerate(labels):
        vals_gb = [(by_group.get(label, {}).get(eps, float("nan")) / 1e9) for eps in eps_order]
        xpos = [xi - group_width / 2 + (i + 0.5) * bar_w for xi in x]
        bars = ax.bar(xpos, vals_gb, width=bar_w * 0.95, color=cmap(i % 10), label=label)
        y_pad = max(0.001, 0.005 * max(v for v in vals_gb if v == v)) if any(v == v for v in vals_gb) else 0.001
        for rect, v in zip(bars, vals_gb):
            if v == v:
                ax.text(rect.get_x() + rect.get_width() / 2,
                        rect.get_height() + y_pad,
                        f"{v:.3f}", ha="center", va="bottom", fontsize=7, rotation=90)

    if dp_ref_candidates:
        dp_ref_gb = (sum(dp_ref_candidates) / len(dp_ref_candidates)) / 1e9
        ax.axhline(dp_ref_gb, color="black", linestyle="--", linewidth=1.2,
                   label=f"DP64 reference: {dp_ref_gb:.3f} GB")

    ax.set_xticks(x)
    ax.set_xticklabels(eps_order)
    ax.set_xlabel("source_epsilon")
    ax.set_ylabel("mixed memory (GB)")
    mode_txt = f", mode={args.mode}" if args.mode else ""
    ax.set_title(f"Memory vs epsilon (NB={args.nb}, {count_col}{mode_txt}, scale={args.scale_granularity})")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)

    if args.report_out:
        rep = Path(args.report_out)
    else:
        rep = out.with_suffix(".tsv")
    rep.parent.mkdir(parents=True, exist_ok=True)
    with rep.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "group", "epsilon", "nb", "count_col",
                "memory_bytes", "memory_MB", "memory_GB",
                "dp64_ref_bytes", "dp64_ref_MB", "dp64_ref_GB",
                "reduction_vs_dp64_pct",
                "scale_groups_per_tile",
            ],
            delimiter='\t'
        )
        w.writeheader()
        for r in sorted(report_rows, key=lambda z: (z["group"], eps_sort_key(z["epsilon"]))):
            w.writerow(r)

    print(str(out))
    print(str(rep))


if __name__ == "__main__":
    main()
