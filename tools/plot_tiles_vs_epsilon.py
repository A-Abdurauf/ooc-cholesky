#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path
import os

import matplotlib.pyplot as plt


def parse_counts(s: str):
    raw = {}
    if not s:
        return raw

    def canon_fmt(fmt: str):
        v = (fmt or "").strip().lower()
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
        return alias.get(v, v)

    for part in s.split(','):
        part = part.strip()
        if not part or '=' not in part:
            continue
        k, v = part.split('=', 1)
        k = canon_fmt(k)
        try:
            raw[k] = raw.get(k, 0) + int(v.strip())
        except ValueError:
            continue

    out = {}
    for k in ["fp64", "fp32", "mx_fp32", "fp16", "mx_fp16", "bf16", "unknown"]:
        if raw.get(k, 0) > 0:
            out[k] = raw[k]

    # Collapse all low formats into one LOW bar with stacked subcomponents.
    low_mxfp8 = sum(raw.get(k, 0) for k in ["mx_e5m2", "mx_e4m3", "fp8_e5m2", "fp8_e4m3"])
    low_mxfp6 = sum(raw.get(k, 0) for k in ["e3m2", "e2m3"])
    low_mxfp4 = raw.get("e2m1", 0)
    low = low_mxfp8 + low_mxfp6 + low_mxfp4
    if low > 0:
        out["low"] = low
        out["low_mxfp8"] = low_mxfp8
        out["low_mxfp6"] = low_mxfp6
        out["low_mxfp4"] = low_mxfp4
    return out


def is_plot_bar_key(fmt: str) -> bool:
    # Only these keys should become standalone grouped bars.
    return fmt in {"fp64", "fp32", "mx_fp32", "fp16", "mx_fp16", "bf16", "low", "unknown"}


def eps_sort_key(eps: str):
    try:
        return float(eps)
    except Exception:
        return float('inf')


def main():
    ap = argparse.ArgumentParser(description="Plot tile-counts vs epsilon from summary TSV")
    ap.add_argument("--summary", required=True, help="Path to summary TSV")
    ap.add_argument("--out", required=True, help="Output PNG path")
    ap.add_argument("--nb", type=int, default=2048, help="Filter NB (default: 2048)")
    ap.add_argument("--use-full", action="store_true", help="Use tile_counts_full instead of tile_counts_lower")
    ap.add_argument("--mode", default="", help="Optional mx_mode filter (e.g. tile, subtile_128)")
    ap.add_argument("--compare-ranges", default="", help='Optional line ranges from summary, e.g. "4-19,20-27,52-59"')
    ap.add_argument("--compare-labels", default="auto,legacy", help='Labels for compare ranges, e.g. "auto,legacy,ieee"')
    ap.add_argument("--out-rel", default="", help="Optional output path for relative-error comparison plot")
    args = ap.parse_args()

    summary = Path(args.summary)
    rows = list(csv.DictReader(summary.open(), delimiter='\t'))

    def filter_rows(base_rows):
        filtered = []
        for r in base_rows:
            eps = (r.get("source_epsilon") or "").strip()
            if not eps or eps.lower() == "na":
                continue
            if args.mode:
                mode_val = (r.get("mx_mode") or "").strip()
                if mode_val != args.mode:
                    continue
            try:
                nb = int((r.get("nb") or "0").strip())
            except ValueError:
                continue
            if args.nb and nb != args.nb:
                continue
            filtered.append(r)
        return filtered

    # --- Optional compare mode: two explicit line ranges ---
    if args.compare_ranges:
        raw_lines = summary.read_text().splitlines()
        if not raw_lines:
            raise SystemExit("Summary is empty.")
        header = raw_lines[0]

        ranges = [s.strip() for s in args.compare_ranges.split(',') if s.strip()]
        if len(ranges) < 2:
            raise SystemExit("--compare-ranges must contain at least two ranges, e.g. 4-19,20-27")
        labels = [s.strip() for s in args.compare_labels.split(',') if s.strip()]
        if len(labels) != len(ranges):
            labels = [f"group{i+1}" for i in range(len(ranges))]

        def parse_range(rng):
            if '-' not in rng:
                raise ValueError
            a, b = rng.split('-', 1)
            return int(a), int(b)

        def average_counts(dicts):
            if not dicts:
                return {}
            acc = defaultdict(float)
            for d in dicts:
                for k, v in d.items():
                    acc[k] += float(v)
            n = float(len(dicts))
            return {k: int(round(v / n)) for k, v in acc.items()}

        groups = []
        for rng, label in zip(ranges, labels):
            selected = []
            # Allow merging multiple ranges into one group: e.g. 44-51+60-67
            sub_ranges = [x.strip() for x in rng.split('+') if x.strip()]
            for sub_rng in sub_ranges:
                a, b = parse_range(sub_rng)
                lo, hi = min(a, b), max(a, b)
                for ln in range(lo, hi + 1):
                    idx = ln - 1
                    if 0 <= idx < len(raw_lines):
                        line = raw_lines[idx].strip()
                        if line:
                            selected.append(line)
            group_rows = list(csv.DictReader([header] + selected, delimiter='\t'))
            group_rows = filter_rows(group_rows)
            groups.append((label, group_rows))

        # Build grouped-bar counts figure (one subplot per group).
        count_col = "tile_counts_full" if args.use_full else "tile_counts_lower"

        fmts = set()
        eps_vals = set()
        by_group = {}
        for label, grows in groups:
            eps_counts = defaultdict(list)
            for r in grows:
                eps = (r.get("source_epsilon") or "").strip()
                counts = parse_counts(r.get(count_col, ""))
                eps_counts[eps].append(counts)
                fmts.update([k for k in counts.keys() if is_plot_bar_key(k)])
                eps_vals.add(eps)
            gdata = {eps: average_counts(cdicts) for eps, cdicts in eps_counts.items()}
            by_group[label] = gdata

        eps_order = sorted(eps_vals, key=eps_sort_key)
        known = ["fp64", "fp32", "mx_fp32", "fp16", "mx_fp16", "bf16", "low", "unknown"]
        fmts_order = [f for f in known if f in fmts] + sorted([f for f in fmts if f not in known])

        n_groups = len(groups)
        fig, axes = plt.subplots(nrows=n_groups, ncols=1,
                                 figsize=(12, max(4.8, 3.2 * n_groups)),
                                 sharex=True)
        if n_groups == 1:
            axes = [axes]
        cmap = plt.get_cmap("tab20")
        color_map = {fmt: cmap(i % 20) for i, fmt in enumerate(fmts_order)}

        display_name = {
            "fp64": "FP64",
            "fp32": "FP32",
            "mx_fp32": "MXFP32",
            "fp16": "FP16",
            "mx_fp16": "MXFP16",
            "bf16": "BF16",
            "low": "LOW",
            "unknown": "UNKNOWN",
        }
        low_comp = [
            ("low_mxfp8", "MXFP8", "#8da0cb"),
            ("low_mxfp6", "MXFP6", "#66c2a5"),
            ("low_mxfp4", "MXFP4", "#fc8d62"),
        ]

        for ax, (label, _) in zip(axes, groups):
            gdata = by_group[label]
            x = list(range(len(eps_order)))
            legacy_label_mode = "legacy" in (label or "").lower()
            active_formats = [fmt for fmt in fmts_order if any(gdata.get(eps, {}).get(fmt, 0) > 0 for eps in eps_order)]
            if active_formats:
                nfmt = len(active_formats)
                group_width = 0.86
                bar_width = group_width / nfmt
                for j, fmt in enumerate(active_formats):
                    xpos = [xi - group_width / 2 + (j + 0.5) * bar_width for xi in x]
                    if fmt == "low":
                        bottoms = [0] * len(eps_order)
                        for ci, (ckey, clabel, ccolor) in enumerate(low_comp):
                            vals = [gdata.get(eps, {}).get(ckey, 0) for eps in eps_order]
                            bars = ax.bar(xpos, vals, bottom=bottoms, color=ccolor,
                                          label=f"LOW:{clabel}", width=bar_width * 0.95)
                            for idx, (rect, v) in enumerate(zip(bars, vals)):
                                if v > 0:
                                    if legacy_label_mode:
                                        y = bottoms[idx] + v + 0.15
                                        va = 'bottom'
                                    else:
                                        y = bottoms[idx] + v / 2.0
                                        va = 'center'
                                    ax.text(rect.get_x() + rect.get_width() / 2,
                                            y,
                                            str(v),
                                            ha='center', va=va, fontsize=7)
                                bottoms[idx] += v
                        # Add LOW total on top of the stacked bar (requested for auto mode).
                        for xi, total_v in enumerate(bottoms):
                            if total_v > 0:
                                ax.text(xpos[xi],
                                        total_v + 0.25,
                                        str(total_v),
                                        ha='center', va='bottom', fontsize=7)
                    else:
                        vals = [gdata.get(eps, {}).get(fmt, 0) for eps in eps_order]
                        bars = ax.bar(xpos, vals, color=color_map[fmt], label=display_name.get(fmt, fmt), width=bar_width * 0.95)
                        for rect, v in zip(bars, vals):
                            if v > 0:
                                ax.text(rect.get_x() + rect.get_width() / 2,
                                        rect.get_height() + 0.25,
                                        str(v),
                                        ha='center', va='bottom', fontsize=7, rotation=90)
            ax.set_title(label)
            ax.set_ylabel("# tiles")
            ax.grid(axis="y", alpha=0.25)

        axes[-1].set_xticks(list(range(len(eps_order))))
        axes[-1].set_xticklabels(eps_order)
        axes[-1].set_xlabel("source_epsilon")
        # Build a unified legend across all axes so no format label is dropped,
        # while preserving the canonical format order.
        legend_map = {}
        for ax in axes:
            hlist, llist = ax.get_legend_handles_labels()
            for h, l in zip(hlist, llist):
                if l and l not in legend_map:
                    legend_map[l] = h
        label_order = [display_name.get(fmt, fmt) for fmt in fmts_order if fmt != "low"]
        if "low" in fmts_order:
            label_order += ["LOW:MXFP8", "LOW:MXFP6", "LOW:MXFP4"]
        labels_for_legend = [lbl for lbl in label_order if lbl in legend_map]
        handles = [legend_map[lbl] for lbl in labels_for_legend]
        fig.legend(handles, labels_for_legend, loc="lower center", bbox_to_anchor=(0.5, -0.01), ncol=min(8, len(labels_for_legend)), frameon=False)
        mode_title = f", mode={args.mode}" if args.mode else ""
        fig.suptitle(f"20k Matrix Tile-format counts comparison by epsilon ({count_col}, NB={args.nb}{mode_title})", y=0.995)
        fig.tight_layout(rect=[0, 0.08, 1, 0.95])
        out = Path(args.out)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=180)

        # Relative error comparison plot (separate file)
        rel_out = Path(args.out_rel) if args.out_rel else out.with_name(out.stem + "_rel_error" + out.suffix)
        fig2, ax2 = plt.subplots(figsize=(10, 4.8))
        for label, grows in groups:
            vals = defaultdict(list)
            for r in grows:
                eps = (r.get("source_epsilon") or "").strip()
                try:
                    v = float((r.get("rel_factor_error") or "nan"))
                    if v == v:
                        vals[eps].append(v)
                except Exception:
                    continue
            y = []
            for eps in eps_order:
                vv = vals.get(eps, [])
                y.append((sum(vv) / len(vv)) if vv else float("nan"))
            ax2.plot(eps_order, y, marker='o', linewidth=2, label=label)
        ax2.set_yscale('log')
        ax2.set_xlabel("source_epsilon")
        ax2.set_ylabel("relative error")
        ax2.set_title(f"Relative error comparison (NB={args.nb}{mode_title})")
        ax2.grid(True, which='both', alpha=0.25)
        ax2.legend(loc='best')
        fig2.tight_layout()
        rel_out.parent.mkdir(parents=True, exist_ok=True)
        fig2.savefig(rel_out, dpi=180)

        print(str(out))
        print(str(rel_out))
        return

    # Keep rows with concrete epsilon and desired NB.
    filtered = filter_rows(rows)

    if not filtered:
        raise SystemExit("No matching rows found in summary for plotting.")

    count_col = "tile_counts_full" if args.use_full else "tile_counts_lower"

    # mode -> eps -> {fmt: count}
    data = defaultdict(lambda: defaultdict(dict))
    fmts = set()
    eps_vals = set()
    modes = set()

    for r in filtered:
        mode = (r.get("mx_mode") or "").strip() or "unknown"
        eps = (r.get("source_epsilon") or "").strip()
        counts = parse_counts(r.get(count_col, ""))
        data[mode][eps] = counts
        fmts.update([k for k in counts.keys() if is_plot_bar_key(k)])
        eps_vals.add(eps)
        modes.add(mode)

    # Stable display order for known formats, then alphabetical.
    known = ["fp64", "fp32", "mx_fp32", "fp16", "mx_fp16", "bf16", "low", "unknown"]
    fmts_order = [f for f in known if f in fmts] + sorted([f for f in fmts if f not in known])

    display_name = {
        "fp64": "FP64",
        "fp32": "FP32",
        "mx_fp32": "MXFP32",
        "fp16": "FP16",
        "mx_fp16": "MXFP16",
        "bf16": "BF16",
        "low": "LOW",
        "unknown": "UNKNOWN",
    }
    low_comp = [
        ("low_mxfp8", "MXFP8", "#8da0cb"),
        ("low_mxfp6", "MXFP6", "#66c2a5"),
        ("low_mxfp4", "MXFP4", "#fc8d62"),
    ]

    eps_order = sorted(eps_vals, key=eps_sort_key)
    mode_order = sorted(modes)

    nrows = len(mode_order)
    fig, axes = plt.subplots(nrows=nrows, ncols=1, figsize=(11, 4.2 * nrows), sharex=True)
    if nrows == 1:
        axes = [axes]

    cmap = plt.get_cmap("tab20")
    color_map = {fmt: cmap(i % 20) for i, fmt in enumerate(fmts_order)}

    for ax, mode in zip(axes, mode_order):
        x = list(range(len(eps_order)))

        # Grouped bars: each bar is one format at a given epsilon.
        active_formats = []
        for fmt in fmts_order:
            if any(data[mode].get(eps, {}).get(fmt, 0) > 0 for eps in eps_order):
                active_formats.append(fmt)

        if not active_formats:
            ax.set_title(f"mx_mode={mode}")
            ax.set_ylabel("# tiles")
            ax.grid(axis="y", alpha=0.25)
            continue

        nfmt = len(active_formats)
        group_width = 0.86
        bar_width = group_width / nfmt

        for j, fmt in enumerate(active_formats):
            xpos = [xi - group_width / 2 + (j + 0.5) * bar_width for xi in x]
            if fmt == "low":
                bottoms = [0] * len(eps_order)
                for ckey, clabel, ccolor in low_comp:
                    vals = [data[mode].get(eps, {}).get(ckey, 0) for eps in eps_order]
                    bars = ax.bar(xpos, vals, bottom=bottoms, color=ccolor,
                                  label=f"LOW:{clabel}", width=bar_width * 0.95)
                    for idx, (rect, v) in enumerate(zip(bars, vals)):
                        if v > 0:
                            y = bottoms[idx] + v / 2.0
                            ax.text(rect.get_x() + rect.get_width() / 2,
                                    y,
                                    str(v),
                                    ha='center', va='center', fontsize=7)
                        bottoms[idx] += v
                # Add LOW total on top of the stacked bar.
                for xi, total_v in enumerate(bottoms):
                    if total_v > 0:
                        ax.text(xpos[xi],
                                total_v + 0.25,
                                str(total_v),
                                ha='center', va='bottom', fontsize=7)
            else:
                vals = [data[mode].get(eps, {}).get(fmt, 0) for eps in eps_order]
                bars = ax.bar(xpos, vals, color=color_map[fmt], label=display_name.get(fmt, fmt), width=bar_width * 0.95)
                for rect, v in zip(bars, vals):
                    if v > 0:
                        ax.text(rect.get_x() + rect.get_width() / 2,
                                rect.get_height() + 0.25,
                                str(v),
                                ha='center', va='bottom', fontsize=7, rotation=90)

        ax.set_title(f"mx_mode={mode}")
        ax.set_ylabel("# tiles")
        ax.grid(axis="y", alpha=0.25)

    axes[-1].set_xticks(list(range(len(eps_order))))
    axes[-1].set_xticklabels(eps_order)
    axes[-1].set_xlabel("source_epsilon")

    # Single legend at the bottom.
    # Build a unified legend across all axes so no format label is dropped,
    # while preserving the canonical format order.
    legend_map = {}
    for ax in axes:
        hlist, llist = ax.get_legend_handles_labels()
        for h, l in zip(hlist, llist):
            if l and l not in legend_map:
                legend_map[l] = h
    label_order = [display_name.get(fmt, fmt) for fmt in fmts_order if fmt != "low"]
    if "low" in fmts_order:
        label_order += ["LOW:MXFP8", "LOW:MXFP6", "LOW:MXFP4"]
    labels = [lbl for lbl in label_order if lbl in legend_map]
    handles = [legend_map[lbl] for lbl in labels]
    fig.legend(handles, labels,
               loc="lower center",
               bbox_to_anchor=(0.5, -0.01),
               ncol=min(8, len(labels)),
               frameon=False)
    mode_title = f", mode={args.mode}" if args.mode else ""
    fig.suptitle(f"20k Matrix Tile-format counts vs epsilon ({count_col}, NB={args.nb}{mode_title})", y=0.995)
    fig.tight_layout(rect=[0, 0.08, 1, 0.95])

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=180)
    print(str(out))


if __name__ == "__main__":
    main()
