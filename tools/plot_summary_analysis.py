#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(description="Plot analysis from summary_rel_error.txt")
    p.add_argument("--summary", default="/home/abduraa/MX_project/logs/mx_ooc_data/summary_rel_error.txt",
                   help="Path to summary_rel_error.txt")
    p.add_argument("--out-dir", default="/home/abduraa/MX_project/logs/mx_ooc_data/plots/analysis",
                   help="Output directory for plots")
    return p.parse_args()


def extract_size(bin_path: str) -> int | None:
    m = re.search(r"(\d+)", Path(bin_path).name)
    return int(m.group(1)) if m else None


def fmt_size(value) -> str:
    try:
        if pd.isna(value):
            return "unknown"
        return str(int(value))
    except Exception:
        return str(value)


def clean_df(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = [c.strip() for c in df.columns]
    # normalize column names if spacing created duplicates
    rename_map = {}
    for c in list(df.columns):
        if c in {"fp64", "fp32", "fp16", "bf16", "mx_fp16", "low"}:
            continue
        if c.endswith("fp64"):
            rename_map[c] = "fp64"
        elif c.endswith("fp32"):
            rename_map[c] = "fp32"
        elif c.endswith("fp16"):
            rename_map[c] = "fp16"
    if rename_map:
        df = df.rename(columns=rename_map)

    df["data_size"] = df["bin"].apply(extract_size)
    # coerce numerics
    for col in ["nb", "n", "cores", "rel_factor_error", "total_tiles",
                "fp64", "fp32", "fp16", "bf16", "mx_fp16", "low", "data_size"]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    if "fp16_bucket" not in df.columns:
        df["fp16_bucket"] = ""
    if "fp32_bucket" not in df.columns:
        df["fp32_bucket"] = ""
    df.loc[(df["fp16_bucket"] == "") & (df["fp16"].fillna(0) > 0) & (df["mx_fp16"].fillna(0) <= 0), "fp16_bucket"] = "fp16"
    df.loc[(df["fp16_bucket"] == "") & (df["mx_fp16"].fillna(0) > 0) & (df["fp16"].fillna(0) <= 0), "fp16_bucket"] = "mx_fp16"
    df.loc[(df["fp32_bucket"] == "") & (df["fp32"].fillna(0) > 0), "fp32_bucket"] = "fp32"
    df = df.dropna(subset=["data_size", "nb", "mx_mode", "format", "rel_factor_error"])
    return df


def read_summary_table(summary_path: Path) -> pd.DataFrame:
    lines = [ln.strip() for ln in summary_path.read_text().splitlines() if ln.strip()]
    if not lines:
        return pd.DataFrame()

    header = lines[0].split()
    # If later rows include bucket columns not in the header, add them.
    if "fp32_bucket" not in header or "fp16_bucket" not in header:
        has_bucket_cols = False
        for line in lines[1:]:
            if line.startswith("-"):
                continue
            parts = line.split()
            if len(parts) >= len(header) + 2:
                has_bucket_cols = True
                break
        if has_bucket_cols:
            if "fp32_bucket" not in header:
                header.append("fp32_bucket")
            if "fp16_bucket" not in header:
                header.append("fp16_bucket")
    rows = []
    n_idx = header.index("n") if "n" in header else None

    for line in lines[1:]:
        if line.startswith("-"):
            continue
        parts = line.split()
        if n_idx is not None and len(parts) == len(header) - 1:
            parts = parts[:n_idx] + [""] + parts[n_idx:]
        if len(parts) < len(header):
            parts += [""] * (len(header) - len(parts))
        elif len(parts) > len(header):
            parts = parts[:len(header)]
        rows.append(parts)

    return pd.DataFrame(rows, columns=header)


def plot_figure1(df: pd.DataFrame, out_dir: Path):
    # Average tile counts per nb/data_size (ignore mx_mode), include zero-count FP16/low cases
    df_use = df.copy()
    group_cols = ["data_size", "nb"]
    agg = df_use.groupby(group_cols)[["fp64", "fp32", "fp16", "low"]].mean().reset_index()

    nb_order = [32, 64, 128, 256]
    categories = ["fp64", "fp32", "fp16", "low"]
    colors = {"fp64": "#d62728", "fp32": "#ff7f0e", "fp16": "#1f77b4", "low": "#2ca02c"}

    for data_size, sub in agg.groupby(["data_size"]):
        fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharey=False)
        axes = axes.flatten()
        row_max = [0.0, 0.0]

        for idx, nb in enumerate(nb_order):
            ax = axes[idx]
            row = sub[sub["nb"] == nb]
            if row.empty:
                ax.set_visible(False)
                continue
            vals = row[categories].iloc[0].fillna(0)
            x = range(len(categories))
            ax.bar(x, vals, color=[colors[c] for c in categories])
            ax.set_xticks(list(x), categories)
            ax.set_title(f"NB={nb}")
            ax.set_ylabel("Avg tile count")
            row_max[idx // 2] = max(row_max[idx // 2], pd.Series(vals).max(skipna=True))

        fig.suptitle(f"Figure 1: Tile counts (size={fmt_size(data_size)})")
        for r in range(2):
            max_y = row_max[r]
            if max_y > 0:
                for c in range(2):
                    ax = axes[r * 2 + c]
                    if ax.get_visible():
                        ax.set_ylim(0, max_y * 1.2)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        out_path = out_dir / f"fig1_counts_size{fmt_size(data_size)}.pdf"
        fig.savefig(out_path, dpi=200)
        plt.close(fig)


def plot_figure2(df: pd.DataFrame, out_dir: Path):
    # Error by format for each data_size + nb; compare tile vs block side-by-side
    df_use = df.copy()
    fmt_map = {
        "fp8e4m3": "fp8_e4m3",
        "fp8_e4m3": "fp8_e4m3",
        "fp8e5m2": "fp8_e5m2",
        "fp8_e5m2": "fp8_e5m2",
        "e4m3": "e4m3",
        "e5m2": "e5m2",
        "e3m2": "e3m2",
        "e2m3": "e2m3",
        "e2m1": "e2m1",
    }
    df_use["format_norm"] = df_use["format"].str.lower().map(fmt_map)
    df_use = df_use[
        (df_use["fp32"].fillna(0) > 0)
        & (df_use["fp16"].fillna(0) > 0)
        & (df_use["low"].fillna(0) > 0)
        & (df_use["mx_fp16"].fillna(0) <= 0)
        & (df_use["format"] != "mx_fp16")
        & (df_use["format_norm"].notna())
    ].copy()
    nb_order = [32, 64, 128, 256]
    format_order = [
        "e4m3", "e5m2", "e3m2", "e2m3", "e2m1",
        "fp8_e4m3", "fp8_e5m2",
    ]

    for data_size, sub_all in df_use.groupby(["data_size"]):
        fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharey=False)
        axes = axes.flatten()
        fig_handles = None
        fig_labels = None
        row_max = [0.0, 0.0]

        for idx, nb in enumerate(nb_order):
            ax = axes[idx]
            sub = sub_all[sub_all["nb"] == nb]
            if sub.empty:
                ax.set_visible(False)
                continue

            tile = sub[sub["mx_mode"] == "tile"]
            block = sub[sub["mx_mode"] == "block"]

            tile2 = tile.drop_duplicates(subset=["format_norm"]).set_index("format_norm")["rel_factor_error"]
            block2 = block.drop_duplicates(subset=["format_norm"]).set_index("format_norm")["rel_factor_error"]

            formats = format_order
            if not formats:
                ax.set_visible(False)
                continue

            x = pd.Series(range(len(formats))).to_numpy()
            width = 0.4
            bars1 = ax.bar(x - width / 2, [tile2.get(f, float("nan")) for f in formats],
                           width=width, label="tile", color="#4c78a8")
            bars2 = ax.bar(x + width / 2, [block2.get(f, float("nan")) for f in formats],
                           width=width, label="block", color="#f58518")

            ax.set_title(f"NB={nb}")
            ax.set_xlabel("Format")
            ax.set_ylabel("rel_factor_error")
            ax.set_xticks(x, formats)
            ax.tick_params(axis="x", rotation=45)
            row_max[idx // 2] = max(row_max[idx // 2], pd.Series([b.get_height() for b in bars1] + [b.get_height() for b in bars2]).max(skipna=True))

            if fig_handles is None:
                fig_handles, fig_labels = ax.get_legend_handles_labels()

            for bars in (bars1, bars2):
                for b in bars:
                    if pd.isna(b.get_height()):
                        continue
                    ax.annotate(f"{b.get_height():.2e}",
                                xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                                xytext=(0, 3), textcoords="offset points",
                                ha="center", va="bottom", fontsize=7, rotation=45)

        fig.suptitle(f"Figure 2: Error by format (size={fmt_size(data_size)})")
        if fig_handles:
            fig.legend(fig_handles, fig_labels, loc="lower center", ncol=2,
                       bbox_to_anchor=(0.5, -0.02), frameon=False)
        for r in range(2):
            max_y = row_max[r]
            if max_y > 0:
                for c in range(2):
                    ax = axes[r * 2 + c]
                    if ax.get_visible():
                        ax.set_ylim(0, max_y * 1.2)
        fig.tight_layout(rect=[0, 0.08, 1, 0.95])
        out_path = out_dir / f"fig2_error_size{fmt_size(data_size)}_tile_vs_block.pdf"
        fig.savefig(out_path, dpi=200)
        plt.close(fig)


def plot_figure3(df: pd.DataFrame, out_dir: Path):
    # Compare FP32+FP16 runs vs FP32+MX_FP16 runs.
    df = df.copy()
    df["fp16_present"] = df["fp16"].fillna(0) > 0
    df["mx_fp16_present"] = df["mx_fp16"].fillna(0) > 0
    df["fp32_present"] = df["fp32"].fillna(0) > 0

    nb_order = [32, 64, 128, 256]
    format_order = ["e4m3", "e5m2", "e3m2", "e2m3", "e2m1"]

    for data_size, sub_all in df.groupby(["data_size"]):
        fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharey=False)
        axes = axes.flatten()
        fig_handles = None
        fig_labels = None
        row_max = [0.0, 0.0]

        for idx, nb in enumerate(nb_order):
            ax = axes[idx]
            sub = sub_all[sub_all["nb"] == nb]

            values = {fmt: (float("nan"), float("nan"), float("nan")) for fmt in format_order}
            for fmt, g in sub.groupby("format"):
                if fmt not in values:
                    continue
                tile_g = g[g["mx_mode"] == "tile"]
                block_g = g[g["mx_mode"] == "block"]

                tile_fp16_rows = tile_g[
                    tile_g["fp32_present"] & tile_g["fp16_present"] & ~tile_g["mx_fp16_present"]
                ]
                tile_mx_rows = tile_g[
                    tile_g["fp32_present"] & tile_g["mx_fp16_present"] & ~tile_g["fp16_present"]
                ]

                if tile_fp16_rows.empty or tile_mx_rows.empty:
                    continue
                if tile_fp16_rows["fp32"].iloc[0] != tile_mx_rows["fp32"].iloc[0]:
                    continue

                tile_fp16 = tile_fp16_rows["rel_factor_error"].iloc[0]
                tile_mx = tile_mx_rows["rel_factor_error"].iloc[0]

                block_mx = block_g[
                    block_g["fp32_present"] & block_g["mx_fp16_present"] & ~block_g["fp16_present"]
                ]["rel_factor_error"]
                block_mx = block_mx.iloc[0] if not block_mx.empty else float("nan")

                values[fmt] = (tile_fp16, tile_mx, block_mx)

            fmts = format_order
            tile_fp16_vals = [values[f][0] for f in fmts]
            tile_mx_vals = [values[f][1] for f in fmts]
            block_mx_vals = [values[f][2] for f in fmts]

            if any(pd.notna(v) for v in tile_fp16_vals + tile_mx_vals + block_mx_vals):
                x = pd.Series(range(len(fmts))).to_numpy()
                width = 0.25
                bars1 = ax.bar(x - width, tile_fp16_vals, width=width, label="tile FP32+FP16", color="#1f77b4")
                bars2 = ax.bar(x, tile_mx_vals, width=width, label="tile FP32+MX_FP16", color="#ff7f0e")
                bars3 = ax.bar(x + width, block_mx_vals, width=width, label="block FP32+MX_FP16", color="#54a24b")
                ax.set_xticks(x, fmts)
                ax.tick_params(axis="x", rotation=45)
                row_max[idx // 2] = max(row_max[idx // 2], pd.Series([b.get_height() for b in bars1]
                                              + [b.get_height() for b in bars2]
                                              + [b.get_height() for b in bars3]).max(skipna=True))
                if fig_handles is None:
                    fig_handles, fig_labels = ax.get_legend_handles_labels()
                for bars in (bars1, bars2, bars3):
                    for b in bars:
                        if pd.isna(b.get_height()):
                            continue
                        ax.annotate(f"{b.get_height():.2e}",
                                    xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                                    xytext=(0, 3), textcoords="offset points",
                                    ha="center", va="bottom", fontsize=7, rotation=45)
            else:
                ax.set_xticks(pd.Series(range(len(fmts))).to_numpy(), fmts)

            ax.set_title(f"NB={nb}")
            ax.set_xlabel("Low precision format")
            ax.set_ylabel("rel_factor_error")

        fig.suptitle(f"Figure 3: FP32+FP16 vs MX_FP16 (tile/block) (size={fmt_size(data_size)})")
        if fig_handles:
            fig.legend(fig_handles, fig_labels, loc="lower center", ncol=3,
                       bbox_to_anchor=(0.5, -0.02), frameon=False)
        for r in range(2):
            max_y = row_max[r]
            if max_y > 0:
                for c in range(2):
                    ax = axes[r * 2 + c]
                    if ax.get_visible():
                        ax.set_ylim(0, max_y * 1.2)
        fig.tight_layout(rect=[0, 0.08, 1, 0.95])
        out_path = out_dir / f"fig3_fp32_fp16_vs_mxfp16_size{fmt_size(data_size)}_combined.pdf"
        fig.savefig(out_path, dpi=200)
        plt.close(fig)


def plot_figure4(df: pd.DataFrame, out_dir: Path):
    # Reference FP32+FP16, compare FP32->MX_FP16 (FP16 kept) and FP32+FP16->MX_FP16
    df = df.copy()
    df["fp16_present"] = df["fp16"].fillna(0) > 0
    df["mx_fp16_present"] = df["mx_fp16"].fillna(0) > 0
    df["fp32_zero"] = df["fp32"].fillna(0) <= 0
    df["fp32_present"] = df["fp32"].fillna(0) > 0

    nb_order = [32, 64, 128, 256]
    format_order = ["e4m3", "e5m2", "e3m2", "e2m3", "e2m1"]

    for data_size, sub_all in df.groupby(["data_size"]):
        fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharey=False)
        axes = axes.flatten()
        fig_handles = None
        fig_labels = None
        row_max = [0.0, 0.0]

        for idx, nb in enumerate(nb_order):
            ax = axes[idx]
            sub = sub_all[sub_all["nb"] == nb]

            values = {
                fmt: (float("nan"), float("nan"), float("nan"), float("nan"), float("nan"), float("nan"))
                for fmt in format_order
            }
            for fmt, g in sub.groupby("format"):
                if fmt not in values:
                    continue
                tile_g = g[g["mx_mode"] == "tile"]
                block_g = g[g["mx_mode"] == "block"]

                tile_ref_rows = tile_g[
                    tile_g["fp32_present"] & tile_g["fp16_present"] & ~tile_g["mx_fp16_present"]
                ]
                tile_fp32_mx_rows = tile_g[
                    tile_g["fp32_zero"] & tile_g["fp16_present"] & tile_g["mx_fp16_present"]
                ]
                tile_both_mx_rows = tile_g[
                    tile_g["fp32_zero"] & ~tile_g["fp16_present"] & tile_g["mx_fp16_present"]
                ]

                if tile_ref_rows.empty or tile_fp32_mx_rows.empty or tile_both_mx_rows.empty:
                    continue

                tile_ref = tile_ref_rows["rel_factor_error"].iloc[0]
                tile_fp32_mx = tile_fp32_mx_rows["rel_factor_error"].iloc[0]
                tile_both_mx = tile_both_mx_rows["rel_factor_error"].iloc[0]

                block_ref_rows = block_g[
                    block_g["fp32_present"] & block_g["fp16_present"] & ~block_g["mx_fp16_present"]
                ]
                block_fp32_mx_rows = block_g[
                    block_g["fp32_zero"] & block_g["fp16_present"] & block_g["mx_fp16_present"]
                ]
                block_both_mx_rows = block_g[
                    block_g["fp32_zero"] & ~block_g["fp16_present"] & block_g["mx_fp16_present"]
                ]

                block_ref = block_ref_rows["rel_factor_error"].iloc[0] if not block_ref_rows.empty else float("nan")
                block_fp32_mx = block_fp32_mx_rows["rel_factor_error"].iloc[0] if not block_fp32_mx_rows.empty else float("nan")
                block_both_mx = block_both_mx_rows["rel_factor_error"].iloc[0] if not block_both_mx_rows.empty else float("nan")

                values[fmt] = (tile_ref, tile_fp32_mx, tile_both_mx, block_ref, block_fp32_mx, block_both_mx)

            fmts = format_order
            tile_ref_vals = [values[f][0] for f in fmts]
            tile_fp32_mx_vals = [values[f][1] for f in fmts]
            tile_both_mx_vals = [values[f][2] for f in fmts]
            block_ref_vals = [values[f][3] for f in fmts]
            block_fp32_mx_vals = [values[f][4] for f in fmts]
            block_both_mx_vals = [values[f][5] for f in fmts]

            if any(pd.notna(v) for v in tile_ref_vals + tile_fp32_mx_vals + tile_both_mx_vals
                    + block_ref_vals + block_fp32_mx_vals + block_both_mx_vals):
                x = pd.Series(range(len(fmts))).to_numpy()
                width = 0.12
                bars1 = ax.bar(x - 2.5 * width, tile_ref_vals, width=width, label="tile FP32+FP16", color="#1f77b4")
                bars2 = ax.bar(x - 1.5 * width, tile_fp32_mx_vals, width=width, label="tile FP32->MX_FP16", color="#ff7f0e")
                bars3 = ax.bar(x - 0.5 * width, tile_both_mx_vals, width=width, label="tile FP32+FP16->MX_FP16", color="#2ca02c")
                bars4 = ax.bar(x + 0.5 * width, block_ref_vals, width=width, label="block FP32+FP16", color="#9467bd")
                bars5 = ax.bar(x + 1.5 * width, block_fp32_mx_vals, width=width, label="block FP32->MX_FP16", color="#8c564b")
                bars6 = ax.bar(x + 2.5 * width, block_both_mx_vals, width=width, label="block FP32+FP16->MX_FP16", color="#e45756")
                ax.set_xticks(x, fmts)
                ax.tick_params(axis="x", rotation=45)
                row_max[idx // 2] = max(row_max[idx // 2], pd.Series([b.get_height() for b in bars1]
                                              + [b.get_height() for b in bars2]
                                              + [b.get_height() for b in bars3]
                                              + [b.get_height() for b in bars4]
                                              + [b.get_height() for b in bars5]
                                              + [b.get_height() for b in bars6]).max(skipna=True))
                if fig_handles is None:
                    fig_handles, fig_labels = ax.get_legend_handles_labels()

                for bars in (bars1, bars2, bars3, bars4, bars5, bars6):
                    for b in bars:
                        if pd.isna(b.get_height()):
                            continue
                        ax.annotate(f"{b.get_height():.2e}",
                                    xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                                    xytext=(0, 3), textcoords="offset points",
                                    ha="center", va="bottom", fontsize=7, rotation=45)
            else:
                ax.set_xticks(pd.Series(range(len(fmts))).to_numpy(), fmts)

            ax.set_title(f"NB={nb}")
            ax.set_xlabel("Low precision format")
            ax.set_ylabel("rel_factor_error")

        fig.suptitle(f"Figure 4: FP32+FP16 reference vs MX_FP16 replacements (size={fmt_size(data_size)})")
        if fig_handles:
            fig.legend(fig_handles, fig_labels, loc="lower center", ncol=3,
                       bbox_to_anchor=(0.5, -0.02), frameon=False)
        for r in range(2):
            max_y = row_max[r]
            if max_y > 0:
                for c in range(2):
                    ax = axes[r * 2 + c]
                    if ax.get_visible():
                        ax.set_ylim(0, max_y * 1.25)
        fig.tight_layout(rect=[0, 0.1, 1, 0.95])
        out_path = out_dir / f"fig4_fp32_replaced_mxfp16_size{fmt_size(data_size)}.pdf"
        fig.savefig(out_path, dpi=200)
        plt.close(fig)


def plot_figure5(df: pd.DataFrame, out_dir: Path):
    # Compare plain FP8 vs MX FP8 counterparts (e4m3/e5m2) for size 512/1024
    df_use = df.copy()
    df_use = df_use[df_use["data_size"].isin([512, 2048])]
    fmt_map = {
        "fp8e4m3": "fp8_e4m3",
        "fp8_e4m3": "fp8_e4m3",
        "fp8e5m2": "fp8_e5m2",
        "fp8_e5m2": "fp8_e5m2",
        "e4m3": "mx_e4m3",
        "mx_e4m3": "mx_e4m3",
        "mx_fp8_e4m3": "mx_e4m3",
        "e5m2": "mx_e5m2",
        "mx_e5m2": "mx_e5m2",
        "mx_fp8_e5m2": "mx_e5m2",
    }
    df_use["format_norm"] = df_use["format"].str.lower().map(fmt_map)
    df_use = df_use[df_use["format_norm"].isin(["fp8_e4m3", "fp8_e5m2", "mx_e4m3", "mx_e5m2"])]

    nb_order = [32, 64, 128, 256]
    pairs = [
        ("fp8_e4m3", ["mx_e4m3"], "E4M3"),
        ("fp8_e5m2", ["mx_e5m2"], "E5M2"),
    ]

    for data_size, sub_all in df_use.groupby(["data_size"]):
        fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharey=False)
        axes = axes.flatten()
        fig_handles = None
        fig_labels = None
        row_max = [0.0, 0.0]

        for idx, nb in enumerate(nb_order):
            ax = axes[idx]
            sub = sub_all[sub_all["nb"] == nb]
            if sub.empty:
                ax.set_visible(False)
                continue

            x = pd.Series(range(len(pairs))).to_numpy()
            width = 0.18
            offsets = [-1.5 * width, -0.5 * width, 0.5 * width, 1.5 * width]

            tile_fp8_vals = []
            tile_mx_vals = []
            block_fp8_vals = []
            block_mx_vals = []

            for fp8_fmt, mx_fmts, _ in pairs:
                tile_fp8 = sub[(sub["format_norm"] == fp8_fmt) & (sub["mx_mode"] == "tile")]
                block_fp8 = sub[(sub["format_norm"] == fp8_fmt) & (sub["mx_mode"] == "block")]
                tile_mx = sub[(sub["format_norm"].isin(mx_fmts)) & (sub["mx_mode"] == "tile")]
                block_mx = sub[(sub["format_norm"].isin(mx_fmts)) & (sub["mx_mode"] == "block")]

                tile_fp8_val = tile_fp8.drop_duplicates(subset=["format_norm", "mx_mode"])["rel_factor_error"]
                block_fp8_val = block_fp8.drop_duplicates(subset=["format_norm", "mx_mode"])["rel_factor_error"]
                tile_mx_val = tile_mx.drop_duplicates(subset=["format_norm", "mx_mode"])["rel_factor_error"]
                block_mx_val = block_mx.drop_duplicates(subset=["format_norm", "mx_mode"])["rel_factor_error"]

                tile_fp8_vals.append(tile_fp8_val.iloc[0] if not tile_fp8_val.empty else float("nan"))
                block_fp8_vals.append(block_fp8_val.iloc[0] if not block_fp8_val.empty else float("nan"))
                tile_mx_vals.append(tile_mx_val.iloc[0] if not tile_mx_val.empty else float("nan"))
                block_mx_vals.append(block_mx_val.iloc[0] if not block_mx_val.empty else float("nan"))

            max_local = pd.Series([
                *tile_fp8_vals,
                *tile_mx_vals,
                *block_fp8_vals,
                *block_mx_vals,
            ]).max(skipna=True)
            if pd.isna(max_local) or max_local <= 0:
                max_local = 1.0

            tile_fp8_plot = [v / max_local if pd.notna(v) else float("nan") for v in tile_fp8_vals]
            tile_mx_plot = [v / max_local if pd.notna(v) else float("nan") for v in tile_mx_vals]
            block_fp8_plot = [v / max_local if pd.notna(v) else float("nan") for v in block_fp8_vals]
            block_mx_plot = [v / max_local if pd.notna(v) else float("nan") for v in block_mx_vals]

            bars1 = ax.bar(x + offsets[0], tile_fp8_plot, width=width, label="tile FP8", color="#4c78a8")
            bars2 = ax.bar(x + offsets[1], tile_mx_plot, width=width, label="tile MX", color="#f58518")
            bars3 = ax.bar(x + offsets[2], block_fp8_plot, width=width, label="block FP8", color="#54a24b")
            bars4 = ax.bar(x + offsets[3], block_mx_plot, width=width, label="block MX", color="#e45756")

            ax.set_xticks(x, [p[2] for p in pairs])
            ax.set_title(f"NB={nb}")
            ax.set_xlabel("Format")
            ax.set_ylabel("rel_factor_error (normalized)")

            row_max[idx // 2] = max(row_max[idx // 2], pd.Series([
                *[b.get_height() for b in bars1],
                *[b.get_height() for b in bars2],
                *[b.get_height() for b in bars3],
                *[b.get_height() for b in bars4],
            ]).max(skipna=True))

            if fig_handles is None:
                fig_handles, fig_labels = ax.get_legend_handles_labels()

            for bars, raw_vals in (
                (bars1, tile_fp8_vals),
                (bars2, tile_mx_vals),
                (bars3, block_fp8_vals),
                (bars4, block_mx_vals),
            ):
                for b, raw in zip(bars, raw_vals):
                    if pd.isna(b.get_height()) or pd.isna(raw):
                        continue
                    ax.annotate(f"{raw:.2e}",
                                xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                                xytext=(0, 3), textcoords="offset points",
                                ha="center", va="bottom", fontsize=7, rotation=45)

        fig.suptitle(f"Figure 5: FP8 vs MX FP8 (size={fmt_size(data_size)})")
        if fig_handles:
            fig.legend(fig_handles, fig_labels, loc="lower center", ncol=4,
                       bbox_to_anchor=(0.5, -0.02), frameon=False)
        for r in range(2):
            max_y = row_max[r]
            if max_y > 0:
                for c in range(2):
                    ax = axes[r * 2 + c]
                    if ax.get_visible():
                        ax.set_ylim(0, max_y * 1.25)
        fig.tight_layout(rect=[0, 0.1, 1, 0.95])
        out_path = out_dir / f"fig5_fp8_vs_mx_fp8_size{fmt_size(data_size)}.pdf"
        fig.savefig(out_path, dpi=200)
        plt.close(fig)


def plot_figure6(df: pd.DataFrame, out_dir: Path):
    # Compare plain FP8 vs MX FP6 counterparts (e3m2/e2m3) for size 512/2048
    df_use = df.copy()
    df_use = df_use[df_use["data_size"].isin([512, 2048])]
    fmt_map = {
        "fp8e4m3": "fp8_e4m3",
        "fp8_e4m3": "fp8_e4m3",
        "fp8e5m2": "fp8_e5m2",
        "fp8_e5m2": "fp8_e5m2",
        "e3m2": "mx_e3m2",
        "e2m3": "mx_e2m3",
    }
    df_use["format_norm"] = df_use["format"].str.lower().map(fmt_map)
    df_use = df_use[df_use["format_norm"].isin([
        "fp8_e4m3", "fp8_e5m2", "mx_e3m2", "mx_e2m3"
    ])]

    nb_order = [32, 64, 128, 256]
    pairs = [
        ("fp8_e4m3", ["mx_e3m2"], "FP8 e4m3 vs MX e3m2"),
        ("fp8_e5m2", ["mx_e2m3"], "FP8 e5m2 vs MX e2m3"),
    ]

    for data_size, sub_all in df_use.groupby(["data_size"]):
        fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharey=False)
        axes = axes.flatten()
        fig_handles = None
        fig_labels = None
        row_max = [0.0, 0.0]

        for idx, nb in enumerate(nb_order):
            ax = axes[idx]
            sub = sub_all[sub_all["nb"] == nb]
            if sub.empty:
                ax.set_visible(False)
                continue

            x = pd.Series(range(len(pairs))).to_numpy()
            width = 0.18
            offsets = [-1.5 * width, -0.5 * width, 0.5 * width, 1.5 * width]

            tile_fp8_vals = []
            tile_mx_vals = []
            block_fp8_vals = []
            block_mx_vals = []

            for fp8_fmt, mx_fmts, _ in pairs:
                tile_fp8 = sub[(sub["format_norm"] == fp8_fmt) & (sub["mx_mode"] == "tile")]
                block_fp8 = sub[(sub["format_norm"] == fp8_fmt) & (sub["mx_mode"] == "block")]
                tile_mx = sub[(sub["format_norm"].isin(mx_fmts)) & (sub["mx_mode"] == "tile")]
                block_mx = sub[(sub["format_norm"].isin(mx_fmts)) & (sub["mx_mode"] == "block")]

                tile_fp8_val = tile_fp8.drop_duplicates(subset=["format_norm", "mx_mode"])["rel_factor_error"]
                block_fp8_val = block_fp8.drop_duplicates(subset=["format_norm", "mx_mode"])["rel_factor_error"]
                tile_mx_val = tile_mx.drop_duplicates(subset=["format_norm", "mx_mode"])["rel_factor_error"]
                block_mx_val = block_mx.drop_duplicates(subset=["format_norm", "mx_mode"])["rel_factor_error"]

                tile_fp8_vals.append(tile_fp8_val.iloc[0] if not tile_fp8_val.empty else float("nan"))
                block_fp8_vals.append(block_fp8_val.iloc[0] if not block_fp8_val.empty else float("nan"))
                tile_mx_vals.append(tile_mx_val.iloc[0] if not tile_mx_val.empty else float("nan"))
                block_mx_vals.append(block_mx_val.iloc[0] if not block_mx_val.empty else float("nan"))

            max_local = pd.Series([
                *tile_fp8_vals,
                *tile_mx_vals,
                *block_fp8_vals,
                *block_mx_vals,
            ]).max(skipna=True)
            if pd.isna(max_local) or max_local <= 0:
                max_local = 1.0

            tile_fp8_plot = [v / max_local if pd.notna(v) else float("nan") for v in tile_fp8_vals]
            tile_mx_plot = [v / max_local if pd.notna(v) else float("nan") for v in tile_mx_vals]
            block_fp8_plot = [v / max_local if pd.notna(v) else float("nan") for v in block_fp8_vals]
            block_mx_plot = [v / max_local if pd.notna(v) else float("nan") for v in block_mx_vals]

            bars1 = ax.bar(x + offsets[0], tile_fp8_plot, width=width, label="tile FP8", color="#4c78a8")
            bars2 = ax.bar(x + offsets[1], tile_mx_plot, width=width, label="tile MX FP6", color="#f58518")
            bars3 = ax.bar(x + offsets[2], block_fp8_plot, width=width, label="block FP8", color="#54a24b")
            bars4 = ax.bar(x + offsets[3], block_mx_plot, width=width, label="block MX FP6", color="#e45756")

            ax.set_xticks(x, [p[2] for p in pairs])
            ax.set_title(f"NB={nb}")
            ax.set_xlabel("Format")
            ax.set_ylabel("rel_factor_error (normalized)")

            row_max[idx // 2] = max(row_max[idx // 2], pd.Series([
                *[b.get_height() for b in bars1],
                *[b.get_height() for b in bars2],
                *[b.get_height() for b in bars3],
                *[b.get_height() for b in bars4],
            ]).max(skipna=True))

            if fig_handles is None:
                fig_handles, fig_labels = ax.get_legend_handles_labels()

            for bars, raw_vals in (
                (bars1, tile_fp8_vals),
                (bars2, tile_mx_vals),
                (bars3, block_fp8_vals),
                (bars4, block_mx_vals),
            ):
                for b, raw in zip(bars, raw_vals):
                    if pd.isna(b.get_height()) or pd.isna(raw):
                        continue
                    ax.annotate(f"{raw:.2e}",
                                xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                                xytext=(0, 3), textcoords="offset points",
                                ha="center", va="bottom", fontsize=7, rotation=45)

        fig.suptitle(f"Figure 6: FP8 vs MX FP6 (size={fmt_size(data_size)})")
        if fig_handles:
            fig.legend(fig_handles, fig_labels, loc="lower center", ncol=4,
                       bbox_to_anchor=(0.5, -0.02), frameon=False)
        for r in range(2):
            max_y = row_max[r]
            if max_y > 0:
                for c in range(2):
                    ax = axes[r * 2 + c]
                    if ax.get_visible():
                        ax.set_ylim(0, max_y * 1.25)
        fig.tight_layout(rect=[0, 0.1, 1, 0.95])
        out_path = out_dir / f"fig6_fp8_vs_mx_fp6_size{fmt_size(data_size)}.pdf"
        fig.savefig(out_path, dpi=200)
        plt.close(fig)


def plot_figure7(df: pd.DataFrame, out_dir: Path):
    # Compare plain FP8 vs MX FP4 (e2m1) for size 512/2048
    df_use = df.copy()
    df_use = df_use[df_use["data_size"].isin([512, 2048])]
    fmt_map = {
        "fp8e4m3": "fp8_e4m3",
        "fp8_e4m3": "fp8_e4m3",
        "fp8e5m2": "fp8_e5m2",
        "fp8_e5m2": "fp8_e5m2",
        "e2m1": "mx_e2m1",
    }
    df_use["format_norm"] = df_use["format"].str.lower().map(fmt_map)
    df_use = df_use[df_use["format_norm"].isin([
        "fp8_e4m3", "fp8_e5m2", "mx_e2m1"
    ])]

    nb_order = [32, 64, 128, 256]
    pairs = [
        ("fp8_e4m3", ["mx_e2m1"], "FP8 e4m3 vs MX e2m1"),
        ("fp8_e5m2", ["mx_e2m1"], "FP8 e5m2 vs MX e2m1"),
    ]

    for data_size, sub_all in df_use.groupby(["data_size"]):
        fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharey=False)
        axes = axes.flatten()
        fig_handles = None
        fig_labels = None
        row_max = [0.0, 0.0]

        for idx, nb in enumerate(nb_order):
            ax = axes[idx]
            sub = sub_all[sub_all["nb"] == nb]
            if sub.empty:
                ax.set_visible(False)
                continue

            x = pd.Series(range(len(pairs))).to_numpy()
            width = 0.18
            offsets = [-1.5 * width, -0.5 * width, 0.5 * width, 1.5 * width]

            tile_fp8_vals = []
            tile_mx_vals = []
            block_fp8_vals = []
            block_mx_vals = []

            for fp8_fmt, mx_fmts, _ in pairs:
                tile_fp8 = sub[(sub["format_norm"] == fp8_fmt) & (sub["mx_mode"] == "tile")]
                block_fp8 = sub[(sub["format_norm"] == fp8_fmt) & (sub["mx_mode"] == "block")]
                tile_mx = sub[(sub["format_norm"].isin(mx_fmts)) & (sub["mx_mode"] == "tile")]
                block_mx = sub[(sub["format_norm"].isin(mx_fmts)) & (sub["mx_mode"] == "block")]

                tile_fp8_val = tile_fp8.drop_duplicates(subset=["format_norm", "mx_mode"])["rel_factor_error"]
                block_fp8_val = block_fp8.drop_duplicates(subset=["format_norm", "mx_mode"])["rel_factor_error"]
                tile_mx_val = tile_mx.drop_duplicates(subset=["format_norm", "mx_mode"])["rel_factor_error"]
                block_mx_val = block_mx.drop_duplicates(subset=["format_norm", "mx_mode"])["rel_factor_error"]

                tile_fp8_vals.append(tile_fp8_val.iloc[0] if not tile_fp8_val.empty else float("nan"))
                block_fp8_vals.append(block_fp8_val.iloc[0] if not block_fp8_val.empty else float("nan"))
                tile_mx_vals.append(tile_mx_val.iloc[0] if not tile_mx_val.empty else float("nan"))
                block_mx_vals.append(block_mx_val.iloc[0] if not block_mx_val.empty else float("nan"))

            max_local = pd.Series([
                *tile_fp8_vals,
                *tile_mx_vals,
                *block_fp8_vals,
                *block_mx_vals,
            ]).max(skipna=True)
            if pd.isna(max_local) or max_local <= 0:
                max_local = 1.0

            tile_fp8_plot = [v / max_local if pd.notna(v) else float("nan") for v in tile_fp8_vals]
            tile_mx_plot = [v / max_local if pd.notna(v) else float("nan") for v in tile_mx_vals]
            block_fp8_plot = [v / max_local if pd.notna(v) else float("nan") for v in block_fp8_vals]
            block_mx_plot = [v / max_local if pd.notna(v) else float("nan") for v in block_mx_vals]

            bars1 = ax.bar(x + offsets[0], tile_fp8_plot, width=width, label="tile FP8", color="#4c78a8")
            bars2 = ax.bar(x + offsets[1], tile_mx_plot, width=width, label="tile MX FP4", color="#f58518")
            bars3 = ax.bar(x + offsets[2], block_fp8_plot, width=width, label="block FP8", color="#54a24b")
            bars4 = ax.bar(x + offsets[3], block_mx_plot, width=width, label="block MX FP4", color="#e45756")

            ax.set_xticks(x, [p[2] for p in pairs])
            ax.set_title(f"NB={nb}")
            ax.set_xlabel("Format")
            ax.set_ylabel("rel_factor_error (normalized)")

            row_max[idx // 2] = max(row_max[idx // 2], pd.Series([
                *[b.get_height() for b in bars1],
                *[b.get_height() for b in bars2],
                *[b.get_height() for b in bars3],
                *[b.get_height() for b in bars4],
            ]).max(skipna=True))

            if fig_handles is None:
                fig_handles, fig_labels = ax.get_legend_handles_labels()

            for bars, raw_vals in (
                (bars1, tile_fp8_vals),
                (bars2, tile_mx_vals),
                (bars3, block_fp8_vals),
                (bars4, block_mx_vals),
            ):
                for b, raw in zip(bars, raw_vals):
                    if pd.isna(b.get_height()) or pd.isna(raw):
                        continue
                    ax.annotate(f"{raw:.2e}",
                                xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                                xytext=(0, 3), textcoords="offset points",
                                ha="center", va="bottom", fontsize=7, rotation=45)

        fig.suptitle(f"Figure 7: FP8 vs MX FP4 (size={fmt_size(data_size)})")
        if fig_handles:
            fig.legend(fig_handles, fig_labels, loc="lower center", ncol=4,
                       bbox_to_anchor=(0.5, -0.02), frameon=False)
        for r in range(2):
            max_y = row_max[r]
            if max_y > 0:
                for c in range(2):
                    ax = axes[r * 2 + c]
                    if ax.get_visible():
                        ax.set_ylim(0, max_y * 1.25)
        fig.tight_layout(rect=[0, 0.1, 1, 0.95])
        out_path = out_dir / f"fig7_fp8_vs_mx_fp4_size{fmt_size(data_size)}.pdf"
        fig.savefig(out_path, dpi=200)
        plt.close(fig)


def plot_figure8(df: pd.DataFrame, out_dir: Path):
    # Compare block FP16 vs block MX-FP8 replacement for sizes 2048/4096
    if "fp16_bucket" not in df.columns:
        return

    df_use = df.copy()
    df_use = df_use[df_use["data_size"].isin([2048, 4096])]
    df_use = df_use[df_use["mx_mode"] == "block"]

    df_use["fp16_bucket_norm"] = df_use["fp16_bucket"].astype(str).str.lower()
    if "fp32_bucket" in df_use.columns:
        df_use["fp32_bucket_norm"] = df_use["fp32_bucket"].astype(str).str.lower()

    nb_order = [32, 64, 128, 256]
    format_order = ["e4m3", "e5m2", "e3m2", "e2m3", "e2m1"]
    mx_fp8_buckets = {"mx_e4m3", "mx_e5m2"}

    for data_size, sub_all in df_use.groupby(["data_size"]):
        fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharey=False)
        axes = axes.flatten()
        fig_handles = None
        fig_labels = None
        row_max = [0.0, 0.0]

        for idx, nb in enumerate(nb_order):
            ax = axes[idx]
            sub = sub_all[sub_all["nb"] == nb]
            if sub.empty:
                ax.set_visible(False)
                continue

            fp16_vals = []
            mx_vals = []
            for fmt in format_order:
                g = sub[sub["format"] == fmt]

                fp16_rows = g[g["fp16_bucket_norm"] == "fp16"]
                mx_rows = g[g["fp16_bucket_norm"].isin(mx_fp8_buckets)]

                fp16_val = fp16_rows.sort_values("fp16_bucket_norm").drop_duplicates(
                    subset=["format", "mx_mode", "fp16_bucket_norm"]
                )
                mx_val = mx_rows.sort_values("fp16_bucket_norm").drop_duplicates(
                    subset=["format", "mx_mode", "fp16_bucket_norm"]
                )

                fp16_vals.append(fp16_val["rel_factor_error"].iloc[0] if not fp16_val.empty else float("nan"))
                mx_vals.append(mx_val["rel_factor_error"].iloc[0] if not mx_val.empty else float("nan"))

            if any(pd.notna(v) for v in fp16_vals + mx_vals):
                x = pd.Series(range(len(format_order))).to_numpy()
                width = 0.3
                bars1 = ax.bar(x - width / 2, fp16_vals, width=width, label="block FP16", color="#1f77b4")
                bars2 = ax.bar(x + width / 2, mx_vals, width=width, label="block MX FP8 (fp16 bucket)", color="#ff7f0e")
                ax.set_xticks(x, format_order)
                ax.tick_params(axis="x", rotation=45)
                row_max[idx // 2] = max(row_max[idx // 2], pd.Series([
                    *[b.get_height() for b in bars1],
                    *[b.get_height() for b in bars2],
                ]).max(skipna=True))
                if fig_handles is None:
                    fig_handles, fig_labels = ax.get_legend_handles_labels()

                for bars in (bars1, bars2):
                    for b in bars:
                        if pd.isna(b.get_height()):
                            continue
                        ax.annotate(f"{b.get_height():.2e}",
                                    xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                                    xytext=(0, 3), textcoords="offset points",
                                    ha="center", va="bottom", fontsize=7, rotation=45)
            else:
                ax.set_xticks(pd.Series(range(len(format_order))).to_numpy(), format_order)

            ax.set_title(f"NB={nb}")
            ax.set_xlabel("Low precision format (MX)")
            ax.set_ylabel("rel_factor_error")

        fig.suptitle(f"Figure 8: Block FP16 vs MX-FP8 replacement (size={fmt_size(data_size)})")
        if fig_handles:
            fig.legend(fig_handles, fig_labels, loc="lower center", ncol=2,
                       bbox_to_anchor=(0.5, -0.02), frameon=False)
        for r in range(2):
            max_y = row_max[r]
            if max_y > 0:
                for c in range(2):
                    ax = axes[r * 2 + c]
                    if ax.get_visible():
                        ax.set_ylim(0, max_y * 1.25)
        fig.tight_layout(rect=[0, 0.1, 1, 0.95])
        out_path = out_dir / f"fig8_block_fp16_vs_mxfp8_size{fmt_size(data_size)}.pdf"
        fig.savefig(out_path, dpi=200)
        plt.close(fig)


def main():
    args = parse_args()
    summary_path = Path(args.summary)
    if not summary_path.exists():
        raise SystemExit(f"Summary file not found: {summary_path}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = read_summary_table(summary_path)
    df = clean_df(df)

    plot_figure1(df, out_dir)
    plot_figure2(df, out_dir)
    plot_figure3(df, out_dir)
    plot_figure4(df, out_dir)
    plot_figure5(df, out_dir)
    plot_figure6(df, out_dir)
    plot_figure7(df, out_dir)
    plot_figure8(df, out_dir)

    print(f"Plots saved to: {out_dir}")


if __name__ == "__main__":
    main()
