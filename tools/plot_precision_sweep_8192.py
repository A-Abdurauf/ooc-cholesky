#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd


def parse_args():
    p = argparse.ArgumentParser(
        description="Plot precision sweep (8k) from summary_precision_sweep_8192.csv"
    )
    p.add_argument(
        "--summary",
        default="/home/abduraa/MX_project/logs/mx_ooc_data/summary_precision_sweep_8192.csv",
        help="Path to summary_precision_sweep_8192.csv",
    )
    p.add_argument(
        "--out-dir",
        default="/home/abduraa/MX_project/logs/mx_ooc_data/plots/precision_sweep_8192",
        help="Output directory for plots",
    )
    return p.parse_args()


def clean_df(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = [c.strip() for c in df.columns]
    for col in [
        "nb",
        "n",
        "cores",
        "rel_factor_error",
        "kl_divergence",
        "total_tiles",
        "fp64",
        "fp32",
        "fp16",
        "bf16",
        "mx_fp16",
        "mx_fp32",
        "mx_e4m3",
        "mx_e5m2",
        "fp8_e4m3",
        "fp8_e5m2",
        "e3m2",
        "e2m3",
        "e2m1",
    ]:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    return df


def norm_format(fmt: str) -> str:
    if not isinstance(fmt, str):
        return ""
    f = fmt.strip().lower()
    aliases = {
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
    return aliases.get(f, f)


def display_format_label(fmt: str) -> str:
    if not isinstance(fmt, str):
        return ""
    f = fmt.strip().lower()
    if f == "e4m3":
        return "mxfp8_e4m3"
    return f


def infer_dataset(bin_path: str) -> str:
    if not isinstance(bin_path, str):
        return "other"
    m = re.search(r"my_cov_(weak|medium)", bin_path.lower())
    if m:
        return m.group(1)
    return "other"


def plot_same_setting_bucket_compare(
    df: pd.DataFrame,
    out_dir: Path,
    target_format: str,
    fig_name: str,
    title: str,
):
    """
    Compare rel_factor_error for same settings, separately for weak/medium:
      - FP32+FP16  (fp32_bucket=fp32, fp16_bucket=fp16)
      - FP32+E4M3  (fp32_bucket=fp32, fp16_bucket=e4m3)
    for a chosen `target_format` (e4m3 or fp8_e4m3).
    """
    if "fp32_bucket" not in df.columns or "fp16_bucket" not in df.columns:
        return

    df_use = df.copy()
    df_use["format_norm"] = df_use["format"].apply(norm_format)
    df_use["dataset"] = df_use["bin"].apply(infer_dataset)
    df_use["fp32_bucket_norm"] = df_use["fp32_bucket"].astype(str).str.lower()
    df_use["fp16_bucket_norm"] = df_use["fp16_bucket"].astype(str).str.lower()

    df_use = df_use[
        (df_use["format_norm"] == target_format)
        & (df_use["dataset"].isin(["weak", "medium"]))
        & (df_use["fp32_bucket_norm"] == "fp32")
        & (df_use["fp16_bucket_norm"].isin(["fp16", "e4m3"]))
    ].copy()

    if df_use.empty:
        return

    fig, axes = plt.subplots(2, 1, figsize=(14, 8), sharey=False)
    datasets = ["weak", "medium"]
    fig_handles = None
    fig_labels = None

    for i, ds in enumerate(datasets):
        ax = axes[i]
        sub = df_use[df_use["dataset"] == ds].copy()
        if sub.empty:
            ax.set_visible(False)
            continue

        sub["setting_key"] = (
            "eps=" + sub["source_epsilon"].astype(str)
            + " | mode=" + sub["mx_mode"].astype(str)
            + " | nb=" + sub["nb"].astype("Int64").astype(str)
        )

        fp16_side = sub[sub["fp16_bucket_norm"] == "fp16"][
            ["setting_key", "rel_factor_error"]
        ].drop_duplicates(subset=["setting_key"])
        e4m3_side = sub[sub["fp16_bucket_norm"] == "e4m3"][
            ["setting_key", "rel_factor_error"]
        ].drop_duplicates(subset=["setting_key"])

        merged = fp16_side.merge(
            e4m3_side,
            on="setting_key",
            how="outer",
            suffixes=("_fp16", "_e4m3"),
        )

        if merged.empty:
            ax.set_visible(False)
            continue

        merged = merged.sort_values("setting_key")
        x = pd.Series(range(len(merged))).to_numpy()
        width = 0.4

        bars1 = ax.bar(
            x - width / 2,
            merged["rel_factor_error_fp16"].to_numpy(),
            width=width,
            label="FP32+FP16",
            color="#4c78a8",
        )
        bars2 = ax.bar(
            x + width / 2,
            merged["rel_factor_error_e4m3"].to_numpy(),
            width=width,
            label="FP32+E4M3",
            color="#f58518",
        )

        ax.set_title(ds)
        ax.set_xlabel("Same settings (epsilon, mode, nb)")
        ax.set_ylabel("rel_factor_error")
        ax.set_xticks(x, merged["setting_key"].tolist())
        ax.tick_params(axis="x", rotation=30)

        for bars in (bars1, bars2):
            for b in bars:
                if pd.isna(b.get_height()):
                    continue
                ax.annotate(
                    f"{b.get_height():.2e}",
                    xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    rotation=45,
                )

        if fig_handles is None:
            fig_handles, fig_labels = ax.get_legend_handles_labels()

    fig.suptitle(title)
    if fig_handles:
        fig.legend(
            fig_handles,
            fig_labels,
            loc="lower center",
            ncol=2,
            bbox_to_anchor=(0.5, -0.02),
            frameon=False,
        )
    fig.tight_layout(rect=[0, 0.08, 1, 0.95])
    fig.savefig(out_dir / fig_name, dpi=200)
    plt.close(fig)


def plot_rel_error_by_format(df: pd.DataFrame, out_dir: Path):
    # Figure A: rel_factor_error by format for each nb, comparing tile vs block
    df_use = df.copy()
    df_use["format_norm"] = df_use["format"].apply(norm_format)

    nb_order = [32, 64, 128, 256]
    format_order = ["fp8_e4m3", "fp8_e5m2", "e4m3", "e5m2", "e3m2", "e2m3", "e2m1"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharey=False)
    axes = axes.flatten()
    row_max = [0.0, 0.0]
    fig_handles = None
    fig_labels = None

    for idx, nb in enumerate(nb_order):
        ax = axes[idx]
        sub = df_use[df_use["nb"] == nb]
        if sub.empty:
            ax.set_visible(False)
            continue

        tile = sub[sub["mx_mode"] == "tile"]
        block = sub[sub["mx_mode"] == "block"]

        tile2 = tile.drop_duplicates(subset=["format_norm"]).set_index("format_norm")[
            "rel_factor_error"
        ]
        block2 = block.drop_duplicates(subset=["format_norm"]).set_index("format_norm")[
            "rel_factor_error"
        ]

        x = pd.Series(range(len(format_order))).to_numpy()
        width = 0.4
        bars1 = ax.bar(
            x - width / 2,
            [tile2.get(f, float("nan")) for f in format_order],
            width=width,
            label="tile",
            color="#4c78a8",
        )
        bars2 = ax.bar(
            x + width / 2,
            [block2.get(f, float("nan")) for f in format_order],
            width=width,
            label="block",
            color="#f58518",
        )

        ax.set_title(f"NB={nb}")
        ax.set_xlabel("Format")
        ax.set_ylabel("rel_factor_error")
        ax.set_xticks(x, [display_format_label(f) for f in format_order])
        ax.tick_params(axis="x", rotation=45)

        row_max[idx // 2] = max(
            row_max[idx // 2],
            pd.Series([b.get_height() for b in bars1] + [b.get_height() for b in bars2]).max(
                skipna=True
            ),
        )

        if fig_handles is None:
            fig_handles, fig_labels = ax.get_legend_handles_labels()

        for bars in (bars1, bars2):
            for b in bars:
                if pd.isna(b.get_height()):
                    continue
                ax.annotate(
                    f"{b.get_height():.2e}",
                    xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    rotation=45,
                )

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

    fig.suptitle("Figure A: rel_factor_error by format (8k)")
    fig.tight_layout(rect=[0, 0.08, 1, 0.95])
    fig.savefig(out_dir / "figA_rel_error_by_format_8k.pdf", dpi=200)
    plt.close(fig)


def plot_bucket_comparison(df: pd.DataFrame, out_dir: Path):
    # Figure B: FP32 bucket vs MX-FP32 bucket (rel_factor_error), grouped by format
    if "fp32_bucket" not in df.columns:
        return

    df_use = df.copy()
    df_use["format_norm"] = df_use["format"].apply(norm_format)
    df_use["fp32_bucket_norm"] = df_use["fp32_bucket"].astype(str).str.lower()

    nb_order = [32, 64, 128, 256]
    format_order = ["fp8_e4m3", "fp8_e5m2", "e4m3", "e5m2", "e3m2", "e2m3", "e2m1"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharey=False)
    axes = axes.flatten()
    row_max = [0.0, 0.0]
    fig_handles = None
    fig_labels = None

    for idx, nb in enumerate(nb_order):
        ax = axes[idx]
        sub = df_use[df_use["nb"] == nb]
        if sub.empty:
            ax.set_visible(False)
            continue

        fp32_vals = []
        mx_vals = []
        for fmt in format_order:
            g = sub[sub["format_norm"] == fmt]
            fp32_vals.append(
                g[g["fp32_bucket_norm"] == "fp32"]["rel_factor_error"].iloc[0]
                if not g[g["fp32_bucket_norm"] == "fp32"].empty
                else float("nan")
            )
            mx_vals.append(
                g[g["fp32_bucket_norm"] == "mx_fp32"]["rel_factor_error"].iloc[0]
                if not g[g["fp32_bucket_norm"] == "mx_fp32"].empty
                else float("nan")
            )

        x = pd.Series(range(len(format_order))).to_numpy()
        width = 0.35
        bars1 = ax.bar(x - width / 2, fp32_vals, width=width, label="fp32 bucket", color="#4c78a8")
        bars2 = ax.bar(x + width / 2, mx_vals, width=width, label="mx_fp32 bucket", color="#f58518")

        ax.set_title(f"NB={nb}")
        ax.set_xlabel("Format")
        ax.set_ylabel("rel_factor_error")
        ax.set_xticks(x, [display_format_label(f) for f in format_order])
        ax.tick_params(axis="x", rotation=45)

        row_max[idx // 2] = max(
            row_max[idx // 2],
            pd.Series([b.get_height() for b in bars1] + [b.get_height() for b in bars2]).max(
                skipna=True
            ),
        )

        if fig_handles is None:
            fig_handles, fig_labels = ax.get_legend_handles_labels()

        for bars in (bars1, bars2):
            for b in bars:
                if pd.isna(b.get_height()):
                    continue
                ax.annotate(
                    f"{b.get_height():.2e}",
                    xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    rotation=45,
                )

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

    fig.suptitle("Figure B: FP32 vs MX-FP32 buckets (8k)")
    fig.tight_layout(rect=[0, 0.08, 1, 0.95])
    fig.savefig(out_dir / "figB_fp32_bucket_vs_mx_fp32_8k.pdf", dpi=200)
    plt.close(fig)


def plot_fp16_bucket_comparison(df: pd.DataFrame, out_dir: Path):
    # Figure C: FP16 bucket vs MX-FP16 bucket (rel_factor_error), grouped by format
    if "fp16_bucket" not in df.columns:
        return

    df_use = df.copy()
    df_use["format_norm"] = df_use["format"].apply(norm_format)
    df_use["fp16_bucket_norm"] = df_use["fp16_bucket"].astype(str).str.lower()

    nb_order = [32, 64, 128, 256]
    format_order = ["fp8_e4m3", "fp8_e5m2", "e4m3", "e5m2", "e3m2", "e2m3", "e2m1"]

    fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharey=False)
    axes = axes.flatten()
    row_max = [0.0, 0.0]
    fig_handles = None
    fig_labels = None

    for idx, nb in enumerate(nb_order):
        ax = axes[idx]
        sub = df_use[df_use["nb"] == nb]
        if sub.empty:
            ax.set_visible(False)
            continue

        fp16_vals = []
        mx_vals = []
        for fmt in format_order:
            g = sub[sub["format_norm"] == fmt]
            fp16_vals.append(
                g[g["fp16_bucket_norm"] == "fp16"]["rel_factor_error"].iloc[0]
                if not g[g["fp16_bucket_norm"] == "fp16"].empty
                else float("nan")
            )
            mx_vals.append(
                g[g["fp16_bucket_norm"] == "mx_fp16"]["rel_factor_error"].iloc[0]
                if not g[g["fp16_bucket_norm"] == "mx_fp16"].empty
                else float("nan")
            )

        x = pd.Series(range(len(format_order))).to_numpy()
        width = 0.35
        bars1 = ax.bar(x - width / 2, fp16_vals, width=width, label="fp16 bucket", color="#4c78a8")
        bars2 = ax.bar(x + width / 2, mx_vals, width=width, label="mx_fp16 bucket", color="#f58518")

        ax.set_title(f"NB={nb}")
        ax.set_xlabel("Format")
        ax.set_ylabel("rel_factor_error")
        ax.set_xticks(x, [display_format_label(f) for f in format_order])
        ax.tick_params(axis="x", rotation=45)

        row_max[idx // 2] = max(
            row_max[idx // 2],
            pd.Series([b.get_height() for b in bars1] + [b.get_height() for b in bars2]).max(
                skipna=True
            ),
        )

        if fig_handles is None:
            fig_handles, fig_labels = ax.get_legend_handles_labels()

        for bars in (bars1, bars2):
            for b in bars:
                if pd.isna(b.get_height()):
                    continue
                ax.annotate(
                    f"{b.get_height():.2e}",
                    xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    rotation=45,
                )

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

    fig.suptitle("Figure C: FP16 vs MX-FP16 buckets (8k)")
    fig.tight_layout(rect=[0, 0.08, 1, 0.95])
    fig.savefig(out_dir / "figC_fp16_bucket_vs_mx_fp16_8k.pdf", dpi=200)
    plt.close(fig)


def plot_single_fp32_fp16_precision_compare(df: pd.DataFrame, out_dir: Path):
    """
    One figure only:
      - Keep same bucket setting: fp32_bucket=fp32 and fp16_bucket=fp16
      - Compare precision format: e4m3 vs fp8_e4m3
      - Separate subplots: weak and medium
      - Mention mode and nb directly in figure title and file name
    """
    if "fp32_bucket" not in df.columns or "fp16_bucket" not in df.columns:
        return

    df_use = df.copy()
    df_use["format_norm"] = df_use["format"].apply(norm_format)
    df_use["dataset"] = df_use["bin"].apply(infer_dataset)
    df_use["fp32_bucket_norm"] = df_use["fp32_bucket"].astype(str).str.lower()
    df_use["fp16_bucket_norm"] = df_use["fp16_bucket"].astype(str).str.lower()

    df_use = df_use[
        (df_use["dataset"].isin(["weak", "medium"]))
        & (df_use["fp32_bucket_norm"] == "fp32")
        & (df_use["fp16_bucket_norm"] == "fp16")
        & (df_use["format_norm"].isin(["e4m3", "fp8_e4m3"]))
    ].copy()

    if df_use.empty:
        return

    # User asked for one figure; pick the dominant mode/nb pair.
    pair_counts = (
        df_use.groupby(["mx_mode", "nb"]).size().sort_values(ascending=False)
    )
    mode, nb = pair_counts.index[0]
    df_use = df_use[(df_use["mx_mode"] == mode) & (df_use["nb"] == nb)].copy()

    n_values = sorted(pd.to_numeric(df_use["n"], errors="coerce").dropna().astype(int).unique().tolist())
    if not n_values:
        n_label = "unknown"
    elif len(n_values) == 1:
        n_label = str(n_values[0])
    else:
        n_label = "/".join(str(v) for v in n_values)

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharey=False)
    datasets = ["weak", "medium"]
    fmt_order = ["e4m3", "fp8_e4m3"]
    fmt_labels = {"e4m3": "mxfp8_e4m3", "fp8_e4m3": "fp8_e4m3"}
    colors = {"e4m3": "#4c78a8", "fp8_e4m3": "#f58518"}

    for i, ds in enumerate(datasets):
        ax = axes[i]
        sub = df_use[df_use["dataset"] == ds].copy()
        if sub.empty:
            ax.set_visible(False)
            continue

        sub["eps_num"] = pd.to_numeric(sub["source_epsilon"], errors="coerce")
        eps_vals = sorted(sub["eps_num"].dropna().unique().tolist())
        if not eps_vals:
            ax.set_visible(False)
            continue

        x = pd.Series(range(len(eps_vals))).to_numpy()
        width = 0.35

        values_by_fmt = {}
        for fmt in fmt_order:
            vals = []
            for eps in eps_vals:
                g = sub[(sub["eps_num"] == eps) & (sub["format_norm"] == fmt)]
                vals.append(g["rel_factor_error"].iloc[0] if not g.empty else float("nan"))
            values_by_fmt[fmt] = vals

        # Log scale requires strictly positive values.
        for fmt in fmt_order:
            values_by_fmt[fmt] = [
                (v if pd.notna(v) and v > 0 else float("nan"))
                for v in values_by_fmt[fmt]
            ]

        bars1 = ax.bar(
            x - width / 2,
            values_by_fmt["e4m3"],
            width=width,
            label=fmt_labels["e4m3"],
            color=colors["e4m3"],
        )
        bars2 = ax.bar(
            x + width / 2,
            values_by_fmt["fp8_e4m3"],
            width=width,
            label=fmt_labels["fp8_e4m3"],
            color=colors["fp8_e4m3"],
        )

        eps_labels = [f"{e:.0e}" for e in eps_vals]
        ax.set_xticks(x, eps_labels)
        ax.set_xlabel("source_epsilon")
        ax.set_ylabel("rel_factor_error (log scale)")
        ax.set_title(ds)
        ax.set_yscale("log")

        positive_vals = [
            v
            for fmt in fmt_order
            for v in values_by_fmt[fmt]
            if pd.notna(v) and v > 0
        ]
        if positive_vals:
            y_min = min(positive_vals)
            y_max = max(positive_vals)
            if y_min == y_max:
                y_min *= 0.8
                y_max *= 1.2
            else:
                y_min *= 0.6
                y_max *= 1.8
            ax.set_ylim(y_min, y_max)

        for bars in (bars1, bars2):
            for b in bars:
                if pd.isna(b.get_height()):
                    continue
                ax.annotate(
                    f"{b.get_height():.2e}",
                    xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                    xytext=(0, 3),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    rotation=45,
                )

        ax.legend(loc="upper left")

    fig.suptitle(
        f"Figure 1 (n={n_label}, mode={mode}, nb={int(nb)}): FP32+FP16 bucket, mxfp8_e4m3 vs fp8_e4m3"
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    mode_tag = re.sub(r"[^a-zA-Z0-9_]+", "_", str(mode))
    out_name = f"fig1_fp32_fp16_e4m3_vs_fp8e4m3_mode_{mode_tag}_nb_{int(nb)}_8k.pdf"
    fig.savefig(out_dir / out_name, dpi=200)
    plt.close(fig)


def plot_tile_count_composition(df: pd.DataFrame, out_dir: Path):
    """
        Figure 2:
            For same setting (fp32_bucket=fp32, fp16_bucket=fp16), show tile counts:
                - FP32
                - FP16
                - Low
            separately for weak and medium.
    """
    needed_cols = {"fp32_bucket", "fp16_bucket", "fp64", "fp32", "fp16", "mx_e4m3", "fp8_e4m3"}
    if not needed_cols.issubset(set(df.columns)):
        return

    df_use = df.copy()
    df_use["format_norm"] = df_use["format"].apply(norm_format)
    df_use["dataset"] = df_use["bin"].apply(infer_dataset)
    df_use["fp32_bucket_norm"] = df_use["fp32_bucket"].astype(str).str.lower()
    df_use["fp16_bucket_norm"] = df_use["fp16_bucket"].astype(str).str.lower()

    df_use = df_use[
        (df_use["dataset"].isin(["weak", "medium"]))
        & (df_use["fp32_bucket_norm"] == "fp32")
        & (df_use["fp16_bucket_norm"] == "fp16")
        & (df_use["format_norm"].isin(["e4m3", "fp8_e4m3"]))
    ].copy()

    if df_use.empty:
        return

    pair_counts = df_use.groupby(["mx_mode", "nb"]).size().sort_values(ascending=False)
    mode, nb = pair_counts.index[0]
    df_use = df_use[(df_use["mx_mode"] == mode) & (df_use["nb"] == nb)].copy()

    n_values = sorted(pd.to_numeric(df_use["n"], errors="coerce").dropna().astype(int).unique().tolist())
    if not n_values:
        n_label = "unknown"
    elif len(n_values) == 1:
        n_label = str(n_values[0])
    else:
        n_label = "/".join(str(v) for v in n_values)

    fig, axes = plt.subplots(2, 1, figsize=(14, 9), sharey=False)
    datasets = ["weak", "medium"]
    for i, ds in enumerate(datasets):
        ax = axes[i]
        sub = df_use[df_use["dataset"] == ds].copy()
        if sub.empty:
            ax.set_visible(False)
            continue

        sub["eps_num"] = pd.to_numeric(sub["source_epsilon"], errors="coerce")
        eps_vals = sorted(sub["eps_num"].dropna().unique().tolist())
        if not eps_vals:
            ax.set_visible(False)
            continue

        x = pd.Series(range(len(eps_vals))).to_numpy()
        width = 0.2

        low_vals, fp64_vals, fp32_vals, fp16_vals = [], [], [], []

        for eps in eps_vals:
            g = sub[sub["eps_num"] == eps]
            if g.empty:
                low_vals.append(float("nan"))
                fp64_vals.append(float("nan"))
                fp32_vals.append(float("nan"))
                fp16_vals.append(float("nan"))
                continue

            # FP32/FP16 counts are shared for this bucket setting.
            fp64_val = pd.to_numeric(g["fp64"], errors="coerce").dropna()
            fp32_val = pd.to_numeric(g["fp32"], errors="coerce").dropna()
            fp16_val = pd.to_numeric(g["fp16"], errors="coerce").dropna()
            fp64_vals.append(fp64_val.iloc[0] if not fp64_val.empty else float("nan"))
            fp32_vals.append(fp32_val.iloc[0] if not fp32_val.empty else float("nan"))
            fp16_vals.append(fp16_val.iloc[0] if not fp16_val.empty else float("nan"))

            # Low count from either representation (e4m3 or fp8_e4m3).
            low_candidates = pd.concat(
                [
                    pd.to_numeric(g["mx_e4m3"], errors="coerce"),
                    pd.to_numeric(g["fp8_e4m3"], errors="coerce"),
                ],
                ignore_index=True,
            ).dropna()
            low_vals.append(low_candidates.max() if not low_candidates.empty else float("nan"))

        bars1 = ax.bar(x - 1.5 * width, fp64_vals, width=width, label="FP64", color="#2f4b7c")
        bars2 = ax.bar(x - 0.5 * width, fp32_vals, width=width, label="FP32", color="#4c78a8")
        bars3 = ax.bar(x + 0.5 * width, fp16_vals, width=width, label="FP16", color="#72b7b2")
        bars4 = ax.bar(x + 1.5 * width, low_vals, width=width, label="Low", color="#f58518")

        eps_labels = [f"{e:.0e}" for e in eps_vals]
        ax.set_xticks(x, eps_labels)
        ax.set_xlabel("source_epsilon")
        ax.set_ylabel("tile count")
        ax.set_title(ds)

        for bars in (bars1, bars2, bars3, bars4):
            for b in bars:
                if pd.isna(b.get_height()):
                    continue
                ax.annotate(
                    f"{b.get_height():.0f}",
                    xy=(b.get_x() + b.get_width() / 2, b.get_height()),
                    xytext=(0, 2),
                    textcoords="offset points",
                    ha="center",
                    va="bottom",
                    fontsize=6,
                    rotation=90,
                )

        ax.legend(loc="upper left", ncol=4, fontsize=8)

    fig.suptitle(
        f"Figure 2 (n={n_label}, mode={mode}, nb={int(nb)}): tile counts for FP32+FP16 bucket (FP64, FP32, FP16, Low)"
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    mode_tag = re.sub(r"[^a-zA-Z0-9_]+", "_", str(mode))
    out_name = f"fig2_tile_counts_fp32_fp16_low_mode_{mode_tag}_nb_{int(nb)}_8k.pdf"
    fig.savefig(out_dir / out_name, dpi=200)
    plt.close(fig)


def main():
    args = parse_args()
    summary_path = Path(args.summary)
    if not summary_path.exists():
        raise SystemExit(f"Summary file not found: {summary_path}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(summary_path)
    df = clean_df(df)

    # One requested figure only.
    plot_single_fp32_fp16_precision_compare(df, out_dir)
    plot_tile_count_composition(df, out_dir)

    print(f"Plots saved to: {out_dir}")


if __name__ == "__main__":
    main()