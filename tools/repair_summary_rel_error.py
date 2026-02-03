#!/usr/bin/env python3
from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path
import re
import shutil
import sys
from datetime import datetime

BASE_DIR = Path("/home/abduraa/MX_project/logs/mx_ooc_data")
SWEEP_DIR = BASE_DIR / "sweep"
SUMMARY_PATH = BASE_DIR / "summary_rel_error.txt"

CONFIG_RE = re.compile(
    r"^=== BIN=(?P<bin>\S+) MX_MODE=(?P<mx_mode>\S+) NB=(?P<nb>\d+) FORMAT=(?P<format>\S+) "
    r"FP32_BUCKET=(?P<fp32_bucket>\S+) FP16_BUCKET=(?P<fp16_bucket>\S+) ===$"
)
ERROR_RE = re.compile(r"^error:\s*([0-9.eE+-]+)")
SUMMARY_COUNTS_RE = re.compile(r"^\[SUMMARY_COUNTS\]\s+([0-9\t ]+)$")


def load_existing(summary_path: Path):
    existing_lines = set()
    header = None
    n_by_bin = Counter()
    cores_by_bin = Counter()

    if not summary_path.exists():
        return header, existing_lines, n_by_bin, cores_by_bin

    lines = summary_path.read_text().splitlines()
    if not lines:
        return header, existing_lines, n_by_bin, cores_by_bin

    header = lines[0]
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        existing_lines.add(line)
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        bin_path = parts[0]
        n_val = parts[4].strip() if len(parts) > 4 else ""
        cores_val = parts[5].strip() if len(parts) > 5 else ""
        if n_val:
            n_by_bin[(bin_path, n_val)] += 1
        if cores_val:
            cores_by_bin[(bin_path, cores_val)] += 1

    return header, existing_lines, n_by_bin, cores_by_bin


def most_common_value(counter: Counter, bin_path: str, default: str = "") -> str:
    candidates = [(k[1], v) for k, v in counter.items() if k[0] == bin_path]
    if not candidates:
        return default
    candidates.sort(key=lambda x: (-x[1], x[0]))
    return candidates[0][0]


def load_configs(sweep_dir: Path):
    configs = []
    files = sorted(sweep_dir.glob("sweep_nohup_*.log"))
    for path in files:
        for line in path.read_text().splitlines():
            line = line.strip()
            m = CONFIG_RE.match(line)
            if m:
                configs.append(m.groupdict())
    return configs


def load_sweep_logs(sweep_dir: Path):
    logs = []
    for path in sweep_dir.glob("sweep_*.log"):
        if path.name.startswith("sweep_nohup_"):
            continue
        if path.name == "sweep_errors.log":
            continue
        logs.append(path)
    return sorted(logs)


def load_run_logs(base_dir: Path):
    return sorted(base_dir.glob("run_*.log"))


def parse_sweep_log(path: Path):
    rel_error = None
    counts = None
    for line in path.read_text().splitlines():
        if rel_error is None:
            m = ERROR_RE.match(line.strip())
            if m:
                rel_error = m.group(1)
        m = SUMMARY_COUNTS_RE.match(line.strip())
        if m:
            counts_text = m.group(1).strip()
            parts = re.split(r"\s+", counts_text)
            if len(parts) >= 7:
                counts = parts[:7]
    return rel_error, counts


def main() -> int:
    header, existing_lines, n_by_bin, cores_by_bin = load_existing(SUMMARY_PATH)

    configs = load_configs(SWEEP_DIR)
    sweep_logs = load_sweep_logs(SWEEP_DIR)
    run_logs = load_run_logs(BASE_DIR)

    if not configs:
        print("No sweep_nohup configs found.")
        return 1
    if not sweep_logs and not run_logs:
        print("No sweep or run logs found.")
        return 1

    log_source = sweep_logs if len(sweep_logs) >= len(run_logs) else run_logs
    source_name = "sweep" if log_source is sweep_logs else "run"

    total_pairs = min(len(configs), len(log_source))
    if len(configs) != len(log_source):
        print(
            f"Warning: configs={len(configs)} {source_name}_logs={len(log_source)}; "
            f"pairing first {total_pairs} entries."
        )

    new_lines = []

    for idx in range(total_pairs):
        cfg = configs[idx]
        log_path = log_source[idx]
        rel_error, counts = parse_sweep_log(log_path)
        if not rel_error or not counts:
            continue

        bin_path = cfg["bin"]
        fmt = cfg["format"]
        mx_mode = cfg["mx_mode"]
        nb = cfg["nb"]
        n_val = most_common_value(n_by_bin, bin_path, default="")
        cores_val = most_common_value(cores_by_bin, bin_path, default="2")

        line = "\t".join([
            bin_path,
            fmt,
            mx_mode,
            nb,
            n_val,
            cores_val,
            rel_error,
            *counts,
        ])

        if line not in existing_lines:
            new_lines.append(line)
            existing_lines.add(line)

    if not new_lines:
        print("No missing entries found to add.")
        return 0

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_path = SUMMARY_PATH.with_suffix(f".txt.bak_{timestamp}")
    if SUMMARY_PATH.exists():
        shutil.copy2(SUMMARY_PATH, backup_path)
        print(f"Backup created: {backup_path}")

    if header is None:
        header = "bin\tformat\tmx_mode\tnb\tn\tcores\trel_factor_error\ttotal_tiles\tfp64\tfp32\tfp16\tbf16\tmx_fp16\tlow"

    with SUMMARY_PATH.open("a") as f:
        if SUMMARY_PATH.stat().st_size == 0:
            f.write(header + "\n")
        for line in new_lines:
            f.write(line + "\n")

    print(f"Added {len(new_lines)} missing entries to {SUMMARY_PATH}.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
