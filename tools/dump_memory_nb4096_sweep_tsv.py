#!/usr/bin/env python3
"""Emit a TSV with the per-(mode, N, eps) tile counts and memory breakdown
used by `plot_paper_memory_nb4096_sweep.py`.

Lower-tri storage with half-diagonal accounting:
  - FP64 diagonals: nb*(nb+1)/2 elements per tile
  - all other lower-tri tiles: nb*nb elements
  - MX-scaled formats add 1 byte per 32-element shared-scale group

Reads tile counts from recreate-ooc-chol/build/runs_all.log and applies the
same FP16 -> MXFP16 reinterpretation as the figure.
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_paper_memory_nb4096_sweep import (  # noqa: E402
    parse_log, memory_breakdown_gb,
    LOG_PATH, EPS_ORDER, N_ORDER, MODE_BARS,
)


COUNT_FMTS = ["FP64", "FP32", "FP16", "MXFP16", "MXFP8", "MXFP4", "FP8_plain"]
GB_FMTS    = COUNT_FMTS + ["Scale"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--log", default=LOG_PATH)
    ap.add_argument("--out",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures/paper_memory_nb4096_sweep.tsv")
    args = ap.parse_args()

    data = parse_log(args.log)

    header = ["mode", "N", "nb", "eps", "M", "tiles_total"]
    header += [f"{f}_tiles" for f in COUNT_FMTS]
    header += [f"{f}_gb" for f in GB_FMTS]
    header += ["total_gb"]

    lines = ["\t".join(header)]
    n_rows = 0
    for mode_lbl, mode_key in MODE_BARS:
        for N in N_ORDER:
            for eps in EPS_ORDER:
                row = data.get((mode_key, N, eps))
                if row is None:
                    continue
                nb = row["nb"]
                M  = N // nb
                tiles_total = M * (M + 1) // 2
                gb = memory_breakdown_gb(row, N)
                tot = sum(gb.get(f, 0.0) for f in GB_FMTS)
                cells = [mode_key, str(N), str(nb), eps, str(M), str(tiles_total)]
                cells += [str(row.get(f, 0)) for f in COUNT_FMTS]
                cells += [f"{gb.get(f, 0.0):.6f}" for f in GB_FMTS]
                cells += [f"{tot:.6f}"]
                lines.append("\t".join(cells))
                n_rows += 1

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text("\n".join(lines) + "\n")
    print(f"{out}  ({n_rows} rows)")


if __name__ == "__main__":
    main()
