#!/usr/bin/env python3
"""Walk every error-data CSV produced by the various sweeps, normalise into
one schema, dedupe, and write a single merged CSV.  Also print a coverage
audit: which (mode, granularity, underflow, N, nb, eps) cells are missing.

Logical schema for the merged CSV:
    mode           : baseline | ladder_ieee | ladder_mx_staircase | ladder_mx_full | dropin_mxfp4
    granularity    : tile | vec1d | n/a
    underflow      : fz | gu | n/a    (n/a where the format is underflow-invariant)
    N              : matrix dimension
    nb             : tile size
    eps            : source epsilon (string: 1e-5/1e-6/1e-7/1e-8)
    rel_factor_error / abs_factor_error / relative_residual : the measured numbers
    tile_breakdown : per-format counts (semicolon-separated lower-tri counts,
                     or comma-separated full-square counts when that is the
                     only thing the source CSV recorded)
    tile_breakdown_kind : "lower" or "full" -- so consumers know what to do
    source_csv     : path of the originating CSV (provenance)

PDF / TSV not touched -- this just consolidates the raw error rows.
"""
import argparse
import csv
from pathlib import Path
from collections import defaultdict


SOURCES = {
    # path -> kind of tile_breakdown column ("lower" = already lower-triangular,
    # "full" = full-square counts that consumers must halve themselves).
    "/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv":         "full",
    "/home/abduraa/MX_project/ooc-cholesky/ladder_rerun_32k/results.csv":          "lower",
    "/home/abduraa/MX_project/ooc-cholesky/ladder_full_gu_missing/results.csv":    "lower",
    "/home/abduraa/MX_project/ooc-cholesky/mx_staircase_vec1d32_gt20k/results.csv":"lower",
    "/home/abduraa/MX_project/ooc-cholesky/mx_staircase_tile_gt20k/results.csv":   "lower",
    "/home/abduraa/MX_project/ooc-cholesky/mx_staircase_40k_nb4096_gu/results.csv":"lower",
    "/home/abduraa/MX_project/ooc-cholesky/fill_4class_matrix/results.csv":        "lower",
    "/home/abduraa/MX_project/ooc-cholesky/true_dropin_mxfp4_32k/results.csv":     "lower",
    "/home/abduraa/MX_project/ooc-cholesky/baseline_ieee_gu/results.csv":           "lower",
}

# (sweep name in source CSV) -> (mode, granularity, underflow)
SWEEP_MAP = {
    # Baseline -- main CSV runs default to FZ (MX_UNDERFLOW_MODE=fz).  The
    # FP8_SUBNORMAL=1 flag does keep subnormals in-format, but the bound
    # selector still respects the underflow mode.  Tag as 'fz' here; a GU
    # run is needed for full coverage of the EXPECTED matrix.
    "requant_baseline_fp8_subnormal_gt20k":          ("baseline",            "n/a",   "fz"),
    "baseline_fp8_subnormal_gu":                     ("baseline",            "n/a",   "gu"),

    # IEEE ladder
    "requant_ladder_ieee_gt20k":                     ("ladder_ieee",         "n/a",   "fz"),
    "ladder_ieee_fz":                                ("ladder_ieee",         "n/a",   "fz"),
    "ladder_ieee_gu":                                ("ladder_ieee",         "n/a",   "gu"),  # baseline_ieee_gu sweep

    # MX staircase (no MXFP16)
    "mx_staircase_tile_fz":                          ("ladder_mx_staircase", "tile",  "fz"),
    "mx_staircase_tile_gu":                          ("ladder_mx_staircase", "tile",  "gu"),
    "mx_staircase_vec1d32_fz":                       ("ladder_mx_staircase", "vec1d", "fz"),
    "mx_staircase_vec1d32_gu":                       ("ladder_mx_staircase", "vec1d", "gu"),

    # Full MX ladder (with MXFP16) -- many alias sweeps end up here
    "requant_ladder_scaled_tile_gt20k":              ("ladder_mx_full",      "tile",  "gu"),  # main CSV: base name == denorm
    "requant_ladder_scaled_tile_FTZ_gt20k":          ("ladder_mx_full",      "tile",  "fz"),
    "requant_ladder_scaled_vec1d32_gt20k":           ("ladder_mx_full",      "vec1d", "gu"),  # main CSV semantics
    "requant_ladder_scaled_vec1d32_FTZ_gt20k":       ("ladder_mx_full",      "vec1d", "fz"),  # archived buggy
    "ladder_full_fz":                                ("ladder_mx_full",      "vec1d", "fz"),  # rerun_32k
    "ladder_full_gu":                                ("ladder_mx_full",      "vec1d", "gu"),  # rerun_32k
    "ladder_full_vec1d_fz":                          ("ladder_mx_full",      "vec1d", "fz"),  # gu_missing / fill_4class
    "ladder_full_vec1d_gu":                          ("ladder_mx_full",      "vec1d", "gu"),
    "ladder_full_tile_fz":                           ("ladder_mx_full",      "tile",  "fz"),
    "ladder_full_tile_gu":                           ("ladder_mx_full",      "tile",  "gu"),

    # Drop-in MXFP4
    "true_mxfp4_dropin_fz":                          ("dropin_mxfp4",        "vec1d", "fz"),
    "true_mxfp4_dropin_gu":                          ("dropin_mxfp4",        "vec1d", "gu"),
}

# Expected cells per logical config (mode, granularity, underflow).  N x NB.
NB_BY_N = {32768: 2048, 40960: 4096, 65536: 4096}
EPS_LIST = ["1e-5", "1e-6", "1e-7", "1e-8"]

EXPECTED = [
    # (mode, granularity, underflow)  --  enumerate cells over NB_BY_N x EPS_LIST
    ("baseline",            "n/a",   "fz"),   # was: "n/a" -- treat FZ as the default we have
    ("baseline",            "n/a",   "gu"),   # NEW: need a GU run for full coverage
    ("ladder_ieee",         "n/a",   "fz"),
    ("ladder_ieee",         "n/a",   "gu"),   # NEW: need a GU run for full coverage
    ("ladder_mx_staircase", "tile",  "fz"),
    ("ladder_mx_staircase", "tile",  "gu"),
    ("ladder_mx_staircase", "vec1d", "fz"),
    ("ladder_mx_staircase", "vec1d", "gu"),
    ("ladder_mx_full",      "tile",  "fz"),
    ("ladder_mx_full",      "tile",  "gu"),
    ("ladder_mx_full",      "vec1d", "fz"),
    ("ladder_mx_full",      "vec1d", "gu"),
]


def _norm_eps(s):
    return (s or "").strip().replace("1e-05","1e-5").replace("1e-06","1e-6") \
                             .replace("1e-07","1e-7").replace("1e-08","1e-8")


# Some source CSVs use slightly different column names.  Normalise.
def _get(row, *candidates, default=""):
    for c in candidates:
        if c in row and row[c] not in (None, ""):
            return row[c]
    return default


def harvest():
    """Return list of normalised rows from every source CSV."""
    out = []
    for path, kind in SOURCES.items():
        p = Path(path)
        if not p.exists():
            print(f"  [WARN] missing: {p}")
            continue
        with p.open() as f:
            for r in csv.DictReader(f):
                sweep = r.get("sweep", "").strip()
                if sweep not in SWEEP_MAP:
                    # silently skip unrecognised sweeps (e.g. legacy lowscale ones)
                    continue
                mode, gran, uflow = SWEEP_MAP[sweep]
                try:
                    N  = int(r["n"])
                    nb = int(r["nb"])
                except (TypeError, ValueError, KeyError):
                    continue
                eps = _norm_eps(_get(r, "source_epsilon", "eps"))
                rel = _get(r, "rel_factor_error", "rel_error", "relative_error")
                abs_ = _get(r, "abs_factor_error", "error", "abs_error")
                res = _get(r, "relative_residual", "rel_residual")
                tb  = _get(r, "tile_breakdown", "tile_counts_full")
                out.append({
                    "mode":               mode,
                    "granularity":        gran,
                    "underflow":          uflow,
                    "N":                  N,
                    "nb":                 nb,
                    "eps":                eps,
                    "rel_factor_error":   rel,
                    "abs_factor_error":   abs_,
                    "relative_residual":  res,
                    "tile_breakdown":     tb.strip().strip('"'),
                    "tile_breakdown_kind": kind,
                    "source_csv":         path,
                    "source_sweep":       sweep,
                })
    return out


def dedupe(rows):
    """Keep one row per (mode, granularity, underflow, N, nb, eps).

    Preference order: prefer rows from fill_4class / staircase_40n4k / rerun
    over the main CSV when both exist (so we use the post-patch numbers).
    """
    PRIORITY = {  # higher = preferred
        "/home/abduraa/MX_project/ooc-cholesky/fill_4class_matrix/results.csv":         50,
        "/home/abduraa/MX_project/ooc-cholesky/ladder_full_gu_missing/results.csv":     45,
        "/home/abduraa/MX_project/ooc-cholesky/ladder_rerun_32k/results.csv":           40,
        "/home/abduraa/MX_project/ooc-cholesky/mx_staircase_40k_nb4096_gu/results.csv": 35,
        "/home/abduraa/MX_project/ooc-cholesky/mx_staircase_vec1d32_gt20k/results.csv": 30,
        "/home/abduraa/MX_project/ooc-cholesky/mx_staircase_tile_gt20k/results.csv":    30,
        "/home/abduraa/MX_project/ooc-cholesky/true_dropin_mxfp4_32k/results.csv":      25,
        "/home/abduraa/MX_project/logs/mx_ooc_data/requant_gt20k_memory.csv":           10,
    }
    by_key = {}
    for r in rows:
        key = (r["mode"], r["granularity"], r["underflow"],
               r["N"], r["nb"], r["eps"])
        prev = by_key.get(key)
        if prev is None or PRIORITY.get(r["source_csv"], 0) > PRIORITY.get(prev["source_csv"], 0):
            by_key[key] = r
    return list(by_key.values())


def coverage_audit(rows):
    """Report missing cells in the EXPECTED matrix."""
    have = defaultdict(set)  # (mode,gran,uflow,N,nb) -> {eps}
    for r in rows:
        key = (r["mode"], r["granularity"], r["underflow"],
               r["N"], r["nb"])
        have[key].add(r["eps"])

    print()
    print("=" * 100)
    print("COVERAGE AUDIT  (expected: 3 N x 4 eps = 12 cells per config)")
    print("=" * 100)
    print(f"{'mode':<22} {'gran':<6} {'uflow':<6} | "
          f"{'32k NB2048':<14} {'40k NB4096':<14} {'65k NB4096':<14} | missing")
    print("-" * 100)
    total_have = total_want = 0
    for (mode, gran, uflow) in EXPECTED:
        cells = []
        missing = []
        for n in (32768, 40960, 65536):
            nb = NB_BY_N[n]
            present = have.get((mode, gran, uflow, n, nb), set())
            total_have += len(present); total_want += 4
            short = "".join({"1e-5":"5","1e-6":"6","1e-7":"7","1e-8":"8"}.get(e,"?")
                            for e in sorted(present))
            cells.append(f"nb{nb}:{short:<8}")
            for e in EPS_LIST:
                if e not in present:
                    missing.append(f"{n//1024}k {e}")
        miss_s = ", ".join(missing) if missing else "—"
        print(f"{mode:<22} {gran:<6} {uflow:<6} | {cells[0]:<14} {cells[1]:<14} {cells[2]:<14} | {miss_s}")
    print("-" * 100)
    print(f"Coverage: {total_have}/{total_want} = {100*total_have/total_want:.1f}%")
    print()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out",
                    default="/home/abduraa/MX_project/ooc-cholesky/figures/all_error_runs.csv")
    args = ap.parse_args()

    raw = harvest()
    print(f"Harvested {len(raw)} rows from {len(SOURCES)} source CSVs")
    rows = dedupe(raw)
    print(f"After dedupe (prefer fill_4class > rerun > main): {len(rows)} rows")

    coverage_audit(rows)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    cols = ["mode", "granularity", "underflow", "N", "nb", "eps",
            "rel_factor_error", "abs_factor_error", "relative_residual",
            "tile_breakdown", "tile_breakdown_kind",
            "source_csv", "source_sweep"]
    # Sort rows for readability.
    rows.sort(key=lambda r: (r["mode"], r["granularity"], r["underflow"],
                              r["N"], r["nb"], r["eps"]))
    with out.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in cols})
    print(f"Wrote {len(rows)} rows -> {out}")


if __name__ == "__main__":
    main()
