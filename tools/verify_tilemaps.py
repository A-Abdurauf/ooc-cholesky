#!/usr/bin/env python3
"""Sanity-check every .tilemap in a directory and rebuild index.tsv.

For each ladder_ftz_N*_nb*_eps*.tilemap:
  - Parse the 3-line header (N, nb, nt, eps).
  - Verify line count == nt*(nt+1)/2.
  - Verify every (k, k) is FP64.
  - Verify each line is "r c FORMAT" with valid FORMAT and r >= c.
  - Count formats.
Writes index.tsv with one row per file: N nb nt eps tiles FP64 FP32 MXFP16
MXFP8 MXFP4 FP16 FP8_plain filename.
"""
import argparse
import re
import sys
from collections import Counter
from pathlib import Path


HEADER_RE = re.compile(r"^# N=(\d+) nb=(\d+) nt=(\d+)$")
EPS_RE = re.compile(r"^# eps=(\S+)$")
ROW_RE = re.compile(r"^(\d+) (\d+) (\S+)$")
FORMATS = ["FP64", "FP32", "MXFP16", "MXFP8", "MXFP4", "FP16", "FP8_plain"]


def check(tilemap: Path) -> dict:
    lines = tilemap.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"{tilemap}: header truncated")
    m = HEADER_RE.match(lines[0])
    if not m:
        raise ValueError(f"{tilemap}: bad header line 0: {lines[0]!r}")
    n, nb, nt = (int(x) for x in m.groups())
    m = EPS_RE.match(lines[1])
    if not m:
        raise ValueError(f"{tilemap}: bad eps line 1: {lines[1]!r}")
    eps = m.group(1)
    if not lines[2].startswith("# mode="):
        raise ValueError(f"{tilemap}: bad mode line 2: {lines[2]!r}")
    body = lines[3:]
    expected = nt * (nt + 1) // 2
    if len(body) != expected:
        raise ValueError(
            f"{tilemap}: {len(body)} body lines, expected nt*(nt+1)/2={expected}"
        )

    counts = Counter()
    diag_seen = {}
    for i, ln in enumerate(body):
        m = ROW_RE.match(ln)
        if not m:
            raise ValueError(f"{tilemap}: bad body line {i+3}: {ln!r}")
        r, c, fmt = int(m.group(1)), int(m.group(2)), m.group(3)
        if r < c:
            raise ValueError(f"{tilemap}: upper-tri tile ({r},{c}) in body")
        if fmt not in FORMATS:
            raise ValueError(f"{tilemap}: unknown FORMAT {fmt!r} at ({r},{c})")
        if r == c:
            diag_seen[r] = fmt
        counts[fmt] += 1

    if len(diag_seen) != nt or any(diag_seen[k] != "FP64" for k in range(nt)):
        bad = [k for k in range(nt) if diag_seen.get(k) != "FP64"]
        raise ValueError(f"{tilemap}: non-FP64 diagonal tiles at {bad}")

    return {
        "n": n, "nb": nb, "nt": nt, "eps": eps,
        "tiles": expected, "counts": counts,
        "file": tilemap.name,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", required=True, type=Path)
    args = ap.parse_args()

    files = sorted(args.dir.glob("ladder_ftz_N*_nb*_eps*.tilemap"))
    if not files:
        print(f"No tilemap files in {args.dir}", file=sys.stderr)
        return 2

    rows = []
    failures = []
    for f in files:
        try:
            rows.append(check(f))
        except ValueError as e:
            failures.append(str(e))

    if failures:
        print("Failures:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)

    # Sort: by n then by numeric eps.
    def keyf(r):
        return (r["n"], float(r["eps"]))

    rows.sort(key=keyf)

    index = args.dir / "index.tsv"
    with index.open("w") as f:
        f.write("n\tnb\tnt\teps\ttiles\t" + "\t".join(FORMATS) + "\tfile\n")
        for r in rows:
            cells = [str(r["n"]), str(r["nb"]), str(r["nt"]), r["eps"],
                     str(r["tiles"])]
            for fmt in FORMATS:
                cells.append(str(r["counts"].get(fmt, 0)))
            cells.append(r["file"])
            f.write("\t".join(cells) + "\n")

    print(f"OK: {len(rows)} tilemap files verified. index.tsv -> {index}")
    if failures:
        print(f"FAIL: {len(failures)} files had errors (see stderr).", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
