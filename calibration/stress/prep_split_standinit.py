#!/usr/bin/env python3
"""Stage 0 (one-time prep): split ENTIRE_FVS_STANDINIT_PLOT.csv by FVS VARIANT.

The full-CONUS stress array runs ~380 tasks. The original conus_100yr_projection.py
rescans the entire 791 MB stand-init on every task just to filter one variant's
batch -- ~380 full-file scans and heavy shared-FS I/O contention. Splitting once
into per-variant files lets each task read only its (small) variant file.

Writes <out_dir>/standinit_<VARIANT>.csv for every variant plus a counts.tsv
(variant, n_stands, n_batches at the given batch size). Rows with a blank VARIANT
are written to standinit__BLANK.csv and excluded from the run manifest.

Usage:
  python prep_split_standinit.py \
      --standinit /fs/scratch/PUOM0008/crsfaaron/FIA/ENTIRE_FVS_STANDINIT_PLOT.csv \
      --out-dir   /fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_by_variant \
      --batch-size 5000
"""
from __future__ import annotations

import argparse
import csv
import math
import os
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--standinit", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--batch-size", type=int, default=5000)
    ap.add_argument("--variant-col", default="VARIANT")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)

    handles: dict[str, object] = {}
    writers: dict[str, object] = {}
    counts: dict[str, int] = {}

    with open(a.standinit, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        try:
            vi = header.index(a.variant_col)
        except ValueError:
            sys.exit(f"variant column {a.variant_col!r} not in header: {header[:6]}...")

        def writer_for(variant):
            key = variant if variant else "_BLANK"
            if key not in writers:
                f = open(os.path.join(a.out_dir, f"standinit_{key}.csv"),
                         "w", newline="")
                w = csv.writer(f)
                w.writerow(header)
                handles[key] = f
                writers[key] = w
                counts[key] = 0
            return key

        for row in reader:
            if len(row) <= vi:
                continue
            variant = row[vi].strip().upper()
            key = writer_for(variant)
            writers[key].writerow(row)
            counts[key] += 1

    for f in handles.values():
        f.close()

    with open(os.path.join(a.out_dir, "counts.tsv"), "w") as out:
        out.write("variant\tn_stands\tn_batches\n")
        total = 0
        for key in sorted(counts, key=lambda k: -counts[k]):
            n = counts[key]
            total += n
            nb = math.ceil(n / a.batch_size)
            out.write(f"{key}\t{n}\t{nb}\n")
        out.write(f"TOTAL\t{total}\t-\n")

    print(f"split {sum(counts.values())} stands into {len(counts)} variant files "
          f"-> {a.out_dir}")
    print("see counts.tsv for per-variant batch counts")


if __name__ == "__main__":
    main()
