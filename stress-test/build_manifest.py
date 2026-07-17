#!/usr/bin/env python3
"""Stage 1: build the SLURM array manifest from counts.tsv.

Emits manifest.tsv with one line per array task:
    array_idx <TAB> variant <TAB> batch_id <TAB> batch_size

Excludes the _BLANK and TOTAL rows. Prints the array index range to use in the
sbatch --array directive (0-based, inclusive).

Usage:
  python build_manifest.py \
      --counts   /fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_by_variant/counts.tsv \
      --manifest /fs/scratch/PUOM0008/crsfaaron/fvs_stress/manifest.tsv \
      --batch-size 5000
"""
from __future__ import annotations

import argparse
import math


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--batch-size", type=int, default=5000)
    a = ap.parse_args()

    idx = 0
    with open(a.counts) as fh, open(a.manifest, "w") as out:
        header = fh.readline()
        for line in fh:
            variant, n_stands, _ = line.rstrip("\n").split("\t")
            if variant in ("TOTAL", "_BLANK"):
                continue
            n = int(n_stands)
            nb = math.ceil(n / a.batch_size)
            for b in range(nb):
                out.write(f"{idx}\t{variant}\t{b}\t{a.batch_size}\n")
                idx += 1

    print(f"wrote {idx} tasks -> {a.manifest}")
    print(f"sbatch --array=0-{idx - 1}%<throttle>   (e.g. %16 to backfill politely)")


if __name__ == "__main__":
    main()
