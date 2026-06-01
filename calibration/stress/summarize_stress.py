#!/usr/bin/env python3
"""Stage 3: aggregate per-task ledgers into one stress-test report.

Reads every ledger_*.json in the output dir and produces:
  - stress_summary.csv : one row per task (variant, batch, counts, failures, sec)
  - stress_failures.csv: every recorded failure (variant, stand_cn, stage, detail)
  - prints overall totals and the failure rate -- the headline stress-test result.

Usage:
  python summarize_stress.py --output-dir /fs/scratch/PUOM0008/crsfaaron/fvs_stress/out
"""
from __future__ import annotations

import argparse
import csv
import glob
import json
import os


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--output-dir", required=True)
    a = ap.parse_args()

    ledgers = sorted(glob.glob(os.path.join(a.output_dir, "ledger_*.json")))
    if not ledgers:
        raise SystemExit(f"no ledger_*.json in {a.output_dir}")

    tot_proj = tot_fail = tot_rows = tot_stands = 0
    per_variant: dict[str, list[int]] = {}
    sum_path = os.path.join(a.output_dir, "stress_summary.csv")
    fail_path = os.path.join(a.output_dir, "stress_failures.csv")

    with open(sum_path, "w", newline="") as sf, open(fail_path, "w", newline="") as ff:
        sw = csv.writer(sf)
        sw.writerow(["task_id", "variant", "batch_id", "n_stands_in_batch",
                     "n_stands_with_trees", "n_stands_projected", "n_output_rows",
                     "n_failures", "elapsed_sec"])
        fw = csv.writer(ff)
        fw.writerow(["variant", "batch_id", "stand_cn", "config", "stage", "detail"])
        for lp in ledgers:
            d = json.load(open(lp))
            v = d["variant"]
            sw.writerow([d["task_id"], v, d["batch_id"], d["n_stands_in_batch"],
                         d["n_stands_with_trees"], d["n_stands_projected"],
                         d["n_output_rows"], d["n_failures"], d["elapsed_sec"]])
            tot_proj += d["n_stands_projected"]
            tot_fail += d["n_failures"]
            tot_rows += d["n_output_rows"]
            tot_stands += d["n_stands_with_trees"]
            agg = per_variant.setdefault(v, [0, 0])
            agg[0] += d["n_stands_projected"]
            agg[1] += d["n_failures"]
            for fdet in d.get("failures", []):
                fw.writerow([v, d["batch_id"], fdet.get("stand_cn"),
                             fdet.get("config"), fdet.get("stage"), fdet.get("detail")])

    print(f"tasks: {len(ledgers)}")
    print(f"stands with trees: {tot_stands}   stands projected: {tot_proj}")
    print(f"output rows: {tot_rows}   failures: {tot_fail}")
    denom = tot_proj * 2 or 1   # 2 configs per projected stand
    print(f"projection failure rate: {100 * tot_fail / denom:.4f}%  "
          f"(failures / (projected stands x 2 configs))")
    print("\nper-variant (projected / failures):")
    for v in sorted(per_variant):
        p, f = per_variant[v]
        print(f"  {v:<6} {p:>8}  {f}")
    print(f"\nwrote {sum_path}\nwrote {fail_path}")


if __name__ == "__main__":
    main()
