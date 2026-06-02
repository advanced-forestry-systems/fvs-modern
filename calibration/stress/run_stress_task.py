#!/usr/bin/env python3
"""Stage 2: one SLURM array task of the full-CONUS stress test.

Looks up (variant, batch_id, batch_size) for this task from the manifest using
SLURM_ARRAY_TASK_ID, reads ONLY that variant's pre-split stand-init, and projects
the batch with the perseus engine (default + calibrated), exactly as
conus_100yr_projection.py does. Adds a per-task ledger so the full run is a real
stress test: it records stand/row counts, per-stand FVS failures (with STAND_CN
and config), and elapsed time to <output-dir>/ledger_<variant>_b<batch>.json.

The engine and projection settings mirror conus_100yr_projection.py so results
are identical; the only changes are (a) read the pre-split file, (b) collect
failures instead of only logging them.

Usage (inside the sbatch array task):
  python run_stress_task.py \
     --manifest    /fs/scratch/PUOM0008/crsfaaron/fvs_stress/manifest.tsv \
     --standinit-dir /fs/scratch/PUOM0008/crsfaaron/fvs_stress/standinit_by_variant \
     --output-dir  /fs/scratch/PUOM0008/crsfaaron/fvs_stress/out \
     --task-id     $SLURM_ARRAY_TASK_ID
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time

import pandas as pd

PROJECT_ROOT = os.environ.get("FVS_PROJECT_ROOT", os.path.expanduser("~/fvs-modern"))
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.join(PROJECT_ROOT, "calibration", "python"))
import perseus_100yr_projection as P  # noqa: E402

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("stress")

FIA_DIR = os.environ.get("FVS_FIA_DATA_DIR", "/fs/scratch/PUOM0008/crsfaaron/FIA")
FIPS = {1: "AL", 2: "AK", 4: "AZ", 5: "AR", 6: "CA", 8: "CO", 9: "CT", 10: "DE",
        12: "FL", 13: "GA", 16: "ID", 17: "IL", 18: "IN", 19: "IA", 20: "KS",
        21: "KY", 22: "LA", 23: "ME", 24: "MD", 25: "MA", 26: "MI", 27: "MN",
        28: "MS", 29: "MO", 30: "MT", 31: "NE", 32: "NV", 33: "NH", 34: "NJ",
        35: "NM", 36: "NY", 37: "NC", 38: "ND", 39: "OH", 40: "OK", 41: "OR",
        42: "PA", 44: "RI", 45: "SC", 46: "SD", 47: "TN", 48: "TX", 49: "UT",
        50: "VT", 51: "VA", 53: "WA", 54: "WV", 55: "WI", 56: "WY"}


def lookup_task(manifest, task_id):
    with open(manifest) as fh:
        for line in fh:
            idx, variant, batch_id, batch_size = line.rstrip("\n").split("\t")
            if int(idx) == task_id:
                return variant, int(batch_id), int(batch_size)
    raise SystemExit(f"task_id {task_id} not found in {manifest}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--standinit-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--task-id", type=int,
                    default=int(os.environ.get("SLURM_ARRAY_TASK_ID", "-1")))
    a = ap.parse_args()
    if a.task_id < 0:
        sys.exit("no task id (set --task-id or SLURM_ARRAY_TASK_ID)")
    os.makedirs(a.output_dir, exist_ok=True)

    variant, batch_id, batch_size = lookup_task(a.manifest, a.task_id)
    t0 = time.time()
    log.info(f"task {a.task_id}: variant={variant} batch={batch_id} size={batch_size}")

    sfile = os.path.join(a.standinit_dir, f"standinit_{variant}.csv")
    si_all = pd.read_csv(sfile, low_memory=False)
    si = si_all.iloc[batch_id * batch_size:(batch_id + 1) * batch_size].reset_index(drop=True)
    log.info(f"  {len(si)} stands in this batch (of {len(si_all)} for {variant})")

    nsbe = P.NSBECalculator(P.NSBE_ROOT)
    rows = []
    failures = []          # the stress-test ledger of FVS failures
    n_stands_attempted = 0
    n_stands_with_trees = 0

    for state_fips, grp in si.groupby("STATE"):
        try:
            state = FIPS[int(float(state_fips))]
        except (KeyError, ValueError, TypeError):
            failures.append({"stand_cn": None, "config": None,
                             "stage": "state_map", "detail": f"unknown STATE={state_fips}",
                             "n_stands": int(len(grp))})
            continue
        plt_cns = [str(int(float(x))) for x in grp["STAND_CN"]]
        try:
            trees = P.load_fia_trees_for_plots(plt_cns, FIA_DIR, state)
        except Exception as e:
            failures.append({"stand_cn": None, "config": None, "stage": "tree_load",
                             "detail": f"{state}: {e}", "n_stands": int(len(grp))})
            continue
        if trees.empty or "PLT_CN" not in trees.columns:
            continue
        tcn = trees["PLT_CN"].apply(lambda x: str(int(float(x))) if pd.notna(x) else "")
        for _, stand in grp.iterrows():
            cn = str(int(float(stand["STAND_CN"])))
            st = trees.loc[tcn == cn]
            if st.empty:
                continue
            n_stands_with_trees += 1
            sid = f"S{cn}"
            inv_year = int(float(stand.get("INV_YEAR") or 2010))
            plot_data = {
                "INVYR": inv_year,
                "LAT": stand.get("LATITUDE"), "LON": stand.get("LONGITUDE"),
                "ELEV": stand.get("ELEVFT") or 500, "SLOPE": stand.get("SLOPE") or 10,
                "ASPECT": stand.get("ASPECT") or 180, "STDAGE": stand.get("AGE") or 50,
            }
            try:
                sdf = P.build_fvs_standinit(plot_data, sid, variant.lower())
                tdf = P.build_fvs_treeinit(st, sid)
            except Exception as e:
                failures.append({"stand_cn": cn, "config": None, "stage": "build",
                                 "detail": str(e)})
                continue
            n_stands_attempted += 1
            for cfg, label in ((None, "default"), ("calibrated", "calibrated")):
                try:
                    fr = P.run_fvs_projection(sdf, tdf, sid, variant.lower(),
                                              config_version=cfg,
                                              num_cycles=20, cycle_length=5)
                    for cy, tl in sorted(fr["treelists"].items()):
                        py = cy - inv_year
                        if py < 0:
                            continue
                        agb = P.compute_plot_agb(tl, nsbe)
                        rows.append({"STAND_CN": cn, "STATE": state, "YEAR": cy,
                                     "PROJ_YEAR": py, "VARIANT": variant.upper(),
                                     "CONFIG": label, "AGB_TONS_AC": round(float(agb), 4)})
                except Exception as e:
                    failures.append({"stand_cn": cn, "config": label,
                                     "stage": "project", "detail": str(e)})

    tag = f"{variant.lower()}_b{batch_id}"
    out_csv = os.path.join(a.output_dir, f"conus_{tag}.csv")
    pd.DataFrame(rows).to_csv(out_csv, index=False)

    ledger = {
        "task_id": a.task_id, "variant": variant.upper(), "batch_id": batch_id,
        "batch_size": batch_size, "n_stands_in_batch": int(len(si)),
        "n_stands_with_trees": n_stands_with_trees,
        "n_stands_projected": n_stands_attempted,
        "n_output_rows": len(rows), "n_failures": len(failures),
        "elapsed_sec": round(time.time() - t0, 1),
        "failures": failures[:500],   # cap stored detail; count is exact above
        "output_csv": out_csv,
    }
    led_path = os.path.join(a.output_dir, f"ledger_{tag}.json")
    with open(led_path, "w") as f:
        json.dump(ledger, f, indent=2)
    log.info(f"task {a.task_id} done: {n_stands_attempted} projected, "
             f"{len(failures)} failures, {len(rows)} rows, "
             f"{ledger['elapsed_sec']}s -> {out_csv}")


if __name__ == "__main__":
    main()
