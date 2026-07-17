#!/usr/bin/env python3
"""
run_conus_perseus.py  thin CONUS wrapper around the working perseus_100yr_projection driver.

The stock driver loads trees for a single state (default ME). The CONUS per-variant plot lists
span many states, so this wrapper loads FIA trees per STATECD (each state's <ST>_TREE.csv) and then
calls the driver's process_plot exactly as the Maine run does. No change to the driver itself.

Usage (one variant, one batch):
  python run_conus_perseus.py --plotlist .../plotlist_sn.csv --variant sn \
      --batch-id 1 --batch-size 5000 --configs default calibrated \
      --output-dir .../out_conus_wo1
"""
import argparse, os, sys
import pandas as pd

FIA_DIR = os.environ.get("FIA_DATA_DIR", "/fs/scratch/PUOM0008/crsfaaron/FIA")
sys.path.insert(0, "/users/PUOM0008/crsfaaron/fvs-conus/python")
from perseus_100yr_projection import (
    load_fia_trees_for_plots, process_plot, NSBECalculator,
)

# FIPS state code -> postal abbreviation (lower 48 + AK)
FIPS = {1:"AL",2:"AK",4:"AZ",5:"AR",6:"CA",8:"CO",9:"CT",10:"DE",12:"FL",13:"GA",16:"ID",
17:"IL",18:"IN",19:"IA",20:"KS",21:"KY",22:"LA",23:"ME",24:"MD",25:"MA",26:"MI",27:"MN",
28:"MS",29:"MO",30:"MT",31:"NE",32:"NV",33:"NH",34:"NJ",35:"NM",36:"NY",37:"NC",38:"ND",
39:"OH",40:"OK",41:"OR",42:"PA",44:"RI",45:"SC",46:"SD",47:"TN",48:"TX",49:"UT",50:"VT",
51:"VA",53:"WA",54:"WV",55:"WI",56:"WY"}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--plotlist", required=True)
    ap.add_argument("--variant", required=True)
    ap.add_argument("--batch-id", type=int, required=True)
    ap.add_argument("--batch-size", type=int, default=5000)
    ap.add_argument("--configs", nargs="+", default=["default", "calibrated"])
    ap.add_argument("--output-dir", required=True)
    a = ap.parse_args()

    configs = [None if c == "default" else c for c in a.configs]
    plots = pd.read_csv(a.plotlist)
    start = (a.batch_id - 1) * a.batch_size
    end = min(start + a.batch_size, len(plots))
    batch = plots.iloc[start:end].copy()
    if batch.empty:
        print(f"batch {a.batch_id}: empty slice"); return
    nsbe = NSBECalculator(os.environ.get("NSBE_ROOT", "/users/PUOM0008/crsfaaron/fvs-modern/data/NSBE"))

    rows = []
    # group by state so each state's TREE.csv is read once
    for statecd, grp in batch.groupby("STATECD"):
        st = FIPS.get(int(statecd))
        if st is None:
            print(f"skip unknown STATECD {statecd}"); continue
        cns = [str(int(float(c))) for c in grp["FIRST_PLTCN"]]
        try:
            trees = load_fia_trees_for_plots(cns, FIA_DIR, state=st)
        except FileNotFoundError as e:
            print(f"state {st}: {e}"); continue
        if trees is None or "PLT_CN" not in getattr(trees, "columns", []):
            print(f"state {st}: no tree records"); continue
        for _, prow in grp.iterrows():
            rows.extend(process_plot(prow.to_dict(), trees, nsbe,
                                     variants=[a.variant], configs=configs))
    if not rows:
        print(f"batch {a.batch_id}: no results"); return
    os.makedirs(a.output_dir, exist_ok=True)
    out = os.path.join(a.output_dir, f"conus_{a.variant}_agb_batch{a.batch_id}.csv")
    pd.DataFrame(rows).to_csv(out, index=False)
    print(f"batch {a.batch_id} variant {a.variant}: wrote {len(rows)} rows to {out}")

if __name__ == "__main__":
    main()
