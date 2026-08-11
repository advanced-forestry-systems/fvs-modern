#!/usr/bin/env python3
"""fvs_managed_scenario.py -- derive the FVS "managed (harvest)" trajectory.

Takes the FVS reserve (no-harvest) per-state density series and applies a
data-driven annual removal regime from the conus_hcs harvest rasters
(state_harvest_rates.csv, sampled at FIA plot locations) plus an annualized
disturbance loss, producing:

  * managed (harvest) density trajectory (agc_live_total, agb_dry)
  * harvest_c_yr carbon flux density (Mg C/ha/yr)

Model (density space, 5-yr steps, transparent and repeatable):
  growth5  = d_reserve[t] - d_reserve[t-5]           (engine's own increment)
  harvest5 = 5 * h * d_managed[t-5]                  (h = expected BA frac/yr)
  dist5    = 5 * a * d_managed[t-5]                  (a = annualized disturbance)
  d_managed[t] = max(0, d_managed[t-5] + growth5 - harvest5 - dist5)
  harvest_flux[t] = h * d_managed[t]
Disturbance annualization: a = 1-(1-p_dist)^(1/W), W default 20 yr (the
p_disturbance raster is a multi-year susceptibility probability, not annual).

Output per engine: <out>/managed_<ST>.csv with the managed densities + flux,
ready for the merge to inject as the "managed (harvest)" bucket + harvest_c_yr.

Usage:
  python3 fvs_managed_scenario.py --series-dir perseus_series_default_v2 \
     --rates state_harvest_rates.csv --start 2030 --window 20 --out managed_default
"""
from __future__ import annotations
import argparse, csv, glob, os
from collections import defaultdict

C_FRACTION = 0.47


def _f(x, default=0.0):
    try:
        v = float(x)
        return v if v == v else default        # NaN -> default
    except (TypeError, ValueError):
        return default


def load_rates(path):
    out = {}
    for r in csv.DictReader(open(path)):
        out[r["ST"]] = (_f(r.get("harvest_frac_yr")),
                        _f(r.get("disturbance_frac_yr")))
    return out


def load_reserve(path, start):
    """metric -> {year: density} for the reserve bucket, years >= start."""
    d = defaultdict(dict)
    for r in csv.DictReader(open(path)):
        if r["mgmt"] != "reserve (no harvest)":
            continue
        y = int(r["year"])
        if y >= start:
            d[r["metric"]][y] = float(r["value"])
    return d


def managed_path(metric, res, h, a):
    """discrete 5-yr harvest+disturbance decrement on a reserve density path."""
    yrs = sorted(res)
    dm = {yrs[0]: res[yrs[0]]}
    flux = {yrs[0]: h * res[yrs[0]]}
    for i in range(1, len(yrs)):
        t, t0 = yrs[i], yrs[i - 1]
        step = t - t0
        growth = res[t] - res[t0]
        removed = step * (h + a) * dm[t0]
        dm[t] = max(0.0, dm[t0] + growth - removed)
        flux[t] = h * dm[t]
    return dm, flux


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--series-dir", required=True)
    ap.add_argument("--rates", required=True)
    ap.add_argument("--start", type=int, default=2030)
    ap.add_argument("--window", type=float, default=20.0,
                    help="years the p_disturbance prob spans (annualize)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rates = load_rates(a.rates)

    summ = []
    for f in sorted(glob.glob(os.path.join(a.series_dir,
                                           "ycx_*_state_series.csv"))):
        st = os.path.basename(f).split("_")[1]
        if st not in rates:
            continue
        h, pdist = rates[st]
        adist = 1.0 - (1.0 - min(pdist, 0.99)) ** (1.0 / a.window)
        res = load_reserve(f, a.start)
        if "agc_live_total" not in res:
            continue
        rows = []
        flux_agc = None
        for metric in ("agc_live_total", "agb_dry"):
            if metric not in res:
                continue
            dm, flux = managed_path(metric, res[metric], h, adist)
            for y in sorted(dm):
                rows.append({"metric": metric, "mgmt": "managed (harvest)",
                             "year": y, "value": round(dm[y], 4)})
            if metric == "agc_live_total":
                flux_agc = flux
        if flux_agc:
            for y in sorted(flux_agc):
                rows.append({"metric": "harvest_c_yr", "mgmt": "managed (harvest)",
                             "year": y, "value": round(flux_agc[y], 5)})
        with open(os.path.join(a.out, f"managed_{st}.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["metric", "mgmt", "year", "value"])
            w.writeheader(); w.writerows(rows)
        r = res["agc_live_total"]
        dm, _ = managed_path("agc", r, h, adist)
        summ.append((st, h, adist, r[max(r)], dm[max(dm)]))

    print(f"{len(summ)} states  (h=harvest/yr, a=disturb/yr, reserve vs managed "
          f"agc density {a.start+95 if False else max(r)})")
    for st, h, adist, rv, mv in summ:
        if st in ("ME", "OR", "GA", "ID", "MN", "WA", "CA", "MS"):
            print(f"  {st}: h={h:.3f} a={adist:.3f}  agc 2125 reserve {rv:.0f} "
                  f"-> managed {mv:.0f} Mg C/ha ({100*mv/rv:.0f}%)")
    print(f"wrote managed_<ST>.csv to {a.out}")


if __name__ == "__main__":
    main()
