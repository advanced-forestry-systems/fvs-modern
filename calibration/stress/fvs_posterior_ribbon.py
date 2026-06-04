#!/usr/bin/env python3
"""fvs_posterior_ribbon.py -- apply the posterior (parameter) CI as the lo/hi
band on the calibrated FVS engine in the dashboard series.

Reads posterior_ci_all.csv (ST, year, rel_lo, rel_hi -- the parameter-uncertainty
multipliers from fvs_posterior_uncertainty.py), and for each anchored state sets
the lo/hi of every point on fvs_national_calibrated_v1's agc_live_total series
(all scenario buckets) to value*rel_lo / value*rel_hi. This makes the calibrated
engine's band reflect Bayesian parameter uncertainty (narrow), distinct from the
structural engine spread. Run AFTER the main merge.

Usage: python3 fvs_posterior_ribbon.py <repo_dir> <posterior_ci_all.csv>
"""
import csv, json, os, sys
from collections import defaultdict

repo, ci_path = sys.argv[1], sys.argv[2]
api = os.path.join(repo, "public", "api")
MODEL = "fvs_national_calibrated_v1"
METRIC = "agc_live_total"

ci = defaultdict(dict)
for r in csv.DictReader(open(ci_path)):
    ci[r["ST"]][int(r["year"])] = (float(r["rel_lo"]), float(r["rel_hi"]))


def band(state, year):
    d = ci[state]
    if year in d:
        return d[year]
    ys = sorted(d)
    if not ys:
        return (1.0, 1.0)
    yy = min(ys, key=lambda x: abs(x - year))   # nearest year
    return d[yy]


touched = 0
for st in ci:
    sp = os.path.join(api, "series", f"{st}.json")
    if not os.path.exists(sp):
        continue
    ser = json.load(open(sp))
    node = ser.get(METRIC, {})
    n = 0
    for bucket in node:
        for s in node[bucket]:
            if s.get("model") != MODEL:
                continue
            for pt in s["pts"]:
                y, v = pt[0], pt[1]
                lo, hi = band(st, y)
                if len(pt) >= 4:
                    pt[2], pt[3] = round(v * lo, 3), round(v * hi, 3)
                else:
                    pt.extend([round(v * lo, 3), round(v * hi, 3)])
                n += 1
    if n:
        json.dump(ser, open(sp, "w"), separators=(",", ":"))
        touched += 1
        print(f"  {st}: parameter band on {n} calibrated points")
print(f"applied posterior CI to {touched} states")
