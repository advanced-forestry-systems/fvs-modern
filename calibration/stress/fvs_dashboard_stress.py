#!/usr/bin/env python3
"""fvs_dashboard_stress.py -- invariant stress test of the live FVS dashboard output.

Scans every state series JSON for the three national FVS engines and checks the
invariants the projections must satisfy. Reports every violation (state, engine,
detail). Read-only.

Checks:
  finite/positive  : all pts values finite and >= 0
  band             : where a point has lo/hi, lo <= value <= hi
  monotone reserve : reserve (no harvest) non-decreasing (no harvest, growth only)
  scenario order   : reserve >= extensive >= harvest >= intensive (each year)
  managed<=reserve : every managed bucket <= reserve
  engine order     : gompit <= default at the final year (density-dep. mortality caps)
  flux sign        : harvest_c_yr >= 0
  anchor           : FIA-anchored states reserve@anchor ~ fia tg_agc (+-2%)

Usage: python3 fvs_dashboard_stress.py <api_dir>
"""
import json, os, sys, glob

api = sys.argv[1]
fia = json.load(open(os.path.join(api, "fia.json")))
MODELS = ["fvs_national_default_v1", "fvs_national_calibrated_v1",
          "fvs_national_gompit_v1"]
ENG = {m: m.split("_")[2] for m in MODELS}
ORDER = ["reserve (no harvest)", "managed (extensive)",
         "managed (harvest)", "managed (intensive)"]
METRIC = "agc_live_total"
TOL = 1e-6
viol = []


def add(st, kind, detail):
    viol.append((st, kind, detail))


def series_pts(node, model):
    for s in node:
        if s.get("model") == model:
            return {p[0]: p for p in s["pts"]}
    return {}


nstates = 0
for f in sorted(glob.glob(os.path.join(api, "series", "*.json"))):
    st = os.path.basename(f)[:-5]
    if st == "US":
        continue
    ser = json.load(open(f))
    agc = ser.get(METRIC, {})
    has_fvs = any(s.get("model") in MODELS
                  for bk in agc for s in agc[bk])
    if not has_fvs:
        continue
    nstates += 1
    # per engine
    for model in MODELS:
        eng = ENG[model]
        buckets = {bk: series_pts(agc.get(bk, []), model) for bk in ORDER}
        # finite/positive + band
        for bk, pts in buckets.items():
            for y, p in pts.items():
                v = p[1]
                if v != v or v in (float("inf"), float("-inf")):
                    add(st, "nonfinite", f"{eng}/{bk}@{y}={v}")
                elif v < -TOL:
                    add(st, "negative", f"{eng}/{bk}@{y}={v}")
                if len(p) >= 4:
                    lo, hi = p[2], p[3]
                    if not (lo - 1e-3 <= v <= hi + 1e-3):
                        add(st, "band", f"{eng}/{bk}@{y} lo{lo}<= {v} <=hi{hi}")
        # monotone reserve (gompit may decline in late succession by design:
        # density-dependent mortality caps over-accumulation)
        r = buckets["reserve (no harvest)"]
        ys = sorted(r)
        for i in range(1, len(ys)):
            if eng == "gompit":
                break
            if r[ys[i]][1] < r[ys[i-1]][1] - 0.5:
                add(st, "reserve_decrease", f"{eng} {ys[i-1]}->{ys[i]}: "
                    f"{r[ys[i-1]][1]:.0f}->{r[ys[i]][1]:.0f}")
                break
        # scenario ordering each year
        for y in ys:
            seq = [buckets[b].get(y, [None, None])[1] for b in ORDER]
            seq = [(b, v) for b, v in zip(ORDER, seq) if v is not None]
            for i in range(1, len(seq)):
                if seq[i][1] > seq[i-1][1] + max(0.5, 0.02*seq[i-1][1]):
                    add(st, "scenario_order",
                        f"{eng}@{y} {seq[i-1][0]}={seq[i-1][1]:.0f} < "
                        f"{seq[i][0]}={seq[i][1]:.0f}")
                    break
    # engine order at final year: gompit <= default
    rd = series_pts(agc.get("reserve (no harvest)", []), "fvs_national_default_v1")
    rg = series_pts(agc.get("reserve (no harvest)", []), "fvs_national_gompit_v1")
    if rd and rg:
        fy = max(set(rd) & set(rg))
        if rg[fy][1] > rd[fy][1] * 1.05:
            add(st, "engine_order", f"gompit {rg[fy][1]:.0f} > default "
                f"{rd[fy][1]:.0f} @{fy}")
    # flux sign
    for bk, node in ser.get("harvest_c_yr", {}).items():
        for s in node:
            if s.get("model") in MODELS:
                for p in s["pts"]:
                    if p[1] < -TOL:
                        add(st, "neg_flux", f"{ENG[s['model']]}/{bk}@{p[0]}={p[1]}")
    # anchor reconciliation
    anc = fia.get(st, {}).get("tg_agc")
    if anc:
        rd = series_pts(agc.get("reserve (no harvest)", []),
                        "fvs_national_default_v1")
        y0 = min(rd) if rd else None
        if y0 and abs(rd[y0][1] - anc) > 0.02 * anc:
            add(st, "anchor", f"reserve@{y0}={rd[y0][1]:.1f} vs FIA {anc}")

print(f"scanned {nstates} states with FVS engines")
print(f"TOTAL violations: {len(viol)}")
from collections import Counter
for kind, n in Counter(v[1] for v in viol).most_common():
    print(f"  {kind}: {n}")
print("\nexamples:")
for v in viol[:25]:
    print(f"  [{v[0]}] {v[1]}: {v[2]}")
