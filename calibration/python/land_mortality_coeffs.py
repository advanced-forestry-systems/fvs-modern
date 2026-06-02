#!/usr/bin/env python3
"""Land Greg's re-fit mortality coefficients into each variant config as a
`categories_conus.mortality` block (per the FVS_CONUS_INTEGRATION_PLAN Option B
schema). SPCD-keyed coefficients are mapped to the FVS species slot via
categories.species_definitions.FIAJSP. Variants/species without a fitted row
fall back to the existing per-variant mortality (handled by config_loader's
categories fallback), so partial coverage is safe.

Model carried in the block (annual-hazard gompit, mortality framing):
    H_annual = exp(b0 + b1*(cr+0.01)^b2 + b3*cch^b4)
    P_survive(T) = exp(-H_annual * T)
cch = crown closure at tree tip.

Usage:
  python land_mortality_coeffs.py --coeffs greg_mortality_coefficients.csv \
      --config-dir config/calibrated [--variants ne,ls,cs]
"""
from __future__ import annotations
import argparse, csv, json, glob, os
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("--coeffs", required=True)
ap.add_argument("--config-dir", required=True)
ap.add_argument("--variants", default="")   # comma list, empty = all
a = ap.parse_args()

# load per-species coefficients keyed by SPCD
coef = {}
with open(a.coeffs) as fh:
    for r in csv.DictReader(fh):
        coef[int(r["SPCD"])] = {k: float(r[k]) for k in ("b0","b1","b2","b3","b4")} | {
            "fallback": str(r.get("fallback","")).upper() == "TRUE"}

want = {v.strip().lower() for v in a.variants.split(",") if v.strip()}
cfgs = [p for p in sorted(glob.glob(os.path.join(a.config_dir, "*.json")))
        if not p.endswith("_draws.json")]

changed = []
for p in cfgs:
    variant = Path(p).stem
    if want and variant not in want:
        continue
    d = json.loads(Path(p).read_text())
    fia = (d.get("categories", {}).get("species_definitions", {}) or {}).get("FIAJSP")
    maxsp = d.get("maxsp")
    if not isinstance(fia, list) or not maxsp:
        continue
    # build length-maxsp coefficient arrays on FVS species slots
    arrs = {k: [None]*maxsp for k in ("b0","b1","b2","b3","b4")}
    n_mapped = 0
    for i in range(min(maxsp, len(fia))):
        try:
            spcd = int(fia[i])
        except (TypeError, ValueError):
            continue
        c = coef.get(spcd)
        if c is None:
            continue
        for k in ("b0","b1","b2","b3","b4"):
            arrs[k][i] = round(c[k], 8)
        n_mapped += 1
    cc = d.setdefault("categories_conus", {})
    cc["mortality"] = {
        "form": "gompit_cr_cch",
        "equation": "H_annual=exp(b0+b1*(cr+0.01)^b2+b3*cch^b4); P_surv_T=exp(-H_annual*T)",
        "covariates": ["cr (crown ratio)", "cch (crown closure at tree tip)"],
        "source": "Greg Johnson CONUS mortality (Johnson/Marshall/Weiskittel 2026-05-26), re-fit per species",
        "n_species_mapped": n_mapped,
        **{k: arrs[k] for k in ("b0","b1","b2","b3","b4")},
    }
    Path(p).write_text(json.dumps(d, indent=2) + "\n")
    changed.append(f"{variant}: mapped {n_mapped}/{maxsp} species")

print(f"landed mortality block into {len(changed)} configs:")
for c in changed:
    print("  ", c)
