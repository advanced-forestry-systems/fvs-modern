#!/usr/bin/env python3
"""fvs_strata_trends.py -- landowner / ecoregion / state carbon trends across
scenarios, with bootstrap uncertainty.

Joins the FVS per-plot reserve trajectories to:
  * ycx_membership_<ST>.csv  -> owner4 (Industrial/NIPF/State/Public-Other),
    prov_name (EPA Level III ecoregion), state
  * plt_area_treemap.csv     -> per-plot CONUS area (TreeMap), for true totals
  * state_harvest_rates.csv  -> per-state harvest h + annualized disturbance a,
    to derive the managed (harvest) per-plot path.

For each (scale in {owner, ecoregion, state}, stratum key, engine, scenario,
year) it computes the area-weighted total live carbon (Tg C) and a bootstrap
95% CI (resampling plots within the stratum, B draws). Density mean is also
reported. Scenarios: reserve (no harvest) and managed (harvest+disturbance).

Output: <out>/strata_trends.csv  (scale,key,engine,scenario,year,
        total_TgC,total_lo,total_hi,mean_MgC_ha,n_plots)

Usage:
  python3 fvs_strata_trends.py --campaign out_fvs_v2 --config default \
     --membership-dir <dir> --areas plt_area_treemap.csv \
     --rates state_harvest_rates.csv --out strata_trends --boot 200
"""
from __future__ import annotations
import argparse, glob, os
import numpy as np, pandas as pd

TONS_AC_TO_MGHA = 2.241702
C_FRACTION = 0.47
YEARS = [2030, 2055, 2080, 2105, 2125]        # PROJ_YEAR = year-2025
PROJ = {y: y - 2025 for y in YEARS}


def load_membership(d):
    fs = glob.glob(os.path.join(d, "ycx_membership_*.csv"))
    m = pd.concat([pd.read_csv(f, dtype={"PLT_CN": str}) for f in fs],
                  ignore_index=True)
    m["PLT_CN"] = m["PLT_CN"].str.replace(r"\.0$", "", regex=True)
    m["owner4"] = m["owner4"].fillna("Unknown")
    m["prov_name"] = m["prov_name"].fillna("Unknown")
    return m[["PLT_CN", "owner4", "prov_name"]].drop_duplicates("PLT_CN")


def managed_factor(reserve_by_year, h, a):
    """per-plot managed density at each year from its reserve path (5-yr steps)."""
    yrs = sorted(reserve_by_year)
    dm = {yrs[0]: reserve_by_year[yrs[0]]}
    for i in range(1, len(yrs)):
        t, t0 = yrs[i], yrs[i - 1]
        g = reserve_by_year[t] - reserve_by_year[t0]
        dm[t] = max(0.0, dm[t0] + g - (t - t0) * (h + a) * dm[t0])
    return dm


def boot_total(dens, area, B, rng):
    """area-weighted total (Tg C) + 95% CI by resampling plots."""
    w = dens * area / 1e6
    tot = float(w.sum())
    if len(w) < 3 or B <= 0:
        return tot, tot, tot
    n = len(w)
    idx = rng.integers(0, n, size=(B, n))
    draws = w.values[idx].sum(axis=1)
    return tot, float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--campaign", required=True)
    ap.add_argument("--config", default="default")
    ap.add_argument("--membership-dir", required=True)
    ap.add_argument("--areas", required=True)
    ap.add_argument("--rates", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--boot", type=int, default=200)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(7)

    mem = load_membership(a.membership_dir)
    area = pd.read_csv(a.areas, dtype={"PLT_CN": str})
    area["PLT_CN"] = area["PLT_CN"].str.replace(r"\.0$", "", regex=True)
    rates = pd.read_csv(a.rates).set_index("ST")

    # per-plot reserve density at all needed years
    keep = set(PROJ.values())
    rows = []
    for f in glob.glob(os.path.join(a.campaign, "conus_*.csv")):
        d = pd.read_csv(f, usecols=["STAND_CN", "STATE", "CONFIG",
                                    "PROJ_YEAR", "AGB_TONS_AC"])
        d = d[(d.CONFIG == a.config) & (d.PROJ_YEAR.isin(keep))]
        if len(d):
            rows.append(d)
    camp = pd.concat(rows, ignore_index=True)
    camp["PLT_CN"] = camp["STAND_CN"].astype(str).str.replace(r"\.0$", "",
                                                              regex=True)
    camp["agc"] = camp.AGB_TONS_AC * TONS_AC_TO_MGHA * C_FRACTION
    wide = camp.pivot_table(index=["PLT_CN", "STATE"], columns="PROJ_YEAR",
                            values="agc", aggfunc="mean")
    wide = wide.dropna().reset_index()
    wide = wide.merge(area[["PLT_CN", "area_ha"]], on="PLT_CN", how="inner")
    wide = wide.merge(mem, on="PLT_CN", how="left")
    wide["owner4"] = wide["owner4"].fillna("Unknown")
    wide["prov_name"] = wide["prov_name"].fillna("Unknown")
    print(f"{len(wide)} plots with full year set + area + strata")

    # managed per-plot densities
    hcol = {st: (rates.loc[st, "harvest_frac_yr"] if st in rates.index else 0.0)
            for st in wide.STATE.unique()}
    pdcol = {st: (rates.loc[st, "disturbance_frac_yr"] if st in rates.index
                  else 0.0) for st in wide.STATE.unique()}
    pys = sorted(keep)
    man = np.zeros((len(wide), len(pys)))
    arr = wide[pys].values
    for i, st in enumerate(wide.STATE.values):
        h = float(hcol.get(st) or 0.0)
        pdz = float(pdcol.get(st) or 0.0)
        adist = 1.0 - (1.0 - min(pdz, 0.99)) ** (1.0 / 20.0)
        rbyy = {py: arr[i, j] for j, py in enumerate(pys)}
        dm = managed_factor(rbyy, h, adist)
        man[i, :] = [dm[py] for py in pys]
    man = pd.DataFrame(man, columns=pys, index=wide.index)

    out = []
    scales = {"owner": "owner4", "ecoregion": "prov_name", "state": "STATE"}
    for scale, col in scales.items():
        for key, g in wide.groupby(col):
            mg = man.loc[g.index]
            for scen, src in (("reserve (no harvest)", g),
                              ("managed (harvest)", mg)):
                for y, py in PROJ.items():
                    dens = src[py] if scen.startswith("reserve") else mg[py]
                    tot, lo, hi = boot_total(dens, g["area_ha"], a.boot, rng)
                    out.append({"scale": scale, "key": str(key),
                                "engine": a.config, "scenario": scen, "year": y,
                                "total_TgC": round(tot, 3),
                                "total_lo": round(lo, 3), "total_hi": round(hi, 3),
                                "mean_MgC_ha": round(float(dens.mean()), 2),
                                "n_plots": int(len(g))})
    res = pd.DataFrame(out)
    res.to_csv(os.path.join(a.out, f"strata_trends_{a.config}.csv"), index=False)
    print(f"wrote strata_trends_{a.config}.csv  "
          f"({res.scale.eq('owner').sum()} owner, "
          f"{res.scale.eq('ecoregion').sum()} ecoregion, "
          f"{res.scale.eq('state').sum()} state rows)")
    own = res[(res.scale == "owner") & (res.year.isin([2030, 2125]))]
    print(own.pivot_table(index=["key", "scenario"], columns="year",
                          values="total_TgC").to_string())


if __name__ == "__main__":
    main()
