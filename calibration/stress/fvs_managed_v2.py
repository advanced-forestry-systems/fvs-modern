#!/usr/bin/env python3
"""fvs_managed_v2.py -- plantation-aware managed scenarios (reserve / extensive /
intensive), per Aaron's refinement: intensive management applies ONLY to
plantation plots (FIADB STDORGCD = 1); natural stands get extensive (partial)
harvest.

Per-plot, per-state. From the empirically-sampled blended harvest rate h (the
conus_hcs expected-BA-removed, which already mixes partial + clearcut over the
landscape) and the plantation AREA fraction p, we split into an extensive rate
h_ext and an intensive rate h_int = k * h_ext (k ~ intensive/extensive ratio,
default 1.9, matching the YC engine's intensive scaling), calibrated so the
plantation-weighted blend reproduces h:
    h = (1-p)*h_ext + p*h_int   ->   h_ext = h / ((1-p) + p*k),  h_int = k*h_ext
This conserves the observed statewide removal while confining the heavier regime
to plantations.

Scenario buckets (live AGC / AGB density, Mg/ha):
  managed (conservation)       : every managed plot gets h_ext (light, all-partial)
  managed (harvest)         : plantation -> h_int, natural -> h_ext (REALISTIC,
                              reproduces the empirical statewide rate)
  managed (intensive)       : every managed plot gets h_int (heavy bound)
plus harvest_c_yr for the realistic bucket. Reserve is the untouched campaign
trajectory (already on the dashboard).

Output <out>/managed_<ST>.csv: metric, mgmt, year, value.

Usage:
  python3 fvs_managed_v2.py --campaign out_fvs_v2 --config default \
     --plantation plt_plantation.csv --rates state_harvest_rates.csv \
     --start 2030 --k 1.9 --window 20 --out managed_default
"""
from __future__ import annotations
import argparse, glob, os
import numpy as np, pandas as pd

TONS_AC_TO_MGHA = 2.241702
C_FRACTION = 0.47


def managed_path(reserve_by_year, h, a):
    yrs = sorted(reserve_by_year)
    dm = {yrs[0]: reserve_by_year[yrs[0]]}
    flux = {yrs[0]: h * reserve_by_year[yrs[0]]}
    for i in range(1, len(yrs)):
        t, t0 = yrs[i], yrs[i - 1]
        g = reserve_by_year[t] - reserve_by_year[t0]
        dm[t] = max(0.0, dm[t0] + g - (t - t0) * (h + a) * dm[t0])
        flux[t] = h * dm[t]
    return dm, flux


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--campaign", required=True)
    ap.add_argument("--config", default="default")
    ap.add_argument("--plantation", required=True)
    ap.add_argument("--rates", required=True)
    ap.add_argument("--start", type=int, default=2030)
    ap.add_argument("--k", type=float, default=1.9, help="intensive/extensive ratio")
    ap.add_argument("--window", type=float, default=20.0)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    plant = pd.read_csv(a.plantation, dtype={"PLT_CN": str})
    plant["PLT_CN"] = plant.PLT_CN.str.replace(r"\.0$", "", regex=True)
    pmap = dict(zip(plant.PLT_CN, plant.plantation))
    rates = pd.read_csv(a.rates).set_index("ST")

    # per-plot reserve densities (agc + agb) at all years
    rows = []
    for f in glob.glob(os.path.join(a.campaign, "conus_*.csv")):
        d = pd.read_csv(f, usecols=["STAND_CN", "STATE", "CONFIG",
                                    "PROJ_YEAR", "AGB_TONS_AC"])
        d = d[d.CONFIG == a.config]
        if len(d):
            rows.append(d)
    camp = pd.concat(rows, ignore_index=True)
    camp["PLT_CN"] = camp.STAND_CN.astype(str).str.replace(r"\.0$", "", regex=True)
    camp["year"] = 2025 + camp.PROJ_YEAR.astype(int)
    camp = camp[camp.year >= a.start]
    camp["plantation"] = camp.PLT_CN.map(pmap).fillna(0).astype(int)

    out_rows = []
    for ST, g in camp.groupby("STATE"):
        if ST not in rates.index:
            continue
        h = float(rates.loc[ST, "harvest_frac_yr"] or 0.0)
        pdz = float(rates.loc[ST, "disturbance_frac_yr"] or 0.0)
        if pdz != pdz:
            pdz = 0.0
        a_dist = 1.0 - (1.0 - min(pdz, 0.99)) ** (1.0 / a.window)
        p = float(g.drop_duplicates("PLT_CN").plantation.mean())   # plantation frac
        denom = (1 - p) + p * a.k
        h_ext = h / denom if denom > 0 else h
        h_int = a.k * h_ext
        # per-plot reserve path (AGB tons/ac)
        for metric, conv in (("agc_live_total", TONS_AC_TO_MGHA * C_FRACTION),
                             ("agb_dry", TONS_AC_TO_MGHA)):
            # build per-plot year->density, aggregate managed per scenario
            ext, real, intn, fluxd = {}, {}, {}, {}
            for cn, gp in g.groupby("PLT_CN"):
                ry = {int(y): float(v) * conv
                      for y, v in zip(gp.year, gp.AGB_TONS_AC)}
                if len(ry) < 2:
                    continue
                is_pl = pmap.get(cn, 0)
                de, _ = managed_path(ry, h_ext, a_dist)
                dr, fr = managed_path(ry, h_int if is_pl else h_ext, a_dist)
                di, _ = managed_path(ry, h_int, a_dist)
                for y in de:
                    ext.setdefault(y, []).append(de[y])
                    real.setdefault(y, []).append(dr[y])
                    intn.setdefault(y, []).append(di[y])
                    if metric == "agc_live_total":
                        fluxd.setdefault(y, []).append(
                            (h_int if is_pl else h_ext) * dr[y])
            for bucket, dd in (("managed (conservation)", ext),
                               ("managed (harvest)", real),
                               ("managed (intensive)", intn)):
                for y in sorted(dd):
                    out_rows.append({"ST": ST, "metric": metric, "mgmt": bucket,
                                     "year": y, "value": round(float(np.mean(dd[y])), 4)})
            if metric == "agc_live_total":
                for y in sorted(fluxd):
                    out_rows.append({"ST": ST, "metric": "harvest_c_yr",
                                     "mgmt": "managed (harvest)", "year": y,
                                     "value": round(float(np.mean(fluxd[y])), 5)})
    res = pd.DataFrame(out_rows)
    for ST, g in res.groupby("ST"):
        g.drop(columns="ST").to_csv(os.path.join(a.out, f"managed_{ST}.csv"),
                                    index=False)
    # summary
    print(f"{res.ST.nunique()} states; plantation-aware buckets written to {a.out}")
    for ST in ["GA", "ME", "OR", "MS"]:
        gg = res[(res.ST == ST) & (res.metric == "agc_live_total")
                 & (res.year == res.year.max())]
        if len(gg):
            d = {r.mgmt.split("(")[1][:3]: r.value for r in gg.itertuples()}
            print(f"  {ST} {res.year.max()} AGC Mg/ha: "
                  f"ext={d.get('ext')} harvest={d.get('har')} int={d.get('int')}")


if __name__ == "__main__":
    main()
