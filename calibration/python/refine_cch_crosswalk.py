#!/usr/bin/env python3
"""Refine the FIA->ORGANON crown group crosswalk and re-fit the affine cch map.

The gompit projection uses cch = crown closure at tree tip, computed at run time
by an ORGANON crown-closure port (cch_organon) and mapped onto the gompit cch
scale by an affine fit cch = A + B*cch_hat. The original crosswalk was a coarse
softwood/hardwood proxy (FIA<300 -> ORGANON group 1 DF, else 16 RA), which fits
PNW conifer well but southern/eastern species poorly.

This script replaces it with a genus / crown-form crosswalk over all 18 ORGANON
SWO groups, recomputes cch_hat on the held cch validation sample (the same
117k-tree, 4000-plot set used to fit the original affine map), and reports
whether Spearman correlation with the panel's stored CCH1 improves. If it does,
the new CCH_A/CCH_B and the GGRP assignment go into cch_organon.py + gompmort.f90.

EXPAN is reconstructed from DBH exactly as validate_cch.R did (6.018 TPA for
DBH>=5 in, 74.965 below), so this reproduces that fit on the identical trees.

Usage:
  python refine_cch_crosswalk.py --sample cch_validation_sample.csv
"""
from __future__ import annotations
import argparse
import numpy as np
import pandas as pd
import cch_organon as CC


# ---------------------------------------------------------------------------
# Refined FIA-SPCD -> ORGANON SWO crown group (1..18).
# Groups: 1 DF, 2 GW(true fir/spruce), 3 PP(hard pine/larch), 4 SP(white pine),
# 5 IC(incense cedar/juniper/cypress), 6 WH(hemlock/redwood/baldcypress),
# 7 RC(Thuja/arborvitae), 8 PY(yew/torreya), 9 MD(madrone), 10 GC(chinquapin),
# 11 TA(tanoak), 12 CL(evergreen live oak), 13 BL(maple), 14 WO(generic broad
# hardwood / white & red oaks / ash / hickory / etc.), 16 RA(alder/birch/aspen),
# 17 PD(dogwood), 18 WI(willow). (15 BO folded into 14; 13 used for maples.)
# ---------------------------------------------------------------------------
WHITE_PINES = {101, 113, 117, 119, 129}                  # soft (white) pine subgenus
LIVE_OAKS   = {801, 805, 807, 808, 810, 843, 846, 838}   # western evergreen oaks


import os
MODE = os.environ.get("GRP_MODE", "refined")  # refined | conifer (hardwoods->16)


def grp(spcd: int) -> int:
    s = int(spcd)
    if MODE == "conifer" and s >= 300:
        return 16   # keep hardwoods on the coarse RA proxy
    # ---- conifers ----
    if s < 300:
        if s == 202 or s == 201:           return 1   # Douglas-fir
        if 10 <= s <= 19:                  return 2   # Abies (true firs)
        if 90 <= s <= 99:                  return 2   # Picea (spruce)
        if 260 <= s <= 269:                return 6   # Tsuga (hemlock)
        if s in (211, 212, 221, 222):      return 6   # redwood/sequoia/baldcypress
        if s in (241, 242):                return 7   # Thuja / arborvitae
        if s in (231, 251):                return 8   # yew / torreya
        if 40 <= s <= 69 or s == 81:       return 5   # cedars / junipers / cypress
        if 70 <= s <= 79:                  return 3   # Larix (larch)
        if 100 <= s <= 140:                            # pines
            return 4 if s in WHITE_PINES else 3
        return 1                                       # other conifer -> DF
    # ---- hardwoods ----
    if 310 <= s <= 329:                    return 13  # Acer (maples)
    if 350 <= s <= 360:                    return 16  # Alnus (alder)
    if s == 361:                           return 9   # Arbutus (madrone)
    if 370 <= s <= 379:                    return 16  # Betula (birch)
    if s == 431:                           return 10  # Chrysolepis (chinquapin)
    if s in (491, 492):                    return 17  # Cornus (dogwood)
    if s == 631:                           return 11  # tanoak
    if 740 <= s <= 749:                    return 16  # Populus (aspen/cottonwood)
    if 920 <= s <= 928:                    return 18  # Salix (willow)
    if 800 <= s <= 849:                                # Quercus (oaks)
        return 12 if s in LIVE_OAKS else 14
    return 14                                          # generic broad hardwood


def cch_hat_for_plot(g: pd.DataFrame) -> np.ndarray:
    trees = []
    for _, r in g.iterrows():
        dbh = float(r["DBH"])
        trees.append(dict(group=int(r["G"]), DBH=dbh, HT=float(r["HT"]),
                          CR=float(r["CR"]),
                          EXPAN=6.018 if dbh >= 5 else 74.965))
    valid = [t for t in trees if np.isfinite(t["HT"]) and t["HT"] > 0
             and 0 < t["CR"] <= 1 and np.isfinite(t["DBH"]) and t["DBH"] > 0]
    if not valid:
        return np.full(len(g), np.nan)
    prof = CC.crown_closure(valid)
    out = []
    for _, r in g.iterrows():
        h = float(r["HT"])
        out.append(CC.tree_cch(h, prof) if np.isfinite(h) and h > 0 else np.nan)
    return np.array(out)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    a = ap.parse_args()
    d = pd.read_csv(a.sample)
    d["G"] = d["SPCD"].apply(grp)
    d["cch_new"] = (d.groupby("PLT_CN", group_keys=False)
                     .apply(lambda g: pd.Series(cch_hat_for_plot(g), index=g.index)))
    m = d.copy()

    def report(col):
        sub = m[np.isfinite(m[col]) & np.isfinite(m["CCH1"])]
        x = sub[col].to_numpy(float)
        y = sub["CCH1"].to_numpy(float)
        pear = np.corrcoef(x, y)[0, 1]
        spear = pd.Series(x).corr(pd.Series(y), method="spearman")
        # ordinary least squares slope/intercept via covariance (robust)
        xm, ym = x.mean(), y.mean()
        b = ((x - xm) * (y - ym)).sum() / ((x - xm) ** 2).sum()
        a0 = ym - b * xm
        return pear, spear, a0, b, len(sub)

    pe_o, sp_o, Ao, Bo, no = report("cch_hat")
    pe_n, sp_n, A, B, nn = report("cch_new")
    print(f"COARSE (sw/hw):  n={no}  Pearson={pe_o:.3f}  Spearman={sp_o:.3f}  (A={Ao:.4f} B={Bo:.6f})")
    print(f"REFINED (genus): n={nn}  Pearson={pe_n:.3f}  Spearman={sp_n:.3f}")
    print(f"refined affine map: CCH_A={A:.4f}  CCH_B={B:.6f}")
    # group usage
    print("group usage (refined):")
    print(d.groupby("G")["SPCD"].nunique().to_string())


if __name__ == "__main__":
    main()
