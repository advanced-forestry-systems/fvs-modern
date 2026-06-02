#!/usr/bin/env python3
"""Projection comparison isolating the mortality effect: base rate vs Greg's
gompit(cr, cch). Same stands, no growth (pure mortality demography), 100 yr in
20x5yr cycles. cch is recomputed each cycle from the live stand (TPA declines),
so each scenario carries its own cch trajectory. Metric: stand basal area.

Inputs:
  --sample  proj_sample.csv  (PLT_CN, SPCD, DBH1 cm, HT1 m, CR1)
  --coeffs  greg_mortality_coefficients.csv
  --base    species_base.csv (SPCD, surv_annual)   base-rate annual survival
Outputs: projection_compare_trajectory.csv (year, BA_base, BA_greg, ...)
"""
from __future__ import annotations
import argparse, sys, math
import numpy as np, pandas as pd

sys.path.insert(0, "/sessions/serene-brave-johnson/mnt/fvs-modern/calibration/python")
import cch_organon as CC
from greg_mortality import GregMortality

CCH_A, CCH_B = 0.062, 0.0036  # validated affine map (35d): stored_cch ~ A + B*cch_hat
def grp(spcd): return 1 if int(spcd) < 300 else 16
def tpa_init(dbh_in): return 6.018 if dbh_in >= 5 else 74.965
def ba(dbh_in, tpa): return math.pi/4.0 * (dbh_in**2) * tpa / 144.0  # ft2/ac (DBH in -> ft2)

def cch_for(trees):
    valid = [t for t in trees if t["TPA"] > 1e-6 and t["HT"] > 0]
    if not valid: return None
    prof = CC.crown_closure([dict(group=t["g"], DBH=t["DBH"], HT=t["HT"],
                                  CR=t["CR"], EXPAN=t["TPA"]) for t in valid])
    return prof

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True); ap.add_argument("--coeffs", required=True)
    ap.add_argument("--base", required=True); ap.add_argument("--out", default="projection_compare_trajectory.csv")
    ap.add_argument("--cycles", type=int, default=20); ap.add_argument("--clen", type=int, default=5)
    a = ap.parse_args()
    greg = GregMortality(a.coeffs)
    base = {int(r.SPCD): float(r.surv_annual) for r in pd.read_csv(a.base).itertuples(index=False)}
    base_default = np.median(list(base.values())) if base else 0.98
    df = pd.read_csv(a.sample)
    df["DBH"] = df.DBH1/2.54; df["HT"] = df.HT1/0.3048; df["CR"] = df.CR1; df["g"] = df.SPCD.map(grp)

    rows = []
    for cn, grp_df in df.groupby("PLT_CN"):
        # two parallel stands (base, greg), each tree carries its own TPA
        st_b = [dict(SPCD=int(r.SPCD), g=r.g, DBH=r.DBH, HT=r.HT, CR=r.CR, TPA=tpa_init(r.DBH)) for r in grp_df.itertuples()]
        st_g = [dict(**t) for t in st_b]
        for cyc in range(a.cycles+1):
            BA_b = sum(ba(t["DBH"], t["TPA"]) for t in st_b)
            BA_g = sum(ba(t["DBH"], t["TPA"]) for t in st_g)
            rows.append(dict(PLT_CN=cn, year=cyc*a.clen, BA_base=BA_b, BA_greg=BA_g))
            if cyc == a.cycles: break
            for st, scen in ((st_b,"base"), (st_g,"greg")):
                prof = cch_for(st)
                for t in st:
                    if t["TPA"] <= 1e-6: continue
                    if scen == "base":
                        s = base.get(t["SPCD"], base_default) ** a.clen
                    else:
                        cch = CCH_A + CCH_B*(CC.tree_cch(t["HT"], prof) if prof is not None else 0.0)
                        s = greg.period_survival(t["SPCD"], t["CR"], cch, a.clen)
                        if s is None: s = base.get(t["SPCD"], base_default) ** a.clen
                    t["TPA"] *= s

    traj = pd.DataFrame(rows).groupby("year").agg(
        BA_base=("BA_base","mean"), BA_greg=("BA_greg","mean"), n=("PLT_CN","nunique")).reset_index()
    traj["greg_minus_base_pct"] = 100*(traj.BA_greg-traj.BA_base)/traj.BA_base
    traj.to_csv(a.out, index=False)
    pd.set_option("display.width", 120)
    print(traj.round(3).to_string(index=False))
    f = traj.iloc[-1]
    print(f"\n100-yr stand BA: base={f.BA_base:.1f}  greg={f.BA_greg:.1f} ft2/ac  "
          f"(Greg {f.greg_minus_base_pct:+.1f}% vs base) over {int(f.n)} stands")

if __name__ == "__main__":
    main()
