#!/usr/bin/env python3
"""Projection-time mortality glue: cch_organon -> affine scale -> GregMortality.

This is the single call perseus.process_plot needs each cycle to apply Greg's
gompit(cr, cch) survival instead of FVS's native mortality. It:
  1. builds the ORGANON crown-closure profile for the cycle treelist,
  2. reads each tree's cch at its tip,
  3. maps it onto the gompit's cch scale via the validated affine fit
     CCH1 ~ 0.062 + 0.0036 * cch_hat  (calibration/R/35d_validate_cch.R;
     Spearman 0.93 vs stored CCH1 over 4000 plots),
  4. returns the per-tree period survival probability (None entries kept native).

Treelist columns expected (rename map configurable): SpeciesFIA, DBH (in), Ht
(ft), CrRatio (0-1), TPA. cch_organon needs an FIA-SPCD -> ORGANON group map;
the coarse default (softwood->1 DF, hardwood->16 RA) matches the validation.
Refine SPP_GROUP for a tighter fit.

Usage in perseus (sketch):
    from project_mortality import cycle_survival
    surv = cycle_survival(treelist_df, greg_mort, years=cycle_length)
    treelist_df["TPA"] *= [s if s is not None else native for s, native in ...]
"""
from __future__ import annotations
import numpy as np
import cch_organon as CC
from greg_mortality import GregMortality

# validated affine map (35d_validate_cch.R): stored-cch ~ A + B * cch_hat
CCH_A, CCH_B = 0.062, 0.0036


def spp_group(spcd: int) -> str:
    """Coarse FIA-SPCD -> ORGANON crown group (matches the cch validation)."""
    return "1" if int(spcd) < 300 else "16"


def cycle_survival(tl, greg: GregMortality, years: float,
                   spcd_col="SpeciesFIA", dbh_col="DBH", ht_col="Ht",
                   cr_col="CrRatio", tpa_col="TPA"):
    """Return a list of period-survival probs (one per tree); None where the
    species has no fitted gompit (caller keeps native mortality)."""
    trees = [dict(group=spp_group(r[spcd_col]), DBH=float(r[dbh_col]),
                  HT=float(r[ht_col]), CR=float(r[cr_col]),
                  EXPAN=float(r[tpa_col])) for _, r in tl.iterrows()]
    valid = [t for t in trees if np.isfinite(t["HT"]) and t["HT"] > 0
             and 0 < t["CR"] <= 1 and np.isfinite(t["DBH"]) and t["DBH"] > 0]
    if not valid:
        return [None] * len(trees)
    profile = CC.crown_closure(valid)
    out = []
    for _, r in tl.iterrows():
        try:
            cch_hat = CC.tree_cch(float(r[ht_col]), profile)
            cch = CCH_A + CCH_B * cch_hat          # onto the gompit's scale
            out.append(greg.period_survival(r[spcd_col], r[cr_col], cch, years))
        except Exception:
            out.append(None)
    return out
