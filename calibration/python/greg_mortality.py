#!/usr/bin/env python3
"""Greg Johnson's CONUS mortality (gompit on crown ratio + crown closure at tree
tip), applied at projection time. This is the engine-facing core that the
perseus projection path calls to override the variant's native survival.

Model (per species), annual hazard with an exposure offset for the cycle length:
    H_annual    = exp(b0 + b1*(cr+0.01)^b2 + b3*cch^b4)
    P_survive(T)= exp(-H_annual * T)          # T = cycle length in years

Coefficients come from the re-fit (greg_mortality_coefficients.csv, columns
SPCD,b0..b4). Species without a fitted row fall back to survival=None so the
caller keeps the variant's native mortality for that tree.

Integration point (perseus_100yr_projection.py): after each cycle's treelist is
produced, compute cr and cch per tree, call period_survival(), and scale TPA by
the survival probability instead of using FVS's own mortality. cch must be the
crown closure at the subject tree's tip; if it is not available from the FVS
output it has to be recomputed from the cycle treelist (crown-width profile) --
that recomputation is the remaining piece before the projection-level
old-vs-new comparison can run.
"""
from __future__ import annotations
import numpy as np
import pandas as pd


class GregMortality:
    def __init__(self, coeff_csv: str):
        df = pd.read_csv(coeff_csv)
        self.coef = {int(r.SPCD): (r.b0, r.b1, r.b2, r.b3, r.b4)
                     for r in df.itertuples(index=False)}

    def has(self, spcd: int) -> bool:
        return int(spcd) in self.coef

    def annual_hazard(self, spcd, cr, cch):
        """Annual mortality hazard H for one tree. Returns None if species unfit."""
        c = self.coef.get(int(spcd))
        if c is None:
            return None
        b0, b1, b2, b3, b4 = c
        cr = min(max(float(cr), 1e-4), 1.0)
        cch = max(float(cch), 0.0)
        eta = b0 + b1 * (cr + 0.01) ** b2 + b3 * (cch ** b4 if cch > 0 else 0.0)
        eta = min(max(eta, -30.0), 30.0)
        return float(np.exp(eta))

    def period_survival(self, spcd, cr, cch, years):
        """P(survive `years`). None if species unfit (caller keeps native rate)."""
        H = self.annual_hazard(spcd, cr, cch)
        if H is None:
            return None
        return float(np.exp(-H * float(years)))

    def apply_to_treelist(self, tl: pd.DataFrame, years: float,
                          spcd_col="SpeciesFIA", cr_col="CrRatio",
                          cch_col="CCH", tpa_col="TPA") -> pd.DataFrame:
        """Scale TPA by period survival per tree. Rows whose species is unfit are
        left unchanged (native mortality assumed already applied upstream)."""
        out = tl.copy()
        surv = []
        for _, r in out.iterrows():
            s = self.period_survival(r[spcd_col], r.get(cr_col, 0.5),
                                     r.get(cch_col, 0.0), years)
            surv.append(1.0 if s is None else s)
        out[tpa_col] = out[tpa_col] * np.array(surv)
        return out
