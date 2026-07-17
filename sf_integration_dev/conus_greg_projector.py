#!/usr/bin/env python3
"""conus_greg_projector.py -- Python CONSUMER for the landed Greg CONUS block.

This is the missing downstream consumer for FvsConfigLoader(variant,
version="conus_greg"). It reads the categories_conus_greg block via the
existing accessors (get_conus_greg_block / get_greg_driver_coefficients) and
evaluates Greg Johnson's deployed fvs_remodeling equation forms per tree:

  DIAMETER GROWTH (in yr^-1), form greg_est_dg (species_params B0..B6):
    z  = B0 + B1*ln((dbh+1)^2 / (cr*ht+1)^B3)
            + B2*bal^B4 / ln(dbh+2.7) + B5*elev + B6*EMT
    dg = exp(clamp(z, -30, 5)),  floored at 0
    units: dbh in inches, cr fraction, ht feet, bal ft^2 ac^-1,
           elev feet, EMT degC (ClimateNA extreme minimum temperature)

  HEIGHT GROWTH (ft yr^-1), form greg_est_hg (species_params B0=max_ht, B1..B8):
    dht = max_ht*b1*b2*cr^b3
          * exp(-b1*ht -b4*ccfl -b8*cch^0.5 -b5*elev +b6*TD^0.5 +b7*EMT)
          * (1 - exp(-b1*ht))^(b2-1),  floored at 0
    units: ht feet, cr fraction, ccfl ft^2 ac^-1, cch fraction (0-1),
           elev feet, TD degC (= MWMT - MCMT), EMT degC

  SURVIVAL (annual probability), form greg_gompit (species_params b0..b4):
    eta    = b0 + b1*(cr+0.01)^b2 + b3*cch^b4     (cch term = 0 when cch <= 0)
    P_surv = 1 - exp(-exp(clamp(eta, -30, 30))),  0 when cr <= 0

These forms and clamps mirror the authoritative reference projector
(/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/conus_eq_projector_greg.R,
functions dg_annual / hg_annual / gomp_surv_annual). Verified 2026-07-17: the
landed species_params equal Greg's raw dg_parms.RDS / hg_parms.RDS /
mort_parm_base_rate_cr_cch.RDS to machine precision.

Species fallback (~9-11% of trees, documented in block.species_fallback):
species outside Greg's fitted set use the softwood / hardwood MEDIAN of the
in-set coefficients (SPCD < 300 = conifer convention, matching the R fallback).
Inside a stand loop the reference projector instead borrows the stand-modal
in-set species; pass fallback_spcd=<modal spcd> to reproduce that exactly.

CLIMATE (EMT/TD): required per stand. Supply a scalar (emt=, td=), a per-stand
CSV lookup (columns STAND_CN,EMT,TD via emt_td_lookup=), or fall back to the
documented ClimateNA 1991-2020 median across greg_emt_td_lookup.rds
(EMT=-28.8 degC, TD=24.8 degC). The scalar default is a placeholder ONLY: a
real per-stand lookup is REQUIRED for production. Whenever the default is used
the projector sets .used_default_climate = True so callers can flag it.

Invocation mirrors the conus_sf consumer so both CONUS options are called the
same way:

    from config.config_loader import FvsConfigLoader
    from conus_greg_projector import GregProjector
    L = FvsConfigLoader("ne", version="conus_greg")
    gp = GregProjector(L)                 # or emt=..., td=..., emt_td_lookup="path.csv"
    dg = gp.dg_annual(97, dbh_in=12, cr=0.5, ht_ft=55, bal=90, elev_ft=1200)
    ps = gp.survival_annual(97, cr=0.5, cch=0.4)   # annual P(survive)

Author: A. Weiskittel + Claude   Date: 2026-07-17
"""
from __future__ import annotations

import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from config.config_loader import FvsConfigLoader  # noqa: E402

# Documented placeholder climate scalars: overall ClimateNA 1991-2020 medians
# across the 1.86M stands in greg_emt_td_lookup.rds. Not a substitute for a
# real per-stand lookup in production.
DEFAULT_EMT = -28.8
DEFAULT_TD = 24.8


def _is_softwood(spcd: int) -> bool:
    """Greg conifer convention used by the reference fallback (SPCD < 300)."""
    return int(spcd) < 300


def _median(xs):
    xs = sorted(x for x in xs if x is not None)
    n = len(xs)
    if n == 0:
        return None
    m = n // 2
    return xs[m] if n % 2 else 0.5 * (xs[m - 1] + xs[m])


def _params_by_spcd(species_params, keys):
    out = {}
    for row in species_params:
        spcd = int(row["SPCD"])
        vals = {}
        for k in keys:
            v = row.get(k, None)
            vals[k] = None if (v is None or v == "NA") else float(v)
        out[spcd] = vals
    return out


class GregProjector:
    """Per-tree Greg-arm predictor built from a version='conus_greg' config."""

    ARM = "conus_greg"

    def __init__(self, loader: FvsConfigLoader, emt=DEFAULT_EMT, td=DEFAULT_TD,
                 emt_td_lookup=None):
        if getattr(loader, "version", None) != "conus_greg":
            raise ValueError(
                "GregProjector requires FvsConfigLoader(..., version='conus_greg'); "
                f"got version={getattr(loader, 'version', None)!r}"
            )
        if not loader.has_conus_greg_block():
            raise ValueError(
                f"variant '{loader.variant}' carries no categories_conus_greg block"
            )
        self.L = loader
        self.variant = loader.variant
        self.default_emt = float(emt)
        self.default_td = float(td)
        self.used_default_climate = False

        dg = loader.get_conus_greg_block("diameter_growth")
        hg = loader.get_conus_greg_block("height_growth")
        su = loader.get_conus_greg_block("survival")
        self.dg = _params_by_spcd(dg["species_params"],
                                  ["B0", "B1", "B2", "B3", "B4", "B5", "B6"])
        self.hg = _params_by_spcd(hg["species_params"],
                                  ["B0", "B1", "B2", "B3", "B4",
                                   "B5", "B6", "B7", "B8"])
        self.su = _params_by_spcd(su["species_params"],
                                  ["b0", "b1", "b2", "b3", "b4"])
        self.dg_set = set(self.dg)
        self.hg_set = set(self.hg)
        self.su_set = set(self.su)
        self.fallback_rule = loader.config.get(
            "categories_conus_greg", {}).get("species_fallback", {})

        self._dg_fb = self._median_fallback(
            self.dg, ["B0", "B1", "B2", "B3", "B4", "B5", "B6"])
        self._hg_fb = self._median_fallback(
            self.hg, ["B0", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8"])
        self._mort_fb = self._median_fallback(
            self.su, ["b0", "b1", "b2", "b3", "b4"])

        self._emt = {}
        self._td = {}
        self.lookup_source = None
        if emt_td_lookup:
            self._load_lookup(emt_td_lookup)

    # ------------------------------------------------------------------ setup
    @staticmethod
    def _median_fallback(pmap, keys):
        sw = {k: _median([v[k] for s, v in pmap.items() if _is_softwood(s)])
              for k in keys}
        hw = {k: _median([v[k] for s, v in pmap.items() if not _is_softwood(s)])
              for k in keys}
        return {"sw": sw, "hw": hw}

    def _load_lookup(self, path):
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh):
                cn = str(row.get("STAND_CN", "")).split(".")[0]
                try:
                    self._emt[cn] = float(row["EMT"])
                    self._td[cn] = float(row["TD"])
                except (KeyError, ValueError, TypeError):
                    continue
        self.lookup_source = path

    def climate(self, stand_cn=None):
        """Return (EMT, TD) for a stand; falls back to the documented scalar."""
        if stand_cn is not None:
            cn = str(stand_cn).split(".")[0]
            if cn in self._emt:
                return self._emt[cn], self._td[cn]
        self.used_default_climate = True
        return self.default_emt, self.default_td

    def _resolve(self, spcd, pmap, fb, fallback_spcd=None):
        spcd = int(spcd)
        if spcd in pmap:
            return pmap[spcd], False
        if fallback_spcd is not None and int(fallback_spcd) in pmap:
            return pmap[int(fallback_spcd)], True
        side = "sw" if _is_softwood(spcd) else "hw"
        return fb[side], True

    # --------------------------------------------------------------- kernels
    def dg_annual(self, spcd, dbh_in, cr, ht_ft, bal, elev_ft,
                  emt=None, fallback_spcd=None):
        """Annual diameter increment (in yr^-1); greg_est_dg form."""
        p, _ = self._resolve(spcd, self.dg, self._dg_fb, fallback_spcd)
        emt = self.default_emt if emt is None else emt
        B0, B1, B2, B3, B4, B5, B6 = (p["B0"], p["B1"], p["B2"], p["B3"],
                                      p["B4"], p["B5"], p["B6"])
        z = (B0
             + B1 * math.log((dbh_in + 1.0) ** 2 / (cr * ht_ft + 1.0) ** B3)
             + B2 * (bal ** B4) / math.log(dbh_in + 2.7)
             + B5 * elev_ft + B6 * emt)
        z = min(max(z, -30.0), 5.0)
        g = math.exp(z)
        return g if (math.isfinite(g) and g > 0.0) else 0.0

    def hg_annual(self, spcd, ht_ft, cr, ccfl, cch, elev_ft,
                  td=None, emt=None, fallback_spcd=None):
        """Annual height increment (ft yr^-1); greg_est_hg form."""
        p, _ = self._resolve(spcd, self.hg, self._hg_fb, fallback_spcd)
        td = self.default_td if td is None else td
        emt = self.default_emt if emt is None else emt
        mx, b1, b2, b3, b4, b5, b6, b7, b8 = (
            p["B0"], p["B1"], p["B2"], p["B3"], p["B4"],
            p["B5"], p["B6"], p["B7"], p["B8"])
        crp = max(cr, 1e-4)
        cchp = max(cch, 0.0)
        dht = (mx * b1 * b2 * crp ** b3
               * math.exp(-b1 * ht_ft - b4 * ccfl - b8 * cchp ** 0.5
                          - b5 * elev_ft + b6 * math.sqrt(max(td, 0.0))
                          + b7 * emt)
               * (1.0 - math.exp(-b1 * ht_ft)) ** (b2 - 1.0))
        return dht if (math.isfinite(dht) and dht > 0.0) else 0.0

    def survival_annual(self, spcd, cr, cch, fallback_spcd=None):
        """Annual survival probability; greg_gompit form."""
        p, _ = self._resolve(spcd, self.su, self._mort_fb, fallback_spcd)
        b0, b1, b2, b3, b4 = (p["b0"], p["b1"], p["b2"], p["b3"], p["b4"])
        crp = max(cr, 1e-4)
        cchp = max(cch, 0.0)
        eta = b0 + b1 * (crp + 0.01) ** b2 + (b3 * cchp ** b4 if cchp > 0 else 0.0)
        eta = min(max(eta, -30.0), 30.0)
        ps = 1.0 - math.exp(-math.exp(eta))
        if cr <= 0:
            ps = 0.0
        if not math.isfinite(ps):
            ps = 1.0
        return min(max(ps, 0.0), 1.0)

    def survival_period(self, spcd, cr, cch, years=5, fallback_spcd=None):
        """Compounded survival probability over `years` (annual^years)."""
        return self.survival_annual(spcd, cr, cch, fallback_spcd) ** years

    def is_fallback(self, spcd, component="diameter_growth"):
        s = {"diameter_growth": self.dg_set, "height_growth": self.hg_set,
             "survival": self.su_set}[component]
        return int(spcd) not in s

    # -------------------------------------------------------- period increment
    def dg_increment(self, spcd, dbh_in, cr, ht_ft, bal, elev_ft,
                     emt=None, years=5, fallback_spcd=None):
        """Total diameter increment (in) over `years`, updating dbh annually.

        Competition (bal), crown ratio and height are held at their entry
        values within the period, matching the reference projector's
        cycle-start competition convention.
        """
        d = float(dbh_in)
        for _ in range(int(years)):
            d += self.dg_annual(spcd, d, cr, ht_ft, bal, elev_ft,
                                emt=emt, fallback_spcd=fallback_spcd)
        return d - float(dbh_in)


def project(loader, arm, trees0, eco=None, drivers=None, years=25, step=5,
            emt=None, td=None, emt_td_lookup=None):
    """Deterministic Greg-arm stand stepper, signature-compatible with the
    conus_sf `project` entry (loader, arm, trees0, eco, drivers, years, step).

    arm must be 'conus_greg'. trees0: list of dicts with keys
    {spcd, dbh_in, cr, ht_ft, tpa} and per-tree competition
    {bal, ccfl, cch} (held constant per step here; the engine supplies live
    competition in production). Optional 'stand_cn' selects EMT/TD from the
    lookup. Returns {'steps': [...], 'trees': [...]} trajectory.
    """
    if arm != GregProjector.ARM:
        raise ValueError(f"project() arm must be 'conus_greg', got {arm!r}")
    gp = GregProjector(loader,
                       emt=DEFAULT_EMT if emt is None else emt,
                       td=DEFAULT_TD if td is None else td,
                       emt_td_lookup=emt_td_lookup)
    trees = [dict(t) for t in trees0]
    steps = []
    for cy in range(0, years + 1, step):
        ba = sum(math.pi / 4 * (t["dbh_in"] ** 2) / 144.0 * t.get("tpa", 1.0)
                 for t in trees)
        steps.append({"proj_year": cy, "ba_ft2ac": ba,
                      "n_trees": len(trees),
                      "tpa": sum(t.get("tpa", 1.0) for t in trees)})
        if cy == years:
            break
        for t in trees:
            e, d = gp.climate(t.get("stand_cn"))
            e = e if emt is None else emt
            d = d if td is None else td
            dinc = gp.dg_increment(t["spcd"], t["dbh_in"], t["cr"], t["ht_ft"],
                                   t.get("bal", 0.0), t.get("elev_ft", 0.0),
                                   emt=e, years=step)
            hinc = step * gp.hg_annual(t["spcd"], t["ht_ft"], t["cr"],
                                       t.get("ccfl", t.get("bal", 0.0)),
                                       t.get("cch", 0.0), t.get("elev_ft", 0.0),
                                       td=d, emt=e)
            surv = gp.survival_period(t["spcd"], t["cr"], t.get("cch", 0.0),
                                      years=step)
            t["dbh_in"] += dinc
            t["ht_ft"] += hinc
            t["tpa"] = t.get("tpa", 1.0) * surv
    return {"steps": steps, "trees": trees, "arm": arm,
            "used_default_climate": gp.used_default_climate}


if __name__ == "__main__":
    L = FvsConfigLoader("ne", version="conus_greg")
    gp = GregProjector(L)
    print("Greg projector loaded for NE: DG=%d HG=%d SURV=%d species"
          % (len(gp.dg_set), len(gp.hg_set), len(gp.su_set)))
    for s in (12, 97, 316, 833):
        print("SPCD %d  dg(12in)=%.4f in/yr  P_surv(ann)=%.4f"
              % (s, gp.dg_annual(s, 12, 0.5, 55, 90, 1200),
                 gp.survival_annual(s, 0.5, 0.4)))
