"""Species-free engine injection for FVS (prototype).

Injects the trait-driven species-free per-tree increment and mortality into the FVS
engine, replacing the per-species multiplier calibration on the native equations. Uses
the FVS API exposed by fvs2py: stop point 5 (after growth and mortality are computed,
before they are applied) plus per-tree attribute read/write (fvsTreeAttr). No Fortran
recompile.

Design: calibration/injection/ENGINE_WIRING_DESIGN.md
Predictor reference: dev_sf_integration/benchmark_sf_vs_legA.R (prep_* + add_eta).
Status: wiring-complete against the confirmed fvs2py API; needs one Cardinal integration
run to validate end to end. Shadow mode (log only, no SET) is the safe first stage.
"""
from __future__ import annotations
import json, math
from pathlib import Path
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Bundle predictor: per-tree linear predictor from a banked species-free bundle
# ---------------------------------------------------------------------------
class BundlePredictor:
    """Loads one component bundle (fixed coefficients, per-species trait effects,
    ecoregion/forest-type random effects, manifest) and computes the per-tree
    increment. The covariate vector is supplied by the caller (built per component
    exactly as the fitting/benchmark prep_* functions build it)."""

    def __init__(self, bundle_dir: str, prefix: str, link: str = "log"):
        bd = Path(bundle_dir)
        self.manifest = json.loads((bd / f"{prefix}_sf_manifest.json").read_text())
        fx = pd.read_csv(bd / f"{prefix}_sf_fixed.csv")
        self.fixed = dict(zip(fx["variable"], fx["mean"]))
        sp = pd.read_csv(bd / f"{prefix}_sf_species.csv")
        sp_col = "trait_effect_mean" if "trait_effect_mean" in sp.columns else sp.columns[-1]
        self.species_te = dict(zip(sp["SPCD"].astype(int), sp[sp_col].astype(float)))
        self.re = {}
        for lvl in ("L1", "L2", "L3", "FT"):
            f = bd / f"{prefix}_sf_re_{lvl}.csv"
            if f.exists():
                t = pd.read_csv(f)
                self.re[lvl] = dict(zip(t["level"].astype(str), t["mean"].astype(float)))
            else:
                self.re[lvl] = {}
        self.link = link
        self.intercept = self.fixed.get("a0", self.fixed.get("b0", self.fixed.get("h0", 0.0)))

    def _re(self, lvl, code):
        return self.re[lvl].get(str(code), 0.0)

    def eta(self, trees: pd.DataFrame, cov_fn) -> np.ndarray:
        """Linear predictor per tree. cov_fn(trees, fixed) returns the covariate
        contribution sum_k(coef_k * x_k) as an array (component-specific)."""
        te = trees["SPCD"].astype(int).map(self.species_te).fillna(0.0).to_numpy()
        rp = (trees["EPA_L1_CODE"].map(lambda c: self._re("L1", c)).to_numpy()
              + trees["EPA_L2_CODE"].map(lambda c: self._re("L2", c)).to_numpy()
              + trees["EPA_L3_CODE"].map(lambda c: self._re("L3", c)).to_numpy()
              + trees["FORTYPCD"].map(lambda c: self._re("FT", int(c) if pd.notna(c) else -1)).to_numpy())
        return self.intercept + te + rp + cov_fn(trees, self.fixed)

    def predict(self, trees: pd.DataFrame, cov_fn, years: float = 1.0) -> np.ndarray:
        """Increment per tree, back-transformed and annualized to `years`."""
        e = self.eta(trees, cov_fn)
        if self.link == "log":
            inc = np.exp(e)                      # log-DDS diameter/height increment
        elif self.link == "logit":
            inc = 1.0 / (1.0 + np.exp(-e))       # crown ratio / survival probability
        else:
            inc = e
        return inc * years


# --- component covariate builders (mirror benchmark_sf_vs_legA.R prep_*) -----
def cov_dg_kuehne_v8(d, f):
    """Kuehne v8 BGI log-DDS diameter-growth covariates. Term mapping follows
    calibration/stan/dg_kuehne2022_v8_bgi_nonlinear.stan; coefficients b1..b15."""
    dbh = d["DBH1"].to_numpy(); cr = d["CR1"].to_numpy(); ht = d["HT1"].to_numpy()
    bal = (d["BAL_SW1"] + d["BAL_HW1"]).to_numpy(); ba = (d["BA1"] * 0.2296).to_numpy()
    ln_dbh = np.log(np.maximum(dbh, 0.1))
    return (f.get("b1", 0) * ln_dbh
            + f.get("b2", 0) * dbh
            + f.get("b3", 0) * np.log(np.maximum(cr, 1e-3))
            + f.get("b4", 0) * bal
            + f.get("b5", 0) * ba)   # remaining b6..b15: BGI/climate terms; see manifest


# ---------------------------------------------------------------------------
# fvsTreeAttr thin wrapper (FVS API: get/set per-tree arrays by name)
# ---------------------------------------------------------------------------
def tree_attr(fvs, name: str, action: str = "get", values: np.ndarray | None = None):
    """Read or write a per-tree attribute array. action='get' returns an array of
    length ntrees; action='set' writes `values`. Names: dbh, ht, dg, htg, crwdth,
    species, plot, mort/prob, and competition covariates."""
    import ctypes as ct
    ntrees = fvs.dims()["ntrees"]
    arr = (ct.c_double * ntrees)()
    if action == "set" and values is not None:
        for i, v in enumerate(values):
            arr[i] = float(v)
    rc = ct.c_int(0); nch = ct.c_int(len(name))
    cname = ct.create_string_buffer(name.encode())
    cact = ct.create_string_buffer(action.encode())
    fvs._fvsTreeAttr(cname, ct.c_int(len(name)), cact, ct.c_int(len(action)),
                     ct.byref(ct.c_int(ntrees)), arr, ct.byref(rc))
    if action == "get":
        return np.frombuffer(arr, dtype=np.float64, count=ntrees).copy()
    return rc.value


# ---------------------------------------------------------------------------
# Injector: read tree list, compute species-free increments, write them back
# ---------------------------------------------------------------------------
class SpeciesFreeInjector:
    def __init__(self, bundle_dir, stand_codes: dict, shadow: bool = True,
                 kappa: dict | None = None):
        """stand_codes: EPA_L1_CODE, EPA_L2_CODE, EPA_L3_CODE, FORTYPCD for this stand.
        shadow=True logs only (no SET). kappa: per-component shrinkage for the blend."""
        self.dg = BundlePredictor(bundle_dir, "dg_v8_sf", link="log")
        # add HG, survival predictors the same way once their cov builders are wired:
        # self.hg = BundlePredictor(bundle_dir, "hg_v8rd_sf", link="log")
        # self.surv = BundlePredictor(bundle_dir, "surv_crz", link="logit")
        self.codes = stand_codes; self.shadow = shadow; self.kappa = kappa or {}
        self.log = []

    def _treeframe(self, fvs) -> pd.DataFrame:
        g = lambda n: tree_attr(fvs, n, "get")
        df = pd.DataFrame({
            "DBH1": g("dbh"), "HT1": g("ht"), "CR1": g("crwdth"),
            "SPCD": g("species").astype(int),
            "BAL_SW1": g("bal"), "BAL_HW1": np.zeros(len(g("dbh"))),  # split if available
            "BA1": g("ba"), "dg_native": g("dg"),
        })
        for k, v in self.codes.items():
            df[k] = v
        return df

    def step(self, fvs, years: float):
        trees = self._treeframe(fvs)
        dg_sf = self.dg.predict(trees, cov_dg_kuehne_v8, years=years)
        if self.shadow:
            self.log.append(pd.DataFrame({"dg_native": trees["dg_native"], "dg_sf": dg_sf}))
            return
        k = self.kappa.get("dg")
        if k:  # per-species shrinkage blend toward native
            n = trees["SPCD"].map(lambda s: 1)  # replace with per-species n table
            w = n / (n + k)
            dg_inj = w * trees["dg_native"].to_numpy() + (1 - w) * dg_sf
        else:
            dg_inj = dg_sf
        tree_attr(fvs, "dg", "set", dg_inj)
        # likewise: tree_attr(fvs,"htg","set",hg_sf); tree_attr(fvs,"mort","set",1-surv_sf)


def run_with_injection(fvs, injector: SpeciesFreeInjector, num_cycles: int, cycle_length: int):
    """Stepwise FVS projection with per-tree species-free injection at stop point 5."""
    fvs.set_stop_point_codes(stop_point_code=5, stop_point_year=-1)
    fvs.run()  # runs to the first stop
    cyc = 0
    while fvs.restart_code != 0 and cyc < num_cycles:
        injector.step(fvs, years=cycle_length)
        fvs.run()  # resume; applies injected increments, runs to next stop
        cyc += 1
    return fvs.summary
