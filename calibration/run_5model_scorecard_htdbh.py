#!/usr/bin/env python3
"""
run_5model_scorecard_cohort.py
==============================
Cohort-matched, forest-type-stratified 5-model standalone FVS scorecard using
FIA remeasurement pairs.  A COPY of run_5model_scorecard_biasfix.py with two
methodological changes; the FVS model systems and keyword logic are unchanged.

CHANGE 1 -- Cohort matching (observed / T2 side):
  Instead of comparing FVS's projection of the initial (T1) cohort against the
  whole observed T2 stand (which includes ingrowth and 5"-recruits), we compare
  against ONLY the surviving members of that same initial cohort.  A T2 live
  tree (STATUSCD==1) is kept iff its PREV_TRE_CN links back to a T1 tree that
  was (a) live at T1 (STATUSCD==1) and (b) in the FVS input list (DIA_T1 >= 5 in).
  Ingrowth (T2 trees with null PREV_TRE_CN) and trees < 5 in at T1 are excluded.
  This makes observed BA/MCuFt/BdFt a like-for-like comparison against FVS's
  projected initial cohort (FVS ingrowth is off).  The DIA >= 5 in floor is kept.
  Method used: PREV_TRE_CN linkage (available in FIA TREE).  A SUBP+TREE fallback
  is provided only if PREV_TRE_CN is entirely missing for a plot.

CHANGE 2 -- Forest-type stratification (per plot, at T1):
  Each plot's softwood share of T1 basal area is computed from the T1 cohort
  (softwood = SPCD < 300, hardwood = SPCD >= 300; the standard FIA convention).
  Classification:  SW if softwood BA share >= 0.75,
                   HW if softwood BA share <= 0.25 (hardwood share >= 0.75),
                   MW otherwise.  A FORTYPE column is added to the output.

Usage:
  python3 run_5model_scorecard_cohort.py --variant ne [--max-pairs 5000] [--seed 42]
"""
from __future__ import annotations

import argparse
import math
import os
import shutil
import sqlite3
import subprocess
import sys
import tempfile

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = "/users/PUOM0008/crsfaaron/fvs-modern"
LIB_DIR      = os.path.join(PROJECT_ROOT, "src-converted", "bin")
CONFIG_DIR   = os.path.join(PROJECT_ROOT, "config")
FIA_DIR      = "/fs/scratch/PUOM0008/crsfaaron/FIA"

_CONF_DG = "/fs/scratch/PUOM0008/crsfaaron/wt-ne-dg/config"

VARIANT_BINS: dict[str, str] = {
    "ne": "/users/PUOM0008/crsfaaron/fvs-modern/src-converted/bin/FVSne",
    "sn": "/users/PUOM0008/crsfaaron/fvs-modern/src-converted/bin/FVSsn",
    "wc": "/users/PUOM0008/crsfaaron/fvs-modern/src-converted/bin/FVSwc",
    "em": "/users/PUOM0008/crsfaaron/fvs-modern/src-converted/bin/FVSem",
    "pn": "/users/PUOM0008/crsfaaron/fvs-modern/src-converted/bin/FVSpn",
}

_HG_BASE: dict[str, str] = {
    "FVS_GREGHG":    "1",
    "FVS_GREG_TD":   "28",
    "FVS_GREG_EMT":  "-30",
    "FVS_GREG_ELEV": "1000",
    "FVS_GREG_MCW_COEF": "/users/PUOM0008/crsfaaron/fvs-modern/config/greg_mcw_coefficients.csv",
}

_DG_ENV: dict[str, str] = {
    "FVS_GREGDG":      "0",
    "FVS_GREGDG_COEF": f"{_CONF_DG}/greg_dg_coefficients.csv",
    "FVS_GREG_DD0":    "1600",
    "FVS_GREG_TD":     "28",
    "FVS_GREG_PPT_SM": "230",
    "FVS_GREG_DD18":   "3500",
}

_DG_KEYS: frozenset[str] = frozenset(_DG_ENV)

_CRW_ENV: dict[str, str] = {
    "FVS_GREGCRW":      "1",
    "FVS_GREGCRW_COEF": f"{_CONF_DG}/greg_crown_change_coefficients.csv",
}

_MORT_ENV: dict[str, str] = {
    "FVS_GOMPIT":      "1",
    "FVS_GOMPIT_COEF": f"{_CONF_DG}/greg_mortality_coefficients_size_bgi.csv",
}

# ---------------------------------------------------------------------------
# Variant metadata
# ---------------------------------------------------------------------------
VARIANTS: dict[str, dict] = {
    "ne": {"states": ["CT", "DE", "MA", "MD", "ME", "NH", "NJ", "NY", "PA", "RI", "VT", "WV"]},
    "sn": {"states": ["AL", "GA", "MS", "SC"]},
    "wc": {"states": ["OR", "WA"]},
    "em": {"states": ["MT", "ND", "SD"]},
    "pn": {"states": ["OR", "WA"]},
}

STATECD_TO_ABBREV: dict[int, str] = {
    1:  "AL",  2:  "AK",  4:  "AZ",  5:  "AR",  6:  "CA",
    8:  "CO",  9:  "CT", 10: "DE",  12: "FL",  13: "GA",
    16: "ID",  17: "IL",  18: "IN",  19: "IA",  20: "KS",
    21: "KY",  22: "LA",  23: "ME",  24: "MD",  25: "MA",
    26: "MI",  27: "MN",  28: "MS",  29: "MO",  30: "MT",
    31: "NE",  32: "NV",  33: "NH",  34: "NJ",  35: "NM",
    36: "NY",  37: "NC",  38: "ND",  39: "OH",  40: "OK",
    41: "OR",  42: "PA",  44: "RI",  45: "SC",  46: "SD",
    47: "TN",  48: "TX",  49: "UT",  50: "VT",  51: "VA",
    53: "WA",  54: "WV",  55: "WI",  56: "WY",
}

# ---------------------------------------------------------------------------
# 5 Standalone Model Systems  (unchanged from biasfix runner)
# ---------------------------------------------------------------------------
STANDALONE_MODELS: dict[str, dict] = {
    "fvs_base": {
        "extra_env":          {},
        "calib_components":   None,
        "variant_compatible": None,
        "states_restrict":    None,
        "states_restrict_dg": None,
    },
    "fvs_regional": {
        "extra_env":          {},
        "calib_components":   {"dg": True, "hg": True, "mort": True},
        "variant_compatible": None,
        "states_restrict":    None,
        "states_restrict_dg": None,
    },
    "organon": {
        "extra_env": {
            **_DG_ENV,
            **_HG_BASE,
            "FVS_GREGHG_COEF": f"{_CONF_DG}/greg_hg_coefficients.csv",
            **_CRW_ENV,
            **_MORT_ENV,
        },
        "calib_components":   {"dg": False, "hg": False, "mort": False},
        "variant_compatible": None,
        "states_restrict":    None,
        "states_restrict_dg": ["CT", "MA", "ME", "NH", "NY", "RI", "VT"],
    },
    "conus_spdep": {
        "extra_env": {
            **_DG_ENV,
            **_HG_BASE,
            "FVS_GREGHG_COEF": f"{_CONF_DG}/greg_hg_coefficients_compound_none.csv",
            **_CRW_ENV,
            **_MORT_ENV,
        },
        "calib_components":   {"dg": False, "hg": False, "mort": False},
        "variant_compatible": None,
        "states_restrict":    None,
        "states_restrict_dg": ["CT", "MA", "ME", "NH", "NY", "RI", "VT"],
    },
    "conus_climate": {
        "extra_env": {
            **_DG_ENV,
            **_HG_BASE,
            "FVS_GREGHG_COEF": f"{_CONF_DG}/greg_hg_coefficients_compound_climate.csv",
            **_CRW_ENV,
            **_MORT_ENV,
        },
        "calib_components":   {"dg": False, "hg": False, "mort": False},
        "variant_compatible": None,
        "states_restrict":    None,
        "states_restrict_dg": ["CT", "MA", "ME", "NH", "NY", "RI", "VT"],
    },
}

MODEL_ORDER = ["fvs_base", "fvs_regional", "organon", "conus_spdep", "conus_climate"]

# ---------------------------------------------------------------------------
# KEYFILE template
# ---------------------------------------------------------------------------
KEYFILE = """\
STDIDENT
{sid}
DATABASE
DSNIN
{db}
DSNOUT
{db}
STANDSQL
SELECT * FROM fvs_standinit WHERE stand_id = '%StandID%'
ENDSQL
TREESQL
SELECT * FROM fvs_treeinit WHERE stand_id = '%StandID%'
ENDSQL
END
DATABASE
SUMMARY            2
END
TIMEINT            0         1
NUMCYCLE          {ncyc}
{calib_kw}
PROCESS
STOP
"""

# ---------------------------------------------------------------------------
# Stand initializer
# ---------------------------------------------------------------------------
def build_standinit(
    sid:      str,
    inv_year: int,
    variant:  str,
    cond_row: dict | None = None,
) -> pd.DataFrame:
    state_cd = 23
    elev     = 1000
    slope    = 5
    site_idx = 42
    age      = 60

    if cond_row is not None:
        def _safe(key: str, default, cast=int):
            try:
                v = cond_row.get(key)
                if v is None or (isinstance(v, float) and math.isnan(v)):
                    return default
                return cast(v)
            except (TypeError, ValueError):
                return default

        state_cd = _safe("STATECD", 23)
        elev     = _safe("ELEV",    1000)
        slope    = _safe("SLOPE",   5)
        site_idx = _safe("SICOND",  42)

    return pd.DataFrame([{
        "stand_id":          sid,
        "variant":           variant.upper(),
        "inv_year":          int(inv_year),
        "latitude":          46.5,
        "longitude":         -68.7,
        "region":            9,
        "forest":            0,
        "district":          0,
        "basal_area_factor": 0.0,
        "inv_plot_size":     1.0,
        "brk_dbh":           999.0,
        "num_plots":         1,
        "age":               age,
        "aspect":            0,
        "slope":             slope,
        "elevft":            elev,
        "site_species":      12,
        "site_index":        site_idx,
        "state":             state_cd,
        "county":            0,
        "forest_type":       121,
        "sam_wt":            1.0,
    }])


# ---------------------------------------------------------------------------
# Wykoff HT-DBH imputation (CONUS preflight, traits0 population posterior)
# ---------------------------------------------------------------------------
# Selected form (ht_dbh_wykoff_simple_lognormal.stan; Stan preflight job
# 12349016 -- the only converged fit):
#   log(HT_m - 1.37) = b0 + b1/(DBH_cm + 1)
#                      + a_bal*BAL + a_ba*sqrt(BA)
#                      + a_cspi*ln(cspi+shift) + a_bard*(BA*rd) + a_blrd*(BAL*rd)
#   HT_m = 1.37 + exp(eta);   HT_ft = HT_m * 3.28084
# Population posterior means from preflight_htdbh_wykoff_traits0.csv.  The fit
# is METRIC (DBH in cm, BA/BAL in m^2/ha).  Species + EPA-L1 random effects and
# trait fixed effects are dropped (population mean = 0) for the population model.
#
# REDUCED FORM: the cohort runner supplies only per-tree DBH, stand BA and
# per-tree BAL.  The site index (cspi) and the two relative-density interaction
# covariates (rd -> ba_x_rd, bal_x_rd) are not available in this harness, so
# a_cspi, a_bard and a_blrd are omitted.  Retained: b0, b1/(DBH+1), a_bal*BAL,
# a_ba*sqrt(BA).  These are the dominant terms (|a_cspi|,|a_bard|,|a_blrd| are
# ~1-2 orders of magnitude smaller in effect, and the two rd terms have opposite
# signs and largely cancel).
WYK_B0    =  2.26219416462
WYK_B1    = -9.3180488111
WYK_A_BAL = -0.0025946649437
WYK_A_BA  =  0.086873289331
_FT2AC_TO_M2HA = 0.2295684        # 1 ft^2 ac^-1 -> m^2 ha^-1
_HT_IMPUTE_STATS = {"imputed": 0, "measured": 0}


def _wykoff_ht_ft(dbh_in: float, bal_ft2ac: float, sqrt_ba_m2ha: float) -> float:
    """Wykoff population HT (feet) for a tree missing a measured height.

    dbh_in       tree DBH in inches (converted to cm for the metric fit)
    bal_ft2ac    basal area (ft^2/ac) in larger trees (converted to m^2/ha)
    sqrt_ba_m2ha sqrt of stand BA already expressed in m^2/ha
    """
    dbh_cm   = dbh_in * 2.54
    bal_m2ha = bal_ft2ac * _FT2AC_TO_M2HA
    eta = (WYK_B0
           + WYK_B1 / (dbh_cm + 1.0)
           + WYK_A_BAL * bal_m2ha
           + WYK_A_BA  * sqrt_ba_m2ha)
    ht_m = 1.37 + math.exp(eta)
    return ht_m * 3.28084


# ---------------------------------------------------------------------------
# Tree list builder
# ---------------------------------------------------------------------------
def build_treeinit(tree_df: pd.DataFrame, sid: str) -> pd.DataFrame:
    """Convert FIA TREE records (live, DIA >= 5 in) to fvs_treeinit.

    HT-DBH CHANGE vs the cohort runner: trees WITHOUT a measured FIA height
    (HT null or <= 0) receive a Wykoff population height imputed from DBH plus
    the stand covariates the runner can compute (stand BA and per-tree BAL).
    Measured heights are never overwritten; every other field is identical.
    """
    # First pass: keep valid DIA>=5 records and their per-tree BA (ft^2/ac).
    recs: list[dict] = []
    for _, r in tree_df.iterrows():
        try:
            spcd = int(r["SPCD"])
            dbh  = float(r["DIA"])
        except (ValueError, TypeError, KeyError):
            continue
        if not (math.isfinite(dbh) and dbh >= 5.0):
            continue
        tpa_raw = r.get("TPA_UNADJ")
        tpa = float(tpa_raw) if (tpa_raw is not None and pd.notna(tpa_raw)) else 6.018
        if tpa <= 0:
            tpa = 6.018
        ht_raw = r.get("HT")
        ht = float(ht_raw) if (ht_raw is not None and pd.notna(ht_raw) and float(ht_raw) > 0) else 0.0
        cr_raw = r.get("CR")
        cr = float(cr_raw) if (cr_raw is not None and pd.notna(cr_raw)) else 50.0
        if cr <= 9:
            cr = cr * 10
        cr = int(max(15, min(99, cr)))
        recs.append({
            "spcd": spcd, "dbh": dbh, "tpa": tpa, "ht": ht, "cr": cr,
            "ba_ft2ac": 0.005454154 * dbh * dbh * tpa,
        })
    if not recs:
        return pd.DataFrame()

    # Stand-level covariate (from the DIA>=5 cohort = the FVS input list).
    stand_ba_ft2ac = sum(rec["ba_ft2ac"] for rec in recs)
    sqrt_ba_m2ha   = math.sqrt(max(stand_ba_ft2ac * _FT2AC_TO_M2HA, 0.0))

    rows = []
    for rec in recs:
        ht = rec["ht"]
        if ht > 0:
            _HT_IMPUTE_STATS["measured"] += 1
        else:
            # BAL = basal area (ft^2/ac) in trees strictly larger in DBH.
            bal_ft2ac = sum(o["ba_ft2ac"] for o in recs if o["dbh"] > rec["dbh"])
            ht = _wykoff_ht_ft(rec["dbh"], bal_ft2ac, sqrt_ba_m2ha)
            _HT_IMPUTE_STATS["imputed"] += 1
        rows.append({
            "stand_id":   sid,
            "plot_id":    1,
            "tree_id":    len(rows) + 1,
            "tree_count": round(rec["tpa"], 4),
            "species":    rec["spcd"],
            "diameter":   round(rec["dbh"], 2),
            "ht":         round(ht, 1),
            "crratio":    rec["cr"],
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Calibration keyword helpers  (unchanged)
# ---------------------------------------------------------------------------
_calib_kw_cache: dict[str, str] = {}


def calibrated_keywords(variant: str) -> str:
    if variant in _calib_kw_cache:
        return _calib_kw_cache[variant]
    try:
        sys.path.insert(0, PROJECT_ROOT)
        from config.config_loader import FvsConfigLoader
        kw = FvsConfigLoader(
            variant.lower(), version="calibrated", config_dir=CONFIG_DIR
        ).generate_keywords(include_comments=False)
        _calib_kw_cache[variant] = kw
        return kw
    except Exception as e:
        sys.stderr.write(f"  calibrated kw unavailable for {variant}: {e}\n")
        kw = "** CALIBRATED CONFIG NOT FOUND"
        _calib_kw_cache[variant] = kw
        return kw


_partial_kw_cache: dict[str, str] = {}


def partial_calibrated_keywords(
    variant:      str,
    include_dg:   bool,
    include_hg:   bool,
    include_mort: bool = True,
) -> str:
    cache_key = (
        f"{variant}_dg{int(include_dg)}_hg{int(include_hg)}_mort{int(include_mort)}"
    )
    if cache_key in _partial_kw_cache:
        return _partial_kw_cache[cache_key]
    try:
        sys.path.insert(0, PROJECT_ROOT)
        from config.config_loader import FvsConfigLoader
        loader = FvsConfigLoader(
            variant.lower(), version="calibrated", config_dir=CONFIG_DIR
        )
        cats        = loader.config.get("categories", {})
        precomputed = loader.config.get("calibration_multipliers") or {}
        parts: list[str] = []

        emit_sdimax = bool(loader.config.get("_emit_sdimax", True))
        if emit_sdimax:
            sdi_values = loader._find_sdi_param(cats)
            if sdi_values is not None:
                parts.append(loader._format_sdimax_keywords(sdi_values, False))

        bamax_values = loader._find_bamax_param(cats)
        if bamax_values is not None:
            parts.append(loader._format_bamax_keywords(bamax_values, False))

        if include_mort:
            mort_mult = loader._array_or_none(precomputed.get("mort_multiplier"))
            if mort_mult is None:
                mort_mult = loader._compute_mortality_multipliers(cats)
            if mort_mult is not None:
                parts.append(loader._format_mortmult_keywords(mort_mult, False))

        if include_dg:
            dds = loader._array_or_none(precomputed.get("dds_multiplier"))
            if dds is not None:
                growth_mult = np.sqrt(np.clip(dds, 0.01, 100.0))
            else:
                growth_mult = loader._compute_growth_multipliers(cats)
            if growth_mult is not None:
                parts.append(loader._format_baimult_keywords(growth_mult, False))

        if include_hg:
            hg_mult = loader._array_or_none(precomputed.get("htg_multiplier"))
            if hg_mult is None:
                hg_mult = loader._compute_height_multipliers(cats)
            if hg_mult is not None:
                parts.append(loader._format_htgmult_keywords(hg_mult, False))

        kw = "\n".join(p for p in parts if p)
        _partial_kw_cache[cache_key] = kw
        return kw
    except Exception as e:
        sys.stderr.write(
            f"  partial calib kw error for {variant} "
            f"(dg={include_dg}, hg={include_hg}, mort={include_mort}): {e}\n"
        )
        kw = "** PARTIAL CALIB NOT FOUND"
        _partial_kw_cache[cache_key] = kw
        return kw


def get_calib_kw(variant: str, calib_components) -> str:
    if calib_components is None:
        return "** DEFAULT PARAMETERS"
    if calib_components == "all":
        return calibrated_keywords(variant)
    return partial_calibrated_keywords(
        variant,
        include_dg=bool(calib_components.get("dg", True)),
        include_hg=bool(calib_components.get("hg", True)),
        include_mort=bool(calib_components.get("mort", True)),
    )


# ---------------------------------------------------------------------------
# FVS runner  (unchanged)
# ---------------------------------------------------------------------------
def run_stand(
    sid:        str,
    stand_df:   pd.DataFrame,
    tree_df:    pd.DataFrame,
    variant:    str,
    model_name: str,
    ncyc:       int,
    *,
    extra_env_override: dict | None = None,
    calib_override=None,
) -> tuple[pd.DataFrame | None, str | None]:
    mdl    = STANDALONE_MODELS[model_name]
    binary = VARIANT_BINS.get(variant.lower())
    if binary is None:
        return None, f"no binary registered for variant {variant}"
    if not os.path.exists(binary):
        return None, f"binary not found: {binary}"
    if len(tree_df) == 0:
        return None, "empty tree list"

    eff_extra_env = extra_env_override if extra_env_override is not None else mdl["extra_env"]
    eff_calib     = calib_override     if calib_override     is not None else mdl["calib_components"]
    calib_kw      = get_calib_kw(variant, eff_calib)

    env = os.environ.copy()
    env.update(eff_extra_env)

    tmp = tempfile.mkdtemp(prefix="fvs5mdl_")
    try:
        db  = os.path.join(tmp, "FVS_Data.db")
        con = sqlite3.connect(db)
        stand_df.to_sql("fvs_standinit", con, if_exists="replace", index=False)
        tree_df.to_sql("fvs_treeinit",   con, if_exists="replace", index=False)
        con.close()

        key = os.path.join(tmp, "run.key")
        with open(key, "w") as fh:
            fh.write(KEYFILE.format(sid=sid, db=db, ncyc=ncyc, calib_kw=calib_kw))

        subprocess.run(
            [binary, f"--keywordfile={key}"],
            cwd=tmp,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            env=env,
            timeout=120,
        )

        con = sqlite3.connect(db)
        df  = pd.read_sql_query("SELECT * FROM FVS_Summary2", con)
        con.close()
        return df, None

    except subprocess.TimeoutExpired:
        return None, "timeout"
    except Exception as exc:
        return None, str(exc)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------------------
# FIA data loaders
# ---------------------------------------------------------------------------
_PLOT_COLS = {"CN", "PREV_PLT_CN", "MEASYEAR", "STATECD", "COUNTYCD", "PLOT"}
# COHORT: added CN and PREV_TRE_CN (needed for tree-level remeasurement linkage)
_TREE_COLS = {"CN", "PLT_CN", "PREV_TRE_CN", "STATECD", "COUNTYCD", "PLOT",
              "SUBP", "TREE", "STATUSCD", "SPCD", "DIA", "HT", "CR",
              "TPA_UNADJ", "VOLCFNET", "VOLBFNET"}
_COND_COLS = {"PLT_CN", "STATECD", "COUNTYCD", "INVYR", "COND_STATUS_CD",
              "SICOND", "BALIVE", "FORTYPCD", "STDAGE", "SLOPE", "ELEV"}

# CN-family columns are read as strings to preserve full 64-bit precision
# (float64 loses precision above 2^53, which would break CN matching).
_TREE_DTYPES = {"CN": "string", "PLT_CN": "string", "PREV_TRE_CN": "string"}


def load_pairs(variant: str, max_pairs: int = 200, seed: int = 42) -> pd.DataFrame:
    states = VARIANTS[variant]["states"]
    frames = []
    for st in states:
        f = os.path.join(FIA_DIR, f"{st}_PLOT.csv")
        if not os.path.exists(f):
            continue
        df = pd.read_csv(f, usecols=lambda c: c in _PLOT_COLS, low_memory=False)
        frames.append(df)
    if not frames:
        return pd.DataFrame()

    plots = pd.concat(frames, ignore_index=True)
    plots["CN"] = plots["CN"].astype("int64")
    yr_map = dict(zip(plots["CN"], plots["MEASYEAR"]))

    rem = plots.dropna(subset=["PREV_PLT_CN"]).copy()
    rem["PREV_PLT_CN"] = rem["PREV_PLT_CN"].astype("int64")
    rem["interval"] = rem.apply(
        lambda r: r["MEASYEAR"] - yr_map.get(r["PREV_PLT_CN"], float("nan")),
        axis=1,
    )
    rem = rem[(rem["interval"] >= 5) & (rem["interval"] <= 15)].copy()
    if len(rem) == 0:
        return pd.DataFrame()

    rng = np.random.default_rng(seed)
    n   = min(max_pairs, len(rem))
    idx = rng.choice(len(rem), size=n, replace=False)
    return rem.iloc[idx].reset_index(drop=True)


def load_trees_for_cns(states: list[str], cns: set) -> pd.DataFrame:
    """Load live FIA TREE records for the given plot CNs.

    CN-family columns are kept as strings for exact matching.  Only STATUSCD==1
    (live) trees are retained -- this covers both the T1 initial cohort and the
    T2 survivors used in cohort matching.
    """
    cns_str = {str(int(c)) for c in cns}
    frames  = []
    for st in states:
        f = os.path.join(FIA_DIR, f"{st}_TREE.csv")
        if not os.path.exists(f):
            continue
        df = pd.read_csv(
            f, usecols=lambda c: c in _TREE_COLS, low_memory=False,
            dtype=_TREE_DTYPES,
        )
        sub = df.loc[df["PLT_CN"].isin(cns_str)].copy()
        if "STATUSCD" in sub.columns:
            sub = sub.loc[sub["STATUSCD"] == 1].copy()
        if not sub.empty:
            frames.append(sub)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def load_cond_for_cns(states: list[str], cns: set) -> pd.DataFrame:
    cns_int = {int(c) for c in cns}
    frames  = []
    for st in states:
        f = os.path.join(FIA_DIR, f"{st}_COND.csv")
        if not os.path.exists(f):
            continue
        df = pd.read_csv(f, usecols=lambda c: c in _COND_COLS, low_memory=False)
        df["PLT_CN"] = df["PLT_CN"].astype("int64")
        sub = df.loc[df["PLT_CN"].isin(cns_int)].copy()
        if not sub.empty:
            frames.append(sub)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


# ---------------------------------------------------------------------------
# COHORT: forest-type classification (per plot, at T1)
# ---------------------------------------------------------------------------
def classify_fortype(t1_df: pd.DataFrame) -> tuple[str, float]:
    """Classify a plot by softwood share of T1 basal area (DIA >= 5 in cohort).

    softwood = SPCD < 300, hardwood = SPCD >= 300 (standard FIA convention).
    Returns (FORTYPE, softwood_BA_share).
      SW if share >= 0.75, HW if share <= 0.25, MW otherwise.
    """
    t = t1_df.loc[
        t1_df["DIA"].notna() & t1_df["TPA_UNADJ"].notna() & (t1_df["DIA"] >= 5.0)
    ].copy()
    if len(t) == 0:
        return "UNK", float("nan")
    ba  = t["TPA_UNADJ"].astype(float) * 0.005454154 * t["DIA"].astype(float) ** 2
    spc = t["SPCD"].astype(float)
    tot = float(ba.sum())
    if tot <= 0:
        return "UNK", float("nan")
    sw_share = float(ba[spc < 300].sum() / tot)
    if sw_share >= 0.75:
        ft = "SW"
    elif sw_share <= 0.25:
        ft = "HW"
    else:
        ft = "MW"
    return ft, sw_share


# ---------------------------------------------------------------------------
# COHORT: observed values restricted to surviving initial cohort
# ---------------------------------------------------------------------------
def obs_ba_ft2ac_cohort(
    tree_df: pd.DataFrame, t2_cn: str, cohort_cns: set, min_dia: float = 5.0
) -> tuple[float, int]:
    """Observed T2 BA (ft2/ac) for surviving members of the T1 cohort.

    Keeps T2 live trees whose PREV_TRE_CN is in cohort_cns (T1 CNs of live,
    DIA_T1 >= 5 in trees).  DIA >= min_dia floor retained.  Returns (BA, n_trees).
    """
    t = tree_df.loc[tree_df["PLT_CN"] == t2_cn].copy()
    if "PREV_TRE_CN" in t.columns:
        t = t.loc[t["PREV_TRE_CN"].notna() & t["PREV_TRE_CN"].isin(cohort_cns)]
    t = t.loc[t["DIA"].notna() & t["TPA_UNADJ"].notna() & (t["DIA"] >= min_dia)]
    if len(t) == 0:
        return float("nan"), 0
    ba = float(
        (t["TPA_UNADJ"].astype(float) * 0.005454154 * t["DIA"].astype(float) ** 2).sum()
    )
    return ba, int(len(t))


def obs_vol_mcuft_bdft_cohort(
    tree_df: pd.DataFrame, t2_cn: str, cohort_cns: set, min_dia: float = 5.0
) -> tuple[float, float]:
    """Observed T2 volume (Mcuft/ac, bdft/ac) for surviving T1 cohort members."""
    t = tree_df.loc[tree_df["PLT_CN"] == t2_cn].copy()
    if "PREV_TRE_CN" in t.columns:
        t = t.loc[t["PREV_TRE_CN"].notna() & t["PREV_TRE_CN"].isin(cohort_cns)]
    t = t.loc[t["DIA"].notna() & (t["DIA"].astype(float) >= min_dia)]
    mcuft: float = float("nan")
    bdft:  float = float("nan")
    if "VOLCFNET" in t.columns:
        cf = t.loc[t["VOLCFNET"].notna() & t["TPA_UNADJ"].notna()]
        if len(cf) > 0:
            mcuft = float((cf["TPA_UNADJ"].astype(float) * cf["VOLCFNET"].astype(float)).sum())
    if "VOLBFNET" in t.columns:
        bf = t.loc[t["VOLBFNET"].notna() & t["TPA_UNADJ"].notna()]
        if len(bf) > 0:
            bdft = float((bf["TPA_UNADJ"].astype(float) * bf["VOLBFNET"].astype(float)).sum())
    return mcuft, bdft


def obs_ba_ft2ac_wholestand(
    tree_df: pd.DataFrame, t2_cn: str, min_dia: float = 5.0
) -> float:
    """Whole-stand observed T2 BA (all live DIA>=5 trees). Diagnostic only."""
    t = tree_df.loc[tree_df["PLT_CN"] == t2_cn].copy()
    t = t.loc[t["DIA"].notna() & t["TPA_UNADJ"].notna() & (t["DIA"] >= min_dia)]
    if len(t) == 0:
        return float("nan")
    return float(
        (t["TPA_UNADJ"].astype(float) * 0.005454154 * t["DIA"].astype(float) ** 2).sum()
    )


# ---------------------------------------------------------------------------
# Bias table printing
# ---------------------------------------------------------------------------
def _bias_stats(grp: pd.DataFrame, pred_col: str, obs_col: str):
    ok = (
        grp[pred_col].notna() & grp[obs_col].notna()
        & (grp[obs_col] > 0) & (grp[pred_col] > 0)
    )
    n = int(ok.sum())
    if n == 0:
        return 0, float("nan"), float("nan")
    p = grp.loc[ok, pred_col].astype(float)
    o = grp.loc[ok, obs_col].astype(float)
    median_bias = float(np.median(100.0 * (p / o - 1.0)))
    agg_bias    = float(100.0 * (p.mean() / o.mean() - 1.0))
    return n, median_bias, agg_bias


def print_stratum_table(df: pd.DataFrame, active_models: list[str], label: str) -> None:
    print("\n" + "=" * 104)
    print(f"  STRATUM = {label}   (n_rows={len(df)})")
    print("=" * 104)
    hdr = (f"  {'model':<15} {'n':>5} | "
           f"{'BA_med%':>8} {'BA_agg%':>8} | "
           f"{'MCuFt_med%':>10} {'MCuFt_agg%':>10} | "
           f"{'BdFt_med%':>10} {'BdFt_agg%':>10}")
    print(hdr)
    print("-" * 104)

    def _f(v):
        return f"{v:+.1f}" if (v == v and abs(v) != float("inf")) else "  n/a"

    for m in active_models:
        grp = df.loc[df["model"] == m]
        if len(grp) == 0:
            print(f"  {m:<15} {'0':>5} |  (no rows)")
            continue
        n_ba, ba_med, ba_agg       = _bias_stats(grp, "BA_PRED", "BA_OBS")
        _,    mc_med, mc_agg       = _bias_stats(grp, "MCuFt_PRED", "MCuFt_OBS")
        _,    bf_med, bf_agg       = _bias_stats(grp, "BdFt_PRED", "BdFt_OBS")
        print(
            f"  {m:<15} {n_ba:>5} | "
            f"{_f(ba_med):>8} {_f(ba_agg):>8} | "
            f"{_f(mc_med):>10} {_f(mc_agg):>10} | "
            f"{_f(bf_med):>10} {_f(bf_agg):>10}"
        )
    print("=" * 104)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    global FIA_DIR
    ap = argparse.ArgumentParser(
        description="Cohort-matched forest-type-stratified 5-model FVS scorecard"
    )
    ap.add_argument("--variant",   required=True, choices=list(VARIANTS))
    ap.add_argument("--max-pairs", type=int, default=200)
    ap.add_argument("--seed",      type=int, default=42)
    ap.add_argument(
        "--outdir",
        default="/users/PUOM0008/crsfaaron/fvs-modern/calibration/output/fia_scorecard",
    )
    ap.add_argument(
        "--fia-dir",
        default=FIA_DIR,
        help="Root directory containing per-state FIA CSV files (default: %(default)s)",
    )
    args = ap.parse_args()
    FIA_DIR = args.fia_dir

    variant = args.variant
    states  = VARIANTS[variant]["states"]

    active_models: list[str] = []
    skipped_models: list[str] = []
    for m in MODEL_ORDER:
        vc = STANDALONE_MODELS[m]["variant_compatible"]
        if vc is None or variant in vc:
            active_models.append(m)
        else:
            skipped_models.append(m)

    if skipped_models:
        print(f"[cohort] NOTE: skipping {skipped_models} "
              f"(not variant-compatible with {variant})")
    print(f"[cohort] variant={variant}  states={states}")
    print(f"[cohort] active models={active_models}")
    sys.stdout.flush()

    pairs = load_pairs(variant, max_pairs=args.max_pairs, seed=args.seed)
    if len(pairs) == 0:
        print(f"[cohort] ERROR: no remeasurement pairs found for {variant}")
        sys.exit(1)
    print(f"[cohort] found {len(pairs)} pairs (interval 5-15 yr), "
          f"sampling up to {args.max_pairs}")
    sys.stdout.flush()

    all_cns = set(pairs["PREV_PLT_CN"].tolist()) | set(pairs["CN"].tolist())
    trees   = load_trees_for_cns(states, all_cns)
    conds   = load_cond_for_cns(states, all_cns)
    print(f"[cohort] loaded {len(trees)} tree recs, {len(conds)} cond recs "
          f"for {len(all_cns)} plot CNs")
    has_prev = "PREV_TRE_CN" in trees.columns and trees["PREV_TRE_CN"].notna().any()
    print(f"[cohort] PREV_TRE_CN linkage available: {has_prev} "
          f"(method = {'PREV_TRE_CN' if has_prev else 'SUBP+TREE fallback'})")
    sys.stdout.flush()

    rows: list[dict] = []
    n_mdl   = len(active_models)
    n_total = len(pairs) * n_mdl
    n_done  = 0
    fortype_counts: dict[str, int] = {"SW": 0, "MW": 0, "HW": 0, "UNK": 0}

    for _, pr in pairs.iterrows():
        t1_cn     = int(pr["PREV_PLT_CN"])
        t2_cn     = int(pr["CN"])
        t1_cn_s   = str(t1_cn)
        t2_cn_s   = str(t2_cn)
        ncyc      = int(pr["interval"])
        state_cd  = int(pr.get("STATECD", 0)) if pd.notna(pr.get("STATECD")) else 0
        measyear2 = int(pr["MEASYEAR"])
        measyear1 = measyear2 - ncyc

        t1 = trees.loc[trees["PLT_CN"] == t1_cn_s].copy()
        if len(t1) < 3:
            n_done += n_mdl
            continue

        # Forest-type classification from T1 cohort (DIA>=5)
        fortype, sw_share = classify_fortype(t1)
        fortype_counts[fortype] = fortype_counts.get(fortype, 0) + 1

        # T1 cohort CN set = live T1 trees with DIA>=5 (the FVS input list)
        t1_cohort = t1.loc[t1["DIA"].notna() & (t1["DIA"].astype(float) >= 5.0)]
        cohort_cns = set(t1_cohort["CN"].dropna().tolist())

        cond_row = None
        c1 = conds.loc[conds["PLT_CN"] == t1_cn]
        if len(c1) > 0:
            cond_row = c1.iloc[0].to_dict()

        # COHORT-MATCHED observed values (T2 survivors of the T1 cohort)
        ba_obs, n_cohort   = obs_ba_ft2ac_cohort(trees, t2_cn_s, cohort_cns)
        mcuft_obs, bdft_obs = obs_vol_mcuft_bdft_cohort(trees, t2_cn_s, cohort_cns)
        # Diagnostic: whole-stand observed (all live DIA>=5 at T2)
        ba_obs_whole       = obs_ba_ft2ac_wholestand(trees, t2_cn_s)

        sid      = f"{variant.upper()}_{t1_cn % 10_000_000:07d}"
        stand_df = build_standinit(sid, measyear1, variant, cond_row)
        tree_df  = build_treeinit(t1, sid)
        if len(tree_df) == 0:
            n_done += n_mdl
            continue

        state_abbrev = STATECD_TO_ABBREV.get(state_cd)

        for model_name in active_models:
            mdl = STANDALONE_MODELS[model_name]

            restrict_dg = mdl.get("states_restrict_dg")
            if restrict_dg and state_abbrev not in restrict_dg:
                eff_extra_env: dict | None = {
                    k: v for k, v in mdl["extra_env"].items()
                    if k not in _DG_KEYS
                }
                base_cc = mdl["calib_components"] or {}
                eff_calib = {**base_cc, "dg": True}
            else:
                eff_extra_env = None
                eff_calib     = None

            result_df, err = run_stand(
                sid, stand_df, tree_df, variant, model_name, ncyc,
                extra_env_override=eff_extra_env,
                calib_override=eff_calib,
            )
            n_done += 1

            if err:
                if n_done % 200 == 0:
                    print(f"  [{n_done}/{n_total}] {sid} {model_name}: ERR {err}")
                    sys.stdout.flush()
                continue
            if result_df is None or len(result_df) == 0:
                continue

            last = result_df.iloc[-1]
            rows.append({
                "PLOT":       t1_cn,
                "STATE":      state_cd,
                "MEASYEAR1":  measyear1,
                "MEASYEAR2":  measyear2,
                "PERIOD_YR":  ncyc,
                "variant":    variant,
                "model":      model_name,
                "FORTYPE":    fortype,
                "SW_BA_SHARE": round(sw_share, 4) if sw_share == sw_share else float("nan"),
                "N_COHORT":   n_cohort,
                "BA_PRED":    float(last.get("BA",     float("nan"))),
                "TPA_PRED":   float(last.get("Tpa",    float("nan"))),
                "QMD_PRED":   float(last.get("QMD",    float("nan"))),
                "MCuFt_PRED": float(last.get("MCuFt",  float("nan"))),
                "BdFt_PRED":  float(last.get("BdFt",   float("nan"))),
                "BA_OBS":       ba_obs,          # cohort-matched
                "BA_OBS_WHOLE": ba_obs_whole,    # diagnostic whole-stand
                "MCuFt_OBS":  mcuft_obs,
                "BdFt_OBS":   bdft_obs,
            })

            if n_done % 500 == 0 or n_done == n_total:
                ba_p = rows[-1]["BA_PRED"] if rows else float("nan")
                print(
                    f"  [{n_done}/{n_total}] {sid} {model_name} [{fortype}]: "
                    f"BA_pred={ba_p:.1f}  BA_obs={ba_obs:.1f}  n_cohort={n_cohort}"
                )
                sys.stdout.flush()

    out_df = pd.DataFrame(rows)
    os.makedirs(args.outdir, exist_ok=True)
    outfile = os.path.join(args.outdir, f"scorecard_{variant}_5model_htdbh.csv")
    out_df.to_csv(outfile, index=False)
    print(f"\n[htdbh] wrote {len(out_df)} rows -> {outfile}")
    print(f"[htdbh] forest-type plot counts (per unique pair): {fortype_counts}")

    _n_imp = _HT_IMPUTE_STATS["imputed"]
    _n_meas = _HT_IMPUTE_STATS["measured"]
    _n_tot = _n_imp + _n_meas
    _frac = (100.0 * _n_imp / _n_tot) if _n_tot else float("nan")
    print(f"[htdbh] Wykoff-imputed heights: {_n_imp} of {_n_tot} input trees "
          f"({_frac:.1f}%); measured FIA heights kept: {_n_meas}")

    # ---------------------------------------------------------------------
    # Bias tables: OVERALL and per forest type
    # ---------------------------------------------------------------------
    print("\n" + "#" * 104)
    print(f"  COHORT-MATCHED 5-MODEL SCORECARD   variant={variant}   "
          f"n_pairs_attempted={len(pairs)}")
    print("  Columns: median (plot-level) and aggregate (mean/mean) percent bias, "
          "PRED vs cohort-matched OBS")
    print("#" * 104)

    if len(out_df) == 0:
        print("  (no usable rows)")
        return

    print_stratum_table(out_df, active_models, "OVERALL")
    for ft in ["SW", "MW", "HW"]:
        sub = out_df.loc[out_df["FORTYPE"] == ft]
        if len(sub) == 0:
            print(f"\n  (no rows for stratum {ft})")
            continue
        print_stratum_table(sub, active_models, ft)


if __name__ == "__main__":
    main()
