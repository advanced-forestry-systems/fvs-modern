#!/usr/bin/env python3
"""
run_fullsystem_fia_scorecard.py
================================
Full-system (DG + HG + CRW) scorecard on FIA remeasurement pairs.

Configurations:
  dg_only                - DG only, no HG, no CRW
  dg_hg_organon          - DG + HG Organon coefficients
  dg_hg_clim             - DG + HG CONUS+climate coefficients
  full_dg_hg_organon_crw - DG + HG Organon + CRW  (full system A)
  full_dg_hg_clim_crw    - DG + HG CONUS+climate + CRW  (full system B)
  hg_organon_crw         - HG Organon + CRW, no DG  (isolate CRW effect)

DG and CRW are NE-fitted; configs with states_restrict are skipped for
non-NE variants. Observed MCuFt/BdFt are joined inline from FIA TREE tables
using the two-step physical-plot lookup from vol_scorecard_fix.py.

Usage:
  python3 run_fullsystem_fia_scorecard.py --variant ne [--max-pairs 200] [--seed 42]

Output:
  {outdir}/scorecard_{variant}_fullsys.csv
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
FIA_DIR = "/fs/scratch/PUOM0008/crsfaaron/FIA"

_BIN  = "/users/PUOM0008/crsfaaron/fvs-modern/src-converted/bin/FVSne"
_CONF = "/fs/scratch/PUOM0008/crsfaaron/wt-ne-dg/config"

_HG_BASE = {
    "FVS_GREGHG":    "1",
    "FVS_GREG_TD":   "28",
    "FVS_GREG_EMT":  "-30",
    "FVS_GREG_ELEV": "1000",
}
_DG_BASE = {
    "FVS_GREGDG":      "1",
    "FVS_GREGDG_COEF": f"{_CONF}/greg_dg_coefficients.csv",
    "FVS_GREG_DD0":    "1600",
    "FVS_GREG_TD":     "28",
    "FVS_GREG_PPT_SM": "230",
    "FVS_GREG_DD18":   "3500",
}
_CRW_BASE = {
    "FVS_GREGCRW": "1",
    "FVS_GREGCRW_COEF": f"{_CONF}/greg_crown_change_coefficients.csv",
}

# ---------------------------------------------------------------------------
# Full-system configurations
# ---------------------------------------------------------------------------
FULLSYS_CONFIGS: dict[str, dict] = {
    # DG only (no HG, no CRW)
    "dg_only": {
        "binary": _BIN,
        "env": {**_DG_BASE},
        "states_restrict": ["CT", "MA", "ME", "NH", "NY", "RI", "VT"],
    },
    # DG + HG Organon (no CRW)
    "dg_hg_organon": {
        "binary": _BIN,
        "env": {
            **_DG_BASE, **_HG_BASE,
            "FVS_GREGHG_COEF": f"{_CONF}/greg_hg_coefficients.csv",
        },
        "states_restrict": ["CT", "MA", "ME", "NH", "NY", "RI", "VT"],
    },
    # DG + HG CONUS+climate (no CRW)
    "dg_hg_clim": {
        "binary": _BIN,
        "env": {
            **_DG_BASE, **_HG_BASE,
            "FVS_GREGHG_COEF": f"{_CONF}/greg_hg_coefficients_compound_climate.csv",
        },
        "states_restrict": ["CT", "MA", "ME", "NH", "NY", "RI", "VT"],
    },
    # DG + HG Organon + CRW  (full system A)
    "full_dg_hg_organon_crw": {
        "binary": _BIN,
        "env": {
            **_DG_BASE, **_HG_BASE,
            "FVS_GREGHG_COEF": f"{_CONF}/greg_hg_coefficients.csv",
            **_CRW_BASE,
        },
        "states_restrict": ["CT", "MA", "ME", "NH", "NY", "RI", "VT"],
    },
    # DG + HG CONUS+climate + CRW  (full system B)
    "full_dg_hg_clim_crw": {
        "binary": _BIN,
        "env": {
            **_DG_BASE, **_HG_BASE,
            "FVS_GREGHG_COEF": f"{_CONF}/greg_hg_coefficients_compound_climate.csv",
            **_CRW_BASE,
        },
        "states_restrict": ["CT", "MA", "ME", "NH", "NY", "RI", "VT"],
    },
    # HG Organon + CRW, no DG  (isolate CRW effect; no state restriction)
    "hg_organon_crw": {
        "binary": _BIN,
        "env": {
            **_HG_BASE,
            "FVS_GREGHG_COEF": f"{_CONF}/greg_hg_coefficients.csv",
            **_CRW_BASE,
        },
    },
}

# ---------------------------------------------------------------------------
# Variant definitions (same as run_7config_fia_scorecard.py)
# ---------------------------------------------------------------------------
VARIANTS: dict[str, dict] = {
    "ne": {"states": ["CT", "DE", "MA", "MD", "ME", "NH", "NJ", "NY", "PA", "RI", "VT", "WV"]},
    "sn": {"states": ["AL", "GA", "MS", "SC"]},
    "wc": {"states": ["OR", "WA"]},
    "em": {"states": ["MT", "ND", "SD"]},
    "pn": {"states": ["OR", "WA"]},
}

# ---------------------------------------------------------------------------
# FIA STATECD (numeric) -> state abbreviation
# ---------------------------------------------------------------------------
STATECD_TO_ABBREV: dict[int, str] = {
    1:  "AL",  2:  "AK",  4:  "AZ",  5:  "AR",  6:  "CA",
    8:  "CO",  9:  "CT", 10: "DE",  12: "FL",  13: "GA",
    16: "ID", 17: "IL",  18: "IN",  19: "IA",  20: "KS",
    21: "KY", 22: "LA",  23: "ME",  24: "MD",  25: "MA",
    26: "MI", 27: "MN",  28: "MS",  29: "MO",  30: "MT",
    31: "NE", 32: "NV",  33: "NH",  34: "NJ",  35: "NM",
    36: "NY", 37: "NC",  38: "ND",  39: "OH",  40: "OK",
    41: "OR", 42: "PA",  44: "RI",  45: "SC",  46: "SD",
    47: "TN", 48: "TX",  49: "UT",  50: "VT",  51: "VA",
    53: "WA", 54: "WV",  55: "WI",  56: "WY",
}

# Reverse map: abbreviation -> set of STATECDs (usually one-to-one)
_ABBREV_TO_STATECDS: dict[str, set[int]] = {}
for _cd, _ab in STATECD_TO_ABBREV.items():
    _ABBREV_TO_STATECDS.setdefault(_ab, set()).add(_cd)

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
** DEFAULT PARAMETERS
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
    """Build a single-row fvs_standinit from FIA COND metadata."""
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
# Tree list builder
# ---------------------------------------------------------------------------
def build_treeinit(tree_df: pd.DataFrame, sid: str) -> pd.DataFrame:
    """Convert FIA TREE records (live, DIA >= 5 in) to fvs_treeinit."""
    rows = []
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
        ht = float(ht_raw) if (
            ht_raw is not None and pd.notna(ht_raw) and float(ht_raw) > 0
        ) else 0.0
        cr_raw = r.get("CR")
        cr = float(cr_raw) if (cr_raw is not None and pd.notna(cr_raw)) else 50.0
        if cr <= 9:
            cr = cr * 10
        cr = int(max(15, min(99, cr)))
        rows.append({
            "stand_id":   sid,
            "plot_id":    1,
            "tree_id":    len(rows) + 1,
            "tree_count": round(tpa, 4),
            "species":    spcd,
            "diameter":   round(dbh, 2),
            "ht":         round(ht, 1),
            "crratio":    cr,
        })
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# FVS runner
# ---------------------------------------------------------------------------
def run_stand(
    sid:      str,
    stand_df: pd.DataFrame,
    tree_df:  pd.DataFrame,
    variant:  str,
    cfg_name: str,
    ncyc:     int,
) -> tuple[pd.DataFrame | None, str | None]:
    """Run FVS for one stand/config; return (Summary2_df, error_str)."""
    cfg    = FULLSYS_CONFIGS[cfg_name]
    binary = cfg["binary"]
    if not os.path.exists(binary):
        return None, f"binary not found: {binary}"
    if len(tree_df) == 0:
        return None, "empty tree list"

    env = os.environ.copy()
    env.update(cfg["env"])

    tmp = tempfile.mkdtemp(prefix="fvsfullsys_")
    try:
        db  = os.path.join(tmp, "FVS_Data.db")
        con = sqlite3.connect(db)
        stand_df.to_sql("fvs_standinit", con, if_exists="replace", index=False)
        tree_df.to_sql("fvs_treeinit",   con, if_exists="replace", index=False)
        con.close()

        key = os.path.join(tmp, "run.key")
        with open(key, "w") as fh:
            fh.write(KEYFILE.format(sid=sid, db=db, ncyc=ncyc))

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
_TREE_COLS = {"PLT_CN", "STATECD", "COUNTYCD", "PLOT", "SUBP", "TREE",
              "STATUSCD", "SPCD", "DIA", "HT", "CR", "TPA_UNADJ"}
_COND_COLS = {"PLT_CN", "STATECD", "COUNTYCD", "INVYR", "COND_STATUS_CD",
              "SICOND", "BALIVE", "FORTYPCD", "STDAGE", "SLOPE", "ELEV"}


def load_pairs(variant: str, max_pairs: int = 200, seed: int = 42) -> pd.DataFrame:
    """Find FIA remeasurement pairs (gap 5-15 yr) for this variant's states."""
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
    cns_int = {int(c) for c in cns}
    frames  = []
    for st in states:
        f = os.path.join(FIA_DIR, f"{st}_TREE.csv")
        if not os.path.exists(f):
            continue
        df = pd.read_csv(f, usecols=lambda c: c in _TREE_COLS, low_memory=False)
        df["PLT_CN"] = df["PLT_CN"].astype("int64")
        sub = df.loc[df["PLT_CN"].isin(cns_int)].copy()
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


def obs_ba_ft2ac(tree_df: pd.DataFrame, plt_cn: int) -> float:
    """Observed BA (ft2/ac) from FIA TREE records at the remeasurement year."""
    t = tree_df.loc[tree_df["PLT_CN"] == plt_cn].copy()
    t = t.loc[t["DIA"].notna() & t["TPA_UNADJ"].notna()]
    if len(t) == 0:
        return float("nan")
    return float((t["TPA_UNADJ"] * 0.005454154 * t["DIA"] ** 2).sum())


# ---------------------------------------------------------------------------
# Observed volume join (two-step FIA PLOT->TREE, same as vol_scorecard_fix.py)
# ---------------------------------------------------------------------------
_VOL_PLOT_COLS = {"CN", "STATECD", "COUNTYCD", "PLOT", "MEASYEAR"}
_VOL_TREE_COLS = {"PLT_CN", "STATUSCD", "VOLCFNET", "VOLBFNET", "TPA_UNADJ"}


def get_obs_volumes(sc_sub: pd.DataFrame, states: list[str]) -> pd.DataFrame:
    """Return DataFrame [_plot_key, MEASYEAR2, MCuFt_OBS, BdFt_OBS].

    _plot_key is the PREV_PLT_CN (measyear1 CN) as string.
    Uses the two-step physical-plot join from vol_scorecard_fix.py.
    """
    obs_rows: list[pd.DataFrame] = []

    for st in states:
        plot_file = os.path.join(FIA_DIR, f"{st}_PLOT.csv")
        tree_file = os.path.join(FIA_DIR, f"{st}_TREE.csv")
        if not os.path.exists(plot_file) or not os.path.exists(tree_file):
            continue

        # Which scorecard rows belong to this state?
        state_cds = _ABBREV_TO_STATECDS.get(st, set())
        state_sc  = sc_sub[sc_sub["STATE"].isin(state_cds)][
            ["_plot_key", "MEASYEAR2"]
        ].drop_duplicates()
        if state_sc.empty:
            continue

        measyear1_cns = set(state_sc["_plot_key"].astype(str))
        needed_myr2   = set(state_sc["MEASYEAR2"].astype(int))

        # Step A: resolve physical (COUNTYCD, PLOT) from the PLOT table
        try:
            plot_df = pd.read_csv(
                plot_file, usecols=lambda c: c in _VOL_PLOT_COLS, low_memory=False,
            )
        except Exception as exc:
            print(f"  WARNING: could not read {plot_file}: {exc}", flush=True)
            continue

        plot_df["CN"]       = plot_df["CN"].astype(str)
        plot_df["MEASYEAR"] = pd.to_numeric(plot_df["MEASYEAR"], errors="coerce")
        plot_df["COUNTYCD"] = pd.to_numeric(
            plot_df["COUNTYCD"], errors="coerce"
        ).astype("Int64")
        plot_df["PLOT"]     = pd.to_numeric(
            plot_df["PLOT"], errors="coerce"
        ).astype("Int64")

        m1_rows    = plot_df[plot_df["CN"].isin(measyear1_cns)][
            ["CN", "COUNTYCD", "PLOT"]
        ].copy()
        cn_to_phys = {row.CN: (row.COUNTYCD, row.PLOT) for _, row in m1_rows.iterrows()}
        if not cn_to_phys:
            continue

        needed_physplots = set(cn_to_phys.values())
        m2_cands = plot_df[plot_df["MEASYEAR"].isin(needed_myr2)].copy()
        m2_cands["phys"] = list(zip(m2_cands["COUNTYCD"], m2_cands["PLOT"]))
        m2_cands = m2_cands[m2_cands["phys"].isin(needed_physplots)]
        phys_yr_to_cn2: dict[tuple, str] = {
            (row.COUNTYCD, row.PLOT, int(row.MEASYEAR)): row.CN
            for _, row in m2_cands.iterrows()
            if pd.notna(row.MEASYEAR)
        }

        sc_to_cn2: dict[tuple, str] = {}
        for _, row in state_sc.iterrows():
            m1_cn = str(row["_plot_key"])
            myr2  = int(row["MEASYEAR2"])
            phys  = cn_to_phys.get(m1_cn)
            if phys is None:
                continue
            cn2 = phys_yr_to_cn2.get((phys[0], phys[1], myr2))
            if cn2 is None:
                continue
            sc_to_cn2[(m1_cn, myr2)] = cn2

        if not sc_to_cn2:
            continue

        needed_cn2s = set(sc_to_cn2.values())

        # Step B: aggregate volume from TREE table
        try:
            chunks = []
            for chunk in pd.read_csv(
                tree_file,
                usecols=lambda c: c in _VOL_TREE_COLS,
                chunksize=400_000,
                low_memory=False,
            ):
                chunk["PLT_CN"] = chunk["PLT_CN"].astype(str)
                filt = chunk[chunk["PLT_CN"].isin(needed_cn2s)]
                if len(filt):
                    chunks.append(filt)
        except Exception as exc:
            print(f"  WARNING: could not read {tree_file}: {exc}", flush=True)
            continue

        if not chunks:
            continue

        trees = pd.concat(chunks, ignore_index=True)
        trees = trees[trees["STATUSCD"] == 1].copy()
        for col in ["VOLCFNET", "VOLBFNET", "TPA_UNADJ"]:
            trees[col] = pd.to_numeric(trees[col], errors="coerce").fillna(0.0)
        trees["cf_contrib"] = trees["VOLCFNET"] * trees["TPA_UNADJ"]
        trees["bf_contrib"] = trees["VOLBFNET"] * trees["TPA_UNADJ"]

        agg = (
            trees.groupby("PLT_CN")
                 .agg(MCuFt_OBS=("cf_contrib", "sum"), BdFt_OBS=("bf_contrib", "sum"))
                 .reset_index()
        )
        agg["PLT_CN"] = agg["PLT_CN"].astype(str)
        cn2_to_sc = {v: k for k, v in sc_to_cn2.items()}
        agg["_plot_key"] = agg["PLT_CN"].map(
            lambda x: cn2_to_sc.get(x, (None, None))[0]
        )
        agg["MEASYEAR2"] = agg["PLT_CN"].map(
            lambda x: cn2_to_sc.get(x, (None, None))[1]
        )
        agg = agg.dropna(subset=["_plot_key"])
        obs_rows.append(agg[["_plot_key", "MEASYEAR2", "MCuFt_OBS", "BdFt_OBS"]])

    if not obs_rows:
        return pd.DataFrame(columns=["_plot_key", "MEASYEAR2", "MCuFt_OBS", "BdFt_OBS"])

    obs = pd.concat(obs_rows, ignore_index=True)
    obs["_plot_key"] = obs["_plot_key"].astype(str)
    obs["MEASYEAR2"] = obs["MEASYEAR2"].astype(int)
    return obs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(
        description="Full-system (DG+HG+CRW) FVS scorecard on FIA remeasurement pairs"
    )
    ap.add_argument("--variant",   default="ne", choices=list(VARIANTS))
    ap.add_argument("--max-pairs", type=int, default=200)
    ap.add_argument("--seed",      type=int, default=42)
    ap.add_argument(
        "--outdir",
        default="/users/PUOM0008/crsfaaron/fvs-modern/calibration/output/fia_scorecard",
    )
    args = ap.parse_args()

    variant = args.variant
    states  = VARIANTS[variant]["states"]

    # Which configs apply to this variant?
    configs = [
        c for c, meta in FULLSYS_CONFIGS.items()
        if not meta.get("states_restrict") or variant == "ne"
    ]
    skipped = [
        c for c, meta in FULLSYS_CONFIGS.items()
        if meta.get("states_restrict") and variant != "ne"
    ]
    if skipped:
        print(
            f"[scorecard] NOTE: skipping {skipped} "
            f"(NE-DG-restricted configs, variant={variant})"
        )

    print(f"[scorecard] variant={variant}  states={states}")
    print(f"[scorecard] configs={configs}")
    sys.stdout.flush()

    # Load remeasurement pairs
    pairs = load_pairs(variant, max_pairs=args.max_pairs, seed=args.seed)
    if len(pairs) == 0:
        print(f"[scorecard] ERROR: no remeasurement pairs found for {variant}")
        sys.exit(1)
    print(
        f"[scorecard] found {len(pairs)} pairs (interval 5-15 yr), "
        f"sampling up to {args.max_pairs}"
    )
    sys.stdout.flush()

    # Load all FIA trees + conds in bulk
    all_cns = set(pairs["PREV_PLT_CN"].tolist()) | set(pairs["CN"].tolist())
    trees   = load_trees_for_cns(states, all_cns)
    conds   = load_cond_for_cns(states, all_cns)
    print(
        f"[scorecard] loaded {len(trees)} tree recs, {len(conds)} cond recs "
        f"for {len(all_cns)} plot CNs"
    )
    sys.stdout.flush()

    rows: list[dict] = []
    n_cfg   = len(configs)
    n_total = len(pairs) * n_cfg
    n_done  = 0

    for _, pr in pairs.iterrows():
        t1_cn     = int(pr["PREV_PLT_CN"])
        t2_cn     = int(pr["CN"])
        ncyc      = int(pr["interval"])
        state_cd  = int(pr.get("STATECD", 0)) if pd.notna(pr.get("STATECD")) else 0
        measyear2 = int(pr["MEASYEAR"])
        measyear1 = measyear2 - ncyc

        t1 = trees.loc[trees["PLT_CN"] == t1_cn].copy()
        if len(t1) < 3:
            n_done += n_cfg
            continue

        cond_row = None
        c1 = conds.loc[conds["PLT_CN"] == t1_cn]
        if len(c1) > 0:
            cond_row = c1.iloc[0].to_dict()

        ba_obs = obs_ba_ft2ac(trees, t2_cn)

        sid      = f"{variant.upper()}_{t1_cn % 10_000_000:07d}"
        stand_df = build_standinit(sid, measyear1, variant, cond_row)
        tree_df  = build_treeinit(t1, sid)
        if len(tree_df) == 0:
            n_done += n_cfg
            continue

        for cfg_name in configs:
            cfg = FULLSYS_CONFIGS[cfg_name]
            # Honour per-config state restriction (DG northern NE only)
            state_restrict = cfg.get("states_restrict")
            if state_restrict:
                state_abbrev = STATECD_TO_ABBREV.get(state_cd)
                if state_abbrev not in state_restrict:
                    n_done += 1
                    continue

            result_df, err = run_stand(sid, stand_df, tree_df, variant, cfg_name, ncyc)
            n_done += 1

            if err:
                if n_done % 20 == 0:
                    print(f"  [{n_done}/{n_total}] {sid} {cfg_name}: ERR {err}")
                    sys.stdout.flush()
                continue
            if result_df is None or len(result_df) == 0:
                continue

            last = result_df.iloc[-1]
            row  = {
                "PLOT":       t1_cn,
                "STATE":      state_cd,
                "MEASYEAR1":  measyear1,
                "MEASYEAR2":  measyear2,
                "PERIOD_YR":  ncyc,
                "variant":    variant,
                "config":     cfg_name,
                "BA_PRED":    float(last.get("BA",    float("nan"))),
                "TPA_PRED":   float(last.get("Tpa",   float("nan"))),
                "QMD_PRED":   float(last.get("QMD",   float("nan"))),
                "MCuFt_PRED": float(last.get("MCuFt", float("nan"))),
                "BdFt_PRED":  float(last.get("BdFt",  float("nan"))),
                "BA_OBS":     ba_obs,
            }
            rows.append(row)

            if n_done % 100 == 0 or n_done == n_total:
                print(
                    f"  [{n_done}/{n_total}] {sid} {cfg_name}: "
                    f"BA_pred={row['BA_PRED']:.1f}  BA_obs={ba_obs:.1f}"
                )
                sys.stdout.flush()

    # Build output DataFrame
    out_df = pd.DataFrame(rows)
    os.makedirs(args.outdir, exist_ok=True)
    outfile = os.path.join(args.outdir, f"scorecard_{variant}_fullsys.csv")

    # Join observed MCuFt/BdFt via two-step FIA PLOT->TREE lookup
    if len(out_df) > 0:
        print("\n[scorecard] joining observed MCuFt/BdFt from FIA ...", flush=True)
        sc_for_vol = out_df[["PLOT", "MEASYEAR2", "STATE"]].copy()
        sc_for_vol["_plot_key"] = sc_for_vol["PLOT"].astype(str)
        obs_vol = get_obs_volumes(sc_for_vol, states)

        out_df["_plot_key"] = out_df["PLOT"].astype(str)
        out_df = out_df.merge(
            obs_vol,
            on=["_plot_key", "MEASYEAR2"],
            how="left",
        ).drop(columns=["_plot_key"])
        n_vol = out_df["MCuFt_OBS"].notna().sum()
        print(
            f"[scorecard] volume join: {n_vol}/{len(out_df)} rows have observed volume "
            f"({100 * n_vol / len(out_df):.1f}%)",
            flush=True,
        )

    out_df.to_csv(outfile, index=False)
    print(f"\n[scorecard] wrote {len(out_df)} rows -> {outfile}")

    # ---------------------------------------------------------------------------
    # Scorecard summary table
    # ---------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print(f"  SCORECARD   variant = {variant}   n_pairs_attempted = {len(pairs)}")
    print("=" * 80)
    hdr = (
        f"  {'config':<24} {'n':>6} {'BA_bias%':>10} "
        f"{'BA_RMSE':>9} {'MCuFt_bias%':>13} {'BdFt_bias%':>12}"
    )
    print(hdr)
    print("-" * 80)

    has_vol = "MCuFt_OBS" in out_df.columns

    for cfg_name in configs:
        grp  = out_df.loc[out_df["config"] == cfg_name].copy()
        ok_ba = (
            grp["BA_PRED"].notna()
            & grp["BA_OBS"].notna()
            & (grp["BA_OBS"] > 0)
            & (grp["BA_PRED"] > 0)
        )
        n = int(ok_ba.sum())
        if n == 0:
            print(f"  {cfg_name:<24} {'0':>6}  (no usable predictions)")
            continue
        bias_pct = 100.0 * (
            grp.loc[ok_ba, "BA_PRED"].mean() / grp.loc[ok_ba, "BA_OBS"].mean() - 1
        )
        rmse = float(
            np.sqrt(((grp.loc[ok_ba, "BA_PRED"] - grp.loc[ok_ba, "BA_OBS"]) ** 2).mean())
        )
        cf_bias = bf_bias = float("nan")
        if has_vol:
            ok_cf = (
                grp["MCuFt_PRED"].notna()
                & grp["MCuFt_OBS"].notna()
                & (grp["MCuFt_OBS"] > 0)
            )
            if ok_cf.sum() > 0:
                cf_bias = 100.0 * (
                    grp.loc[ok_cf, "MCuFt_PRED"].mean()
                    / grp.loc[ok_cf, "MCuFt_OBS"].mean()
                    - 1
                )
            ok_bf = (
                grp["BdFt_PRED"].notna()
                & grp["BdFt_OBS"].notna()
                & (grp["BdFt_OBS"] > 0)
            )
            if ok_bf.sum() > 0:
                bf_bias = 100.0 * (
                    grp.loc[ok_bf, "BdFt_PRED"].mean()
                    / grp.loc[ok_bf, "BdFt_OBS"].mean()
                    - 1
                )
        cf_str = f"{cf_bias:+.1f}" if math.isfinite(cf_bias) else "NA"
        bf_str = f"{bf_bias:+.1f}" if math.isfinite(bf_bias) else "NA"
        print(
            f"  {cfg_name:<24} {n:>6} {bias_pct:>+10.2f} {rmse:>9.2f} "
            f"{cf_str:>13} {bf_str:>12}"
        )
    print("=" * 80)


if __name__ == "__main__":
    main()
