#!/usr/bin/env python3
"""
run_7config_fia_scorecard.py
============================
7-configuration x 5-variant FVS scorecard using FIA remeasurement pairs.

Configurations tested per variant:
  default        - stock FVS parameters
  calibrated     - Bayesian posterior keyword injection via FvsConfigLoader
  greg_dg        - Greg DG modifier binary (NE northern states only)
  hg_organon     - Greg HG Organon-calibrated coefficients
  hg_conus_sp    - Greg HG CONUS species-free
  hg_conus_clim  - Greg HG CONUS + climate

NOTE: hg_bayes removed 2026-07-10. No genuine Bayesian HG coefficient set
exists on disk; greg_hg_coefficients_bayes.csv was byte-identical to
greg_hg_coefficients.csv (same B0-B8 for all 95 species, different
formatting only). Re-add when Bayesian HG fitting is complete.

Usage:
  python3 run_7config_fia_scorecard.py --variant ne [--max-pairs 200] [--seed 42]

Output:
  {outdir}/scorecard_{variant}.csv
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
_DG_BIN  = "/fs/scratch/PUOM0008/crsfaaron/wt-ne-dg/src-converted/bin/FVSne"
_HG_BASE = {
    "FVS_GREGHG":    "1",
    "FVS_GREG_TD":   "28",
    "FVS_GREG_EMT":  "-30",
    "FVS_GREG_ELEV": "1000",
}

# ---------------------------------------------------------------------------
# Variant to primary FIA state abbreviations
# ---------------------------------------------------------------------------
VARIANTS: dict[str, dict] = {
    "ne": {"states": ["CT", "DE", "MA", "MD", "ME", "NH", "NJ", "NY", "PA", "RI", "VT", "WV"]},
    "sn": {"states": ["AL", "GA", "MS", "SC"]},          # Southern variant: AL GA MS SC
    "wc": {"states": ["OR", "WA"]},                       # Westside Cascades: OR WA
    "em": {"states": ["MT", "ND", "SD"]},                 # Eastern Montana: MT ND SD
    "pn": {"states": ["OR", "WA"]},                       # Pacific Northwest: OR WA
}

# ---------------------------------------------------------------------------
# FIA STATECD (numeric) -> state abbreviation
# ---------------------------------------------------------------------------
STATECD_TO_ABBREV: dict[int, str] = {
    1:  "AL", 2:  "AK", 4:  "AZ", 5:  "AR", 6:  "CA",
    8:  "CO", 9:  "CT", 10: "DE", 12: "FL", 13: "GA",
    16: "ID", 17: "IL", 18: "IN", 19: "IA", 20: "KS",
    21: "KY", 22: "LA", 23: "ME", 24: "MD", 25: "MA",
    26: "MI", 27: "MN", 28: "MS", 29: "MO", 30: "MT",
    31: "NE", 32: "NV", 33: "NH", 34: "NJ", 35: "NM",
    36: "NY", 37: "NC", 38: "ND", 39: "OH", 40: "OK",
    41: "OR", 42: "PA", 44: "RI", 45: "SC", 46: "SD",
    47: "TN", 48: "TX", 49: "UT", 50: "VT", 51: "VA",
    53: "WA", 54: "WV", 55: "WI", 56: "WY",
}

# ---------------------------------------------------------------------------
# 7 Configurations
# ---------------------------------------------------------------------------
CONFIGS: dict[str, dict] = {
    "default": {
        "extra_env": {},
        "use_calib": False,
    },
    "calibrated": {
        "extra_env": {},
        "use_calib": True,
    },
    "greg_dg": {
        "extra_env": {
            "FVS_GREGDG":      "1",
            "FVS_GREGDG_COEF": f"{_CONF_DG}/greg_dg_coefficients.csv",
            "FVS_GREG_DD0":    "1600",
            "FVS_GREG_TD":     "28",
            "FVS_GREG_PPT_SM": "230",
            "FVS_GREG_DD18":   "3500",
        },
        "use_calib":       False,
        "binary_override": _DG_BIN,
        "ne_only":         True,
        # DG coefficients fit primarily on northern plots; restrict to states
        # where the model reduces bias (ME/NH/VT/NY/MA/CT/RI). Systematically
        # overpredicts BA in MD, WV, PA — excluded until southern re-fit.
        "states_restrict": ["CT", "MA", "ME", "NH", "NY", "RI", "VT"],
    },
    "hg_organon": {
        "extra_env": {
            **_HG_BASE,
            "FVS_GREGHG_COEF": f"{_CONF_DG}/greg_hg_coefficients.csv",
        },
        "use_calib": False,
    },
    "hg_conus_sp": {
        "extra_env": {
            **_HG_BASE,
            "FVS_GREGHG_COEF": f"{_CONF_DG}/greg_hg_coefficients_compound_none.csv",
        },
        "use_calib": False,
    },
    "hg_conus_clim": {
        "extra_env": {
            **_HG_BASE,
            "FVS_GREGHG_COEF": f"{_CONF_DG}/greg_hg_coefficients_compound_climate.csv",
        },
        "use_calib": False,
    },
}

# ---------------------------------------------------------------------------
# KEYFILE template (DATABASE/SUMMARY pattern from run_fvs_on_cfi.py)
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
        ht = float(ht_raw) if (ht_raw is not None and pd.notna(ht_raw) and float(ht_raw) > 0) else 0.0
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
# Calibrated keyword cache
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
        kw = "** DEFAULT (calibrated config not found)"
        _calib_kw_cache[variant] = kw
        return kw


# ---------------------------------------------------------------------------
# FVS runner (subprocess + SQLite, same pattern as run_fvs_on_cfi.py)
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
    cfg    = CONFIGS[cfg_name]
    binary = cfg.get("binary_override") or os.path.join(LIB_DIR, f"FVS{variant.lower()}")
    if not os.path.exists(binary):
        return None, f"binary not found: {binary}"
    if len(tree_df) == 0:
        return None, "empty tree list"

    calib_kw = (calibrated_keywords(variant)
                if cfg["use_calib"] else "** DEFAULT PARAMETERS")

    env = os.environ.copy()
    env.update(cfg["extra_env"])

    tmp = tempfile.mkdtemp(prefix="fvs7cfg_")
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
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser(description="7-config FVS scorecard on FIA pairs")
    ap.add_argument("--variant",   required=True, choices=list(VARIANTS))
    ap.add_argument("--max-pairs", type=int, default=200)
    ap.add_argument("--seed",      type=int, default=42)
    ap.add_argument(
        "--outdir",
        default="/users/PUOM0008/crsfaaron/fvs-modern/calibration/output/fia_scorecard",
    )
    ap.add_argument(
        "--out-suffix",
        default="",
        help="Suffix inserted before .csv in output filename (e.g. _fix)",
    )
    args = ap.parse_args()

    variant = args.variant
    states  = VARIANTS[variant]["states"]

    # Which configs apply to this variant?
    configs = [
        c for c, meta in CONFIGS.items()
        if not meta.get("ne_only") or variant == "ne"
    ]
    skipped = [c for c, meta in CONFIGS.items() if meta.get("ne_only") and variant != "ne"]
    if skipped:
        print(f"[scorecard] NOTE: skipping {skipped} (NE-only configs, variant={variant})")

    print(f"[scorecard] variant={variant}  states={states}")
    print(f"[scorecard] configs={configs}")
    sys.stdout.flush()

    # Load remeasurement pairs
    pairs = load_pairs(variant, max_pairs=args.max_pairs, seed=args.seed)
    if len(pairs) == 0:
        print(f"[scorecard] ERROR: no remeasurement pairs found for {variant}")
        sys.exit(1)
    print(f"[scorecard] found {len(pairs)} pairs (interval 5-15 yr), "
          f"sampling up to {args.max_pairs}")
    sys.stdout.flush()

    # Load all FIA trees + conds in bulk
    all_cns = set(pairs["PREV_PLT_CN"].tolist()) | set(pairs["CN"].tolist())
    trees   = load_trees_for_cns(states, all_cns)
    conds   = load_cond_for_cns(states, all_cns)
    print(f"[scorecard] loaded {len(trees)} tree recs, {len(conds)} cond recs "
          f"for {len(all_cns)} plot CNs")
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

        # Cond row for stand metadata
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
            cfg = CONFIGS[cfg_name]
            # Honour per-config state restriction (e.g. greg_dg northern NE only)
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
                "BA_PRED":    float(last.get("BA",     float("nan"))),
                "TPA_PRED":   float(last.get("Tpa",    float("nan"))),
                "QMD_PRED":   float(last.get("QMD",    float("nan"))),
                "MCuFt_PRED": float(last.get("MCuFt",  float("nan"))),
                "BdFt_PRED":  float(last.get("BdFt",   float("nan"))),
                "BA_OBS":     ba_obs,
            }
            rows.append(row)

            if n_done % 100 == 0 or n_done == n_total:
                ba_p = row["BA_PRED"]
                print(
                    f"  [{n_done}/{n_total}] {sid} {cfg_name}: "
                    f"BA_pred={ba_p:.1f}  BA_obs={ba_obs:.1f}"
                )
                sys.stdout.flush()

    # Write results
    out_df = pd.DataFrame(rows)
    os.makedirs(args.outdir, exist_ok=True)
    out_suffix = getattr(args, "out_suffix", "")
    outfile = os.path.join(args.outdir, f"scorecard_{variant}{out_suffix}.csv")
    out_df.to_csv(outfile, index=False)
    print(f"\n[scorecard] wrote {len(out_df)} rows -> {outfile}")

    # Print scorecard summary table
    print("\n" + "=" * 76)
    print(f"  SCORECARD   variant = {variant}   n_pairs_attempted = {len(pairs)}")
    print("=" * 76)
    hdr = (f"  {'config':<18} {'n':>6} {'BA_bias%':>10} "
           f"{'BA_RMSE':>10} {'MCuFt_mean':>12}")
    print(hdr)
    print("-" * 76)
    for cfg_name in configs:
        grp = out_df.loc[out_df["config"] == cfg_name].copy()
        ok  = (
            grp["BA_PRED"].notna()
            & grp["BA_OBS"].notna()
            & (grp["BA_OBS"] > 0)
            & (grp["BA_PRED"] > 0)
        )
        n = int(ok.sum())
        if n == 0:
            print(f"  {cfg_name:<18} {'0':>6}  (no usable predictions)")
            continue
        bias_pct = 100.0 * (
            grp.loc[ok, "BA_PRED"].mean() / grp.loc[ok, "BA_OBS"].mean() - 1
        )
        rmse = float(
            np.sqrt(((grp.loc[ok, "BA_PRED"] - grp.loc[ok, "BA_OBS"]) ** 2).mean())
        )
        mcuft_ok = grp["MCuFt_PRED"].notna() & (grp["MCuFt_PRED"] > 0)
        mcuft    = grp.loc[mcuft_ok, "MCuFt_PRED"].mean() if mcuft_ok.any() else float("nan")
        print(
            f"  {cfg_name:<18} {n:>6} {bias_pct:>+10.2f} {rmse:>10.2f} {mcuft:>12.1f}"
        )
    print("=" * 76)


if __name__ == "__main__":
    main()
