#!/usr/bin/env python3
"""
run_ba_decomp.py
================
BA decomposition diagnostic: isolates Greg DG, Greg HG+CRW, and GOMPIT
mortality contributions to the systematic -27% BA underprediction seen in
organon / conus_spdep / conus_climate across NE and SN variants.

Six configs (3 x diagnostic + 2 x bracketing + 1 x interaction):
  1. fvs_base      - native FVS, no Greg hooks (upper bound / control)
  2. organon       - full Greg system: DG + HG_ORGANON + CRW + GOMPIT (-27% case)
  3. diag_dg_only  - Greg DG env only; native HG, native varmort, no CRW
  4. diag_hg_only  - Greg HG_ORGANON + CRW env only; native DG, native varmort
  5. diag_mort_only- GOMPIT mortality env only; native DG, native HG
  6. diag_dg_mort  - Greg DG + GOMPIT; no HG/CRW (DG x mortality interaction)

Diagnostic logic:
  If diag_dg_only  BA_bias ~ -27%  -> Greg DG suppresses individual tree DG
  If diag_mort_only BA_bias ~ -27% -> GOMPIT is over-killing trees (TPA driven)
  If diag_hg_only  BA_bias ~   0%  -> HG/CRW has negligible BA effect (expected)

Reuses all data-loading, build, and FVS-run functions from run_5model_scorecard.py
via importlib (avoids code duplication and keeps DG-swap logic in one place).

Usage:
  python3 run_ba_decomp.py --variant ne [--max-pairs 200] [--seed 42]

Output:
  calibration/output/fia_scorecard/ba_decomp_{variant}.csv
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import sys

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Load parent module dynamically (avoids sys.argv side-effects from its main)
# ---------------------------------------------------------------------------
_PARENT = "/users/PUOM0008/crsfaaron/fvs-modern/calibration/run_5model_scorecard.py"
_spec = importlib.util.spec_from_file_location("sc", _PARENT)
sc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sc)

# Pull constants and functions we need
_CONF_DG    = sc._CONF_DG
_DG_ENV     = sc._DG_ENV
_HG_BASE    = sc._HG_BASE
_CRW_ENV    = sc._CRW_ENV
_MORT_ENV   = sc._MORT_ENV
_DG_KEYS    = sc._DG_KEYS
VARIANTS    = sc.VARIANTS
STATECD_TO_ABBREV = sc.STATECD_TO_ABBREV

load_pairs             = sc.load_pairs
load_trees_for_cns     = sc.load_trees_for_cns
load_cond_for_cns      = sc.load_cond_for_cns
obs_ba_ft2ac           = sc.obs_ba_ft2ac
obs_vol_mcuft_bdft     = sc.obs_vol_mcuft_bdft
build_standinit        = sc.build_standinit
build_treeinit         = sc.build_treeinit
run_stand              = sc.run_stand

# ---------------------------------------------------------------------------
# Six diagnostic model configs
# ---------------------------------------------------------------------------
_HG_ORGANON_COEF = f"{_CONF_DG}/greg_hg_coefficients.csv"
_NO_CALIB        = {"dg": False, "hg": False, "mort": False}  # SDIMAX+BAMAX only
_NE_DG_STATES   = ["CT", "MA", "ME", "NH", "NY", "RI", "VT"]

DECOMP_MODELS: dict[str, dict] = {
    # Control: pure native FVS
    "fvs_base": {
        "extra_env":          {},
        "calib_components":   None,
        "states_restrict_dg": None,
    },
    # Full Greg system: the -27% case being diagnosed
    "organon": {
        "extra_env": {
            **_DG_ENV,
            **_HG_BASE,
            "FVS_GREGHG_COEF": _HG_ORGANON_COEF,
            **_CRW_ENV,
            **_MORT_ENV,
        },
        "calib_components":   _NO_CALIB,
        "states_restrict_dg": _NE_DG_STATES,
    },
    # Greg DG env only; everything else native
    "diag_dg_only": {
        "extra_env":          {**_DG_ENV},
        "calib_components":   _NO_CALIB,
        "states_restrict_dg": _NE_DG_STATES,
    },
    # Greg HG ORGANON + CRW only; native DG and varmort
    "diag_hg_only": {
        "extra_env": {
            **_HG_BASE,
            "FVS_GREGHG_COEF": _HG_ORGANON_COEF,
            **_CRW_ENV,
        },
        "calib_components":   _NO_CALIB,
        "states_restrict_dg": None,  # no DG to restrict
    },
    # GOMPIT mortality only; native DG and HG
    "diag_mort_only": {
        "extra_env":          {**_MORT_ENV},
        "calib_components":   _NO_CALIB,
        "states_restrict_dg": None,
    },
    # Greg DG + GOMPIT interaction (no HG/CRW)
    "diag_dg_mort": {
        "extra_env":          {**_DG_ENV, **_MORT_ENV},
        "calib_components":   _NO_CALIB,
        "states_restrict_dg": _NE_DG_STATES,
    },
}

DECOMP_ORDER = [
    "fvs_base", "organon",
    "diag_dg_only", "diag_hg_only", "diag_mort_only", "diag_dg_mort",
]


def _resolve_env_calib(model_name: str, state_abbrev: str | None):
    """Apply per-plot DG swap (same logic as run_5model_scorecard main loop)."""
    mdl = DECOMP_MODELS[model_name]
    eff_extra_env = dict(mdl["extra_env"])
    eff_calib     = mdl["calib_components"]

    restrict_dg = mdl.get("states_restrict_dg")
    if restrict_dg and state_abbrev not in restrict_dg:
        eff_extra_env = {k: v for k, v in eff_extra_env.items() if k not in _DG_KEYS}
        base_cc       = eff_calib or {}
        eff_calib     = {**base_cc, "dg": True}

    return eff_extra_env, eff_calib


def main() -> None:
    ap = argparse.ArgumentParser(
        description="BA decomposition diagnostic for Greg DG / HG / GOMPIT"
    )
    ap.add_argument("--variant",   required=True, choices=list(VARIANTS))
    ap.add_argument("--max-pairs", type=int, default=200)
    ap.add_argument("--seed",      type=int, default=42)
    ap.add_argument(
        "--outdir",
        default="/users/PUOM0008/crsfaaron/fvs-modern/calibration/output/fia_scorecard",
    )
    args = ap.parse_args()

    variant = args.variant
    states  = VARIANTS[variant]["states"]
    print(f"[ba_decomp] variant={variant}  states={states}")
    print(f"[ba_decomp] models={DECOMP_ORDER}")
    sys.stdout.flush()

    pairs = load_pairs(variant, max_pairs=args.max_pairs, seed=args.seed)
    if len(pairs) == 0:
        print(f"[ba_decomp] ERROR: no remeasurement pairs found for {variant}")
        sys.exit(1)
    print(f"[ba_decomp] {len(pairs)} remeasurement pairs (5-15 yr interval)")
    sys.stdout.flush()

    all_cns = set(pairs["PREV_PLT_CN"].tolist()) | set(pairs["CN"].tolist())
    trees   = load_trees_for_cns(states, all_cns)
    conds   = load_cond_for_cns(states, all_cns)
    print(f"[ba_decomp] {len(trees)} tree recs, {len(conds)} cond recs, "
          f"{len(all_cns)} plot CNs")
    sys.stdout.flush()

    rows: list[dict] = []
    n_total = len(pairs) * len(DECOMP_ORDER)
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
            n_done += len(DECOMP_ORDER)
            continue

        cond_row = None
        c1 = conds.loc[conds["PLT_CN"] == t1_cn]
        if len(c1) > 0:
            cond_row = c1.iloc[0].to_dict()

        ba_obs              = obs_ba_ft2ac(trees, t2_cn)
        mcuft_obs, bdft_obs = obs_vol_mcuft_bdft(trees, t2_cn)

        sid      = f"{variant.upper()}_{t1_cn % 10_000_000:07d}"
        stand_df = build_standinit(sid, measyear1, variant, cond_row)
        tree_df  = build_treeinit(t1, sid)
        if len(tree_df) == 0:
            n_done += len(DECOMP_ORDER)
            continue

        state_abbrev = STATECD_TO_ABBREV.get(state_cd)

        for model_name in DECOMP_ORDER:
            eff_extra_env, eff_calib = _resolve_env_calib(model_name, state_abbrev)

            # run_stand() looks up binary by variant; model_name "fvs_base" is a dummy
            # so we can use extra_env_override/calib_override to control the run fully.
            result_df, err = run_stand(
                sid, stand_df, tree_df, variant, "fvs_base", ncyc,
                extra_env_override=eff_extra_env,
                calib_override=eff_calib,
            )
            n_done += 1

            if err or result_df is None or len(result_df) == 0:
                if n_done % 50 == 0:
                    print(f"  [{n_done}/{n_total}] {sid} {model_name}: ERR {err}")
                    sys.stdout.flush()
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
                "BA_PRED":    float(last.get("BA",     float("nan"))),
                "TPA_PRED":   float(last.get("Tpa",    float("nan"))),
                "QMD_PRED":   float(last.get("QMD",    float("nan"))),
                "MCuFt_PRED": float(last.get("MCuFt",  float("nan"))),
                "BdFt_PRED":  float(last.get("BdFt",   float("nan"))),
                "BA_OBS":     ba_obs,
                "MCuFt_OBS":  mcuft_obs,
                "BdFt_OBS":   bdft_obs,
            })

            if n_done % 100 == 0 or n_done == n_total:
                ba_p = rows[-1]["BA_PRED"] if rows else float("nan")
                print(f"  [{n_done}/{n_total}] {sid} {model_name}: "
                      f"BA_pred={ba_p:.1f}  BA_obs={ba_obs:.1f}")
                sys.stdout.flush()

    out_df = pd.DataFrame(rows)
    os.makedirs(args.outdir, exist_ok=True)
    outfile = os.path.join(args.outdir, f"ba_decomp_{variant}.csv")
    out_df.to_csv(outfile, index=False)
    print(f"\n[ba_decomp] wrote {len(out_df)} rows -> {outfile}")

    # Summary table
    print("\n" + "=" * 95)
    print(f"  BA DECOMP SCORECARD   variant={variant}   n_pairs_attempted={len(pairs)}")
    print("=" * 95)
    hdr = (f"  {'model':<18} {'n':>5} {'BA_bias%':>10} {'BA_RMSE':>9} "
           f"{'TPA_mean':>10} {'QMD_mean':>9}  interpretation")
    print(hdr)
    print("-" * 95)
    interp = {
        "fvs_base":       "control (native FVS varmort)",
        "organon":        "full Greg system [-27% target]",
        "diag_dg_only":   "Greg DG only  -> ~-27% = DG culprit",
        "diag_hg_only":   "Greg HG+CRW only -> ~0% expected",
        "diag_mort_only":  "GOMPIT only   -> ~-27% = mort culprit",
        "diag_dg_mort":   "Greg DG+GOMPIT -> additive or synergistic?",
    }
    for model_name in DECOMP_ORDER:
        grp = out_df.loc[out_df["model"] == model_name]
        ok  = (grp["BA_PRED"].notna() & grp["BA_OBS"].notna()
               & (grp["BA_OBS"] > 0) & (grp["BA_PRED"] > 0))
        n   = int(ok.sum())
        if n == 0:
            print(f"  {model_name:<18} {'0':>5}  (no usable predictions)")
            continue
        ba_bias_pct = 100.0 * (grp.loc[ok,"BA_PRED"].mean() / grp.loc[ok,"BA_OBS"].mean() - 1)
        ba_rmse     = float(np.sqrt(((grp.loc[ok,"BA_PRED"] - grp.loc[ok,"BA_OBS"])**2).mean()))
        tpa_mean    = grp.loc[ok,"TPA_PRED"].mean()
        qmd_mean    = grp.loc[ok,"QMD_PRED"].mean()
        note        = interp.get(model_name, "")
        print(f"  {model_name:<18} {n:>5} {ba_bias_pct:>+10.2f} {ba_rmse:>9.2f} "
              f"{tpa_mean:>10.1f} {qmd_mean:>9.2f}  {note}")
    print("=" * 95)


if __name__ == "__main__":
    main()
