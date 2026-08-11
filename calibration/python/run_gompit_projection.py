#!/usr/bin/env python3
"""Gompit-mortality projection runner (Greg Johnson CONUS survival in the loop).

This wires Greg Johnson's gompit survival -- annual hazard
    H = exp(b0 + b1*(cr+0.01)^b2 + b3*cch^b4),  P(survive T) = exp(-H*T)
-- into the perseus FVS projection path so a stand can be projected with gompit
mortality substituted for FVS's native mortality, then compared against the
default and the (over-reducing) old-calibrated arms.

The mortality science is Aaron's: the gompit coefficients, the ORGANON crown
closure proxy (cch_organon), and the validated affine map onto the gompit cch
scale (project_mortality.cycle_survival, CCH_A/CCH_B from 35d_validate_cch.R,
Spearman 0.93). This script is the *wiring + validation* layer only.

Two gompit arms are implemented because they trade rigor against cost:

  iterative  -- cycle-by-cycle re-entry. Run FVS one cycle with calibrated
                growth and MORTMULT 0 (mortality off so no trees are dropped),
                read the grown treelist, apply gompit survival to TPA, rebuild
                the treeinit from the thinned/grown stand, and feed it into the
                next cycle. Growth therefore sees the gompit-reduced density each
                cycle. Correct, but costs num_cycles FVS calls per stand.

  posthoc    -- single-shot. Run FVS once for all cycles with calibrated growth
                and MORTMULT 0, then walk the returned per-cycle treelists and
                reduce each tree's TPA by the cumulative product of its gompit
                survivals (tracked by the stable FVS TreeId). One FVS call per
                stand; growth does NOT see the thinning (density slightly high).

Validation goal: confirm posthoc ~= iterative on a sample so the cheap posthoc
arm can be used for the full CONUS campaign, and confirm the gompit AGB sits
between default and the over-reduced old-calibrated arm.

FVS treelist note: the live-tree table carries PctCr (percent), not a 0-1 crown
ratio, and SpeciesFIA / DBH / Ht / TPA / TreeId. We map CrRatio = PctCr / 100.

Usage:
  python run_gompit_projection.py --variant ne --n 20 \
      --coeff conus_mort/full_out/greg_mortality_coefficients.csv \
      --standinit-dir .../standinit_by_variant \
      --treeinit-dir  .../FIA_fresh/treeinit \
      --output out_gompit/ne_validation.csv [--arms default,calibrated,iterative,posthoc]
"""
from __future__ import annotations

import argparse
import copy
import json
import logging
import os
import sys

import numpy as np
import pandas as pd

PROJECT_ROOT = os.environ.get("FVS_PROJECT_ROOT", os.path.expanduser("~/fvs-modern"))
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.join(PROJECT_ROOT, "calibration", "python"))
import perseus_100yr_projection as P  # noqa: E402

# gompit glue (Aaron's validated pieces)
from greg_mortality import GregMortality  # noqa: E402
from project_mortality import cycle_survival  # noqa: E402

logging.basicConfig(level=logging.INFO,
                    format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("gompit")

def set_regeneration(enabled: bool) -> None:
    """Toggle FVS automatic establishment in the shared keyfile template.

    The default eastern keyfiles fire the establishment model every cycle, so
    ingrowth keeps adding stems. That is fine for a regen-on reference run, but
    it confounds a *mortality-model* comparison: with FVS mortality disabled
    (gompit arm), unchecked regeneration makes biomass run away. For a clean
    closed-cohort comparison (same trees, three mortality models) we inject
    NOAUTOES so no new cohorts are established in any arm."""
    base = getattr(P, "_KEYFILE_TEMPLATE_ORIG", None)
    if base is None:
        P._KEYFILE_TEMPLATE_ORIG = P.KEYFILE_TEMPLATE
        base = P.KEYFILE_TEMPLATE
    if enabled:
        P.KEYFILE_TEMPLATE = base
    else:
        # NOAUTOES = turn off automatic establishment (closed cohort)
        P.KEYFILE_TEMPLATE = base.replace("ECHOSUM\n", "ECHOSUM\nNOAUTOES\n", 1)


FIPS = {1: "AL", 2: "AK", 4: "AZ", 5: "AR", 6: "CA", 8: "CO", 9: "CT", 10: "DE",
        12: "FL", 13: "GA", 16: "ID", 17: "IL", 18: "IN", 19: "IA", 20: "KS",
        21: "KY", 22: "LA", 23: "ME", 24: "MD", 25: "MA", 26: "MI", 27: "MN",
        28: "MS", 29: "MO", 30: "MT", 31: "NE", 32: "NV", 33: "NH", 34: "NJ",
        35: "NM", 36: "NY", 37: "NC", 38: "ND", 39: "OH", 40: "OK", 41: "OR",
        42: "PA", 44: "RI", 45: "SC", 46: "SD", 47: "TN", 48: "TX", 49: "UT",
        50: "VT", 51: "VA", 53: "WA", 54: "WV", 55: "WI", 56: "WY"}


# ---------------------------------------------------------------------------
# gompit config (calibrated growth, mortality disabled so FVS keeps every tree)
# ---------------------------------------------------------------------------
def build_gompit_config_dir(base_config_dir: str, variant: str, dest_dir: str) -> str:
    """Write a calibrated config with mort_multiplier zeroed (MORTMULT 0).

    Keeps the calibrated DG/HI/HD/CR multipliers (growth calibration) but sets
    mortality to zero so FVS does not remove trees -- gompit owns mortality.
    Returns the config dir to hand to perseus (P.CONFIG_DIR)."""
    src = os.path.join(base_config_dir, "calibrated", f"{variant.lower()}.json")
    with open(src) as fh:
        cfg = json.load(fh)
    maxsp = int(cfg.get("maxsp"))
    cm = cfg.setdefault("calibration_multipliers", {})
    cm["mort_multiplier"] = [0.0] * maxsp  # MORTMULT 0 -> mortality off
    prov = cm.get("provenance", {})
    if isinstance(prov, dict):
        prov["mort_mapping"] = "zeroed (gompit arm; FVS mortality disabled)"
    os.makedirs(os.path.join(dest_dir, "calibrated"), exist_ok=True)
    out = os.path.join(dest_dir, "calibrated", f"{variant.lower()}.json")
    with open(out, "w") as fh:
        json.dump(cfg, fh)
    return dest_dir


# ---------------------------------------------------------------------------
# treelist (grown) -> fvs_treeinit table (for cycle re-entry)
# ---------------------------------------------------------------------------
def treelist_to_treeinit(tl: pd.DataFrame, stand_id: str,
                         tpa_override: dict | None = None) -> pd.DataFrame:
    """Convert a grown FVS treelist back into the fvs_treeinit input schema.

    tpa_override: optional {TreeId(str): tpa} to inject gompit-thinned expansion.
    Trees whose override TPA is <= 0 are dropped."""
    recs = []
    for _, r in tl.iterrows():
        tid = str(r["TreeId"])
        tpa = float(r["TPA"])
        if tpa_override is not None:
            tpa = float(tpa_override.get(tid, 0.0))
        if not np.isfinite(tpa) or tpa <= 0:
            continue
        dbh = float(r.get("DBH", 0) or 0)
        if not np.isfinite(dbh) or dbh < 1.0:
            continue
        ht = r.get("Ht", 0)
        pctcr = r.get("PctCr", 0)
        recs.append({
            "stand_id": stand_id,
            "plot_id": int(float(r.get("PtIndex", 1) or 1)),
            "tree_id": int(float(r.get("TreeIndex", len(recs) + 1))),
            "tree_count": tpa,
            "species": int(float(r.get("SpeciesFIA", 0) or 0)),
            "diameter": round(dbh, 1),
            "ht": round(float(ht), 0) if pd.notna(ht) and float(ht) > 0 else 0,
            "crratio": int(float(pctcr)) if pd.notna(pctcr) and float(pctcr) > 0 else 0,
        })
    return pd.DataFrame(recs)


def _treelist_survival(tl: pd.DataFrame, greg: GregMortality, years: float) -> dict:
    """Per-TreeId gompit survival for one cycle. CrRatio = PctCr/100.

    Returns {TreeId(str): survival in (0,1]}. Species without a fitted gompit
    get survival = 1.0 (no native fallback because FVS mortality is off)."""
    t = tl.copy()
    t["CrRatio"] = pd.to_numeric(t["PctCr"], errors="coerce").fillna(0) / 100.0
    surv = cycle_survival(t, greg, years=years,
                          spcd_col="SpeciesFIA", dbh_col="DBH", ht_col="Ht",
                          cr_col="CrRatio", tpa_col="TPA")
    out = {}
    for (_, r), s in zip(t.iterrows(), surv):
        out[str(r["TreeId"])] = 1.0 if s is None else float(s)
    return out


# ---------------------------------------------------------------------------
# the four arms
# ---------------------------------------------------------------------------
def arm_native(sdf, tdf, sid, variant, nsbe, inv_year, config_version,
               num_cycles, cycle_length):
    """default (None) or calibrated -- single FVS call, FVS owns mortality."""
    fr = P.run_fvs_projection(sdf, tdf, sid, variant, config_version=config_version,
                              num_cycles=num_cycles, cycle_length=cycle_length)
    out = []
    for cy, tl in sorted(fr["treelists"].items()):
        if cy - inv_year < 0:
            continue
        out.append((cy, cy - inv_year, P.compute_plot_agb(tl, nsbe)))
    return out


def arm_gompit_posthoc(sdf, tdf, sid, variant, nsbe, inv_year, gompit_dir,
                       greg, num_cycles, cycle_length):
    """Single FVS call (calibrated growth, MORTMULT 0); cumulative gompit TPA
    override tracked by stable TreeId. Growth does NOT see the thinning."""
    saved = P.CONFIG_DIR
    P.CONFIG_DIR = gompit_dir
    try:
        fr = P.run_fvs_projection(sdf, tdf, sid, variant, config_version="calibrated",
                                  num_cycles=num_cycles, cycle_length=cycle_length)
    finally:
        P.CONFIG_DIR = saved
    years = sorted(fr["treelists"].keys())
    out = []
    surv_tpa: dict[str, float] = {}
    for i, cy in enumerate(years):
        tl = fr["treelists"][cy].copy()
        # initialise / carry surviving expansion by TreeId (FVS TPA == initial,
        # constant, because mortality is off)
        for _, r in tl.iterrows():
            tid = str(r["TreeId"])
            if tid not in surv_tpa:
                surv_tpa[tid] = float(r["TPA"])
        tl["TPA"] = [surv_tpa.get(str(r["TreeId"]), 0.0) for _, r in tl.iterrows()]
        if cy - inv_year >= 0:
            out.append((cy, cy - inv_year, P.compute_plot_agb(tl, nsbe)))
        # apply this cycle's gompit survival to carry TPA into the next cycle
        if i < len(years) - 1:
            yr_len = years[i + 1] - cy
            s = _treelist_survival(fr["treelists"][cy], greg, yr_len)
            for tid in list(surv_tpa):
                surv_tpa[tid] *= s.get(tid, 1.0)
    return out


def arm_gompit_iterative(plot_data, st_rows, sid, variant, nsbe, inv_year,
                         gompit_dir, greg, num_cycles, cycle_length,
                         build_treeinit):
    """Cycle-by-cycle re-entry: one FVS cycle at a time, gompit thinning fed back
    so growth sees the reduced density. Correct but num_cycles FVS calls."""
    saved = P.CONFIG_DIR
    P.CONFIG_DIR = gompit_dir
    out = []
    try:
        cur_tdf = build_treeinit(st_rows, sid)
        if cur_tdf.empty:
            return out
        cur_year = inv_year
        elapsed = 0
        for c in range(num_cycles):
            # rebuild the standinit each re-entry so FVS projects from cur_year,
            # not the original inventory year (otherwise every cycle re-projects
            # the same inv_year..inv_year+cycle window).
            pdat = dict(plot_data)
            pdat["INVYR"] = cur_year
            sdf = P.build_fvs_standinit(pdat, sid, variant)
            fr = P.run_fvs_projection(sdf, cur_tdf, sid, variant,
                                      config_version="calibrated",
                                      num_cycles=1, cycle_length=cycle_length)
            yrs = sorted(fr["treelists"].keys())
            if not yrs:
                break
            y0 = yrs[0]
            if c == 0:  # record the starting state once (proj_year 0)
                out.append((inv_year, 0, P.compute_plot_agb(fr["treelists"][y0], nsbe)))
            if len(yrs) < 2:
                break  # no projected end state -> stand stopped
            y_end = yrs[-1]
            grown = fr["treelists"][y_end]
            # gompit survival over this cycle, from the START crown state
            s = _treelist_survival(fr["treelists"][y0], greg, y_end - y0)
            thinned_tpa = {str(r["TreeId"]): float(r["TPA"]) * s.get(str(r["TreeId"]), 1.0)
                           for _, r in grown.iterrows()}
            tl = grown.copy()
            tl["TPA"] = [thinned_tpa.get(str(r["TreeId"]), 0.0) for _, r in tl.iterrows()]
            elapsed += (y_end - y0)
            out.append((inv_year + elapsed, elapsed, P.compute_plot_agb(tl, nsbe)))
            # feed the grown + gompit-thinned stand into the next cycle
            cur_tdf = treelist_to_treeinit(grown, sid, tpa_override=thinned_tpa)
            cur_year = y_end
            if cur_tdf.empty:
                break
    finally:
        P.CONFIG_DIR = saved
    return out


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--variant", required=True)
    ap.add_argument("--n", type=int, default=20, help="stands to sample")
    ap.add_argument("--coeff", required=True, help="greg_mortality_coefficients.csv")
    ap.add_argument("--standinit-dir", required=True)
    ap.add_argument("--treeinit-dir", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--config-dir", default=P.CONFIG_DIR)
    ap.add_argument("--gompit-config-dir", default="/tmp/gompit_config")
    ap.add_argument("--num-cycles", type=int, default=20)
    ap.add_argument("--cycle-length", type=int, default=5)
    ap.add_argument("--arms", default="default,calibrated,posthoc,iterative")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--regen", choices=["on", "off"], default="off",
                    help="automatic establishment; 'off' = closed-cohort "
                         "comparison (recommended for mortality-model isolation)")
    a = ap.parse_args()

    arms = [x.strip() for x in a.arms.split(",") if x.strip()]
    variant = a.variant.lower()
    set_regeneration(a.regen == "on")
    log.info(f"regeneration: {a.regen} (closed cohort)" if a.regen == "off"
             else "regeneration: on")

    # the treeinit builder lives in the canonical CONUS runner; import it
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    try:
        import run_conus_task_fvstreeinit as RC  # noqa: E402
        build_treeinit = RC.treeinit_for_stand
    except Exception:
        # fallback: identical near-passthrough builder
        def build_treeinit(fvs_rows, stand_id):
            recs = []
            for i, t in enumerate(fvs_rows.itertuples(index=False)):
                d = getattr(t, "DIAMETER", np.nan)
                if pd.isna(d) or float(d) < 1.0:
                    continue
                ht = getattr(t, "HT", 0)
                cr = getattr(t, "CRRATIO", 0)
                recs.append({"stand_id": stand_id,
                             "plot_id": int(float(getattr(t, "PLOT_ID", 1) or 1)),
                             "tree_id": i + 1,
                             "tree_count": float(getattr(t, "TREE_COUNT", 1.0) or 1.0),
                             "species": int(float(getattr(t, "SPECIES", 0) or 0)),
                             "diameter": round(float(d), 1),
                             "ht": round(float(ht), 0) if pd.notna(ht) and ht not in ("", None) and float(ht) > 0 else 0,
                             "crratio": int(float(cr)) if pd.notna(cr) and cr not in ("", None) and float(cr) > 0 else 0})
            return pd.DataFrame(recs)

    greg = GregMortality(a.coeff)
    gompit_dir = build_gompit_config_dir(a.config_dir, variant, a.gompit_config_dir)
    log.info(f"gompit config -> {gompit_dir} (mort zeroed); fitted species: {len(greg.coef)}")

    nsbe = P.NSBECalculator(P.NSBE_ROOT)

    si = pd.read_csv(os.path.join(a.standinit_dir, f"standinit_{variant.upper()}.csv"),
                     low_memory=False)
    si["STAND_CN"] = si["STAND_CN"].apply(lambda x: str(int(float(x))) if pd.notna(x) else "")

    # collect stands that have matching trees, up to --n
    selected = []
    cache = {}
    for state_fips, grp in si.groupby("STATE"):
        if len(selected) >= a.n:
            break
        try:
            state = FIPS[int(float(state_fips))]
        except (KeyError, ValueError, TypeError):
            continue
        tfile = os.path.join(a.treeinit_dir, f"{state}_FVS_TREEINIT_PLOT.csv")
        if not os.path.exists(tfile):
            continue
        if state not in cache:
            tt = pd.read_csv(tfile, low_memory=False)
            tt["STAND_CN"] = tt["STAND_CN"].apply(lambda x: str(int(float(x))) if pd.notna(x) else "")
            cache[state] = {k: v for k, v in tt.groupby("STAND_CN")}
        by_cn = cache[state]
        for _, stand in grp.iterrows():
            if len(selected) >= a.n:
                break
            cn = stand["STAND_CN"]
            rows = by_cn.get(cn)
            if rows is None or rows.empty:
                continue
            selected.append((state, cn, stand, rows))

    log.info(f"{variant}: {len(selected)} stands with trees selected")
    out_rows = []
    for state, cn, stand, rows in selected:
        sid = f"S{cn}"
        inv_year = int(float(stand.get("INV_YEAR") or 2010))
        plot_data = {"INVYR": inv_year, "LAT": stand.get("LATITUDE"),
                     "LON": stand.get("LONGITUDE"),
                     "ELEV": stand.get("ELEVFT") or 500, "SLOPE": stand.get("SLOPE") or 10,
                     "ASPECT": stand.get("ASPECT") or 180, "STDAGE": stand.get("AGE") or 50}
        try:
            sdf = P.build_fvs_standinit(plot_data, sid, variant)
            tdf = build_treeinit(rows, sid)
        except Exception as e:
            log.warning(f"{cn}: build failed: {e}")
            continue
        if tdf.empty:
            continue

        results = {}
        try:
            if "default" in arms:
                results["default"] = arm_native(sdf, tdf, sid, variant, nsbe, inv_year,
                                                 None, a.num_cycles, a.cycle_length)
            if "calibrated" in arms:
                results["calibrated"] = arm_native(sdf, tdf, sid, variant, nsbe, inv_year,
                                                    "calibrated", a.num_cycles, a.cycle_length)
            if "posthoc" in arms:
                results["gompit_posthoc"] = arm_gompit_posthoc(
                    sdf, tdf, sid, variant, nsbe, inv_year, gompit_dir, greg,
                    a.num_cycles, a.cycle_length)
            if "iterative" in arms:
                results["gompit_iterative"] = arm_gompit_iterative(
                    plot_data, rows, sid, variant, nsbe, inv_year, gompit_dir, greg,
                    a.num_cycles, a.cycle_length, build_treeinit)
        except Exception as e:
            log.warning(f"{cn}: projection failed: {e}")
            continue

        for arm, series in results.items():
            for cy, py, agb in series:
                out_rows.append({"STAND_CN": cn, "STATE": state, "VARIANT": variant.upper(),
                                 "ARM": arm, "YEAR": cy, "PROJ_YEAR": py,
                                 "AGB_TONS_AC": round(float(agb), 4)})
        log.info(f"  {cn}: arms={list(results.keys())}")

    os.makedirs(os.path.dirname(os.path.abspath(a.output)), exist_ok=True)
    df = pd.DataFrame(out_rows)
    df.to_csv(a.output, index=False)
    log.info(f"wrote {len(df)} rows -> {a.output}")

    # quick validation summary at the terminal year
    if not df.empty:
        last = df["PROJ_YEAR"].max()
        piv = (df[df["PROJ_YEAR"] == last]
               .groupby("ARM")["AGB_TONS_AC"].mean().round(2))
        log.info(f"mean AGB at proj_year {last} by arm:\n{piv.to_string()}")


if __name__ == "__main__":
    main()
