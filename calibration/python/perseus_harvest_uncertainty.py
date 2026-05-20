#!/usr/bin/env python3
from __future__ import annotations

"""
PERSEUS 100-Year FVS Projection: Harvest x Climate x Posterior Uncertainty
==========================================================================

Combines the three existing PERSEUS drivers into one scenario-aware
uncertainty driver.

  * From perseus_uncertainty_projection.py
      Posterior parameter draws per variant (config/calibrated/{v}_draws.json)
      Default and calibrated MAP point estimates
      UncertaintyEngine integration (get_draw, generate_keywords_for_draw)

  * From perseus_harvest_projection.py
      50 percent BA harvest via ThinBBA on scheduled cycles (4/27 absolute
      residual patch)
      HARVEST_KEYFILE_TEMPLATE with calibration and harvest plug-in slots

  * From perseus_climate_projection.py
      Site index multiplier by plot x RCP applied to stand_init via
      modify_standinit_si_delta

Output: per-cycle AGB rows (PLOT, YEAR, PROJ_YEAR, VARIANT, CONFIG, DRAW,
AGB_TONS_AC) tagged with SCENARIO, RCP, HARVEST. Schema is compatible with
perseus_uncertainty_aggregate.py and aggregate_with_uncertainty when grouped
on SCENARIO.

Per scenario per plot:
  2 variants x (1 default + 1 calibrated MAP + N posterior draws) FVS runs

Usage (single scenario, single batch):
  python perseus_harvest_uncertainty.py \
      --batch-id 1 --batch-size 100 --n-draws 500 \
      --rcp 4.5 --harvest \
      --variants ne acd \
      --perseus-csv ... --si-ratio-csv ... \
      --output-dir <base>/perseus

Outputs (per scenario subdir):
  <output-dir>/<scenario>/perseus_<scenario>_unc_point_batch<N>.csv
  <output-dir>/<scenario>/perseus_<scenario>_unc_draws_batch<N>.csv
  <output-dir>/<scenario>/perseus_<scenario>_unc_meta_batch<N>.json

Author: A. Weiskittel (and Claude, picking up the 2026-04-24 DRAFT)
Date: 2026-05-06
"""

import argparse
import json
import logging
import os
import sqlite3
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Environment and path setup (mirror perseus_uncertainty_projection.py)
# ---------------------------------------------------------------------------

PROJECT_ROOT = os.environ.get(
    "FVS_PROJECT_ROOT",
    str(Path(__file__).resolve().parents[2])
)
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.join(PROJECT_ROOT, "deployment", "fvs2py"))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "deployment", "microfvs"))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "calibration", "python"))

FIA_DIR = os.environ.get("FVS_FIA_DATA_DIR", os.path.expanduser("~/fia_data"))
FVS_LIB_DIR = os.environ.get("FVS_LIB_DIR", os.path.join(PROJECT_ROOT, "lib"))
NSBE_ROOT = os.environ.get("NSBE_ROOT", os.path.join(PROJECT_ROOT, "data", "NSBE"))
CONFIG_DIR = os.environ.get("FVS_CONFIG_DIR", os.path.join(PROJECT_ROOT, "config"))
OUTPUT_DIR = os.environ.get(
    "PERSEUS_OUTPUT_DIR",
    os.path.join(PROJECT_ROOT, "calibration", "output", "perseus")
)
PERSEUS_CSV = os.environ.get(
    "PERSEUS_CSV",
    os.path.join(PROJECT_ROOT, "calibration", "data", "perseus_plots.csv")
)
SI_RATIO_CSV_DEFAULT = os.path.expanduser(
    "~/SiteIndex/future_SI/perseus_si_ratios.csv"
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Imports from the three existing drivers
# ---------------------------------------------------------------------------

from perseus_100yr_projection import (
    NSBECalculator,
    load_fia_trees_for_plots,
    build_fvs_standinit,
    build_fvs_treeinit,
    compute_plot_agb,
    compute_initial_agb,
    _run_via_subprocess,
    _run_via_fvs2py,
)

from perseus_harvest_projection import (
    HARVEST_KEYFILE_TEMPLATE,
    NUM_CYCLES,
    CYCLE_LENGTH,
    determine_harvest_eligibility,
    build_harvest_schedule,
    generate_harvest_keywords,
    load_fia_cond_for_plots,
)

from perseus_climate_projection import (
    load_si_ratios,
    get_plot_si_delta_ft,
    modify_standinit_si_delta,
)

# Match canonical pipeline constants
CARBON_FRACTION = 0.5

# ---------------------------------------------------------------------------
# FVS runner: harvest keyfile combined with calibration / draw keywords
# ---------------------------------------------------------------------------

def run_fvs_scenario(
    stand_init_df: pd.DataFrame,
    tree_init_df: pd.DataFrame,
    stand_id: str,
    variant: str,
    calibration_keywords: str,
    harvest_cycles: list[int],
    inv_year: int,
    current_ba: float | None = None,  # kept for forward compatibility
    num_cycles: int = NUM_CYCLES,
    cycle_length: int = CYCLE_LENGTH,
) -> dict:
    """Run FVS using HARVEST_KEYFILE_TEMPLATE with both keyword slots filled.

    The calibration_keywords slot holds either:
      - ""                        -> default parameters
      - calibrated MAP keywords   -> from FvsConfigLoader
      - draw keywords             -> from UncertaintyEngine.generate_keywords_for_draw

    The harvest slot is generated by perseus_harvest_projection.generate_harvest_keywords
    and degrades to a comment when harvest_cycles is empty.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "FVS_Data.db")
        conn = sqlite3.connect(db_path)
        stand_init_df.to_sql("fvs_standinit", conn, if_exists="replace", index=False)
        tree_init_df.to_sql("fvs_treeinit", conn, if_exists="replace", index=False)
        conn.close()

        # generate_harvest_keywords signature on Cardinal (post-rollback of
        # 4/27 absolute-residual patch) does not accept current_ba. Match the
        # live signature; the existing -0.50 retain-50pct behavior is what
        # the existing factorial in the workbook is already running with, so
        # this is the right thing for adding posterior bands on top.
        harvest_kw = generate_harvest_keywords(
            harvest_cycles, inv_year, cycle_length
        )

        keyfile_content = HARVEST_KEYFILE_TEMPLATE.format(
            stand_id=stand_id,
            db_path=db_path,
            calibration_keywords=(
                calibration_keywords if calibration_keywords
                else "** DEFAULT PARAMETERS"
            ),
            harvest_keywords=harvest_kw,
            num_cycles=num_cycles,
            cycle_length=cycle_length,
        )

        keyfile_path = os.path.join(tmpdir, f"{variant}_{stand_id}_hu.key")
        with open(keyfile_path, "w") as f:
            f.write(keyfile_content)

        exe_candidates = [
            os.path.join(FVS_LIB_DIR, f"FVS{variant.lower()}"),
            os.path.join(FVS_LIB_DIR, f"FVS{variant.upper()}"),
            f"/usr/local/bin/FVS{variant.lower()}",
        ]
        for exe_path in exe_candidates:
            if os.path.exists(exe_path) and os.access(exe_path, os.X_OK):
                return _run_via_subprocess(exe_path, keyfile_path,
                                           db_path, tmpdir)

        lib_path = os.path.join(FVS_LIB_DIR, f"FVS{variant.lower()}.so")
        if os.path.exists(lib_path):
            return _run_via_fvs2py(lib_path, keyfile_path, db_path,
                                   variant, None)

        raise FileNotFoundError(
            f"Cannot find FVS executable or library for {variant}"
        )


# ---------------------------------------------------------------------------
# Per plot processing (scenario aware, posterior draws aware)
# ---------------------------------------------------------------------------

def _scenario_tag(rcp: str, do_harvest: bool) -> str:
    """Match the canonical perseus_harvest_projection scenario directory names."""
    h = "harvest" if do_harvest else "noharvest"
    if rcp == "none":
        return f"{h}_noclimate"
    return f"{h}_rcp{rcp.replace('.', '')}"


def _initial_ba_from_trees(plot_trees: pd.DataFrame) -> float:
    """Approximate stand basal area (sq ft / ac) from FIA tree DBH and TPA."""
    if plot_trees.empty:
        return 0.0
    dbh = plot_trees.get("DIA", pd.Series(dtype=float))
    tpa = plot_trees.get("TPA_UNADJ", pd.Series(dtype=float))
    if dbh.empty or tpa.empty:
        return 0.0
    ba_per_tree = 0.005454154 * dbh ** 2  # sq ft per tree
    return float((ba_per_tree * tpa).sum())


def process_plot_scenario(
    plot_row,
    fia_trees: pd.DataFrame,
    nsbe: NSBECalculator,
    uncertainty_engines: dict,
    default_configs: dict,
    n_draws: int,
    variants: list[str],
    rcp: str,
    do_harvest: bool,
    si_ratios: pd.DataFrame | None,
    harvest_schedule: dict[int, set],
    seed: int = 42,
) -> tuple[list[dict], list[dict]]:
    """Process one plot under a single scenario (rcp x harvest) with draws.

    Mirrors perseus_uncertainty_projection.process_plot_with_uncertainty but
    splices in the SI delta and harvest schedule. Returns (point_records,
    draw_records) that include SCENARIO / RCP / HARVEST tags on every row.
    """
    plot_id = plot_row["PLOT"]
    raw_cn = plot_row["FIRST_PLTCN"]
    plt_cn = str(int(float(raw_cn))) if pd.notna(raw_cn) else ""
    stand_id = f"P{int(float(plot_id))}"

    if fia_trees.empty or "PLT_CN" not in fia_trees.columns:
        logger.warning(f"Plot {plot_id}: no tree data available")
        return [], []

    fia_cn_str = fia_trees["PLT_CN"].apply(
        lambda x: str(int(float(x))) if pd.notna(x) else ""
    )
    plot_trees = fia_trees.loc[fia_cn_str == plt_cn]
    if plot_trees.empty:
        logger.warning(f"Plot {plot_id}: no trees for PLT_CN={plt_cn}")
        return [], []

    plot_data = plot_trees.iloc[0].to_dict()
    plot_data["COUNTYCD"] = plot_row.get("COUNTYCD", 19)
    inv_year = int(plot_data.get("INVYR", plot_row.get("FIRST_INVYR", 2000)))
    countycd = int(plot_row.get("COUNTYCD", 0))

    # Climate: SI delta (ft) for this plot under this RCP
    si_delta_ft = 0.0
    if rcp != "none" and si_ratios is not None:
        try:
            si_delta_ft = get_plot_si_delta_ft(
                si_ratios, countycd, int(float(plot_id)), rcp
            )
        except Exception as e:
            logger.warning(f"SI delta lookup failed plot {plot_id}: {e}")

    # Harvest schedule: which cycles (1-based) does this plot get a ThinBBA event
    harvest_cycles: list[int] = []
    if do_harvest:
        for cycle, plot_set in harvest_schedule.items():
            if int(float(plot_id)) in plot_set:
                harvest_cycles.append(cycle)

    # Approximate current BA for the 4/27 ThinBBA absolute residual patch
    current_ba = _initial_ba_from_trees(plot_trees)

    initial_agb = compute_initial_agb(plot_trees, nsbe)
    scenario = _scenario_tag(rcp, do_harvest)
    point_records: list[dict] = []
    draw_records: list[dict] = []

    for variant in variants:
        stand_df = build_fvs_standinit(plot_data, stand_id, variant)
        tree_df = build_fvs_treeinit(plot_trees, stand_id)
        if tree_df.empty:
            logger.warning(f"Plot {plot_id}/{variant}: no valid trees")
            continue
        # Apply SI delta for the climate scenario (in place on standinit)
        if abs(si_delta_ft) > 0.01:
            stand_df = modify_standinit_si_delta(stand_df, si_delta_ft)

        # ------------------- Default and calibrated point -------------------
        for cfg_label, cfg_kw in [
            ("default", ""),
            ("calibrated", _calibrated_keywords(variant, default_configs)),
        ]:
            point_records.append({
                "PLOT": plot_id, "FIRST_PLTCN": plt_cn,
                "YEAR": inv_year, "PROJ_YEAR": 0,
                "VARIANT": variant.upper(), "CONFIG": cfg_label,
                "AGB_TONS_AC": round(initial_agb, 4),
                "SCENARIO": scenario, "RCP": rcp, "HARVEST": do_harvest,
                "SI_DELTA_FT": round(si_delta_ft, 2),
            })
            try:
                res = run_fvs_scenario(
                    stand_df, tree_df, stand_id, variant,
                    calibration_keywords=cfg_kw,
                    harvest_cycles=harvest_cycles, inv_year=inv_year,
                )
                for cycle_year, treelist in sorted(res["treelists"].items()):
                    proj_year = cycle_year - inv_year
                    agb = compute_plot_agb(treelist, nsbe)
                    point_records.append({
                        "PLOT": plot_id, "FIRST_PLTCN": plt_cn,
                        "YEAR": cycle_year, "PROJ_YEAR": proj_year,
                        "VARIANT": variant.upper(), "CONFIG": cfg_label,
                        "AGB_TONS_AC": round(agb, 4),
                        "SCENARIO": scenario, "RCP": rcp, "HARVEST": do_harvest,
                        "SI_DELTA_FT": round(si_delta_ft, 2),
                    })
            except Exception as e:
                logger.error(f"Plot {plot_id}/{variant}/{cfg_label}: {e}")

        # ------------------- Posterior draws (calibrated) -------------------
        engine = uncertainty_engines.get(variant.lower())
        default_cfg = default_configs.get(variant.lower())
        if engine is None or default_cfg is None:
            logger.warning(f"No uncertainty engine for {variant}; skipping draws")
            continue

        rng = np.random.default_rng(seed + int(float(plot_id)))
        n_avail = engine.n_draws
        n_use = min(n_draws, n_avail)
        draw_indices = rng.choice(n_avail, size=n_use, replace=False)

        for di, draw_idx in enumerate(draw_indices):
            draw = engine.get_draw(int(draw_idx))
            draw_keywords = engine.generate_keywords_for_draw(
                draw, default_cfg, draw_idx=int(draw_idx)
            )
            # Engine now emits BAIMULT (a real top-level FVS keyword, per
            # base/keywds.f90 TABLE position 58) instead of GROWMULT (which is
            # a CLIMATE-extension sub-keyword that triggers STOP 20 at the
            # outer level). The legacy GROWMULT filter is no longer needed.
            # Defensive: still strip any stray GROWMULT lines if present.
            draw_keywords = "\n".join(
                line for line in draw_keywords.splitlines()
                if not line.startswith("GROWMULT")
            )

            draw_records.append({
                "PLOT": plot_id, "FIRST_PLTCN": plt_cn,
                "YEAR": inv_year, "PROJ_YEAR": 0,
                "VARIANT": variant.upper(), "DRAW": int(draw_idx),
                "AGB_TONS_AC": round(initial_agb, 4),
                "SCENARIO": scenario, "RCP": rcp, "HARVEST": do_harvest,
                "SI_DELTA_FT": round(si_delta_ft, 2),
            })

            try:
                res = run_fvs_scenario(
                    stand_df, tree_df, stand_id, variant,
                    calibration_keywords=draw_keywords,
                    harvest_cycles=harvest_cycles, inv_year=inv_year,
                )
                for cycle_year, treelist in sorted(res["treelists"].items()):
                    proj_year = cycle_year - inv_year
                    agb = compute_plot_agb(treelist, nsbe)
                    draw_records.append({
                        "PLOT": plot_id, "FIRST_PLTCN": plt_cn,
                        "YEAR": cycle_year, "PROJ_YEAR": proj_year,
                        "VARIANT": variant.upper(), "DRAW": int(draw_idx),
                        "AGB_TONS_AC": round(agb, 4),
                        "SCENARIO": scenario, "RCP": rcp, "HARVEST": do_harvest,
                        "SI_DELTA_FT": round(si_delta_ft, 2),
                    })
            except Exception as e:
                logger.error(
                    f"Plot {plot_id}/{variant}/draw{int(draw_idx)}: {e}"
                )

    return point_records, draw_records


def _calibrated_keywords(variant: str, default_configs: dict) -> str:
    """Get calibrated MAP keywords for this variant.

    Loads via FvsConfigLoader.generate_keywords. Returns "" on failure to
    fall back to default parameters in the same scenario row (logged).
    """
    try:
        from config.config_loader import FvsConfigLoader
        loader = FvsConfigLoader(
            variant.lower(), version="calibrated", config_dir=CONFIG_DIR
        )
        return loader.generate_keywords(include_comments=False)
    except Exception as e:
        logger.warning(f"{variant} calibrated keyword load failed: {e}")
        return ""


# ---------------------------------------------------------------------------
# Engine and config priming (called once, reused across all plots)
# ---------------------------------------------------------------------------

def _build_engines(variants: list[str], seed: int):
    from config.uncertainty import UncertaintyEngine
    from config.config_loader import FvsConfigLoader
    engines = {}
    default_cfgs = {}
    for v in variants:
        try:
            eng = UncertaintyEngine(
                v.lower(), config_dir=CONFIG_DIR, seed=seed
            )
            if eng.draws_available:
                engines[v.lower()] = eng
                logger.info(
                    f"{v.upper()} engine: {eng.n_draws} draws available"
                )
        except Exception as e:
            logger.warning(f"{v} engine init: {e}")
        try:
            # Canonical pattern: store loader.config (the dict) so the engine's
            # generate_keywords_for_draw(draw, default_cfg, draw_idx) can call
            # default_cfg.get(...). Passing the FvsConfigLoader instance breaks
            # with: 'FvsConfigLoader' object has no attribute 'get'.
            loader = FvsConfigLoader(
                v.lower(), version="default", config_dir=CONFIG_DIR
            )
            default_cfgs[v.lower()] = loader.config
        except Exception as e:
            logger.warning(f"{v} default loader init: {e}")
    return engines, default_cfgs


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(
        description="PERSEUS 100yr FVS scenario factorial with posterior uncertainty"
    )
    p.add_argument("--batch-id", type=int, default=None)
    p.add_argument("--batch-size", type=int, default=100)
    p.add_argument("--all", action="store_true")
    p.add_argument("--perseus-csv", type=str, default=PERSEUS_CSV)
    p.add_argument("--si-ratio-csv", type=str, default=SI_RATIO_CSV_DEFAULT)
    p.add_argument("--variants", type=str, nargs="+", default=["acd", "ne"])
    p.add_argument("--n-draws", type=int, default=500)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument("--start-year", type=int, default=1999)
    p.add_argument("--end-year", type=int, default=2004)
    p.add_argument("--output-dir", type=str, default=OUTPUT_DIR)
    p.add_argument("--rcp", type=str, default="none",
                   choices=["none", "4.5", "8.5"], help="Climate scenario")
    p.add_argument("--harvest", action="store_true",
                   help="Apply 50pct BA harvest on scheduled cycles")
    p.add_argument("--exclude-plots", type=str, default="",
                   help="Comma list of plot IDs to skip, or path to a file")
    args = p.parse_args()

    # ------------------------- Load PERSEUS plot list -----------------------
    perseus = pd.read_csv(args.perseus_csv)
    mask = (
        (perseus["FIRST_INVYR"] >= args.start_year) &
        (perseus["FIRST_INVYR"] <= args.end_year)
    )
    perseus = perseus.loc[mask].copy()
    logger.info(f"PERSEUS rows after year filter: {len(perseus)}")

    if args.exclude_plots:
        raw = args.exclude_plots
        if os.path.exists(raw):
            with open(raw) as f:
                excludes = {line.strip() for line in f if line.strip()}
        else:
            excludes = {s.strip() for s in raw.split(",") if s.strip()}
        before = len(perseus)
        perseus = perseus.loc[
            ~perseus["PLOT"].astype(str).isin(excludes)
        ].copy()
        logger.info(f"Excluded {before - len(perseus)} plots: {sorted(excludes)}")

    # ------------------------- Build harvest schedule -----------------------
    # Schedule must be deterministic across batches for reproducibility, so
    # we build it from the FULL filtered perseus list (not the batch slice).
    harvest_schedule: dict[int, set] = {}
    if args.harvest:
        plt_cns_full = [
            str(int(float(x)))
            for x in perseus["FIRST_PLTCN"].dropna().unique()
        ]
        fia_cond = load_fia_cond_for_plots(plt_cns_full, FIA_DIR)
        eligibility = determine_harvest_eligibility(perseus, fia_cond)
        harvest_schedule = build_harvest_schedule(
            eligibility, NUM_CYCLES, seed=args.seed
        )
        n_events = sum(len(v) for v in harvest_schedule.values())
        logger.info(
            f"Harvest schedule: {n_events} plot-events across "
            f"{len(harvest_schedule)} cycles"
        )

    # ------------------------- Slice the batch ------------------------------
    if args.batch_id is not None:
        s = (args.batch_id - 1) * args.batch_size
        e = min(s + args.batch_size, len(perseus))
        perseus_batch = perseus.iloc[s:e].copy()
        logger.info(f"Batch {args.batch_id}: plots {s+1} to {e}")
    elif args.all:
        perseus_batch = perseus.copy()
    else:
        logger.error("Specify --batch-id or --all")
        sys.exit(1)

    if perseus_batch.empty:
        logger.info("No plots in batch. Exiting.")
        return

    # ------------------------- Prime engines and configs --------------------
    nsbe = NSBECalculator(NSBE_ROOT)
    engines, default_cfgs = _build_engines(args.variants, args.seed)

    # ------------------------- SI ratios ------------------------------------
    si_ratios = None
    if args.rcp != "none":
        if not os.path.exists(args.si_ratio_csv):
            logger.error(f"SI ratio file not found: {args.si_ratio_csv}")
            sys.exit(1)
        si_ratios = load_si_ratios(args.si_ratio_csv)
        logger.info(f"Loaded SI ratios for RCP {args.rcp}")

    # ------------------------- Load FIA tree data ---------------------------
    plt_cns_batch = [
        str(int(float(x)))
        for x in perseus_batch["FIRST_PLTCN"].dropna().unique()
    ]
    fia_trees = load_fia_trees_for_plots(plt_cns_batch, FIA_DIR)
    logger.info(
        f"Loaded {len(fia_trees)} tree rows for {len(plt_cns_batch)} PLT_CNs"
    )

    # ------------------------- Process plots --------------------------------
    all_point: list[dict] = []
    all_draws: list[dict] = []
    t0 = time.time()
    n_total = len(perseus_batch)
    for idx, (_, plot_row) in enumerate(perseus_batch.iterrows()):
        try:
            pr, dr = process_plot_scenario(
                plot_row, fia_trees, nsbe, engines, default_cfgs,
                args.n_draws, args.variants,
                args.rcp, args.harvest, si_ratios, harvest_schedule,
                seed=args.seed,
            )
            all_point.extend(pr)
            all_draws.extend(dr)
        except Exception as exc:
            logger.error(f"Plot {plot_row.get('PLOT')}: {exc}")
        if (idx + 1) % 5 == 0:
            elapsed = max(time.time() - t0, 1)
            rate = (idx + 1) / elapsed * 60
            logger.info(
                f"Processed {idx+1}/{n_total} plots ({rate:.1f}/min, "
                f"{len(all_point)} point + {len(all_draws)} draw rows)"
            )

    # ------------------------- Write outputs --------------------------------
    scenario = _scenario_tag(args.rcp, args.harvest)
    scen_dir = os.path.join(args.output_dir, scenario)
    os.makedirs(scen_dir, exist_ok=True)
    batch_tag = f"batch{args.batch_id}" if args.batch_id else "all"

    point_path = os.path.join(
        scen_dir, f"perseus_{scenario}_unc_point_{batch_tag}.csv"
    )
    draws_path = os.path.join(
        scen_dir, f"perseus_{scenario}_unc_draws_{batch_tag}.csv"
    )
    meta_path = os.path.join(
        scen_dir, f"perseus_{scenario}_unc_meta_{batch_tag}.json"
    )

    pd.DataFrame(all_point).to_csv(point_path, index=False)
    pd.DataFrame(all_draws).to_csv(draws_path, index=False)
    with open(meta_path, "w") as f:
        json.dump({
            "scenario": scenario, "rcp": args.rcp, "harvest": args.harvest,
            "batch_id": args.batch_id, "batch_size": args.batch_size,
            "n_plots_in_batch": int(n_total),
            "n_draws_requested": args.n_draws,
            "n_point_records": len(all_point),
            "n_draw_records": len(all_draws),
            "wall_clock_s": round(time.time() - t0, 1),
            "variants": args.variants,
        }, f, indent=2)
    logger.info(f"Wrote point records: {point_path}")
    logger.info(f"Wrote draw records:  {draws_path}")
    logger.info(f"Complete in {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
