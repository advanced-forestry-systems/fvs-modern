#!/usr/bin/env python3
from __future__ import annotations

"""
PERSEUS 100-Year FVS Projection with Uncertainty Estimation
============================================================

Extends the PERSEUS projection pipeline to include Bayesian uncertainty
quantification. For each plot x variant combination, runs FVS with:
  1. Default parameters (single run)
  2. Calibrated MAP parameters (single run)
  3. N posterior draws from the calibrated distribution (ensemble)

The posterior draws propagate parametric uncertainty through the FVS
projection, producing credible intervals on state-level MMT estimates.

This script is designed for FIA plots with initial measurements in
1999-2004 (Maine FIA annual panel design, cycles 5-6), enabling
100-year projections from the earliest available data.

Pipeline:
  1. Read PERSEUS CSV, filter to FIRST_INVYR in [1999, 2004]
  2. Load FIA tree data for those PLT_CNs
  3. For each plot:
     a. Run FVS-ACD default + calibrated
     b. Run FVS-NE default + calibrated
     c. For each variant, run n_draws posterior draws
  4. Aggregate to state-level MMT with credible intervals
  5. Output: plot-level CSV + state-level MMT with uncertainty bands

Usage:
  python perseus_uncertainty_projection.py --batch-id 1 --batch-size 100
  python perseus_uncertainty_projection.py --all --n-draws 100
  python perseus_uncertainty_projection.py --all --n-draws 50 --seed 42

Author: A. Weiskittel
Date: 2026-04-11
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
# Environment and path setup
# ---------------------------------------------------------------------------

PROJECT_ROOT = os.environ.get(
    "FVS_PROJECT_ROOT",
    str(Path(__file__).resolve().parents[2])
)
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.join(PROJECT_ROOT, "deployment", "fvs2py"))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "deployment", "microfvs"))

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

# Import shared components from the base projection script
from perseus_100yr_projection import (
    NSBECalculator,
    load_fia_trees_for_plots,
    build_fvs_standinit,
    build_fvs_treeinit,
    run_fvs_projection,
    compute_plot_agb,
    compute_initial_agb,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ===========================================================================
# Constants
# ===========================================================================

MAX_FVS_FAILURES_PER_PLOT = 4  # Skip remaining combos after this many failures

# Carbon fraction: convert aboveground biomass (AGB) to aboveground carbon (AGC)
# Standard value 0.5; species-specific values range 0.47-0.52 (FIA CRM docs)
CARBON_FRACTION = 0.5


# ===========================================================================
# Uncertainty-Aware FVS Runner
# ===========================================================================

def run_fvs_with_draw(
    stand_init_df: pd.DataFrame,
    tree_init_df: pd.DataFrame,
    stand_id: str,
    variant: str,
    draw_keywords: str,
    num_cycles: int = 20,
    cycle_length: int = 5,
) -> dict:
    """Run FVS with a specific posterior draw injected as keywords.

    This function is identical to run_fvs_projection() except that it
    accepts pre-formatted keyword strings from the UncertaintyEngine
    rather than loading from the calibrated config.

    Args:
        stand_init_df: standinit table (single row)
        tree_init_df: treeinit table
        stand_id: stand identifier
        variant: FVS variant code
        draw_keywords: FVS keyword block from UncertaintyEngine
        num_cycles: number of projection cycles
        cycle_length: years per cycle

    Returns:
        dict with summary, treelists, exit_code
    """
    from perseus_100yr_projection import (
        KEYFILE_TEMPLATE,
        _run_via_subprocess,
        _run_via_fvs2py,
        _extract_treelists,
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "FVS_Data.db")
        conn = sqlite3.connect(db_path)
        stand_init_df.to_sql("fvs_standinit", conn, if_exists="replace", index=False)
        tree_init_df.to_sql("fvs_treeinit", conn, if_exists="replace", index=False)
        conn.close()

        keyfile_content = KEYFILE_TEMPLATE.format(
            stand_id=stand_id,
            db_path=db_path,
            calibration_keywords=draw_keywords,
            num_cycles=num_cycles,
            cycle_length=cycle_length,
        )

        keyfile_path = os.path.join(tmpdir, f"{variant}_{stand_id}_draw.key")
        with open(keyfile_path, "w") as f:
            f.write(keyfile_content)

        # Try subprocess first, then fvs2py
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


# ===========================================================================
# Process Plot with Uncertainty Draws
# ===========================================================================

def process_plot_with_uncertainty(
    plot_row: dict,
    fia_trees: pd.DataFrame,
    nsbe: NSBECalculator,
    uncertainty_engines: dict,
    default_configs: dict,
    n_draws: int = 50,
    variants: list[str] = None,
    seed: int = 42,
) -> tuple[list[dict], list[dict]]:
    """Process one PERSEUS plot with point estimates and uncertainty ensemble.

    Returns:
        Tuple of (point_results, draw_results)
        point_results: list of dicts for default + calibrated runs
        draw_results: list of dicts for posterior draw ensemble
    """
    if variants is None:
        variants = ["acd", "ne"]

    plot_id = plot_row["PLOT"]
    raw_cn = plot_row["FIRST_PLTCN"]
    plt_cn = str(int(float(raw_cn))) if pd.notna(raw_cn) else ""
    stand_id = f"P{int(float(plot_id))}"

    # Guard against empty tree data
    if fia_trees.empty or "PLT_CN" not in fia_trees.columns:
        logger.warning(f"Plot {plot_id}: no tree data available")
        return [], []

    # Filter trees for this plot
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

    initial_agb = compute_initial_agb(plot_trees, nsbe)

    point_results = []
    draw_results = []

    for variant in variants:
        stand_df = build_fvs_standinit(plot_data, stand_id, variant)
        tree_df = build_fvs_treeinit(plot_trees, stand_id)

        if tree_df.empty:
            logger.warning(f"Plot {plot_id}/{variant}: no valid trees")
            continue

        # --- Point estimates: default and calibrated ---
        for config in [None, "calibrated"]:
            config_label = config if config else "default"

            # Add initial condition (year 0)
            point_results.append({
                "PLOT": plot_id,
                "FIRST_PLTCN": plt_cn,
                "YEAR": inv_year,
                "PROJ_YEAR": 0,
                "VARIANT": variant.upper(),
                "CONFIG": config_label,
                "AGB_TONS_AC": round(initial_agb, 4),
            })

            try:
                fvs_result = run_fvs_projection(
                    stand_df, tree_df, stand_id, variant,
                    config_version=config,
                    num_cycles=20,
                    cycle_length=5,
                )
                for cycle_year, treelist in sorted(fvs_result["treelists"].items()):
                    proj_year = cycle_year - inv_year
                    agb = compute_plot_agb(treelist, nsbe)
                    point_results.append({
                        "PLOT": plot_id,
                        "FIRST_PLTCN": plt_cn,
                        "YEAR": cycle_year,
                        "PROJ_YEAR": proj_year,
                        "VARIANT": variant.upper(),
                        "CONFIG": config_label,
                        "AGB_TONS_AC": round(agb, 4),
                    })
            except Exception as e:
                logger.error(f"Plot {plot_id}/{variant}/{config_label}: {e}")

        # --- Posterior draws for calibrated uncertainty ---
        engine = uncertainty_engines.get(variant.lower())
        default_cfg = default_configs.get(variant.lower())

        if engine is None or default_cfg is None:
            logger.warning(f"No uncertainty engine for {variant}; skipping draws")
            continue

        rng = np.random.default_rng(seed + int(float(plot_id)))
        draw_indices = rng.choice(engine.n_draws, size=n_draws, replace=False)

        for di, draw_idx in enumerate(draw_indices):
            draw = engine.get_draw(int(draw_idx))
            draw_keywords = engine.generate_keywords_for_draw(
                draw, default_cfg, draw_idx=int(draw_idx)
            )
            # GROWMULT keywords trigger FVS STOP 20; filter them out
            draw_keywords = "\n".join(
                line for line in draw_keywords.splitlines()
                if not line.startswith("GROWMULT")
            )

            # Add initial condition for this draw
            draw_results.append({
                "PLOT": plot_id,
                "FIRST_PLTCN": plt_cn,
                "YEAR": inv_year,
                "PROJ_YEAR": 0,
                "VARIANT": variant.upper(),
                "DRAW": int(draw_idx),
                "AGB_TONS_AC": round(initial_agb, 4),
            })

            try:
                fvs_result = run_fvs_with_draw(
                    stand_df, tree_df, stand_id, variant,
                    draw_keywords=draw_keywords,
                    num_cycles=20,
                    cycle_length=5,
                )
                for cycle_year, treelist in sorted(fvs_result["treelists"].items()):
                    proj_year = cycle_year - inv_year
                    agb = compute_plot_agb(treelist, nsbe)
                    draw_results.append({
                        "PLOT": plot_id,
                        "FIRST_PLTCN": plt_cn,
                        "YEAR": cycle_year,
                        "PROJ_YEAR": proj_year,
                        "VARIANT": variant.upper(),
                        "DRAW": int(draw_idx),
                        "AGB_TONS_AC": round(agb, 4),
                    })
            except Exception as e:
                logger.error(
                    f"Plot {plot_id}/{variant}/draw{draw_idx}: {e}"
                )

    return point_results, draw_results


# ===========================================================================
# State-Level Aggregation with Uncertainty
# ===========================================================================

def aggregate_with_uncertainty(
    point_df: pd.DataFrame,
    draw_df: pd.DataFrame,
    expns_df: pd.DataFrame | None = None,
    n_plots_total: int = 3586,
    quantiles: list[float] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Aggregate plot-level results to state-level MMT with credible intervals.

    Args:
        point_df: Point estimate results (default + calibrated)
        draw_df: Posterior draw ensemble results
        expns_df: Optional DataFrame with PLOT and EXPNS columns
        n_plots_total: Total plots for approximate expansion
        quantiles: Quantile levels for credible intervals

    Returns:
        Tuple of (point_mmt, uncertainty_mmt)
    """
    if quantiles is None:
        quantiles = [0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975]

    # Deduplicate: if both a manual initial record and an FVS cycle-0
    # treelist exist at PROJ_YEAR(_BIN)=0, keep only the first (FIA-based)
    yr_col_dedup = "PROJ_YEAR_BIN" if "PROJ_YEAR_BIN" in point_df.columns else "PROJ_YEAR"
    for df in [point_df, draw_df]:
        if df.empty:
            continue
        dedup_cols = ["PLOT", "VARIANT"]
        if "CONFIG" in df.columns:
            dedup_cols.append("CONFIG")
        if "DRAW" in df.columns:
            dedup_cols.append("DRAW")
        dedup_cols.append(yr_col_dedup)
        before = len(df)
        df.drop_duplicates(subset=dedup_cols, keep="first", inplace=True)
        after = len(df)
        if before != after:
            logger.info(
                f"Deduplicated {before - after} duplicate records "
                f"(cycle-0 treelist overlap)"
            )

    # Add expansion factors
    maine_forest_acres = 17_600_000

    for df in [point_df, draw_df]:
        if df.empty:
            continue
        if expns_df is not None:
            df_merged = df.merge(expns_df[["PLOT", "EXPNS"]], on="PLOT", how="left")
            df_merged["EXPNS"] = df_merged["EXPNS"].fillna(
                maine_forest_acres / n_plots_total
            )
            df.loc[:, "EXPNS"] = df_merged["EXPNS"].values
        else:
            approx_expns = maine_forest_acres / n_plots_total
            df.loc[:, "EXPNS"] = approx_expns

    # --- Point estimates ---
    # Use PROJ_YEAR_BIN if available (cycle-aligned bins), else PROJ_YEAR
    point_mmt = pd.DataFrame()
    if not point_df.empty:
        yr_col = "PROJ_YEAR_BIN" if "PROJ_YEAR_BIN" in point_df.columns else "PROJ_YEAR"
        point_mmt = (
            point_df
            .groupby(["VARIANT", "CONFIG", yr_col])
            .apply(
                lambda g: pd.Series({
                    "MMT_AGB": (g["AGB_TONS_AC"] * g["EXPNS"]).sum() * 0.907185 / 1e6,
                    "MMT": (g["AGB_TONS_AC"] * g["EXPNS"]).sum() * 0.907185 / 1e6 * CARBON_FRACTION,
                    "MEAN_AGB_TONS_AC": g["AGB_TONS_AC"].mean(),
                    "N_PLOTS": g["PLOT"].nunique(),
                }),
                include_groups=False,
            )
            .reset_index()
        )
        # Normalize column name to PROJ_YEAR for downstream compatibility
        if yr_col != "PROJ_YEAR":
            point_mmt = point_mmt.rename(columns={yr_col: "PROJ_YEAR"})

    # --- Uncertainty ensemble ---
    uncertainty_mmt = pd.DataFrame()
    if not draw_df.empty:
        # Compute MMT per draw (AGC = AGB * carbon fraction)
        draw_mmt = (
            draw_df
            .groupby(["VARIANT", "DRAW", "PROJ_YEAR"])
            .apply(
                lambda g: pd.Series({
                    "MMT": (g["AGB_TONS_AC"] * g["EXPNS"]).sum() * 0.907185 / 1e6 * CARBON_FRACTION,
                    "N_PLOTS": g["PLOT"].nunique(),
                }),
                include_groups=False,
            )
            .reset_index()
        )

        # Summarize across draws
        uncertainty_mmt = (
            draw_mmt
            .groupby(["VARIANT", "PROJ_YEAR"])
            .agg(
                MMT_MEAN=("MMT", "mean"),
                MMT_SD=("MMT", "std"),
                N_DRAWS=("DRAW", "nunique"),
                **{
                    f"MMT_Q{int(q*1000):04d}": ("MMT", lambda x, q=q: x.quantile(q))
                    for q in quantiles
                },
            )
            .reset_index()
        )

    return point_mmt, uncertainty_mmt


# ===========================================================================
# Main
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description="PERSEUS 100-year FVS projection with uncertainty"
    )
    parser.add_argument("--batch-id", type=int, default=None)
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--perseus-csv", type=str, default=PERSEUS_CSV)
    parser.add_argument(
        "--variants", type=str, nargs="+", default=["acd", "ne"],
    )
    parser.add_argument(
        "--n-draws", type=int, default=50,
        help="Number of posterior draws per variant per plot"
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--start-year", type=int, default=1999,
        help="Earliest FIRST_INVYR to include"
    )
    parser.add_argument(
        "--end-year", type=int, default=2004,
        help="Latest FIRST_INVYR to include"
    )
    parser.add_argument(
        "--expansion-csv", type=str, default=None,
        help="CSV with PLOT and EXPNS columns"
    )
    parser.add_argument("--output-dir", type=str, default=OUTPUT_DIR)
    parser.add_argument(
        "--exclude-plots", type=str, default="",
        help="Comma-separated plot IDs to skip (e.g., 2869.0,3104.5), or path to "
             "a text file with one plot ID per line. Plots are matched against "
             "PLOT column in perseus_plots.csv as string."
    )
    args = parser.parse_args()

    # Load PERSEUS plots and filter to requested year range
    perseus = pd.read_csv(args.perseus_csv)
    logger.info(f"Loaded {len(perseus)} total PERSEUS plots")

    mask = (
        (perseus["FIRST_INVYR"] >= args.start_year) &
        (perseus["FIRST_INVYR"] <= args.end_year)
    )
    perseus = perseus.loc[mask].copy()
    logger.info(
        f"Filtered to {len(perseus)} plots with FIRST_INVYR in "
        f"[{args.start_year}, {args.end_year}]"
    )

    year_counts = perseus["FIRST_INVYR"].value_counts().sort_index()
    for yr, cnt in year_counts.items():
        logger.info(f"  FIRST_INVYR={yr}: {cnt} plots")

    # Apply plot exclusion filter (for pathological plots that hang FVS)
    if args.exclude_plots:
        raw = args.exclude_plots
        if os.path.exists(raw):
            with open(raw) as f:
                excludes = {line.strip() for line in f if line.strip()}
        else:
            excludes = {s.strip() for s in raw.split(",") if s.strip()}
        if excludes:
            # Type-robust normalization: "2869", "2869.0", 2869, 2869.0 all match
            def _norm(x):
                try:
                    return str(int(float(x)))
                except Exception:
                    return str(x).strip()
            norm_excludes = {_norm(x) for x in excludes}
            before = len(perseus)
            perseus = perseus.loc[
                ~perseus["PLOT"].apply(_norm).isin(norm_excludes)
            ].copy()
            dropped = before - len(perseus)
            logger.info(
                f"Excluded {dropped} plot(s) matching "
                f"{sorted(excludes)} (requested {len(excludes)})"
            )

    # Batch slicing
    if args.batch_id is not None:
        start = (args.batch_id - 1) * args.batch_size
        end = min(start + args.batch_size, len(perseus))
        perseus = perseus.iloc[start:end]
        logger.info(f"Batch {args.batch_id}: plots {start+1} to {end}")
    elif not args.all:
        logger.error("Specify --batch-id or --all")
        sys.exit(1)

    if perseus.empty:
        logger.info("No plots in this batch. Exiting.")
        return

    # Initialize NSBE calculator
    nsbe = NSBECalculator(NSBE_ROOT)

    # Initialize uncertainty engines for each variant
    from config.uncertainty import UncertaintyEngine
    from config.config_loader import FvsConfigLoader

    uncertainty_engines = {}
    default_configs = {}

    for variant in args.variants:
        # Uncertainty engine for posterior draws
        try:
            engine = UncertaintyEngine(variant, config_dir=CONFIG_DIR, seed=args.seed)
            if engine.draws_available:
                uncertainty_engines[variant] = engine
                logger.info(
                    f"{variant.upper()}: {engine.n_draws} posterior draws available"
                )
            else:
                logger.warning(f"{variant.upper()}: no draws file found")
        except Exception as e:
            logger.warning(f"{variant.upper()} uncertainty engine failed: {e}")

        # Default config for multiplier computation
        try:
            loader = FvsConfigLoader(variant, version="default", config_dir=CONFIG_DIR)
            default_configs[variant] = loader.config
        except Exception as e:
            logger.warning(f"{variant.upper()} default config failed: {e}")

    # Load FIA trees
    plt_cns = [
        str(int(float(x)))
        for x in perseus["FIRST_PLTCN"].dropna().unique()
    ]
    logger.info(f"Loading FIA trees for {len(plt_cns)} unique PLT_CNs...")
    fia_trees = load_fia_trees_for_plots(plt_cns, FIA_DIR)
    logger.info(f"Loaded {len(fia_trees)} tree records")

    # Load expansion factors if provided
    expns_df = None
    if args.expansion_csv and os.path.exists(args.expansion_csv):
        expns_df = pd.read_csv(args.expansion_csv)
        logger.info(f"Loaded expansion factors for {len(expns_df)} plots")

    # Process plots
    all_point = []
    all_draws = []
    t0 = time.time()

    for idx, (_, plot_row) in enumerate(perseus.iterrows()):
        try:
            point, draws = process_plot_with_uncertainty(
                plot_row.to_dict(),
                fia_trees,
                nsbe,
                uncertainty_engines,
                default_configs,
                n_draws=args.n_draws,
                variants=args.variants,
                seed=args.seed,
            )
            all_point.extend(point)
            all_draws.extend(draws)
        except Exception as exc:
            logger.warning(
                f"Plot {plot_row.get('PLOT', '?')} failed: {exc}"
            )

        if (idx + 1) % 25 == 0:
            elapsed = time.time() - t0
            rate = (idx + 1) / elapsed * 60
            logger.info(
                f"Processed {idx+1}/{len(perseus)} plots "
                f"({rate:.1f} plots/min, {elapsed:.0f}s, "
                f"{len(all_point)} point + {len(all_draws)} draw records)"
            )

    # Save results
    os.makedirs(args.output_dir, exist_ok=True)
    batch_suffix = f"_batch{args.batch_id}" if args.batch_id else ""

    point_df = pd.DataFrame(all_point)
    draws_df = pd.DataFrame(all_draws)

    point_path = os.path.join(
        args.output_dir,
        f"perseus_uncertainty_point{batch_suffix}.csv"
    )
    draws_path = os.path.join(
        args.output_dir,
        f"perseus_uncertainty_draws{batch_suffix}.csv"
    )

    if not point_df.empty:
        point_df.to_csv(point_path, index=False)
        logger.info(f"Point estimates: {len(point_df)} records -> {point_path}")

    if not draws_df.empty:
        draws_df.to_csv(draws_path, index=False)
        logger.info(f"Draw results: {len(draws_df)} records -> {draws_path}")

    # Aggregate to state level with uncertainty
    if not point_df.empty or not draws_df.empty:
        point_mmt, uncertainty_mmt = aggregate_with_uncertainty(
            point_df, draws_df,
            expns_df=expns_df,
            n_plots_total=len(perseus),
        )

        if not point_mmt.empty:
            mmt_path = os.path.join(
                args.output_dir,
                f"perseus_uncertainty_mmt_point{batch_suffix}.csv"
            )
            point_mmt.to_csv(mmt_path, index=False)
            logger.info(f"Point MMT: {mmt_path}")

            print("\n" + "=" * 80)
            print("Point Estimates: State-Level AGB (MMT)")
            print("=" * 80)
            pivot = point_mmt.pivot_table(
                index="PROJ_YEAR",
                columns=["VARIANT", "CONFIG"],
                values="MMT",
            )
            print(pivot.to_string(float_format="{:.2f}".format))

        if not uncertainty_mmt.empty:
            unc_path = os.path.join(
                args.output_dir,
                f"perseus_uncertainty_mmt_bands{batch_suffix}.csv"
            )
            uncertainty_mmt.to_csv(unc_path, index=False)
            logger.info(f"Uncertainty bands: {unc_path}")

            print("\n" + "=" * 80)
            print("Calibrated Uncertainty Bands (MMT): Median [90% CI]")
            print("=" * 80)
            for variant in args.variants:
                vdf = uncertainty_mmt[
                    uncertainty_mmt["VARIANT"] == variant.upper()
                ].sort_values("PROJ_YEAR")
                if vdf.empty:
                    continue
                print(f"\n{variant.upper()} calibrated:")
                for _, row in vdf.iterrows():
                    pyr = int(row["PROJ_YEAR"])
                    med = row.get("MMT_Q0500", row.get("MMT_MEAN", 0))
                    lo = row.get("MMT_Q0100", row.get("MMT_Q0025", 0))
                    hi = row.get("MMT_Q0900", row.get("MMT_Q0975", 0))
                    print(f"  Year +{pyr:>3d}: {med:8.1f} [{lo:8.1f}, {hi:8.1f}]")

    # Save metadata
    meta = {
        "description": "PERSEUS 100yr projection with Bayesian uncertainty",
        "year_range": [args.start_year, args.end_year],
        "n_plots": len(perseus),
        "n_draws": args.n_draws,
        "seed": args.seed,
        "variants": args.variants,
        "year_distribution": year_counts.to_dict(),
        "elapsed_seconds": round(time.time() - t0, 1),
    }
    meta_path = os.path.join(
        args.output_dir,
        f"perseus_uncertainty_meta{batch_suffix}.json"
    )
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=2, default=str)

    elapsed = time.time() - t0
    logger.info(f"Complete in {elapsed:.0f}s")


if __name__ == "__main__":
    main()
