#!/usr/bin/env python3
"""
Bakuzis Matrix with Parametric Uncertainty: Default vs. Calibrated
===================================================================

Extends bakuzis_100yr_comparison.py to include posterior draw ensembles
for long term (100 year) projection uncertainty quantification.

For each of 36 scenarios (4 species groups x 3 sites x 3 densities) and
each variant, runs:
  - 1 default parameter projection
  - 1 calibrated MAP projection (posterior median)
  - N posterior draw projections (credible bands)

The posterior draws propagate parametric uncertainty from the Bayesian
calibration through to stand level BA, TPA, volume, and mortality
trajectories over 100 years. The ensemble spread represents the range
of plausible outcomes given the posterior.

This is complementary to perseus_uncertainty_projection.py: Perseus runs
plot level projections from actual FIA remeasurements, whereas this
script uses a synthetic factorial design to expose how uncertainty
scales with species, site quality, and initial density.

Output (long format CSV, one row per year x scenario x config x draw):
  year, scenario, species_group, site_class, density_class,
  variant, config, draw_id, age, tpa, atba, tcuft, mcuft, mort, ...

  config is one of: default, calibrated_map, posterior
  draw_id is NaN for default/calibrated_map, 0..N-1 for posterior

Usage:
  python bakuzis_uncertainty_comparison.py --variant ne --n-draws 50
  python bakuzis_uncertainty_comparison.py --all-variants --n-draws 100 --seed 42
  python bakuzis_uncertainty_comparison.py --variant ne --n-draws 50 \
      --batch-id 0 --batch-size 6   # for SLURM array jobs

Author: A. Weiskittel
Date: 2026-04-23
"""

from __future__ import annotations

import argparse
import json
import logging
import os
import sqlite3
import sys
import tempfile
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# Import shared scaffolding from the point estimate version
SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from bakuzis_100yr_comparison import (  # noqa: E402
    DENSITY_CLASSES,
    KEYFILE_TEMPLATE,
    SITE_CLASSES,
    SPECIES_GROUPS as SPECIES_GROUPS_EAST,
    generate_synthetic_stand,
)

# Variant specific species groups. The eastern set (NE, ACD) uses the
# original Spruce Fir / Northern Hardwood / Pine / Oak Pine mixes from
# the point estimate Bakuzis runner. PN, SN, and IE use locally
# appropriate species drawn from each variant's MAXSP table.

SPECIES_GROUPS_PN = {
    "Douglas-fir-Hemlock": {
        "label": "Douglas-fir / Hemlock",
        "species": [
            {"spcd": 202, "pct": 55, "name": "Douglas-fir"},
            {"spcd": 263, "pct": 25, "name": "Western hemlock"},
            {"spcd": 17,  "pct": 12, "name": "Grand fir"},
            {"spcd": 15,  "pct": 8,  "name": "White fir"},
        ],
        "site_species": 202,
        "fortypcd": 201,
    },
    "Mixed-Conifer": {
        "label": "Mixed Conifer",
        "species": [
            {"spcd": 15,  "pct": 30, "name": "White fir"},
            {"spcd": 17,  "pct": 25, "name": "Grand fir"},
            {"spcd": 202, "pct": 25, "name": "Douglas-fir"},
            {"spcd": 263, "pct": 20, "name": "Western hemlock"},
        ],
        "site_species": 15,
        "fortypcd": 270,
    },
    "Pine": {
        "label": "Pacific Northwest Pine",
        "species": [
            {"spcd": 122, "pct": 55, "name": "Ponderosa pine"},
            {"spcd": 108, "pct": 20, "name": "Lodgepole pine"},
            {"spcd": 119, "pct": 15, "name": "Western white pine"},
            {"spcd": 116, "pct": 10, "name": "Jeffrey pine"},
        ],
        "site_species": 122,
        "fortypcd": 221,
    },
    "Spruce-Fir": {
        "label": "Spruce-Fir (Pacific)",
        "species": [
            {"spcd": 98,  "pct": 40, "name": "Sitka spruce"},
            {"spcd": 263, "pct": 30, "name": "Western hemlock"},
            {"spcd": 19,  "pct": 20, "name": "Subalpine fir"},
            {"spcd": 15,  "pct": 10, "name": "White fir"},
        ],
        "site_species": 98,
        "fortypcd": 264,
    },
}

SPECIES_GROUPS_SN = {
    "Loblolly-Shortleaf": {
        "label": "Loblolly-Shortleaf Pine",
        "species": [
            {"spcd": 131, "pct": 55, "name": "Loblolly pine"},
            {"spcd": 110, "pct": 25, "name": "Shortleaf pine"},
            {"spcd": 132, "pct": 15, "name": "Virginia pine"},
            {"spcd": 316, "pct": 5,  "name": "Red maple"},
        ],
        "site_species": 131,
        "fortypcd": 161,
    },
    "Longleaf-Slash": {
        "label": "Longleaf-Slash Pine",
        "species": [
            {"spcd": 121, "pct": 50, "name": "Longleaf pine"},
            {"spcd": 111, "pct": 30, "name": "Slash pine"},
            {"spcd": 110, "pct": 15, "name": "Shortleaf pine"},
            {"spcd": 131, "pct": 5,  "name": "Loblolly pine"},
        ],
        "site_species": 121,
        "fortypcd": 141,
    },
    "Mixed-Pine-Hardwood": {
        "label": "Mixed Pine-Hardwood",
        "species": [
            {"spcd": 131, "pct": 40, "name": "Loblolly pine"},
            {"spcd": 316, "pct": 30, "name": "Red maple"},
            {"spcd": 129, "pct": 20, "name": "Eastern white pine"},
            {"spcd": 110, "pct": 10, "name": "Shortleaf pine"},
        ],
        "site_species": 131,
        "fortypcd": 401,
    },
    "Bottomland-Hardwood": {
        "label": "Bottomland Hardwood",
        "species": [
            {"spcd": 316, "pct": 35, "name": "Red maple"},
            {"spcd": 311, "pct": 30, "name": "Florida maple"},
            {"spcd": 313, "pct": 20, "name": "Boxelder"},
            {"spcd": 221, "pct": 15, "name": "Baldcypress"},
        ],
        "site_species": 316,
        "fortypcd": 601,
    },
}

SPECIES_GROUPS_IE = {
    "Douglas-fir-Mix": {
        "label": "Douglas-fir Mix",
        "species": [
            {"spcd": 202, "pct": 50, "name": "Douglas-fir"},
            {"spcd": 122, "pct": 25, "name": "Ponderosa pine"},
            {"spcd": 119, "pct": 15, "name": "Western white pine"},
            {"spcd": 17,  "pct": 10, "name": "Grand fir"},
        ],
        "site_species": 202,
        "fortypcd": 201,
    },
    "Lodgepole-Subalpine": {
        "label": "Lodgepole-Subalpine",
        "species": [
            {"spcd": 108, "pct": 50, "name": "Lodgepole pine"},
            {"spcd": 19,  "pct": 30, "name": "Subalpine fir"},
            {"spcd": 93,  "pct": 15, "name": "Engelmann spruce"},
            {"spcd": 73,  "pct": 5,  "name": "Western larch"},
        ],
        "site_species": 108,
        "fortypcd": 281,
    },
    "White-Pine-Larch": {
        "label": "White Pine-Larch",
        "species": [
            {"spcd": 119, "pct": 40, "name": "Western white pine"},
            {"spcd": 73,  "pct": 30, "name": "Western larch"},
            {"spcd": 202, "pct": 20, "name": "Douglas-fir"},
            {"spcd": 122, "pct": 10, "name": "Ponderosa pine"},
        ],
        "site_species": 119,
        "fortypcd": 211,
    },
    "Ponderosa-Mixed": {
        "label": "Ponderosa-Mixed",
        "species": [
            {"spcd": 122, "pct": 50, "name": "Ponderosa pine"},
            {"spcd": 202, "pct": 25, "name": "Douglas-fir"},
            {"spcd": 119, "pct": 15, "name": "Western white pine"},
            {"spcd": 108, "pct": 10, "name": "Lodgepole pine"},
        ],
        "site_species": 122,
        "fortypcd": 221,
    },
}


def get_species_groups(variant: str) -> dict:
    """Return the species group dictionary appropriate for the given variant."""
    v = variant.lower()
    if v in ("pn",):
        return SPECIES_GROUPS_PN
    if v in ("sn",):
        return SPECIES_GROUPS_SN
    if v in ("ie",):
        return SPECIES_GROUPS_IE
    # Default: eastern groups for NE, ACD, and any future eastern variants
    return SPECIES_GROUPS_EAST


# Variant specific stand location defaults. The synthetic stand
# generator in bakuzis_100yr_comparison.py hardcodes Maine values
# (region 9, state 23, county 19, lat 45, lon -69). Western and
# Southern variants need region/state/lat/lon appropriate to their
# coverage area or FVS rejects the stand on input.
LOCATION_DEFAULTS = {
    "ne":  dict(region=9, state=23, county=19, latitude=45.0,  longitude=-69.0),
    "acd": dict(region=9, state=23, county=19, latitude=45.0,  longitude=-69.0),
    "pn":  dict(region=6, state=41, county=51, latitude=44.5,  longitude=-122.5),  # Linn Co OR
    "sn":  dict(region=8, state=1,  county=125, latitude=33.5,  longitude=-86.8),   # Tuscaloosa AL
    "ie":  dict(region=1, state=16, county=49, latitude=46.5,  longitude=-115.5),  # Idaho County ID
}


def patch_stand_location(stand_df: pd.DataFrame, variant: str) -> pd.DataFrame:
    """Override the hardcoded Maine location values with variant defaults."""
    v = variant.lower()
    loc = LOCATION_DEFAULTS.get(v, LOCATION_DEFAULTS["ne"])
    df = stand_df.copy()
    for k, val in loc.items():
        if k in df.columns:
            df[k] = val
    return df


# Backwards compatible alias so existing imports continue to work
SPECIES_GROUPS = SPECIES_GROUPS_EAST

# ---------------------------------------------------------------------------
# Environment and path setup (no hardcoded Cardinal paths; all via env)
# ---------------------------------------------------------------------------

PROJECT_ROOT = os.environ.get(
    "FVS_PROJECT_ROOT",
    str(Path(__file__).resolve().parents[2]),
)
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.join(PROJECT_ROOT, "deployment", "fvs2py"))
sys.path.insert(0, os.path.join(PROJECT_ROOT, "deployment", "microfvs"))

FVS_LIB_DIR = os.environ.get("FVS_LIB_DIR", os.path.join(PROJECT_ROOT, "lib"))
CONFIG_DIR = os.environ.get("FVS_CONFIG_DIR", os.path.join(PROJECT_ROOT, "config"))
OUTPUT_DIR = os.environ.get(
    "BAKUZIS_OUTPUT_DIR",
    os.path.join(PROJECT_ROOT, "calibration", "output", "bakuzis"),
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Scenario enumeration
# ---------------------------------------------------------------------------

def enumerate_scenarios(species_groups: dict | None = None) -> list[tuple[str, str, str, int]]:
    """Return the full factorial of scenarios as (species, site, density, id).

    Scenario ids are 1 based and stable across batches so that SLURM
    array jobs can address specific scenarios by id. Pass a variant
    specific species_groups dict to enumerate scenarios for variants
    other than NE/ACD; default uses the eastern groups.
    """
    if species_groups is None:
        species_groups = SPECIES_GROUPS_EAST
    scenarios = []
    scenario_id = 0
    for sp_key in species_groups:
        for site_key in SITE_CLASSES:
            for dens_key in DENSITY_CLASSES:
                scenario_id += 1
                scenarios.append((sp_key, site_key, dens_key, scenario_id))
    return scenarios


def slice_batch(
    scenarios: list[tuple[str, str, str, int]],
    batch_id: int,
    batch_size: int,
) -> list[tuple[str, str, str, int]]:
    """Return the batch_id'th contiguous slice of size batch_size."""
    if batch_size <= 0:
        return scenarios
    start = batch_id * batch_size
    end = start + batch_size
    return scenarios[start:end]


# ---------------------------------------------------------------------------
# FVS runner (extended for posterior draws)
# ---------------------------------------------------------------------------

def run_scenario(
    stand_df: pd.DataFrame,
    tree_df: pd.DataFrame,
    stand_id: str,
    variant: str,
    config_version: str,
    draw_keywords: Optional[str] = None,
) -> pd.DataFrame:
    """Run one FVS scenario and return the summary DataFrame.

    Args:
        stand_df, tree_df: synthetic stand tables.
        stand_id: unique scenario identifier (used in STDIDENT).
        variant: FVS variant code (lowercased for library name).
        config_version: 'default' | 'calibrated' | 'posterior'.
        draw_keywords: If config_version == 'posterior', pass the keyword
            block from UncertaintyEngine.generate_keywords_for_draw().
            Ignored for 'default' and 'calibrated'.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        # Build the calibration keyword block
        cal_keywords = ""
        if config_version == "calibrated":
            try:
                from config.config_loader import FvsConfigLoader
                loader = FvsConfigLoader(
                    variant.lower(),
                    version="calibrated",
                    config_dir=CONFIG_DIR,
                )
                cal_keywords = loader.generate_keywords(include_comments=False)
            except Exception as e:
                logger.warning(f"Could not load calibrated config: {e}")
        elif config_version == "posterior" and draw_keywords:
            cal_keywords = draw_keywords

        # Eastern variants (NE, ACD) use the original DATABASE/STANDSQL
        # path that has been verified for them. Non-eastern variants
        # use INVENTORY mode (STDINFO + TREEFMT + TREEDATA inline) to
        # avoid the FVS_StandInit schema mismatch.
        eastern_variants = ("ne", "acd")
        if variant.lower() in eastern_variants:
            db_path = os.path.join(tmpdir, "FVS_Data.db")
            conn = sqlite3.connect(db_path)
            stand_df.to_sql("fvs_standinit", conn, if_exists="replace", index=False)
            tree_df.to_sql("fvs_treeinit", conn, if_exists="replace", index=False)
            conn.close()

            keyfile_content = KEYFILE_TEMPLATE.format(
                stand_id=stand_id,
                db_path=db_path,
                calibration_keywords=cal_keywords or "** DEFAULT PARAMETERS",
            )
        else:
            # INVENTORY mode for PN, SN, IE, and other non-eastern variants.
            from inventory_keyfile import make_inventory_keyfile
            keyfile_content = make_inventory_keyfile(
                stand_df=stand_df,
                tree_df=tree_df,
                stand_id=stand_id,
                variant=variant,
                inv_year=int(stand_df.iloc[0].get("inv_year", 2000)),
                num_cycles=20,
                calibration_keywords=cal_keywords,
            )

        keyfile_path = os.path.join(tmpdir, f"{variant}_{stand_id}.key")
        with open(keyfile_path, "w") as f:
            f.write(keyfile_content)

        # Shared library pathway (preferred)
        lib_path = os.path.join(FVS_LIB_DIR, f"FVS{variant.lower()}.so")
        if os.path.exists(lib_path):
            from fvs2py import FVS
            fvs = FVS(
                lib_path=lib_path,
                config_version=(
                    "calibrated" if config_version in ("calibrated", "posterior") else "default"
                ),
                config_dir=CONFIG_DIR,
            )
            fvs.load_keyfile(keyfile_path)
            fvs.run()
            summary = fvs.summary
            if summary is not None:
                return summary.copy()
            return pd.DataFrame()

        # Subprocess fallback
        import subprocess
        exe_path = f"/usr/local/bin/FVS{variant.lower()}"
        if os.path.exists(exe_path):
            subprocess.run(
                [exe_path, f"--keywordfile={os.path.basename(keyfile_path)}"],
                capture_output=True, cwd=tmpdir, timeout=180,
            )
            if os.path.exists(db_path):
                conn = sqlite3.connect(db_path)
                try:
                    return pd.read_sql_query("SELECT * FROM FVS_Summary2", conn)
                except Exception:
                    return pd.DataFrame()
                finally:
                    conn.close()

        logger.error(f"FVS not available for variant {variant}; expected {lib_path}")
        return pd.DataFrame()


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def run_all(
    variant: str,
    n_draws: int,
    seed: int,
    batch_id: int,
    batch_size: int,
    output_dir: str,
    use_fia: bool = False,
    fia_n_plots: int = 5,
) -> pd.DataFrame:
    """Run the full scenario x config x draw matrix for one variant.

    When use_fia is True, real FIA stands are sampled from
    fia_stand_generator instead of the synthetic generator. Each
    Bakuzis cell yields up to fia_n_plots replicates so the BA bin
    is genuinely represented rather than approximated by one
    synthetic stand.

    Returns a long DataFrame with columns identifying scenario,
    config, draw_id, replicate, and all FVS_Summary2 outputs.
    """
    np.random.seed(seed)

    # Initialize UncertaintyEngine once per variant
    engine = None
    default_config = None
    draw_indices: list[int] = []
    if n_draws > 0:
        try:
            from config.uncertainty import UncertaintyEngine
            from config.config_loader import FvsConfigLoader
            engine = UncertaintyEngine(
                variant.lower(), config_dir=CONFIG_DIR, seed=seed,
            )
            if not engine.draws_available:
                logger.warning(
                    f"No posterior draws file for variant {variant}; "
                    "running default + calibrated MAP only"
                )
                engine = None
            else:
                avail = engine.n_draws
                logger.info(f"Variant {variant}: {avail} posterior draws available")
                # Sample without replacement if n_draws <= avail, else with replacement
                if n_draws <= avail:
                    draw_indices = list(
                        np.random.default_rng(seed).choice(avail, size=n_draws, replace=False)
                    )
                else:
                    draw_indices = list(
                        np.random.default_rng(seed).choice(avail, size=n_draws, replace=True)
                    )
                # Load default config for the keyword generator
                default_loader = FvsConfigLoader(
                    variant.lower(), version="default", config_dir=CONFIG_DIR,
                )
                default_config = default_loader.config
        except Exception as e:
            logger.warning(f"UncertaintyEngine setup failed for {variant}: {e}")
            engine = None

    species_groups = get_species_groups(variant)
    scenarios = slice_batch(
        enumerate_scenarios(species_groups), batch_id, batch_size
    )
    logger.info(
        f"Variant {variant}: {len(scenarios)} scenarios "
        f"(batch {batch_id} of size {batch_size}) using "
        f"{len(species_groups)} species groups"
    )

    all_rows = []

    # Optional FIA generator setup
    fia_generator = None
    if use_fia:
        try:
            from fia_stand_generator import generate_real_stand
            fia_generator = generate_real_stand
            logger.info(
                f"FIA mode enabled: sampling up to {fia_n_plots} real "
                f"stands per cell for variant {variant}"
            )
        except Exception as e:
            logger.error(
                f"--use-fia requested but fia_stand_generator import failed: {e}. "
                "Falling back to synthetic generator."
            )
            fia_generator = None

    for sp_key, site_key, dens_key, scenario_num in scenarios:
        stand_id = f"BK{scenario_num:02d}"
        logger.info(
            f"  Scenario {scenario_num}: {sp_key} / {site_key} / {dens_key}"
        )

        # Choose generator: FIA (real plots) or synthetic
        scenario_stands: list[tuple] = []
        if fia_generator is not None:
            try:
                fia_stands = fia_generator(
                    variant=variant,
                    site_class=site_key,
                    density_class=dens_key,
                    n_plots=fia_n_plots,
                    seed=seed + scenario_num,
                )
                # Returns list of (stand_df, tree_df, cond_row); attach BK id
                for rep_idx, (sdf, tdf, _row) in enumerate(fia_stands, start=1):
                    rep_stand_id = f"{stand_id}R{rep_idx}"
                    sdf = sdf.copy()
                    sdf["stand_id"] = rep_stand_id
                    tdf = tdf.copy()
                    tdf["stand_id"] = rep_stand_id
                    scenario_stands.append((rep_stand_id, sdf, tdf, rep_idx))
                if not scenario_stands:
                    logger.warning(
                        f"    No FIA plots match {sp_key}/{site_key}/{dens_key}; "
                        "skipping this scenario."
                    )
                    continue
                logger.info(f"    sampled {len(scenario_stands)} real FIA stands")
            except Exception as e:
                logger.error(f"    FIA sampling failed: {e}; using synthetic")
                fia_generator = None

        if not scenario_stands:
            # Synthetic single-stand fallback (no FIA mode or no FIA matches)
            stand_df, tree_df = generate_synthetic_stand(
                species_groups[sp_key],
                SITE_CLASSES[site_key],
                DENSITY_CLASSES[dens_key],
                stand_id,
                variant,
            )
            # Override hardcoded Maine location values with variant defaults
            # so FVS-PN/SN/IE accepts the stand instead of rejecting it.
            stand_df = patch_stand_location(stand_df, variant)
            scenario_stands = [(stand_id, stand_df, tree_df, 1)]

        # Iterate over each stand replicate (1 for synthetic, up to fia_n_plots
        # for FIA mode) and run default + calibrated MAP + posterior draws.
        for rep_stand_id, stand_df, tree_df, rep_idx in scenario_stands:
            # -- Default and calibrated MAP ----------------------------------
            for config_label, config_version in [
                ("default", "default"),
                ("calibrated_map", "calibrated"),
            ]:
                try:
                    summary = run_scenario(
                        stand_df, tree_df, rep_stand_id, variant, config_version,
                    )
                    if summary is not None and not summary.empty:
                        summary = summary.copy()
                        summary["species_group"] = sp_key
                        summary["site_class"] = site_key
                        summary["density_class"] = dens_key
                        summary["variant"] = variant.upper()
                        summary["config"] = config_label
                        summary["draw_id"] = np.nan
                        summary["scenario"] = scenario_num
                        summary["replicate"] = rep_idx
                        all_rows.append(summary)
                except Exception as e:
                    logger.error(f"    {config_label} rep {rep_idx} failed: {e}")

            # -- Posterior draws --------------------------------------------
            if engine is not None and default_config is not None:
                for i, draw_idx in enumerate(draw_indices):
                    try:
                        draw = engine.get_draw(int(draw_idx))
                        draw_kw = engine.generate_keywords_for_draw(
                            draw, default_config, draw_idx=int(draw_idx),
                        )
                        summary = run_scenario(
                            stand_df, tree_df, rep_stand_id, variant,
                            "posterior", draw_keywords=draw_kw,
                        )
                        if summary is not None and not summary.empty:
                            summary = summary.copy()
                            summary["species_group"] = sp_key
                            summary["site_class"] = site_key
                            summary["density_class"] = dens_key
                            summary["variant"] = variant.upper()
                            summary["config"] = "posterior"
                            summary["draw_id"] = int(draw_idx)
                            summary["scenario"] = scenario_num
                            summary["replicate"] = rep_idx
                            all_rows.append(summary)
                    except Exception as e:
                        logger.error(f"    draw {draw_idx} rep {rep_idx} failed: {e}")

                    # Progress log every 10 draws
                    if (i + 1) % 10 == 0:
                        logger.info(
                            f"    rep {rep_idx} draws done: "
                            f"{i + 1}/{len(draw_indices)}"
                        )

    if not all_rows:
        return pd.DataFrame()

    df = pd.concat(all_rows, ignore_index=True)

    # Long form output file name
    suffix = f"_batch{batch_id:02d}" if batch_size > 0 else ""
    outfile = os.path.join(
        output_dir,
        f"bakuzis_uncertainty_{variant.lower()}_n{n_draws}{suffix}.csv",
    )
    os.makedirs(output_dir, exist_ok=True)
    df.to_csv(outfile, index=False)
    logger.info(f"Wrote {len(df)} rows to {outfile}")
    return df


def main():
    parser = argparse.ArgumentParser(
        description="Bakuzis matrix 100 year uncertainty comparison",
    )
    parser.add_argument("--variant", type=str, default="ne")
    parser.add_argument(
        "--all-variants", action="store_true",
        help="Run both NE and ACD (extend with more as needed)",
    )
    parser.add_argument(
        "--n-draws", type=int, default=50,
        help="Number of posterior draws per scenario (0 to skip posterior)",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--batch-id", type=int, default=0,
        help="SLURM array index; 0 for single run",
    )
    parser.add_argument(
        "--batch-size", type=int, default=0,
        help="Scenarios per batch (0 means all 36)",
    )
    parser.add_argument("--output-dir", type=str, default=OUTPUT_DIR)
    parser.add_argument(
        "--use-fia", action="store_true",
        help="Use real FIA stands (fia_stand_generator) instead of synthetic"
    )
    parser.add_argument(
        "--fia-n-plots", type=int, default=5,
        help="Number of real FIA plots to sample per Bakuzis cell"
    )
    args = parser.parse_args()

    variants = ["ne", "acd"] if args.all_variants else [args.variant]

    t0 = time.time()
    total_rows = 0
    for variant in variants:
        logger.info(f"\n{'=' * 60}")
        logger.info(f"Variant: {variant.upper()}")
        logger.info(f"{'=' * 60}")
        df = run_all(
            variant,
            n_draws=args.n_draws,
            seed=args.seed,
            batch_id=args.batch_id,
            batch_size=args.batch_size,
            output_dir=args.output_dir,
            use_fia=args.use_fia,
            fia_n_plots=args.fia_n_plots,
        )
        total_rows += len(df)

    elapsed = time.time() - t0
    logger.info(f"\nTotal rows: {total_rows} in {elapsed:.0f}s")


if __name__ == "__main__":
    main()
