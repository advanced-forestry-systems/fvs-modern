#!/usr/bin/env python3
from __future__ import annotations

"""
PERSEUS Harvest + Climate Factorial FVS Projection
====================================================

Implements the harvest arm of the PERSEUS Round 2 factorial design:
  2 (Harvest: Yes/No) x 3 (Climate: None/RCP4.5/RCP8.5) = 6 scenarios

Harvest design:
  - ~2% of ALL FIA plots harvested annually (~10% per 5-year FVS cycle)
  - 50% basal area removal (ThinBBA from below, proportional)
  - Targeting criteria (must meet ALL):
      1. Privately owned (OWNGRPCD 40-46)
      2. Not in reserve status (RESERVCD = 0)
      3. Not ecologically sensitive (COND_STATUS_CD = 1, i.e. accessible forest)
      4. High merchantable volume (ranked by VOLCFNET or BALIVE)
  - Plots can be re-entered on ~50-year rotation
  - Non-eligible plots are never harvested

Harvest scheduling:
  With N_total plots and 10% harvested per cycle, we need
  N_harvest_per_cycle = round(0.10 * N_total) plots harvested each cycle.
  Among eligible plots, we rank by basal area / merchantable volume and
  assign harvest cycles using a round-robin on the ranked list, with
  re-entry allowed after MIN_REENTRY_CYCLES (10 cycles = 50 years).

FVS keyword approach:
  ThinBBA with field 2 = -0.50 (negative = proportion of BA to RETAIN).
  This removes 50% of basal area from below (smallest trees first).

  Keyword line format (10-char fields):
    ThinBBA           0      -0.5       1.0         0       999         0

  Scheduled via TIMEINT + date field in ThinBBA to fire at specific cycles.

Usage:
  # All 6 scenarios, batch mode for SLURM
  python perseus_harvest_projection.py --batch-id 1 --batch-size 100

  # Single scenario
  python perseus_harvest_projection.py --rcp 4.5 --harvest --batch-id 1

  # No-harvest baseline (equivalent to climate script)
  python perseus_harvest_projection.py --rcp none --no-harvest --all

Author: A. Weiskittel
Date: 2026-04-18
"""

import argparse
import json
import logging
import os
import signal
import sqlite3
import sys
import tempfile
import time

import numpy as np
import pandas as pd
from pathlib import Path

# ---------------------------------------------------------------------------
# Environment and path setup
# ---------------------------------------------------------------------------

PROJECT_ROOT = os.environ.get(
    "FVS_PROJECT_ROOT",
    str(Path(__file__).resolve().parents[2])
)
sys.path.insert(0, PROJECT_ROOT)
sys.path.insert(0, os.path.join(PROJECT_ROOT, "calibration", "python"))
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
SI_RATIO_CSV = os.environ.get(
    "SI_RATIO_CSV",
    os.path.expanduser("~/SiteIndex/future_SI/perseus_si_ratios.csv")
)

# Import shared components
from perseus_100yr_projection import (
    NSBECalculator,
    load_fia_trees_for_plots,
    build_fvs_standinit,
    build_fvs_treeinit,
    run_fvs_projection,
    compute_plot_agb,
    compute_initial_agb,
)
from perseus_uncertainty_projection import (
    aggregate_with_uncertainty,
)
from config.uncertainty import UncertaintyEngine

from perseus_uncertainty_projection import run_fvs_with_draw  # added for draw support
from perseus_climate_projection import (
    load_si_ratios,
    get_plot_si_multiplier,
    get_plot_si_delta_ft,
    modify_standinit_si,
    modify_standinit_si_delta,
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FVS_CALL_TIMEOUT = 300  # 5 minutes per FVS call

# Harvest parameters
HARVEST_ANNUAL_RATE = 0.02       # 2% of all plots per year
CYCLE_LENGTH = 5                  # years per FVS cycle
NUM_CYCLES = 20                   # 100-year projection
HARVEST_CYCLE_RATE = HARVEST_ANNUAL_RATE * CYCLE_LENGTH  # 10% per cycle
MIN_REENTRY_CYCLES = 10           # 50 years before a plot can be re-harvested
BA_REMOVAL_PROPORTION = 0.50      # Remove 50% of basal area
BA_RETAIN_PROPORTION = 1.0 - BA_REMOVAL_PROPORTION  # Retain 50%

# FIA ownership codes for private land
PRIVATE_OWNGRPCD = set(range(40, 47))  # 40-46 inclusive

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Signal-based timeout (same as climate script)
# ---------------------------------------------------------------------------

class _TimeoutError(Exception):
    pass


def _timeout_handler(signum, frame):
    raise _TimeoutError(f"FVS call exceeded {FVS_CALL_TIMEOUT}s timeout")


def run_fvs_safe(func, *args, **kwargs):
    """Run an FVS function with signal-based timeout."""
    old_handler = signal.signal(signal.SIGALRM, _timeout_handler)
    signal.alarm(FVS_CALL_TIMEOUT)
    try:
        result = func(*args, **kwargs)
        signal.alarm(0)
        return result
    except _TimeoutError:
        raise TimeoutError(f"FVS call exceeded {FVS_CALL_TIMEOUT}s timeout")
    finally:
        signal.signal(signal.SIGALRM, old_handler)
        signal.alarm(0)


# ===========================================================================
# Harvest Eligibility and Scheduling
# ===========================================================================

def load_fia_cond_for_plots(plt_cns: list[str], fia_dir: str,
                            state: str = "ME") -> pd.DataFrame:
    """Load FIA COND table fields needed for harvest eligibility.

    Fetches OWNGRPCD, RESERVCD, COND_STATUS_CD, BALIVE, VOLCFNET
    for condition 1 of each plot.
    """
    # Try SQLite first
    db_candidates = [
        os.path.join(fia_dir, f"{state}_FIA.db"),
        os.path.join(fia_dir, state, f"{state}_FIA.db"),
        os.path.join(fia_dir, f"FIADB_{state}.db"),
        os.path.join(fia_dir, state, f"FIADB_{state}.db"),
    ]

    for db_path in db_candidates:
        if os.path.exists(db_path):
            return _load_cond_sqlite(db_path, plt_cns)

    # Fallback to CSV
    csv_candidates = [
        os.path.join(fia_dir, state, f"{state}_COND.csv"),
        os.path.join(fia_dir, f"{state}_COND.csv"),
    ]
    for csv_path in csv_candidates:
        if os.path.exists(csv_path):
            return _load_cond_csv(csv_path, plt_cns)

    raise FileNotFoundError(
        f"Cannot find FIA COND data for {state} in {fia_dir}"
    )


def _load_cond_sqlite(db_path: str, plt_cns: list[str]) -> pd.DataFrame:
    """Load COND data from SQLite."""
    logger.info(f"Loading FIA COND from SQLite: {db_path}")
    conn = sqlite3.connect(db_path)
    placeholders = ",".join(["?"] * len(plt_cns))

    sql = f"""
    SELECT c.PLT_CN, c.CONDID,
           c.OWNGRPCD, c.RESERVCD, c.COND_STATUS_CD,
           c.BALIVE, c.VOLCFNET, c.FORTYPCD,
           c.TPAGROW_UNADJ, c.TPAMORT_UNADJ,
           c.SICOND, c.SISP
    FROM COND c
    WHERE c.PLT_CN IN ({placeholders})
      AND c.CONDID = 1
    """

    df = pd.read_sql_query(sql, conn, params=plt_cns)
    conn.close()
    logger.info(f"Loaded COND records for {len(df)} plots")
    return df


def _load_cond_csv(csv_path: str, plt_cns: list[str]) -> pd.DataFrame:
    """Load COND data from CSV."""
    logger.info(f"Loading FIA COND from CSV: {csv_path}")
    plt_cn_set = set(str(int(float(cn))) for cn in plt_cns)
    plt_cn_int_set = set(int(float(cn)) for cn in plt_cns)

    chunks = []
    cols_needed = [
        "PLT_CN", "CONDID", "OWNGRPCD", "RESERVCD",
        "COND_STATUS_CD", "BALIVE", "VOLCFNET", "FORTYPCD",
        "TPAGROW_UNADJ", "TPAMORT_UNADJ",
        "SICOND", "SISP",
    ]

    for chunk in pd.read_csv(csv_path, chunksize=50000, low_memory=False):
        # Keep only columns that exist
        available = [c for c in cols_needed if c in chunk.columns]
        mask = (
            chunk["PLT_CN"].isin(plt_cn_int_set)
            | chunk["PLT_CN"].astype(str).isin(plt_cn_set)
        )
        condid_mask = chunk.get("CONDID", pd.Series([1]*len(chunk))) == 1
        subset = chunk.loc[mask & condid_mask, available]
        if not subset.empty:
            chunks.append(subset)

    if not chunks:
        logger.warning("No COND records found for requested plots")
        return pd.DataFrame()

    df = pd.concat(chunks, ignore_index=True)
    logger.info(f"Loaded COND records for {len(df)} plots from CSV")
    return df


def determine_harvest_eligibility(
    perseus_df: pd.DataFrame,
    cond_df: pd.DataFrame,
) -> pd.DataFrame:
    """Determine which plots are eligible for harvest.

    Criteria:
      1. Private ownership (OWNGRPCD in 40-46)
      2. Not reserved (RESERVCD = 0)
      3. Accessible forest land (COND_STATUS_CD = 1)
      4. Has valid BA or volume data for ranking

    Returns:
        DataFrame with PLOT, FIRST_PLTCN, ELIGIBLE (bool), BALIVE, VOLCFNET,
        and a HARVEST_PRIORITY rank (lower = harvest sooner).
    """
    # Normalize PLT_CN types for merging
    perseus = perseus_df.copy()
    perseus["PLT_CN_STR"] = perseus["FIRST_PLTCN"].apply(
        lambda x: str(int(float(x))) if pd.notna(x) else ""
    )

    cond = cond_df.copy()
    cond["PLT_CN_STR"] = cond["PLT_CN"].apply(
        lambda x: str(int(float(x))) if pd.notna(x) else ""
    )

    # Merge condition data onto PERSEUS plots (only columns that exist)
    merge_cols = ["PLT_CN_STR", "OWNGRPCD", "RESERVCD", "COND_STATUS_CD", "BALIVE"]
    if "VOLCFNET" in cond.columns:
        merge_cols.append("VOLCFNET")
    available_merge = [c for c in merge_cols if c in cond.columns]
    merged = perseus.merge(
        cond[available_merge],
        on="PLT_CN_STR",
        how="left",
    )

    # Ensure expected columns exist (fill missing with NaN)
    for col in ["OWNGRPCD", "RESERVCD", "COND_STATUS_CD", "BALIVE", "VOLCFNET"]:
        if col not in merged.columns:
            merged[col] = np.nan

    # Apply eligibility criteria
    private = merged["OWNGRPCD"].isin(PRIVATE_OWNGRPCD)
    not_reserve = merged["RESERVCD"].fillna(1) == 0
    accessible = merged["COND_STATUS_CD"].fillna(2) == 1
    has_volume = merged["BALIVE"].fillna(0) > 0

    merged["ELIGIBLE"] = private & not_reserve & accessible & has_volume

    # For volume ranking, use VOLCFNET if available, fall back to BALIVE
    merged["RANK_VOLUME"] = merged["VOLCFNET"].fillna(0)
    low_vol = merged["RANK_VOLUME"] <= 0
    merged.loc[low_vol, "RANK_VOLUME"] = merged.loc[low_vol, "BALIVE"].fillna(0)

    # Rank eligible plots by volume (highest first)
    eligible_mask = merged["ELIGIBLE"]
    merged["HARVEST_PRIORITY"] = np.nan
    if eligible_mask.any():
        # Higher volume = lower priority number (harvested first)
        merged.loc[eligible_mask, "HARVEST_PRIORITY"] = (
            merged.loc[eligible_mask, "RANK_VOLUME"]
            .rank(ascending=False, method="first")
        )

    n_eligible = eligible_mask.sum()
    n_total = len(merged)
    logger.info(
        f"Harvest eligibility: {n_eligible}/{n_total} plots eligible "
        f"({100*n_eligible/n_total:.1f}%)"
    )
    logger.info(
        f"  Private: {private.sum()}, Not-reserve: {not_reserve.sum()}, "
        f"Accessible: {accessible.sum()}, Has-volume: {has_volume.sum()}"
    )

    if n_eligible > 0:
        elig_vol = merged.loc[eligible_mask, "RANK_VOLUME"]
        logger.info(
            f"  Eligible plot volume: mean={elig_vol.mean():.1f}, "
            f"median={elig_vol.median():.1f}, "
            f"range=[{elig_vol.min():.1f}, {elig_vol.max():.1f}]"
        )

    return merged


def build_harvest_schedule(
    eligibility_df: pd.DataFrame,
    n_total_plots: int,
    num_cycles: int = NUM_CYCLES,
    seed: int = 42,
) -> dict[int, set[int]]:
    """Build a deterministic harvest schedule across all projection cycles.

    Returns:
        dict mapping cycle_number (1-based) -> set of PLOT IDs to harvest
        in that cycle.

    Logic:
      - n_harvest_per_cycle = round(HARVEST_CYCLE_RATE * n_total_plots)
      - Sort eligible plots by HARVEST_PRIORITY (highest volume first)
      - Assign plots to cycles in round-robin order
      - A plot can be re-entered MIN_REENTRY_CYCLES after its last harvest
      - If not enough eligible plots, harvest fewer in that cycle

    The schedule is deterministic given the seed and eligibility ranking,
    ensuring reproducibility across SLURM batches.
    """
    n_per_cycle = max(1, round(HARVEST_CYCLE_RATE * n_total_plots))

    eligible = eligibility_df.loc[
        eligibility_df["ELIGIBLE"]
    ].sort_values("HARVEST_PRIORITY")

    eligible_plots = eligible["PLOT"].tolist()
    n_eligible = len(eligible_plots)

    if n_eligible == 0:
        logger.warning("No eligible plots for harvest!")
        return {c: set() for c in range(1, num_cycles + 1)}

    logger.info(
        f"Harvest schedule: {n_per_cycle} plots/cycle from "
        f"{n_eligible} eligible ({n_per_cycle/n_eligible*100:.1f}% of eligible pool)"
    )

    # Build schedule with re-entry tracking
    schedule = {c: set() for c in range(1, num_cycles + 1)}
    last_harvest = {}  # plot_id -> last cycle harvested

    rng = np.random.default_rng(seed)

    for cycle in range(1, num_cycles + 1):
        # Find plots available this cycle (not recently harvested)
        available = []
        for pid in eligible_plots:
            last = last_harvest.get(pid, -MIN_REENTRY_CYCLES)
            if (cycle - last) >= MIN_REENTRY_CYCLES:
                available.append(pid)

        # Select up to n_per_cycle from available, prioritizing high volume
        # (available list is already sorted by priority from the eligible sort)
        n_to_harvest = min(n_per_cycle, len(available))
        if n_to_harvest > 0:
            selected = available[:n_to_harvest]
            schedule[cycle] = set(selected)
            for pid in selected:
                last_harvest[pid] = cycle

    # Summary
    total_harvested = sum(len(s) for s in schedule.values())
    unique_plots = len(set().union(*schedule.values()))
    logger.info(
        f"Harvest schedule summary: {total_harvested} total harvest events "
        f"across {num_cycles} cycles, {unique_plots} unique plots"
    )
    for c in [1, 5, 10, 15, 20]:
        if c <= num_cycles:
            logger.info(f"  Cycle {c:2d}: {len(schedule[c]):4d} plots")

    return schedule


# ===========================================================================
# FVS Keyword Template with Harvest
# ===========================================================================

HARVEST_KEYFILE_TEMPLATE = """STDIDENT
{stand_id}

DATABASE
DSNIN
{db_path}

DSNOUT
{db_path}

STANDSQL
SELECT *
FROM fvs_standinit
WHERE stand_id = '%StandID%'
ENDSQL

TREESQL
SELECT *
FROM fvs_treeinit
WHERE stand_id = '%StandID%'
ENDSQL

END

DATABASE
SUMMARY            2
TREELIDB           2         2
CARBREDB           2
COMPUTDB
END

CALBSTAT
TREELIST           0
ECHOSUM

TIMEINT            0{cycle_length:>10d}
NUMCYCLE{num_cycles:>10d}

{calibration_keywords}

{harvest_keywords}

PROCESS
STOP
"""



def _compute_current_ba(tree_init_df) -> float:
    """Compute stand basal area (sq ft/ac) from a fvs_treeinit DataFrame.

    BA per tree (sq ft) = 0.005454 * diameter^2 (DBH in inches).
    BA per acre = sum over trees of (tree_count * 0.005454 * diameter^2).
    """
    if tree_init_df is None or len(tree_init_df) == 0:
        return 0.0
    if "diameter" not in tree_init_df.columns or "tree_count" not in tree_init_df.columns:
        return 0.0
    dbh = tree_init_df["diameter"].astype(float)
    tpa = tree_init_df["tree_count"].astype(float)
    ba_per_tree = 0.005454 * (dbh ** 2)
    return float((tpa * ba_per_tree).sum())


def generate_harvest_keywords(
    harvest_cycles: list[int],
    inv_year: int,
    cycle_length: int = CYCLE_LENGTH,
    current_ba: float | None = None,
) -> str:
    """Generate FVS ThinBBA keywords for scheduled harvest cycles.

    For each harvest cycle, inserts a ThinBBA keyword that removes 50%
    of the standing basal area from below (smallest trees first).

    FVS ThinBBA field layout (10-char fields):
      Field 1: Date (calendar year) when thinning occurs
      Field 2: Residual BA target. Negative = proportion of current BA
               to RETAIN. So -0.50 = retain 50%, remove 50%.
      Field 3: Efficiency (1.0 = all marked trees removed)
      Field 4: Min DBH to cut (0 = no minimum)
      Field 5: Max DBH to cut (999 = no maximum)
      Field 6: Species code (0 = all species)

    Args:
        harvest_cycles: list of cycle numbers (1-based) when this plot
                       is harvested
        inv_year: inventory year (to compute calendar year for each cycle)
        cycle_length: years per FVS cycle

    Returns:
        String of FVS keywords (one ThinBBA per harvest event)
    """
    if not harvest_cycles:
        return "** NO HARVEST TREATMENT"

    # PATCH 27 April 2026: switch to absolute residual target. Original code
    # used negative target (-0.5) interpreted by FVS-NE as a clearcut sentinel
    # and silently cancelled by FVS-ACD. Absolute target = 0.5 * current_ba
    # produces intended ~50% retention in NE. ACD cut routine still broken
    # at executable level (see PATCH_fvs_acd_harvest.md).
    if current_ba is None or current_ba <= 0:
        target_ba = 50.0
    else:
        target_ba = max(10.0, 0.5 * current_ba)

    lines = [f"** HARVEST: thin to absolute residual {target_ba:.1f} sq ft/ac (50% retention)"]

    for cycle in sorted(harvest_cycles):
        harvest_year = inv_year + (cycle - 1) * cycle_length
        line = (
            f"ThinBBA   "
            f"{harvest_year:>10d}"
            f"{target_ba:>10.1f}"
            f"       1.0"
            f"         0"
            f"       999"
            f"         0"
        )
        lines.append(line)

    return "\n".join(lines)


# ===========================================================================
# Harvest-Aware FVS Runner
# ===========================================================================

def run_fvs_harvest_projection(
    stand_init_df: pd.DataFrame,
    tree_init_df: pd.DataFrame,
    stand_id: str,
    variant: str,
    harvest_cycles: list[int],
    inv_year: int,
    config_version: str = "calibrated",
    num_cycles: int = NUM_CYCLES,
    cycle_length: int = CYCLE_LENGTH,
    draw_keywords: str | None = None,
) -> dict:
    """Run FVS projection with harvest keywords injected.

    Similar to run_fvs_projection but uses the harvest keyfile template
    and injects ThinBBA keywords for scheduled harvest cycles.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path = os.path.join(tmpdir, "FVS_Data.db")
        conn = sqlite3.connect(db_path)
        stand_init_df.to_sql("fvs_standinit", conn, if_exists="replace", index=False)
        tree_init_df.to_sql("fvs_treeinit", conn, if_exists="replace", index=False)
        conn.close()

        # Calibration keywords
        cal_keywords = ""
        if draw_keywords is not None:
            # Use posterior draw keywords directly (bypasses config loading)
            cal_keywords = draw_keywords
        elif config_version == "calibrated":
            try:
                from config.config_loader import FvsConfigLoader
                loader = FvsConfigLoader(
                    variant.lower(),
                    version="calibrated",
                    config_dir=CONFIG_DIR,
                )
                cal_keywords = loader.generate_keywords(include_comments=False)
            except Exception as e:
                logger.warning(f"Could not load calibrated config for {variant}: {e}")

        # Harvest keywords (compute current_ba so target = 0.5 * BA)
        current_ba = _compute_current_ba(tree_init_df)
        harvest_kw = generate_harvest_keywords(
            harvest_cycles, inv_year, cycle_length, current_ba=current_ba
        )

        keyfile_content = HARVEST_KEYFILE_TEMPLATE.format(
            stand_id=stand_id,
            db_path=db_path,
            calibration_keywords=cal_keywords if cal_keywords else "** DEFAULT PARAMETERS",
            harvest_keywords=harvest_kw,
            num_cycles=num_cycles,
            cycle_length=cycle_length,
        )

        keyfile_path = os.path.join(tmpdir, f"{variant}_{stand_id}.key")
        with open(keyfile_path, "w") as f:
            f.write(keyfile_content)

        # Prefer subprocess (executable) over ctypes
        exe_candidates = [
            os.path.join(FVS_LIB_DIR, f"FVS{variant.lower()}"),
            os.path.join(FVS_LIB_DIR, f"FVS{variant.upper()}"),
            f"/usr/local/bin/FVS{variant.lower()}",
        ]

        # Import the subprocess runner from the base projection
        from perseus_100yr_projection import _run_via_subprocess, _run_via_fvs2py

        for exe_path in exe_candidates:
            if os.path.exists(exe_path) and os.access(exe_path, os.X_OK):
                return _run_via_subprocess(exe_path, keyfile_path, db_path, tmpdir)

        # Fallback to fvs2py
        lib_path = os.path.join(FVS_LIB_DIR, f"FVS{variant.lower()}.so")
        if os.path.exists(lib_path):
            return _run_via_fvs2py(
                lib_path, keyfile_path, db_path, variant, config_version
            )

        raise FileNotFoundError(
            f"Cannot find FVS executable or library for {variant}. "
            f"Tried: {exe_candidates + [lib_path]}"
        )


# ===========================================================================
# Process a Single Plot
# ===========================================================================

def process_plot_harvest(
    plot_row: dict,
    fia_trees: pd.DataFrame,
    nsbe: NSBECalculator,
    si_ratios: pd.DataFrame | None,
    rcp: str,
    harvest_schedule: dict[int, set[int]],
    do_harvest: bool,
    cond_df: pd.DataFrame | None = None,
    variants: list[str] = None,
    configs: list[str] = None,
    seed: int = 42,
    engines: dict | None = None,
    n_draws: int = 0,
    default_configs: dict | None = None,
) -> list[dict]:
    """Process one PERSEUS plot under a harvest + climate scenario.

    Runs point estimates for each variant x config combination.
    No posterior draws for harvest scenarios to keep computation manageable.

    Args:
        plot_row: dict with PLOT, FIRST_PLTCN, COUNTYCD, etc.
        fia_trees: all FIA tree records (will be filtered to this plot)
        nsbe: NSBE biomass calculator
        si_ratios: climate SI ratio table (None if rcp='none')
        rcp: climate scenario ("4.5", "8.5", or "none")
        harvest_schedule: {cycle -> set of plot IDs} schedule
        do_harvest: whether to apply harvest to this scenario
        cond_df: FIA COND table with SICOND, SISP (for site index)
        variants: list of FVS variants to run
        configs: list of config versions ("default", "calibrated")
        seed: random seed

    Returns:
        list of result dicts with per-cycle AGB values
    """
    if variants is None:
        variants = ["ne", "acd"]
    if configs is None:
        configs = ["default", "calibrated"]

    plot_id = plot_row["PLOT"]
    raw_cn = plot_row["FIRST_PLTCN"]
    plt_cn = str(int(float(raw_cn))) if pd.notna(raw_cn) else ""
    stand_id = f"P{int(float(plot_id))}"
    countycd = int(plot_row.get("COUNTYCD", 19))

    # Determine climate SI delta (absolute change in feet)
    si_delta_ft = 0.0
    si_mult = 1.0  # kept for output metadata
    if rcp != "none" and si_ratios is not None:
        si_delta_ft = get_plot_si_delta_ft(
            si_ratios, countycd, int(float(plot_id)), rcp
        )
        si_mult = get_plot_si_multiplier(
            si_ratios, countycd, int(float(plot_id)), rcp
        )

    # Find which cycles this plot is harvested
    harvest_cycles = []
    if do_harvest:
        for cycle, plot_set in harvest_schedule.items():
            if int(float(plot_id)) in plot_set:
                harvest_cycles.append(cycle)

    # Filter trees for this plot
    if fia_trees.empty or "PLT_CN" not in fia_trees.columns:
        logger.warning(f"Plot {plot_id}: no tree data available")
        return []

    fia_cn_str = fia_trees["PLT_CN"].apply(
        lambda x: str(int(float(x))) if pd.notna(x) else ""
    )
    plot_trees = fia_trees.loc[fia_cn_str == plt_cn]

    if plot_trees.empty:
        logger.warning(f"Plot {plot_id}: no trees for PLT_CN={plt_cn}")
        return []

    plot_data = plot_trees.iloc[0].to_dict()
    plot_data["COUNTYCD"] = countycd

    # Inject SICOND / SISP from COND table (not present in TREE table)
    if cond_df is not None and not cond_df.empty:
        cond_match = cond_df.loc[
            cond_df["PLT_CN"].astype(str).str.strip()
            == str(int(float(plt_cn)))
        ]
        if not cond_match.empty:
            cond_row = cond_match.iloc[0]
            for col in ("SICOND", "SISP"):
                val = cond_row.get(col)
                if pd.notna(val):
                    plot_data[col] = val

    inv_year = int(plot_data.get("INVYR", plot_row.get("FIRST_INVYR", 2000)))

    initial_agb = compute_initial_agb(plot_trees, nsbe)

    # Scenario label
    rcp_label = f"rcp{rcp.replace('.', '')}" if rcp != "none" else "noclimate"
    harvest_label = "harvest" if do_harvest else "noharvest"

    results = []

    for variant in variants:
        stand_df_base = build_fvs_standinit(plot_data, stand_id, variant)
        tree_df = build_fvs_treeinit(plot_trees, stand_id)

        if tree_df.empty:
            logger.warning(f"Plot {plot_id}/{variant}: no valid trees")
            continue

        # Apply climate SI modification (absolute delta in feet)
        if abs(si_delta_ft) > 0.01:
            stand_df_base = modify_standinit_si_delta(stand_df_base, si_delta_ft)

        variant_fail_count = 0

        for config in configs:
            scenario = f"{harvest_label}_{rcp_label}_{config}"

            # Initial condition (year 0, same for both configs)
            results.append({
                "PLOT": plot_id,
                "FIRST_PLTCN": plt_cn,
                "YEAR": inv_year,
                "PROJ_YEAR": 0,
                "PROJ_YEAR_BIN": 0,
                "VARIANT": variant.upper(),
                "CONFIG": config,
                "SCENARIO": scenario,
                "RCP": rcp,
                "HARVEST": do_harvest,
                "N_HARVESTS": len(harvest_cycles),
                "HARVEST_CYCLES": ",".join(str(c) for c in harvest_cycles) if harvest_cycles else "",
                "SI_MULT": round(si_mult, 4),
                "SI_DELTA_FT": round(si_delta_ft, 2),
                "AGB_TONS_AC": round(initial_agb, 4),
            })

            try:
                if harvest_cycles:
                    fvs_result = run_fvs_safe(
                        run_fvs_harvest_projection,
                        stand_df_base, tree_df, stand_id, variant,
                        harvest_cycles=harvest_cycles,
                        inv_year=inv_year,
                        config_version=config,
                        num_cycles=NUM_CYCLES,
                        cycle_length=CYCLE_LENGTH,
                    )
                else:
                    # No harvest: use standard projection
                    fvs_result = run_fvs_safe(
                        run_fvs_projection,
                        stand_df_base, tree_df, stand_id, variant,
                        config_version=config,
                        num_cycles=NUM_CYCLES,
                        cycle_length=CYCLE_LENGTH,
                    )

                ec = fvs_result.get("exit_code", 0)
                if ec not in (0, 10):
                    variant_fail_count += 1
                    logger.warning(
                        f"Plot {plot_id}/{variant}/{scenario}: exit code {ec}"
                    )
                else:
                    variant_fail_count = 0
                    for cycle_year, treelist in sorted(
                        fvs_result["treelists"].items()
                    ):
                        proj_year = cycle_year - inv_year
                        # Skip cycle-0 treelist — we already have the
                        # manual initial record from FIA tree data above
                        if proj_year <= 0:
                            continue
                        # Snap to nearest CYCLE_LENGTH-year bin so plots
                        # with different INVYRs align to the same time steps
                        proj_year_bin = round(proj_year / CYCLE_LENGTH) * CYCLE_LENGTH
                        agb = compute_plot_agb(treelist, nsbe)
                        results.append({
                            "PLOT": plot_id,
                            "FIRST_PLTCN": plt_cn,
                            "YEAR": cycle_year,
                            "PROJ_YEAR": proj_year,
                            "PROJ_YEAR_BIN": proj_year_bin,
                            "VARIANT": variant.upper(),
                            "CONFIG": config,
                            "SCENARIO": scenario,
                            "RCP": rcp,
                            "HARVEST": do_harvest,
                            "N_HARVESTS": len(harvest_cycles),
                            "HARVEST_CYCLES": ",".join(str(c) for c in harvest_cycles) if harvest_cycles else "",
                            "SI_MULT": round(si_mult, 4),
                            "SI_DELTA_FT": round(si_delta_ft, 2),
                            "AGB_TONS_AC": round(agb, 4),
                        })

            except Exception as e:
                variant_fail_count += 1
                logger.error(f"Plot {plot_id}/{variant}/{scenario}: {e}")

            # If both configs fail for this variant, skip remaining
            if variant_fail_count >= 2:
                logger.warning(
                    f"Plot {plot_id}/{variant}: both configs failed, "
                    f"skipping variant"
                )
                break

        # --- Posterior draws for this variant --------------------------------
        engine = (engines or {}).get(variant)
        if engine is not None and engine.draws_available and n_draws > 0:
            draw_limit = min(n_draws, engine.n_draws)
            for draw_idx in range(draw_limit):
                try:
                    draw_obj = engine.get_draw(int(draw_idx))
                    draw_kw = engine.generate_keywords_for_draw(
                        draw_obj,
                        (default_configs or {}).get(variant, None),
                        draw_idx=int(draw_idx),
                    )
                    # PATCH 28 April 2026: keyword renamed GROWMULT -> BAIMULT
                    # in config/uncertainty.py + config_loader.py to match FVS-NE
                    # spec. The old filter is now a no-op since BAIMULT is valid;
                    # left commented for reference.
                    # draw_kw = "\n".join(line for line in draw_kw.splitlines() if not line.startswith("GROWMULT"))
                except Exception as e:
                    logger.warning(
                        f"Plot {plot_id}/{variant}/draw{draw_idx} generate_keywords_for_draw: {e}"
                    )
                    continue
                scenario = f"{harvest_label}_{rcp_label}_draw{draw_idx:04d}"
                # Initial condition row (year 0)
                results.append({
                    "PLOT": plot_id,
                    "FIRST_PLTCN": plt_cn,
                    "YEAR": inv_year,
                    "PROJ_YEAR": 0,
                    "PROJ_YEAR_BIN": 0,
                    "VARIANT": variant.upper(),
                    "CONFIG": "draw",
                    "DRAW": int(draw_idx),
                    "SCENARIO": scenario,
                    "RCP": rcp,
                    "HARVEST": do_harvest,
                    "N_HARVESTS": len(harvest_cycles),
                    "HARVEST_CYCLES": ",".join(str(c) for c in harvest_cycles) if harvest_cycles else "",
                    "SI_MULT": round(si_mult, 4),
                    "SI_DELTA_FT": round(si_delta_ft, 2),
                    "AGB_TONS_AC": round(initial_agb, 4),
                })
                try:
                    if harvest_cycles:
                        fvs_result = run_fvs_safe(
                            run_fvs_harvest_projection,
                            stand_df_base, tree_df, stand_id, variant,
                            harvest_cycles=harvest_cycles,
                            inv_year=inv_year,
                            draw_keywords=draw_kw,
                            num_cycles=NUM_CYCLES,
                            cycle_length=CYCLE_LENGTH,
                        )
                    else:
                        # No-harvest draws route through KEYFILE_TEMPLATE
                        # (HARVEST_KEYFILE_TEMPLATE with draw_keywords returns
                        # FVS STOP 20; baseline uncertainty pipeline proves
                        # run_fvs_with_draw works for draw-only runs.)
                        fvs_result = run_fvs_safe(
                            run_fvs_with_draw,
                            stand_df_base, tree_df, stand_id, variant,
                            draw_keywords=draw_kw,
                            num_cycles=NUM_CYCLES,
                            cycle_length=CYCLE_LENGTH,
                        )
                    ec = fvs_result.get("exit_code", 0)
                    if ec not in (0, 10):
                        logger.warning(
                            f"Plot {plot_id}/{variant}/draw{draw_idx}: exit code {ec}"
                        )
                        continue
                    for cycle_year, treelist in sorted(
                        fvs_result["treelists"].items()
                    ):
                        proj_year = cycle_year - inv_year
                        if proj_year <= 0:
                            continue
                        proj_year_bin = round(proj_year / CYCLE_LENGTH) * CYCLE_LENGTH
                        agb = compute_plot_agb(treelist, nsbe)
                        results.append({
                            "PLOT": plot_id,
                            "FIRST_PLTCN": plt_cn,
                            "YEAR": cycle_year,
                            "PROJ_YEAR": proj_year,
                            "PROJ_YEAR_BIN": proj_year_bin,
                            "VARIANT": variant.upper(),
                            "CONFIG": "draw",
                            "DRAW": int(draw_idx),
                            "SCENARIO": scenario,
                            "RCP": rcp,
                            "HARVEST": do_harvest,
                            "N_HARVESTS": len(harvest_cycles),
                            "HARVEST_CYCLES": ",".join(str(c) for c in harvest_cycles) if harvest_cycles else "",
                            "SI_MULT": round(si_mult, 4),
                            "SI_DELTA_FT": round(si_delta_ft, 2),
                            "AGB_TONS_AC": round(agb, 4),
                        })
                except Exception as e:
                    logger.error(
                        f"Plot {plot_id}/{variant}/draw{draw_idx}: {e}"
                    )

    return results


# ===========================================================================
# Main
# ===========================================================================

def main():
    parser = argparse.ArgumentParser(
        description="PERSEUS harvest + climate factorial FVS projection"
    )
    parser.add_argument(
        "--rcp", type=str, default="all",
        choices=["4.5", "8.5", "none", "all"],
        help="Climate scenario (default: all = run none + 4.5 + 8.5)"
    )
    parser.add_argument(
        "--harvest", action="store_true", default=False,
        help="Run harvest scenario only"
    )
    parser.add_argument(
        "--no-harvest", action="store_true", default=False,
        help="Run no-harvest scenario only"
    )
    parser.add_argument("--batch-id", type=int, default=None)
    parser.add_argument("--batch-size", type=int, default=100)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--perseus-csv", type=str, default=PERSEUS_CSV)
    parser.add_argument("--si-ratio-csv", type=str, default=SI_RATIO_CSV)
    parser.add_argument(
        "--variants", type=str, nargs="+", default=["ne", "acd"],
        help="FVS variants to run (default: ne acd)"
    )
    parser.add_argument(
        "--configs", type=str, nargs="+", default=["default", "calibrated"],
        choices=["default", "calibrated"],
        help="Config versions to run (default: both default and calibrated)"
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--start-year", type=int, default=1999)
    parser.add_argument("--end-year", type=int, default=2004)
    parser.add_argument("--expansion-csv", type=str, default=None)
    parser.add_argument("--output-dir", type=str, default=OUTPUT_DIR)
    parser.add_argument("--n-draws", type=int, default=0,
                        help="Posterior draws per variant per plot (0 = off)")
    parser.add_argument("--seed-draws", type=int, default=42,
                        help="Seed for UncertaintyEngine")
    parser.add_argument("--exclude-plots", type=str, default="",
                        help="Comma-separated plot IDs to skip, or path to file")
    args = parser.parse_args()

    # Determine which harvest/climate combos to run
    if args.harvest and args.no_harvest:
        harvest_modes = [True, False]
    elif args.harvest:
        harvest_modes = [True]
    elif args.no_harvest:
        harvest_modes = [False]
    else:
        harvest_modes = [True, False]  # default: both

    if args.rcp == "all":
        rcp_scenarios = ["none", "4.5", "8.5"]
    else:
        rcp_scenarios = [args.rcp]

    scenarios = [
        (h, r) for h in harvest_modes for r in rcp_scenarios
    ]
    logger.info(
        f"Running {len(scenarios)} harvest x climate combos: "
        + ", ".join(
            f"{'harvest' if h else 'noharvest'}_rcp{r.replace('.','')}"
            if r != 'none' else f"{'harvest' if h else 'noharvest'}_noclimate"
            for h, r in scenarios
        )
    )
    logger.info(
        f"Variants: {args.variants}, Configs: {args.configs} "
        f"({len(args.variants) * len(args.configs)} FVS runs per plot per scenario)"
    )

    # Load PERSEUS plots (FULL set needed for harvest scheduling)
    perseus_full = pd.read_csv(args.perseus_csv)
    logger.info(f"Loaded {len(perseus_full)} total PERSEUS plots")

    mask = (
        (perseus_full["FIRST_INVYR"] >= args.start_year)
        & (perseus_full["FIRST_INVYR"] <= args.end_year)
    )
    perseus_full = perseus_full.loc[mask].copy()
    n_total = len(perseus_full)
    logger.info(
        f"Filtered to {n_total} plots with FIRST_INVYR in "
        f"[{args.start_year}, {args.end_year}]"
    )

    # Apply plot exclusion filter (e.g., for pathological plots that hang FVS)
    if args.exclude_plots:
        _raw = args.exclude_plots
        if os.path.exists(_raw):
            with open(_raw) as f:
                _excludes = {line.strip() for line in f if line.strip()}
        else:
            _excludes = {s.strip() for s in _raw.split(",") if s.strip()}
        if _excludes:
            # Normalize both sides: accept "2869", "2869.0", "2869.00" etc.
            def _norm(x):
                try:
                    return str(int(float(x)))
                except Exception:
                    return str(x).strip()
            _norm_excludes = {_norm(x) for x in _excludes}
            _before = len(perseus_full)
            perseus_full = perseus_full.loc[
                ~perseus_full["PLOT"].apply(_norm).isin(_norm_excludes)
            ].copy()
            logger.info(
                f"Excluded {_before - len(perseus_full)} plot(s): "
                f"{sorted(_excludes)}"
            )

    # Load FIA tree data for all plots
    plt_cns = [
        str(int(float(x)))
        for x in perseus_full["FIRST_PLTCN"].dropna().unique()
    ]
    logger.info(f"Loading FIA trees for {len(plt_cns)} unique PLT_CNs...")
    fia_trees = load_fia_trees_for_plots(plt_cns, FIA_DIR)
    logger.info(f"Loaded {len(fia_trees)} tree records")

    # Load COND data for harvest eligibility
    logger.info("Loading FIA COND data for harvest eligibility...")
    cond_df = load_fia_cond_for_plots(plt_cns, FIA_DIR)

    # Determine harvest eligibility (uses FULL dataset for consistent schedule)
    eligibility = determine_harvest_eligibility(perseus_full, cond_df)

    # Build harvest schedule (deterministic, same for all batches)
    harvest_schedule = build_harvest_schedule(
        eligibility, n_total_plots=n_total, seed=args.seed
    )

    # Load SI ratios for climate scenarios
    si_ratios = None
    if any(r != "none" for _, r in scenarios):
        si_ratios = load_si_ratios(args.si_ratio_csv)

    # Slice to this batch
    if args.batch_id is not None:
        start = (args.batch_id - 1) * args.batch_size
        end = min(start + args.batch_size, n_total)
        perseus_batch = perseus_full.iloc[start:end]
        logger.info(f"Batch {args.batch_id}: plots {start+1} to {end}")
    elif args.all:
        perseus_batch = perseus_full
    else:
        logger.error("Specify --batch-id or --all")
        sys.exit(1)

    if perseus_batch.empty:
        logger.info("No plots in this batch. Exiting.")
        return

    # Load expansion factors if provided
    expns_df = None
    if args.expansion_csv and os.path.exists(args.expansion_csv):
        expns_df = pd.read_csv(args.expansion_csv)
        logger.info(f"Loaded expansion factors for {len(expns_df)} plots")

    # Initialize NSBE calculator
    nsbe = NSBECalculator(NSBE_ROOT)
    # Initialize uncertainty engines + default configs once (if --n-draws > 0)
    engines = {}
    default_configs = {}
    if args.n_draws > 0:
        for _v in args.variants:
            try:
                _eng = UncertaintyEngine(
                    _v, config_dir=CONFIG_DIR, seed=args.seed_draws
                )
                if _eng.draws_available:
                    engines[_v] = _eng
                    logger.info(
                        f"{_v.upper()} has {_eng.n_draws} posterior draws available"
                    )
                else:
                    logger.warning(f"{_v.upper()} no draws available")
            except Exception as _e:
                logger.warning(f"{_v} UncertaintyEngine init failed: {_e}")
            # Load default config (needed for generate_keywords_for_draw)
            try:
                from config.config_loader import FvsConfigLoader
                _loader = FvsConfigLoader(
                    _v, version="default", config_dir=CONFIG_DIR
                )
                default_configs[_v] = _loader.config
            except Exception as _e:
                logger.warning(f"{_v} default config load failed: {_e}")


    # Process each scenario
    for do_harvest, rcp in scenarios:
        rcp_tag = f"rcp{rcp.replace('.', '')}" if rcp != "none" else "noclimate"
        harvest_tag = "harvest" if do_harvest else "noharvest"
        scenario_tag = f"{harvest_tag}_{rcp_tag}"

        scenario_output_dir = os.path.join(args.output_dir, scenario_tag)
        os.makedirs(scenario_output_dir, exist_ok=True)

        logger.info(f"\n{'='*70}")
        logger.info(f"SCENARIO: {scenario_tag}")
        logger.info(f"{'='*70}")

        all_results = []
        t0 = time.time()

        for idx, (_, plot_row) in enumerate(perseus_batch.iterrows()):
            try:
                results = process_plot_harvest(
                    plot_row.to_dict(),
                    fia_trees,
                    nsbe,
                    si_ratios,
                    rcp,
                    harvest_schedule,
                    do_harvest=do_harvest,
                    cond_df=cond_df,
                    variants=args.variants,
                    configs=args.configs,
                    seed=args.seed,
                    engines=engines,
                    n_draws=args.n_draws,
                    default_configs=default_configs,
                )
                all_results.extend(results)
            except Exception as exc:
                logger.warning(
                    f"Plot {plot_row.get('PLOT', '?')} failed: {exc}"
                )

            if (idx + 1) % 25 == 0:
                elapsed = time.time() - t0
                rate = (idx + 1) / elapsed * 60
                logger.info(
                    f"  [{scenario_tag}] {idx+1}/{len(perseus_batch)} plots "
                    f"({rate:.1f} plots/min, {len(all_results)} records)"
                )

        # Save results
        batch_suffix = f"_batch{args.batch_id}" if args.batch_id else ""

        result_df = pd.DataFrame(all_results)
        if not result_df.empty:
            out_path = os.path.join(
                scenario_output_dir,
                f"perseus_{scenario_tag}{batch_suffix}.csv"
            )
            result_df.to_csv(out_path, index=False)
            logger.info(
                f"  {scenario_tag}: {len(result_df)} records -> {out_path}"
            )

            # Quick aggregation to state level
            point_mmt, _ = aggregate_with_uncertainty(
                result_df, pd.DataFrame(),
                expns_df=expns_df,
                n_plots_total=n_total,
            )

            if not point_mmt.empty:
                mmt_path = os.path.join(
                    scenario_output_dir,
                    f"perseus_{scenario_tag}_mmt{batch_suffix}.csv"
                )
                point_mmt.to_csv(mmt_path, index=False)

                print(f"\n{'='*70}")
                print(f"Scenario: {scenario_tag} — State-Level AGC (MMT)")
                print(f"{'='*70}")
                for _, row in point_mmt.sort_values("PROJ_YEAR").iterrows():
                    yr = int(row["PROJ_YEAR"])
                    mmt = row["MMT"]
                    n = int(row.get("N_PLOTS", 0))
                    print(f"  Year +{yr:>3d}: {mmt:>8.2f} MMT ({n} plots)")

        elapsed = time.time() - t0
        logger.info(f"  {scenario_tag} complete in {elapsed:.0f}s")

    # Save harvest schedule metadata
    schedule_meta = {
        "description": "PERSEUS Round 2 harvest schedule",
        "harvest_rate_annual": HARVEST_ANNUAL_RATE,
        "harvest_rate_per_cycle": HARVEST_CYCLE_RATE,
        "ba_removal_proportion": BA_REMOVAL_PROPORTION,
        "min_reentry_years": MIN_REENTRY_CYCLES * CYCLE_LENGTH,
        "n_total_plots": n_total,
        "eligibility_criteria": {
            "ownership": "Private (OWNGRPCD 40-46)",
            "reserve": "Not reserved (RESERVCD = 0)",
            "condition": "Accessible forest (COND_STATUS_CD = 1)",
            "volume": "BALIVE > 0",
        },
        "n_eligible": int(eligibility["ELIGIBLE"].sum()),
        "pct_eligible": round(
            100 * eligibility["ELIGIBLE"].sum() / n_total, 1
        ),
        "schedule_summary": {
            str(c): len(s) for c, s in harvest_schedule.items()
        },
        "total_harvest_events": sum(
            len(s) for s in harvest_schedule.values()
        ),
        "unique_plots_harvested": len(
            set().union(*harvest_schedule.values())
        ),
        "seed": args.seed,
        "scenarios_run": [
            f"{'harvest' if h else 'noharvest'}_rcp{r.replace('.','')}"
            if r != "none"
            else f"{'harvest' if h else 'noharvest'}_noclimate"
            for h, r in scenarios
        ],
    }

    meta_path = os.path.join(
        args.output_dir,
        f"perseus_harvest_meta{('_batch' + str(args.batch_id)) if args.batch_id else ''}.json"
    )
    with open(meta_path, "w") as f:
        json.dump(schedule_meta, f, indent=2, default=str)
    logger.info(f"Harvest metadata: {meta_path}")


if __name__ == "__main__":
    main()
