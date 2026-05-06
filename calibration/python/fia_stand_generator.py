"""
FIA-based stand generator for Bakuzis matrix scenarios.

Replaces the synthetic stand generator in bakuzis_100yr_comparison.py
with a stratified random sample of real FIA condition pairs. For each
(variant, site_class, density_class) cell, returns one or more real
stands whose initial conditions match the target SI band and BA band,
drawn from FIA plots within the variant's state coverage.

This is the Option A path from manuscript/C2_FIA_GENERATOR_OUTLINE.md.
The stand and tree dataframes returned have the same schema the
bakuzis runner already uses (fvs_standinit and fvs_treeinit columns),
so the runner only needs to swap which generator it calls.

Required inputs (Cardinal paths shown; override with env vars):
  FIA_DATA_DIR   per-state CSVs: STATE_PLOT.csv, STATE_COND.csv,
                  STATE_TREE.csv. Default: ~/fia_data
  CONUS_RDS      consolidated remeasurement pairs (optional; if
                  available, overrides per-state CSV reads)

Variant -> state coverage map. Conservative subset of states with
substantial forested area in each FVS variant's geographic domain.
The runner can pass `n_plots=5` to draw five stratified samples per
cell so each (species_group, site, density) gets replicated coverage.

Author: A. Weiskittel
Date: 2026-04-25
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# Variant -> tuple of FIPS state codes covered. Conservative;
# extend as needed for additional regional fidelity.
VARIANT_STATES = {
    "ne":  (9, 23, 25, 33, 36, 44, 50),    # CT ME MA NH NY RI VT
    "acd": (23, 33, 50),                    # ME NH VT (Acadian core)
    "pn":  (41, 53),                        # OR WA
    "sn":  (1, 12, 13, 28, 45, 47),         # AL FL GA MS SC TN
    "ie":  (16, 30, 53),                    # ID MT WA Inland Empire portion
}

# Site index bins (in feet at base age 50, FIA SICOND).
SI_BINS = {
    "Low":    (30, 50),
    "Medium": (50, 70),
    "High":   (70, 100),
}

# Basal area bins (ft2/ac, FIA BALIVE).
BA_BINS = {
    "Low":    (50, 90),
    "Medium": (90, 150),
    "High":   (150, 240),
}


def _state_csv(prefix: str, table: str, fia_dir: Path) -> Optional[Path]:
    """Return the path to a per-state FIA CSV if it exists."""
    p = fia_dir / f"{prefix}_{table}.csv"
    return p if p.exists() else None


def _state_abbrev(state_code: int) -> str:
    """Map FIPS state code to USPS abbreviation. Subset matching VARIANT_STATES."""
    table = {
        1: "AL", 9: "CT", 12: "FL", 13: "GA", 16: "ID",
        23: "ME", 25: "MA", 28: "MS", 30: "MT", 33: "NH",
        36: "NY", 41: "OR", 44: "RI", 45: "SC", 47: "TN",
        50: "VT", 53: "WA",
    }
    return table.get(state_code, "")


def load_fia_conditions(
    variant: str,
    fia_dir: Path,
) -> pd.DataFrame:
    """Load and concatenate per-state COND tables for one variant.

    Returns columns: STATECD, COUNTYCD, PLOT, INVYR, COND_STATUS_CD,
    SICOND, SIBASE, BALIVE, FORTYPCD, STDAGE, ASPECT, SLOPE, ELEV,
    PLT_CN. Missing files are skipped.
    """
    states = VARIANT_STATES.get(variant.lower(), ())
    frames = []
    for code in states:
        ab = _state_abbrev(code)
        cond_csv = _state_csv(ab, "COND", fia_dir)
        if cond_csv is None:
            continue
        df = pd.read_csv(
            cond_csv,
            usecols=lambda c: c
            in (
                "STATECD", "COUNTYCD", "PLOT", "INVYR",
                "COND_STATUS_CD", "SICOND", "SIBASE", "BALIVE",
                "FORTYPCD", "STDAGE", "ASPECT", "SLOPE", "ELEV", "PLT_CN",
            ),
            low_memory=False,
        )
        # Filter to forested conditions only
        if "COND_STATUS_CD" in df.columns:
            df = df.loc[df["COND_STATUS_CD"] == 1].copy()
        frames.append(df)
    if not frames:
        raise FileNotFoundError(
            f"No FIA COND CSVs found for variant {variant} in {fia_dir}. "
            "Convert from the consolidated RDS or download per-state CSVs."
        )
    return pd.concat(frames, ignore_index=True)


def load_fia_trees(
    variant: str,
    plt_cns: list[int],
    fia_dir: Path,
) -> pd.DataFrame:
    """Load TREE records for the given PLT_CN list, from each variant state."""
    states = VARIANT_STATES.get(variant.lower(), ())
    frames = []
    for code in states:
        ab = _state_abbrev(code)
        tree_csv = _state_csv(ab, "TREE", fia_dir)
        if tree_csv is None:
            continue
        df = pd.read_csv(
            tree_csv,
            usecols=lambda c: c
            in (
                "PLT_CN", "STATECD", "COUNTYCD", "PLOT", "SUBP", "TREE",
                "STATUSCD", "SPCD", "DIA", "HT", "CR", "TPA_UNADJ",
                "DAMTYP1", "DECAYCD", "STANDING_DEAD_CD",
            ),
            low_memory=False,
        )
        df = df.loc[df["PLT_CN"].isin(plt_cns)].copy()
        if "STATUSCD" in df.columns:
            # Live trees only
            df = df.loc[df["STATUSCD"] == 1].copy()
        if not df.empty:
            frames.append(df)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True)


def stratified_sample_conditions(
    cond_df: pd.DataFrame,
    site_class: str,
    density_class: str,
    n_plots: int = 5,
    seed: int = 42,
) -> pd.DataFrame:
    """Filter conditions to the target SI band and BA band; sample n_plots."""
    si_lo, si_hi = SI_BINS[site_class]
    ba_lo, ba_hi = BA_BINS[density_class]

    df = cond_df.copy()
    if "SICOND" in df.columns:
        df = df.loc[df["SICOND"].between(si_lo, si_hi)].copy()
    if "BALIVE" in df.columns:
        df = df.loc[df["BALIVE"].between(ba_lo, ba_hi)].copy()

    # Conservative requirements: forested, has site index, has BA
    df = df.dropna(subset=["SICOND", "BALIVE"])

    if len(df) == 0:
        return df
    if len(df) <= n_plots:
        return df.copy()
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(df), size=n_plots, replace=False)
    return df.iloc[idx].reset_index(drop=True)


def build_stand_init(
    cond_row: pd.Series,
    stand_id: str,
    variant: str,
) -> pd.DataFrame:
    """Convert one COND row to an fvs_standinit single-row DataFrame.

    Schema matches the existing synthetic generator so the bakuzis
    runner can swap generators without changing downstream code.
    """
    return pd.DataFrame([{
        "stand_id": stand_id,
        "variant": variant.upper(),
        "inv_year": int(cond_row.get("INVYR", 2000)),
        "latitude": 0.0,            # FIA coordinates are perturbed; leave 0
        "longitude": 0.0,
        "region": _variant_region(variant),
        "forest": 0,
        "district": 0,
        "basal_area_factor": 0.0,
        "inv_plot_size": 6.0,
        "brk_dbh": 5.0,
        "num_plots": 1,
        "age": int(cond_row.get("STDAGE", 40) or 40),
        "aspect": float(cond_row.get("ASPECT", 0) or 0),
        "slope": float(cond_row.get("SLOPE", 15) or 15),
        "elevft": float(cond_row.get("ELEV", 1000) or 1000),
        "site_species": _variant_default_site_species(variant),
        "site_index": float(cond_row.get("SICOND", 60) or 60),
        "state": int(cond_row.get("STATECD", 0) or 0),
        "county": int(cond_row.get("COUNTYCD", 0) or 0),
        "forest_type": int(cond_row.get("FORTYPCD", 0) or 0),
        "sam_wt": 1.0,
    }])


def build_tree_init(
    tree_df: pd.DataFrame,
    stand_id: str,
) -> pd.DataFrame:
    """Convert TREE records to fvs_treeinit format."""
    if tree_df.empty:
        return pd.DataFrame(columns=[
            "stand_id", "plot_id", "tree_id", "tree_count",
            "species", "diameter", "ht", "crratio",
        ])
    out = pd.DataFrame({
        "stand_id":   stand_id,
        "plot_id":    1,
        "tree_id":    range(1, len(tree_df) + 1),
        "tree_count": tree_df["TPA_UNADJ"].fillna(6.0).round(2),
        "species":    tree_df["SPCD"].astype(int),
        "diameter":   tree_df["DIA"].fillna(5.0).round(1),
        "ht":         tree_df["HT"].fillna(40).round(0),
        "crratio":    (tree_df["CR"] * 10 if tree_df["CR"].max() <= 9 else tree_df["CR"])
                          .fillna(50).clip(15, 99).astype(int),
    })
    return out


def _variant_region(variant: str) -> int:
    table = {"ne": 9, "acd": 9, "pn": 6, "sn": 8, "ie": 1}
    return table.get(variant.lower(), 9)


def _variant_default_site_species(variant: str) -> int:
    """Default FIA SPCD used as the site index reference species."""
    table = {
        "ne":  12,    # Balsam fir
        "acd": 12,    # Balsam fir
        "pn":  202,   # Douglas-fir
        "sn":  131,   # Loblolly pine
        "ie":  202,   # Douglas-fir
    }
    return table.get(variant.lower(), 12)


def generate_real_stand(
    variant: str,
    site_class: str,
    density_class: str,
    fia_dir: Optional[Path] = None,
    n_plots: int = 5,
    seed: int = 42,
) -> list[tuple[pd.DataFrame, pd.DataFrame, pd.Series]]:
    """Return up to n_plots real (stand_df, tree_df, cond_row) tuples for a cell.

    Each tuple is one stand. Empty list if no plots match the bin.
    Raises FileNotFoundError if the FIA per-state CSVs for the variant
    are not available (e.g., trying PN without OR/WA CSVs).
    """
    fia_dir = fia_dir or Path(os.environ.get("FIA_DATA_DIR", str(Path.home() / "fia_data")))
    cond = load_fia_conditions(variant, fia_dir)
    sample = stratified_sample_conditions(
        cond, site_class, density_class, n_plots=n_plots, seed=seed,
    )
    if sample.empty:
        return []

    plt_cns = sample["PLT_CN"].astype(int).tolist()
    trees = load_fia_trees(variant, plt_cns, fia_dir)

    out = []
    for i, row in sample.iterrows():
        stand_id = f"FIA{int(row['PLT_CN']) % 1_000_000:06d}"
        this_trees = trees.loc[trees["PLT_CN"] == row["PLT_CN"]].copy()
        stand_df = build_stand_init(row, stand_id, variant)
        tree_df = build_tree_init(this_trees, stand_id)
        out.append((stand_df, tree_df, row))
    return out


__all__ = [
    "VARIANT_STATES",
    "SI_BINS",
    "BA_BINS",
    "generate_real_stand",
    "load_fia_conditions",
    "load_fia_trees",
    "build_stand_init",
    "build_tree_init",
]
