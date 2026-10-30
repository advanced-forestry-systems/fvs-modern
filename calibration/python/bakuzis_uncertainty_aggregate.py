#!/usr/bin/env python3
"""
Bakuzis Uncertainty Aggregation and Benchmarking
=================================================

Consumes the raw ensemble CSV from bakuzis_uncertainty_comparison.py and
produces three derived products:

  1. Long ensemble summary
     (scenario, variant, config, year, variable, median, q025, q975, mean, sd)
     config is default, calibrated_map, or posterior_band.

  2. Wide benchmarking table: one row per scenario with BA, TPA, and
     volume summaries at year 50 and year 100, band width, and the
     default vs calibrated divergence expressed as percent of default.

  3. Bakuzis biological law compliance at year 100:
        - Sukachev effect (competition reduces per tree growth)
        - Eichhorn's rule (better site increases stand volume)
        - Crown recession under density
        - Mortality size pattern (small trees heavier mortality)

The four laws are evaluated for each variant and config so you can
compare how the default, calibrated MAP, and posterior band perform.

Usage:
  python bakuzis_uncertainty_aggregate.py \
      --input calibration/output/bakuzis/bakuzis_uncertainty_ne_n50.csv
  python bakuzis_uncertainty_aggregate.py --input-dir calibration/output/bakuzis
  python bakuzis_uncertainty_aggregate.py --input-dir ... --quantiles 0.025,0.10,0.5,0.90,0.975

Author: A. Weiskittel
Date: 2026-04-23
"""

from __future__ import annotations

import argparse
import glob
import logging
import os
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# FVS summary variables of interest for the Bakuzis benchmark
CORE_VARS = ["tpa", "atba", "tcuft", "mcuft", "bdft", "attopht", "mort", "age"]
DEFAULT_QUANTILES = [0.025, 0.25, 0.5, 0.75, 0.975]


# ---------------------------------------------------------------------------
# Load helpers
# ---------------------------------------------------------------------------

def load_ensemble(paths: Iterable[str]) -> pd.DataFrame:
    """Read one or more ensemble CSVs and return a single concatenated frame."""
    frames = []
    for p in paths:
        if not os.path.exists(p):
            logger.warning(f"Missing file: {p}")
            continue
        df = pd.read_csv(p)
        frames.append(df)
    if not frames:
        return pd.DataFrame()
    df = pd.concat(frames, ignore_index=True)
    # Keep only variables we understand. The "replicate" column is
    # present when the runner used FIA mode (one row per real plot per
    # config-draw); it is absent in synthetic mode. Either way we
    # preserve it so summarize_ensemble can pool across replicates.
    vars_present = [v for v in CORE_VARS if v in df.columns]
    keep = [
        "year", "scenario", "species_group", "site_class", "density_class",
        "variant", "config", "draw_id", "replicate",
    ] + vars_present
    keep = [c for c in keep if c in df.columns]
    return df[keep].copy()


def resolve_inputs(args) -> list[str]:
    if args.input:
        return [args.input]
    if args.input_dir:
        return sorted(glob.glob(os.path.join(args.input_dir, "bakuzis_uncertainty_*.csv")))
    return []


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def _add_horizon(df: pd.DataFrame) -> pd.DataFrame:
    """Years since each run's own start, on the 5-year grid, as 2000 + horizon
    (the same alignment summarize_ensemble applies)."""
    df = df.copy()
    keys = ["variant", "scenario", "replicate", "config", "draw_id"]
    h = df["year"] - df.groupby(keys, dropna=False)["year"].transform("min")
    df["horizon_yr"] = (h / 5).round().astype(int) * 5
    df["year"] = 2000 + df["horizon_yr"]
    return df


def summarize_stand_bands(df: pd.DataFrame, lo: float = 0.025, hi: float = 0.975):
    """Parameter-only (within-stand) bands, kept separate from stand spread.

    Returns (stand_bands, cell_bands).

    stand_bands: one row per variant x cell x FIA stand (replicate) x year x
    variable, with the 2.5/50/97.5 percent quantiles of that stand's posterior
    draws (parameter uncertainty only), its default and calibrated MAP values,
    and band_pct = 100 * (q975 - q025) / q500.

    cell_bands: one row per cell x year x variable with
      within_*      the within-stand band summarised across the cell's stands
                    (median of the stand quantiles and of band_pct; min and max
                    of band_pct)
      map_*         between-stand spread of the stands' calibrated MAP values
                    (median, min, max, range_pct = 100 * (max - min) / median, sd)
      default_median, n_stands
    The two are never pooled.
    """
    if "replicate" not in df.columns:
        df = df.assign(replicate=1)
    vars_present = [v for v in CORE_VARS if v in df.columns]
    d = _add_horizon(df)
    sk = ["variant", "scenario", "species_group", "site_class", "density_class", "replicate", "year"]
    rows = []
    for var in vars_present:
        post = d[d["config"] == "posterior"].groupby(sk, dropna=False)[var]
        b = pd.DataFrame({
            "n_draws": post.count(),
            "post_q025": post.quantile(lo),
            "post_q500": post.median(),
            "post_q975": post.quantile(hi),
        })
        for cfg, col in (("calibrated_map", "map"), ("default", "default")):
            b[col] = d[d["config"] == cfg].groupby(sk, dropna=False)[var].median()
        b = b.reset_index()
        b["variable"] = var
        rows.append(b)
    stand = pd.concat(rows, ignore_index=True)
    stand["band_width"] = stand["post_q975"] - stand["post_q025"]
    stand["band_pct"] = 100 * stand["band_width"] / stand["post_q500"].replace(0, np.nan)

    ck = ["variant", "scenario", "species_group", "site_class", "density_class", "year", "variable"]
    g = stand.groupby(ck, dropna=False)
    cell = pd.DataFrame({
        "n_stands": g["replicate"].nunique(),
        "within_q025_med": g["post_q025"].median(),
        "within_q500_med": g["post_q500"].median(),
        "within_q975_med": g["post_q975"].median(),
        "within_band_pct_med": g["band_pct"].median(),
        "within_band_pct_min": g["band_pct"].min(),
        "within_band_pct_max": g["band_pct"].max(),
        "map_median": g["map"].median(),
        "map_min": g["map"].min(),
        "map_max": g["map"].max(),
        "map_sd": g["map"].std(),
        "default_median": g["default"].median(),
    }).reset_index()
    cell["map_range_pct"] = 100 * (cell["map_max"] - cell["map_min"]) / cell["map_median"].replace(0, np.nan)
    return stand, cell


def summarize_ensemble(df: pd.DataFrame, quantiles: list[float]) -> pd.DataFrame:
    """Collapse posterior draws (and FIA replicates) into quantile bands.

    Default and calibrated_map rows are aggregated across FIA replicates
    so each cell yields a single median (and SD across replicates,
    capturing within-cell stand-to-stand variability). Posterior rows
    are summarized across the cross product of draw_id and replicate
    so the credible band reflects both parameter and stand uncertainty.

    In synthetic-stand mode (no replicate column) behavior is unchanged
    because there is exactly one stand per cell.
    """
    vars_present = [v for v in CORE_VARS if v in df.columns]
    if not vars_present:
        return pd.DataFrame()

    has_replicate = "replicate" in df.columns
    n_replicates = (
        int(df["replicate"].nunique()) if has_replicate else 1
    )

    # When replicate is present, FIA-sampled plots may carry different
    # FIA INVYRs, so the calendar `year` column is misaligned across
    # replicates. Compute a horizon (years since each replicate's own
    # start) and aggregate on that, then convert horizon back to a
    # canonical calendar year (2000 + horizon) so downstream figures
    # and benchmark tables continue to work without modification.
    df = df.copy()
    if has_replicate:
        plot_keys = ["variant", "scenario", "replicate", "config", "draw_id"]
        df["horizon_yr"] = (
            df["year"]
            - df.groupby(plot_keys, dropna=False)["year"].transform("min")
        )
        # Quantize horizon to the nearest 5-year cycle so small
        # within-cycle differences across plots merge cleanly.
        df["horizon_yr"] = (df["horizon_yr"] / 5).round().astype(int) * 5
        df["year"] = 2000 + df["horizon_yr"]

    group_cols = [
        "year", "scenario", "species_group", "site_class", "density_class",
        "variant",
    ]

    # Default and calibrated MAP: aggregate across replicates so each
    # cell-year gets a single point estimate plus an SD across plots.
    point_mask = df["config"].isin(["default", "calibrated_map"])
    point = df.loc[point_mask].copy()
    pieces_point = []
    for var in vars_present:
        for cfg, sub in point.groupby("config"):
            tmp = sub.groupby(group_cols, dropna=False, as_index=False).agg(
                median=(var, lambda s: float(np.nanmedian(s))),
                mean=(var, "mean"),
                sd=(var, "std"),
                n_draws=(var, "count"),
            )
            for q in quantiles:
                tmp[f"q{int(round(q * 1000)):03d}"] = (
                    sub.groupby(group_cols, dropna=False)[var]
                    .quantile(q)
                    .reset_index(drop=True)
                    .values
                )
            tmp["variable"] = var
            tmp["config"] = cfg
            # Preserve old contract: SD of point configs is across
            # replicates only (not across draws). Replace NaN with 0
            # for synthetic mode where n=1.
            tmp["sd"] = tmp["sd"].fillna(0.0)
            pieces_point.append(tmp)
    point_long = (
        pd.concat(pieces_point, ignore_index=True)
        if pieces_point
        else pd.DataFrame()
    )

    # Posterior band: pool across draws AND replicates so the band
    # captures both parameter and stand-to-stand uncertainty.
    post_mask = df["config"] == "posterior"
    post = df.loc[post_mask].copy()
    if post.empty:
        return point_long

    pieces = []
    for var in vars_present:
        tmp = post.groupby(group_cols, dropna=False, as_index=False).agg(
            median=(var, lambda s: float(np.nanmedian(s))),
            mean=(var, "mean"),
            sd=(var, "std"),
            n_draws=(var, "count"),
        )
        for q in quantiles:
            tmp[f"q{int(round(q * 1000)):03d}"] = (
                post.groupby(group_cols, dropna=False)[var]
                .quantile(q)
                .reset_index(drop=True)
                .values
            )
        tmp["variable"] = var
        tmp["config"] = "posterior_band"
        # Kept for backward compatibility. This band pools the 100 draws of
        # every FIA stand in the cell, so it mixes parameter uncertainty with
        # stand-to-stand differences. For the parameter-only band use
        # bakuzis_uncertainty_stand_bands.csv / _cell_bands.csv.
        tmp["band_basis"] = "pooled_draws_x_stands"
        pieces.append(tmp)

    post_long = pd.concat(pieces, ignore_index=True)
    out = pd.concat([point_long, post_long], ignore_index=True, sort=False)
    out.attrs["n_replicates"] = n_replicates
    return out


def build_benchmark_table(summary: pd.DataFrame) -> pd.DataFrame:
    """Wide per scenario table at year 50 and year 100.

    Columns like: ba_default_y050, ba_calmap_y050, ba_postmed_y050,
                  ba_post_q025_y050, ba_post_q975_y050, ba_diff_pct_y050.
    """
    wanted_years = [2025 + 25, 2025 + 50, 2025 + 75, 2025 + 100]  # 2050, 2075, 2100, 2125
    # Actual inv_year is 2000 in the synthetic generator, so projection years:
    # 2000, 2050, 2075, 2100 are year 0, 50, 75, 100 of projection.
    # Emit both calendar year and horizon so the user can pivot either way.
    # Horizon from first year of each scenario:
    summary = summary.copy()
    first_year_per_scenario = (
        summary.groupby(["variant", "scenario"])["year"].transform("min")
    )
    summary["horizon_yr"] = summary["year"] - first_year_per_scenario

    rows = []
    for (var, scen, variant, sp, site, dens), group in summary.groupby(
        ["variable", "scenario", "variant", "species_group", "site_class", "density_class"],
        dropna=False,
    ):
        row = {
            "variable": var,
            "scenario": scen,
            "variant": variant,
            "species_group": sp,
            "site_class": site,
            "density_class": dens,
        }
        for horizon in [25, 50, 75, 100]:
            sub = group.loc[group["horizon_yr"] == horizon]
            if sub.empty:
                continue
            for config in ["default", "calibrated_map", "posterior_band"]:
                rec = sub.loc[sub["config"] == config]
                if rec.empty:
                    continue
                tag = {
                    "default": "def",
                    "calibrated_map": "cal",
                    "posterior_band": "post",
                }[config]
                row[f"{tag}_median_y{horizon:03d}"] = float(rec["median"].iloc[0])
                if config == "posterior_band":
                    row[f"{tag}_q025_y{horizon:03d}"] = float(rec["q025"].iloc[0])
                    row[f"{tag}_q975_y{horizon:03d}"] = float(rec["q975"].iloc[0])
                    row[f"{tag}_band_y{horizon:03d}"] = float(
                        rec["q975"].iloc[0] - rec["q025"].iloc[0]
                    )
            # Percent divergence: (cal median minus def) / def x 100
            if (
                f"def_median_y{horizon:03d}" in row
                and f"cal_median_y{horizon:03d}" in row
                and row[f"def_median_y{horizon:03d}"] not in (0, None)
                and not np.isnan(row[f"def_median_y{horizon:03d}"])
            ):
                row[f"diff_pct_y{horizon:03d}"] = (
                    100.0
                    * (row[f"cal_median_y{horizon:03d}"] - row[f"def_median_y{horizon:03d}"])
                    / row[f"def_median_y{horizon:03d}"]
                )
        rows.append(row)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Bakuzis biological law compliance
# ---------------------------------------------------------------------------

def evaluate_bakuzis_laws(summary: pd.DataFrame) -> pd.DataFrame:
    """Evaluate the four Bakuzis biological laws at the 100 yr horizon.

    Returns a long DataFrame indexed by (variant, config, law) with a
    Boolean 'satisfied' column and a descriptive 'evidence' string.
    """
    summary = summary.copy()
    summary["horizon_yr"] = summary.groupby(["variant", "scenario"])["year"].transform(
        lambda s: s - s.min()
    )
    terminal = summary.loc[summary["horizon_yr"] == 100].copy()
    if terminal.empty:
        # fallback: last available year per scenario
        last_year = summary.groupby(["variant", "scenario"])["horizon_yr"].transform("max")
        terminal = summary.loc[summary["horizon_yr"] == last_year].copy()

    results = []
    for (variant, config), group in terminal.groupby(["variant", "config"]):
        # Pivot to wide: one row per scenario, columns per variable
        wide = group.pivot_table(
            index=["scenario", "species_group", "site_class", "density_class"],
            columns="variable",
            values="median",
        ).reset_index()

        # Law 1: Sukachev effect. At the same site and species, per tree
        # stand volume (tcuft / tpa) should be lower in the high density
        # class than in the low density class, because competition
        # compresses individual growth.
        if set(["tpa", "tcuft"]).issubset(wide.columns):
            wide["vol_per_tree"] = wide["tcuft"] / wide["tpa"].replace(0, np.nan)
            pivots = []
            for (sp, site), g in wide.groupby(["species_group", "site_class"]):
                gd = g.set_index("density_class")["vol_per_tree"]
                if "Low" in gd.index and "High" in gd.index:
                    pivots.append(gd["High"] < gd["Low"])
            sukachev = (
                float(np.mean(pivots)) if pivots else float("nan")
            )
        else:
            sukachev = float("nan")

        # Law 2: Eichhorn's rule. At the same species and density, BA
        # should increase with site class (Low < Medium < High).
        if "atba" in wide.columns:
            hits = []
            for (sp, dens), g in wide.groupby(["species_group", "density_class"]):
                gd = g.set_index("site_class")["atba"]
                if {"Low", "Medium", "High"}.issubset(gd.index):
                    hits.append(gd["Low"] <= gd["Medium"] <= gd["High"])
            eichhorn = float(np.mean(hits)) if hits else float("nan")
        else:
            eichhorn = float("nan")

        # Law 3: Crown recession. We don't have CR in the summary; use
        # mortality as a competition proxy: higher density should yield
        # greater cumulative mortality.
        if "mort" in wide.columns:
            hits = []
            for (sp, site), g in wide.groupby(["species_group", "site_class"]):
                gd = g.set_index("density_class")["mort"]
                if "Low" in gd.index and "High" in gd.index:
                    hits.append(gd["High"] >= gd["Low"])
            crown_proxy = float(np.mean(hits)) if hits else float("nan")
        else:
            crown_proxy = float("nan")

        # Law 4: Mortality size pattern. Within a scenario the mean
        # tree size (qmd proxy from BA and TPA) at year 100 should be
        # larger than at year 0 (stand self thins). Using aggregate
        # wide only gives terminal, so approximate by checking whether
        # QMD_100 > sqrt(BA_0 * k / TPA_0) stays positive-trending:
        # without the initial state we fall back to reporting NA.
        mort_size = float("nan")

        results.append({
            "variant": variant,
            "config": config,
            "law_sukachev_fraction": sukachev,
            "law_eichhorn_fraction": eichhorn,
            "law_density_mortality_fraction": crown_proxy,
            "law_mortality_size_fraction": mort_size,
            "n_scenarios": int(wide["scenario"].nunique()),
        })
    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Aggregate Bakuzis uncertainty ensemble into summaries",
    )
    parser.add_argument("--input", type=str, default=None,
                        help="Single ensemble CSV to process")
    parser.add_argument("--input-dir", type=str, default=None,
                        help="Directory containing bakuzis_uncertainty_*.csv files")
    parser.add_argument("--output-dir", type=str,
                        default="calibration/output/bakuzis")
    parser.add_argument("--quantiles", type=str,
                        default=",".join(str(q) for q in DEFAULT_QUANTILES))
    args = parser.parse_args()

    paths = resolve_inputs(args)
    if not paths:
        logger.error("No input CSV found. Use --input or --input-dir.")
        sys.exit(1)

    quantiles = [float(q) for q in args.quantiles.split(",")]
    df = load_ensemble(paths)
    logger.info(f"Loaded {len(df)} rows from {len(paths)} file(s)")

    if df.empty:
        logger.error("Ensemble frame is empty; nothing to aggregate.")
        sys.exit(1)

    os.makedirs(args.output_dir, exist_ok=True)

    summary = summarize_ensemble(df, quantiles=quantiles)
    summary_path = os.path.join(args.output_dir, "bakuzis_uncertainty_summary_long.csv")
    summary.to_csv(summary_path, index=False)
    logger.info(f"Wrote ensemble summary: {summary_path} ({len(summary)} rows)")

    benchmark = build_benchmark_table(summary)
    bench_path = os.path.join(args.output_dir, "bakuzis_benchmark_wide.csv")
    benchmark.to_csv(bench_path, index=False)
    logger.info(f"Wrote benchmarking table: {bench_path} ({len(benchmark)} rows)")

    if "posterior" in set(df["config"].unique()):
        stand_bands, cell_bands = summarize_stand_bands(df)
        sp = os.path.join(args.output_dir, "bakuzis_uncertainty_stand_bands.csv")
        cp = os.path.join(args.output_dir, "bakuzis_uncertainty_cell_bands.csv")
        stand_bands.to_csv(sp, index=False)
        cell_bands.to_csv(cp, index=False)
        logger.info(f"Wrote per-stand (parameter-only) bands: {sp} ({len(stand_bands)} rows)")
        logger.info(f"Wrote cell within-stand band and between-stand MAP spread: {cp} ({len(cell_bands)} rows)")

    laws = evaluate_bakuzis_laws(summary)
    laws_path = os.path.join(args.output_dir, "bakuzis_laws_compliance.csv")
    laws.to_csv(laws_path, index=False)
    logger.info(f"Wrote Bakuzis law compliance: {laws_path}")

    # Console summary for the user
    print("\nBakuzis law compliance (fraction of scenario pairs that satisfy each law):")
    print(laws.to_string(index=False, float_format="{:.2f}".format))


if __name__ == "__main__":
    main()
