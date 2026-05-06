#!/usr/bin/env python3
"""
Bakuzis preview figures from the existing point estimate CSV
============================================================

Reads calibration/output/bakuzis/bakuzis_100yr_results.csv (produced by
bakuzis_100yr_comparison.py, point estimates only) and renders three
preview figures using matplotlib:

  1. bakuzis_trajectories_preview.png
     3 x 4 grid of BA trajectories, one row per site class, one column
     per species group, lines colored by density class and linetype
     by default vs calibrated.

  2. bakuzis_year100_divergence_preview.png
     Scatter of calibrated minus default BA at year 100 vs default BA
     at year 100, colored by species group and shaped by site class.

  3. bakuzis_laws_preview.png
     Bar chart of the three evaluable Bakuzis law fractions (Sukachev,
     Eichhorn, density driven mortality) for default vs calibrated.

These are intentionally preview quality for review. The publication
grade figures live in the R script (calibration/R/17_bakuzis_uncertainty_figures.R)
which is meant to be run after the uncertainty pipeline produces the
summary long frame with posterior bands.

Usage:
  python bakuzis_preview_figures.py
  python bakuzis_preview_figures.py --input path/to/results.csv --outdir /tmp/figs

Author: A. Weiskittel
Date: 2026-04-23
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

# Colors by density class (colorblind friendly)
DENSITY_COLORS = {
    "Low": "#1b9e77",
    "Medium": "#d95f02",
    "High": "#7570b3",
}
SITE_ORDER = ["Low", "Medium", "High"]
SPECIES_ORDER = ["Spruce-Fir", "Northern-Hardwood", "Pine", "Oak-Pine"]
DENSITY_ORDER = ["Low", "Medium", "High"]


def load_results(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    # Normalize column names (handle 'year' vs 'Year')
    if "Year" in df.columns and "year" not in df.columns:
        df = df.rename(columns={"Year": "year"})
    # Compute projection horizon year since the scenario starts at inv_year
    df["horizon_yr"] = df.groupby(["variant", "scenario"])["year"].transform(
        lambda s: s - s.min()
    )
    return df


def plot_trajectory_grid(df: pd.DataFrame, outpath: str) -> None:
    """3 x 4 grid of BA vs horizon year."""
    fig, axes = plt.subplots(
        len(SITE_ORDER), len(SPECIES_ORDER),
        figsize=(13, 9), sharex=True, sharey="row",
    )
    for i, site in enumerate(SITE_ORDER):
        for j, sp in enumerate(SPECIES_ORDER):
            ax = axes[i, j]
            for dens in DENSITY_ORDER:
                sub = df.loc[
                    (df["site_class"] == site)
                    & (df["species_group"] == sp)
                    & (df["density_class"] == dens)
                ]
                if sub.empty:
                    continue
                for config, ls in [("default", "--"), ("calibrated", "-")]:
                    g = sub.loc[sub["config"] == config]
                    if g.empty:
                        continue
                    g = g.sort_values("horizon_yr")
                    ax.plot(
                        g["horizon_yr"],
                        g["atba"],
                        ls=ls,
                        color=DENSITY_COLORS[dens],
                        alpha=0.85,
                        lw=1.5,
                    )
            if i == 0:
                ax.set_title(sp, fontsize=10)
            if j == 0:
                ax.set_ylabel(f"{site} site\nBA (ft²/ac)", fontsize=9)
            if i == len(SITE_ORDER) - 1:
                ax.set_xlabel("Projection year", fontsize=9)
            ax.grid(alpha=0.3)

    # Combined legend: density colors + config linestyles
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], color=DENSITY_COLORS[d], lw=2, label=f"{d} density")
        for d in DENSITY_ORDER
    ] + [
        Line2D([0], [0], color="black", lw=1.5, ls="--", label="Default"),
        Line2D([0], [0], color="black", lw=1.5, ls="-", label="Calibrated"),
    ]
    fig.legend(
        handles=handles, loc="lower center", ncol=len(handles),
        frameon=False, fontsize=9, bbox_to_anchor=(0.5, -0.02),
    )
    fig.suptitle(
        "Bakuzis matrix: 100 year BA trajectories "
        "(default vs calibrated, point estimates)",
        fontsize=12, y=0.99,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.97])
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Wrote {outpath}")


def plot_year100_divergence(df: pd.DataFrame, outpath: str) -> None:
    """Scatter of (calibrated minus default) BA at year 100 against default BA."""
    # Take terminal horizon year per scenario
    terminal = df.loc[df["horizon_yr"] == df["horizon_yr"].max()].copy()
    if terminal.empty:
        logger.warning("No terminal horizon rows; skipping divergence figure")
        return
    wide = terminal.pivot_table(
        index=["species_group", "site_class", "density_class", "scenario"],
        columns="config",
        values="atba",
    ).reset_index()
    if "default" not in wide.columns or "calibrated" not in wide.columns:
        logger.warning("Need both default and calibrated to plot divergence")
        return
    wide["diff_pct"] = 100 * (wide["calibrated"] - wide["default"]) / wide["default"].replace(0, np.nan)

    fig, ax = plt.subplots(figsize=(8, 6))
    markers = {"Low": "o", "Medium": "s", "High": "^"}
    for sp, g_sp in wide.groupby("species_group"):
        for site, g_site in g_sp.groupby("site_class"):
            ax.scatter(
                g_site["default"],
                g_site["diff_pct"],
                marker=markers.get(site, "x"),
                label=f"{sp} / {site}" if site == "Medium" else None,
                s=55, alpha=0.75, edgecolor="black", linewidth=0.5,
            )
    ax.axhline(0, color="gray", lw=0.8, ls=":")
    ax.set_xlabel("Default BA at year 100 (ft²/ac)")
    ax.set_ylabel("Calibrated minus default (percent of default)")
    ax.set_title("Year 100 BA divergence: calibrated vs default")
    ax.grid(alpha=0.3)
    # Legend: species groups (color) + site classes (marker)
    from matplotlib.lines import Line2D
    color_handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=8,
               color="C" + str(i), label=sp)
        for i, sp in enumerate(sorted(wide["species_group"].unique()))
    ]
    shape_handles = [
        Line2D([0], [0], marker=markers[s], linestyle="", markersize=8,
               color="gray", label=f"{s} site")
        for s in SITE_ORDER
    ]
    ax.legend(handles=color_handles + shape_handles, loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Wrote {outpath}")


def evaluate_preview_laws(df: pd.DataFrame) -> pd.DataFrame:
    """Lightweight Bakuzis law compliance on point estimates."""
    terminal = df.loc[df["horizon_yr"] == df["horizon_yr"].max()].copy()
    results = []
    for (variant, config), g in terminal.groupby(["variant", "config"]):
        wide = g.pivot_table(
            index=["scenario", "species_group", "site_class", "density_class"],
            values=["atba", "tpa", "tcuft", "mort"],
        ).reset_index()

        # Sukachev: volume per tree (tcuft / tpa) higher under low density
        wide["vol_per_tree"] = wide["tcuft"] / wide["tpa"].replace(0, np.nan)
        hits = []
        for (sp, site), sg in wide.groupby(["species_group", "site_class"]):
            d = sg.set_index("density_class")["vol_per_tree"]
            if "Low" in d.index and "High" in d.index and d.notna().all():
                hits.append(float(d["High"] < d["Low"]))
        sukachev = float(np.mean(hits)) if hits else float("nan")

        # Eichhorn: BA increases with site quality
        hits = []
        for (sp, dens), sg in wide.groupby(["species_group", "density_class"]):
            d = sg.set_index("site_class")["atba"]
            if {"Low", "Medium", "High"}.issubset(d.index):
                hits.append(float((d["Low"] <= d["Medium"]) and (d["Medium"] <= d["High"])))
        eichhorn = float(np.mean(hits)) if hits else float("nan")

        # Density driven mortality
        hits = []
        for (sp, site), sg in wide.groupby(["species_group", "site_class"]):
            d = sg.set_index("density_class")["mort"]
            if "Low" in d.index and "High" in d.index:
                hits.append(float(d["High"] >= d["Low"]))
        dens_mort = float(np.mean(hits)) if hits else float("nan")

        results.append({
            "variant": variant,
            "config": config,
            "Sukachev (competition reduces per tree volume)": sukachev,
            "Eichhorn (better site, more BA)": eichhorn,
            "Density drives mortality": dens_mort,
        })
    return pd.DataFrame(results)


def plot_laws(laws: pd.DataFrame, outpath: str) -> None:
    if laws.empty:
        logger.warning("No law results to plot")
        return
    law_cols = [c for c in laws.columns if c not in ("variant", "config")]
    # One panel per variant, grouped bars by config
    variants = laws["variant"].unique()
    fig, axes = plt.subplots(1, len(variants), figsize=(5 * len(variants), 4.5), sharey=True)
    if len(variants) == 1:
        axes = [axes]
    for ax, v in zip(axes, variants):
        sub = laws.loc[laws["variant"] == v]
        x = np.arange(len(law_cols))
        width = 0.35
        for i, config in enumerate(sub["config"].unique()):
            vals = sub.loc[sub["config"] == config, law_cols].values.flatten()
            ax.bar(x + i * width, vals, width, label=config)
        ax.set_xticks(x + width / 2)
        ax.set_xticklabels(
            [c.split(" (")[0] for c in law_cols], rotation=20, ha="right",
        )
        ax.set_ylim(0, 1.05)
        ax.set_title(v)
        ax.grid(axis="y", alpha=0.3)
        ax.set_ylabel("Fraction of scenario pairs satisfying law")
    axes[-1].legend(loc="lower right", fontsize=9)
    fig.suptitle("Bakuzis biological law compliance at year 100 (preview)", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Wrote {outpath}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input", type=str,
        default="calibration/output/bakuzis/bakuzis_100yr_results.csv",
    )
    parser.add_argument("--outdir", type=str, default="calibration/output/bakuzis")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        raise SystemExit(f"Input CSV not found: {args.input}")
    os.makedirs(args.outdir, exist_ok=True)

    df = load_results(args.input)
    logger.info(
        f"Loaded {len(df)} rows, variants={df['variant'].unique().tolist()}, "
        f"configs={df['config'].unique().tolist()}"
    )

    plot_trajectory_grid(df, os.path.join(args.outdir, "bakuzis_trajectories_preview.png"))
    plot_year100_divergence(df, os.path.join(args.outdir, "bakuzis_year100_divergence_preview.png"))

    laws = evaluate_preview_laws(df)
    laws_csv = os.path.join(args.outdir, "bakuzis_laws_preview.csv")
    laws.to_csv(laws_csv, index=False)
    logger.info(f"Wrote {laws_csv}")
    plot_laws(laws, os.path.join(args.outdir, "bakuzis_laws_preview.png"))

    print("\nPreview Bakuzis law compliance (fraction satisfying each law):")
    print(laws.to_string(index=False, float_format="{:.2f}".format))


if __name__ == "__main__":
    main()
