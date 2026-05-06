#!/usr/bin/env python3
"""
Uncertainty Bakuzis figures via matplotlib (sandbox friendly).

Mirrors 17_bakuzis_uncertainty_figures.R for environments without R.
Consumes the long summary frame from bakuzis_uncertainty_aggregate.py
and renders the four publication panels as PNG.

Usage:
  python bakuzis_uncertainty_figures.py
  python bakuzis_uncertainty_figures.py --summary path/to/summary.csv \
      --bench path/to/benchmark.csv --laws path/to/laws.csv \
      --outdir figures/

Author: A. Weiskittel
Date: 2026-04-23
"""

from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger(__name__)

CONFIG_COLORS = {
    "default": "#666666",
    "calibrated_map": "#d95f02",
    "posterior_band": "#1b9e77",
}
CONFIG_LABELS = {
    "default": "Default parameters",
    "calibrated_map": "Calibrated MAP",
    "posterior_band": "Calibrated posterior (95% band)",
}
SITE_ORDER = ["Low", "Medium", "High"]
SPECIES_ORDER = ["Spruce-Fir", "Northern-Hardwood", "Pine", "Oak-Pine"]
DENSITY_ORDER = ["Low", "Medium", "High"]


def add_horizon(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["horizon_yr"] = df.groupby(["variant", "scenario"])["year"].transform(
        lambda s: s - s.min()
    )
    return df


def plot_trajectory_bands(summary: pd.DataFrame, outpath: str,
                           variant: str = "NE") -> None:
    ba = summary.loc[
        (summary["variable"] == "atba") & (summary["variant"] == variant)
    ].copy()
    if ba.empty:
        logger.warning(f"No atba rows for variant {variant}; skipping trajectory figure")
        return

    fig, axes = plt.subplots(
        len(SITE_ORDER), len(SPECIES_ORDER),
        figsize=(14, 10), sharex=True, sharey="row",
    )
    for i, site in enumerate(SITE_ORDER):
        for j, sp in enumerate(SPECIES_ORDER):
            ax = axes[i, j]
            for dens in DENSITY_ORDER:
                sub = ba.loc[
                    (ba["site_class"] == site)
                    & (ba["species_group"] == sp)
                    & (ba["density_class"] == dens)
                ].sort_values(["config", "horizon_yr"])
                if sub.empty:
                    continue
                # Ribbon for posterior band
                post = sub.loc[sub["config"] == "posterior_band"]
                if not post.empty and "q025" in post.columns and "q975" in post.columns:
                    ax.fill_between(
                        post["horizon_yr"], post["q025"], post["q975"],
                        color=CONFIG_COLORS["posterior_band"], alpha=0.22,
                    )
                # Lines per config
                for config, ls in [
                    ("default", "--"),
                    ("calibrated_map", "-"),
                    ("posterior_band", "-"),
                ]:
                    g = sub.loc[sub["config"] == config]
                    if g.empty:
                        continue
                    lw = 0.9 if config == "posterior_band" else 1.5
                    ax.plot(
                        g["horizon_yr"], g["median"],
                        color=CONFIG_COLORS[config], ls=ls, lw=lw,
                        alpha=0.9,
                    )
            if i == 0:
                ax.set_title(sp, fontsize=11)
            if j == 0:
                ax.set_ylabel(f"{site} site\nBA (ft²/ac)", fontsize=10)
            if i == len(SITE_ORDER) - 1:
                ax.set_xlabel("Projection year", fontsize=10)
            ax.grid(alpha=0.3)

    # Facet annotation: density class band in faint gray
    handles = [
        Line2D([0], [0], color=CONFIG_COLORS["default"], lw=1.5, ls="--",
               label=CONFIG_LABELS["default"]),
        Line2D([0], [0], color=CONFIG_COLORS["calibrated_map"], lw=1.8,
               label=CONFIG_LABELS["calibrated_map"]),
        Line2D([0], [0], color=CONFIG_COLORS["posterior_band"], lw=1.8,
               label=CONFIG_LABELS["posterior_band"]),
    ]
    fig.legend(
        handles=handles, loc="lower center", ncol=len(handles),
        frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.01),
    )
    fig.suptitle(
        f"Bakuzis matrix ({variant}): 100 year BA trajectories with parametric uncertainty\n"
        "Ribbon is 95 percent posterior credible band around the calibrated median",
        fontsize=12, y=0.99,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.96])
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Wrote {outpath}")


def plot_divergence(benchmark: pd.DataFrame, outpath: str,
                     variant: str = "NE") -> None:
    ba = benchmark.loc[
        (benchmark["variable"] == "atba") & (benchmark["variant"] == variant)
    ].copy()
    if ba.empty or "def_median_y100" not in ba.columns:
        logger.warning(f"Benchmark table missing y100 for {variant}; skipping divergence")
        return
    ba["band_pct"] = (
        100
        * (ba.get("post_band_y100", 0))
        / ba["post_median_y100"].replace(0, np.nan)
    )
    ba["diff_pct"] = ba.get("diff_pct_y100")

    fig, ax = plt.subplots(figsize=(9, 6.5))
    markers = {"Low": "o", "Medium": "s", "High": "^"}
    species_cmap = {sp: f"C{i}" for i, sp in enumerate(SPECIES_ORDER)}

    for sp in SPECIES_ORDER:
        sub_sp = ba[ba["species_group"] == sp]
        for site in SITE_ORDER:
            g = sub_sp[sub_sp["site_class"] == site]
            sizes = g["band_pct"].fillna(4.0)
            sizes = np.clip(sizes * 2, 20, 220)
            ax.scatter(
                g["def_median_y100"], g["diff_pct"],
                marker=markers.get(site, "x"), s=sizes,
                alpha=0.75, edgecolor="black", linewidth=0.5,
                color=species_cmap.get(sp, "gray"),
            )

    ax.axhline(0, color="gray", lw=0.8, ls=":")
    ax.set_xlabel("Default BA at year 100 (ft²/ac)")
    ax.set_ylabel("Calibrated posterior median minus default (% of default)")
    ax.set_title(
        f"Year 100 BA divergence ({variant}) with posterior band width as marker size"
    )
    ax.grid(alpha=0.3)

    species_handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=8,
               color=species_cmap[sp], label=sp, markeredgecolor="black",
               markeredgewidth=0.4)
        for sp in SPECIES_ORDER
    ]
    site_handles = [
        Line2D([0], [0], marker=markers[s], linestyle="", markersize=8,
               color="gray", label=f"{s} site", markeredgecolor="black",
               markeredgewidth=0.4)
        for s in SITE_ORDER
    ]
    leg1 = ax.legend(
        handles=species_handles, loc="upper left", fontsize=9,
        title="Species group", frameon=False,
    )
    ax.add_artist(leg1)
    ax.legend(
        handles=site_handles, loc="lower right", fontsize=9,
        title="Site class", frameon=False,
    )

    fig.tight_layout()
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Wrote {outpath}")


def plot_laws(laws: pd.DataFrame, outpath: str) -> None:
    if laws.empty:
        logger.warning("Empty laws frame; skipping")
        return
    law_cols = [c for c in laws.columns if c.startswith("law_")]
    pretty = {
        "law_sukachev_fraction": "Sukachev\n(competition reduces per tree)",
        "law_eichhorn_fraction": "Eichhorn\n(site quality increases BA)",
        "law_density_mortality_fraction": "Density drives mortality",
        "law_mortality_size_fraction": "Mortality size pattern",
    }
    variants = laws["variant"].unique()
    fig, axes = plt.subplots(1, len(variants), figsize=(5 * len(variants), 5),
                             sharey=True, squeeze=False)
    configs = ["default", "calibrated_map", "posterior_band"]
    width = 0.27
    for ax_idx, v in enumerate(variants):
        ax = axes[0, ax_idx]
        sub = laws[laws["variant"] == v]
        x = np.arange(len(law_cols))
        for i, cfg in enumerate(configs):
            row = sub[sub["config"] == cfg]
            if row.empty:
                continue
            vals = row[law_cols].values.flatten()
            ax.bar(
                x + i * width - width, vals, width,
                color=CONFIG_COLORS.get(cfg, "gray"),
                label=CONFIG_LABELS.get(cfg, cfg),
                edgecolor="black", linewidth=0.4,
            )
        ax.set_xticks(x)
        ax.set_xticklabels([pretty.get(c, c) for c in law_cols],
                           fontsize=9)
        ax.set_ylim(0, 1.05)
        ax.set_title(v)
        ax.set_ylabel("Fraction of scenario pairs satisfying law")
        ax.grid(axis="y", alpha=0.3)
    axes[0, -1].legend(loc="lower right", fontsize=9, frameon=False)
    fig.suptitle(
        "Bakuzis biological law compliance at year 100",
        fontsize=13, y=0.98,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Wrote {outpath}")


def plot_band_growth(summary: pd.DataFrame, outpath: str,
                      variant: str = "NE") -> None:
    ba = summary.loc[
        (summary["variable"] == "atba")
        & (summary["config"] == "posterior_band")
        & (summary["variant"] == variant)
    ].copy()
    if ba.empty:
        logger.warning(f"No posterior rows for variant {variant}; skipping band growth")
        return
    ba = add_horizon(ba)
    ba["band_pct"] = 100 * (ba["q975"] - ba["q025"]) / ba["median"].replace(0, np.nan)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex=True, sharey=True)
    palette = {"Low": "#1b9e77", "Medium": "#d95f02", "High": "#7570b3"}
    for ax, sp in zip(axes.flat, SPECIES_ORDER):
        sub = ba[ba["species_group"] == sp]
        for site in SITE_ORDER:
            subsub = sub[sub["site_class"] == site]
            for scen, g in subsub.groupby("scenario"):
                ax.plot(
                    g["horizon_yr"], g["band_pct"],
                    color=palette[site], alpha=0.55, lw=0.9,
                )
        ax.set_title(sp, fontsize=11)
        ax.set_xlabel("Projection year")
        ax.set_ylabel("BA 95% band width (% of median)")
        ax.grid(alpha=0.3)
    site_handles = [
        Line2D([0], [0], color=palette[s], lw=2, label=f"{s} site")
        for s in SITE_ORDER
    ]
    axes[0, 1].legend(handles=site_handles, loc="upper left",
                      fontsize=9, frameon=False)
    fig.suptitle(
        f"Posterior uncertainty growth with projection horizon ({variant})",
        fontsize=13, y=0.99,
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(outpath, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Wrote {outpath}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--summary", type=str,
        default="calibration/output/bakuzis/bakuzis_uncertainty_summary_long.csv",
    )
    parser.add_argument(
        "--bench", type=str,
        default="calibration/output/bakuzis/bakuzis_benchmark_wide.csv",
    )
    parser.add_argument(
        "--laws", type=str,
        default="calibration/output/bakuzis/bakuzis_laws_compliance.csv",
    )
    parser.add_argument(
        "--outdir", type=str,
        default="calibration/output/comparisons/manuscript_figures",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    summary = pd.read_csv(args.summary)
    summary = add_horizon(summary)

    # Render per-variant trajectory and band-growth figures for every
    # variant present in the summary. Filenames carry the variant
    # suffix; the historical NE-only filename is preserved by also
    # writing the NE outputs without a suffix.
    variants_present = sorted(summary["variant"].dropna().unique())
    for v in variants_present:
        suffix = "" if v.upper() == "NE" else f"_{v.lower()}"
        plot_trajectory_bands(
            summary,
            os.path.join(args.outdir, f"fig_bakuzis_trajectories{suffix}.png"),
            variant=v,
        )
        plot_band_growth(
            summary,
            os.path.join(args.outdir, f"fig_bakuzis_band_growth{suffix}.png"),
            variant=v,
        )

    if os.path.exists(args.bench):
        bench = pd.read_csv(args.bench)
        for v in variants_present:
            suffix = "" if v.upper() == "NE" else f"_{v.lower()}"
            plot_divergence(
                bench,
                os.path.join(args.outdir, f"fig_bakuzis_divergence{suffix}.png"),
                variant=v,
            )
    if os.path.exists(args.laws):
        laws = pd.read_csv(args.laws)
        plot_laws(laws, os.path.join(args.outdir, "fig_bakuzis_laws.png"))

    logger.info("Done.")


if __name__ == "__main__":
    main()
