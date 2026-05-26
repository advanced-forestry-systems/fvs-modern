#!/usr/bin/env python3
"""
Figure 4: Candidate-covariate importance for the species-free diameter
growth model (DG-Kuehne, B1 species-free), distilled from the 2026-05-15
covariate-exploration job 9572122 (100K-row subsample, 90s wall).

The figure renders two side-by-side panels:

  A. Univariate R^2 on annual DG (cm/yr) for candidate site variables.
     The headline result: SICOND (FIA Site Index, base age 50, +sign)
     beats climate_si (the variable the early manuscript candidate used)
     by approximately 6x and carries the correct positive biological
     direction.

  B. Conditional R^2 added on top of the core covariates
     (ln_dbh + ln_cr_adj + ln_bal_sw_adj) for the 8 strongest candidate
     predictors. STDAGE leads at 0.121 R^2 added, followed closely by
     is_plantation at 0.107 and ELEV at 0.070.

The figure documents the data-driven covariate selection that shifted
the v2 manuscript candidate architecture from v4_full + climate_si to
v6 + SICOND + is_plantation + ELEV.

Input data is hard-coded from the strategic-close handoff (the
exploration script's printed output is the source of truth and the
table is small enough to embed directly).

Output:
  manuscript/v2_figures/fig4_covariate_importance.png
  manuscript/v2_figures/fig4_covariate_importance.pdf
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
OUT_DIR = REPO_ROOT / "manuscript" / "v2_figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Data (from job 9572122, strategic_close.md)
#
# Univariate R^2 against annual DG (cm/yr) for the candidate site
# variables. Sign indicates whether the coefficient is positive or
# negative on DG when fit alone.
# ---------------------------------------------------------------------------

SITE_VARS = [
    {"name": "SICOND",       "r2": 0.044, "sign": "+", "is_winner": True,
     "note": "FIA Site Index, base age 50"},
    {"name": "ELEV",         "r2": 0.039, "sign": "-", "is_winner": False,
     "note": "elevation (topographic position)"},
    {"name": "bgi",          "r2": 0.032, "sign": "-", "is_winner": False,
     "note": "biomass growth index (proxy artifact?)"},
    {"name": "cspi",         "r2": 0.024, "sign": "-", "is_winner": False,
     "note": "composite site productivity index"},
    {"name": "max_biomass",  "r2": 0.017, "sign": "-", "is_winner": False,
     "note": ""},
    {"name": "climate_si",   "r2": 0.007, "sign": "+", "is_winner": False,
     "note": "early manuscript candidate; near-useless univariately"},
    # byi, clim_pca1, clim_pca2 ~ 0; data issue, near-constant. Skipped from chart.
]

# Conditional R^2 added on top of (ln_dbh + ln_cr_adj + ln_bal_sw_adj)
COND_VARS = [
    {"name": "STDAGE",         "r2_added": 0.121, "category": "stand",
     "note": "stand age; reconsider despite -999 issues"},
    {"name": "is_plantation",  "r2_added": 0.107, "category": "disturbance",
     "note": "main-effect upgrade from prior modifier-only role"},
    {"name": "ELEV",           "r2_added": 0.070, "category": "topography",
     "note": "adds beyond climate"},
    {"name": "SICOND",         "r2_added": 0.069, "category": "site",
     "note": "strongest site variable"},
    {"name": "BAL1",           "r2_added": 0.053, "category": "competition",
     "note": "Curtis BAL total; already in models"},
    {"name": "CCFL1",          "r2_added": 0.042, "category": "competition",
     "note": "crown competition factor for larger trees"},
    {"name": "sdi_complexity", "r2_added": 0.010, "category": "structure",
     "note": "sdi_additive1 / SDI_Reineke; modest signal"},
    {"name": "rd_additive",    "r2_added": 0.00008, "category": "structure",
     "note": "no signal beyond Reineke RD"},
]

# Color palette
COL_PRIMARY = "#2E5C8A"
COL_WINNER  = "#1B7837"
COL_LOSER   = "#7A0019"
COL_NEUTRAL = "#888888"

CATEGORY_COLORS = {
    "stand":        "#5DA1D8",
    "disturbance":  "#D77949",
    "topography":   "#A6CB66",
    "site":         "#1B7837",
    "competition":  "#7B7B7B",
    "structure":    "#C84B31",
}

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

fig, (axA, axB) = plt.subplots(1, 2, figsize=(11, 6))

# ---- Panel A: univariate R^2 across candidate site variables ----
site_df = pd.DataFrame(SITE_VARS)
site_df = site_df.sort_values("r2", ascending=True).reset_index(drop=True)

bar_colors = [
    COL_WINNER if row["is_winner"]
    else (COL_LOSER if row["name"] == "climate_si" else COL_PRIMARY)
    for _, row in site_df.iterrows()
]

bars = axA.barh(site_df["name"], site_df["r2"], color=bar_colors, edgecolor="white", height=0.7)
for bar, sign, r2 in zip(bars, site_df["sign"], site_df["r2"]):
    axA.text(
        bar.get_width() + 0.0008,
        bar.get_y() + bar.get_height() / 2,
        f"{r2:.3f} ({sign})",
        va="center", ha="left",
        fontsize=9, color="#333333",
    )

axA.set_xlabel("Univariate R² on annual DG (cm/yr)", fontsize=11)
axA.set_title(
    "A. Site variable head-to-head\nSICOND beats climate_si by ~6× with correct + sign",
    fontsize=11, fontweight="bold", loc="left", pad=10,
)
axA.set_xlim(0, max(site_df["r2"]) * 1.35)
axA.spines["top"].set_visible(False)
axA.spines["right"].set_visible(False)
axA.grid(axis="x", linewidth=0.3, color="#dddddd", zorder=-3)
axA.set_axisbelow(True)
axA.tick_params(axis="y", length=0)

# ---- Panel B: conditional R^2 added ----
cond_df = pd.DataFrame(COND_VARS)
cond_df = cond_df.sort_values("r2_added", ascending=True).reset_index(drop=True)

bar_colors_b = [CATEGORY_COLORS[c] for c in cond_df["category"]]

bars = axB.barh(cond_df["name"], cond_df["r2_added"], color=bar_colors_b, edgecolor="white", height=0.7)
for bar, r2 in zip(bars, cond_df["r2_added"]):
    text = "≈0" if r2 < 0.001 else f"{r2:.3f}"
    axB.text(
        bar.get_width() + 0.002,
        bar.get_y() + bar.get_height() / 2,
        text,
        va="center", ha="left",
        fontsize=9, color="#333333",
    )

axB.set_xlabel("Conditional R² added\non top of ln_dbh + ln_cr_adj + ln_bal_sw_adj", fontsize=11)
axB.set_title(
    "B. Top 8 candidate covariates after core\nSTDAGE leads; is_plantation deserves main-effect status",
    fontsize=11, fontweight="bold", loc="left", pad=10,
)
axB.set_xlim(0, max(cond_df["r2_added"]) * 1.25)
axB.spines["top"].set_visible(False)
axB.spines["right"].set_visible(False)
axB.grid(axis="x", linewidth=0.3, color="#dddddd", zorder=-3)
axB.set_axisbelow(True)
axB.tick_params(axis="y", length=0)

# Category legend for panel B
from matplotlib.patches import Patch
cat_legend = [Patch(facecolor=col, label=cat) for cat, col in CATEGORY_COLORS.items()]
axB.legend(
    handles=cat_legend,
    loc="lower right",
    fontsize=8.5,
    frameon=True,
    framealpha=0.9,
    edgecolor="#cccccc",
    title="Covariate category",
    title_fontsize=9,
)

fig.suptitle(
    "Figure 4. Candidate covariate importance for species-free DG (job 9572122, 100K subsample)",
    fontsize=13, fontweight="bold", y=0.98,
)
fig.text(
    0.5, 0.92,
    "The 90-second exploration that shifted the manuscript candidate from v4_full + climate_si to v6 + SICOND + is_plantation + ELEV",
    ha="center", fontsize=10, color="#555555",
)

plt.tight_layout(rect=[0, 0, 1, 0.91])

png_path = OUT_DIR / "fig4_covariate_importance.png"
pdf_path = OUT_DIR / "fig4_covariate_importance.pdf"
fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(pdf_path, bbox_inches="tight")

# Numerical summary for the Discussion section
print("\n=== Panel A summary (univariate site R^2) ===")
for _, row in site_df.iloc[::-1].iterrows():
    flag = "[winner]" if row["is_winner"] else ""
    print(f"  {row['name']:14s} R²={row['r2']:.3f} ({row['sign']})  {flag}")

print("\n=== Panel B summary (conditional R^2 added) ===")
for _, row in cond_df.iloc[::-1].iterrows():
    print(f"  {row['name']:18s} R²_added={row['r2_added']:.4f}  [{row['category']}]")

print(f"\nWrote:\n  {png_path}\n  {pdf_path}")
