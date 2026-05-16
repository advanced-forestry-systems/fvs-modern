#!/usr/bin/env python3
"""
Figure 3: Modifier coefficient comparison across the 6 CONUS Phase 4
components at the production lambda=10 regularization choice.

Each row is a modifier coefficient (alpha_plant, alpha_cutting,
beta_plant_dia, etc.). Each column-grouped marker is one of the 6
component fits: HG, HT_DBH, CR, HCB, Mortality, DG (Kuehne). 90% CIs
are rendered as horizontal error bars.

The figure makes the cross-disturbance, cross-component pattern legible:
  1. Planted and cutting alphas have strong, sign-consistent effects on
     most components — these are the well-identified disturbance signals.
  2. Insect, disease, and fire alphas have substantial mortality and
     crown-base effects but weak height-DBH effects.
  3. Site prep, wind, and harvest alphas span zero in most components,
     indicating they aren't sharply identified at this lambda.

Input:
  manuscript/v2_inputs/manuscript_table1_lambda10_20260510.csv

Output:
  manuscript/v2_figures/fig3_modifier_coefficients.png
  manuscript/v2_figures/fig3_modifier_coefficients.pdf

Usage:
  python3 manuscript/v2_figures/fig3_modifier_coefficients.py
"""

from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent.parent
INPUT_DIR = REPO_ROOT / "manuscript" / "v2_inputs"
OUT_DIR = REPO_ROOT / "manuscript" / "v2_figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)

INPUT_PATH = INPUT_DIR / "manuscript_table1_lambda10_20260510.csv"

# Display labels for the disturbance-class alpha coefficients
# (intercept terms in the modifier model per disturbance category)
ALPHA_LABELS = {
    "alpha_plant":     "Planted (afforestation / replanting)",
    "alpha_cutting":   "Cutting (silvicultural treatment)",
    "alpha_siteprep":  "Site preparation",
    "alpha_insect":    "Insect damage",
    "alpha_disease":   "Disease damage",
    "alpha_fire":      "Fire",
    "alpha_wind":      "Wind",
    "alpha_harvest":   "Harvest",
}

# Component column labels in the input plus display names
COMPONENT_DISPLAY = [
    ("HG",     "Height growth (HG)"),
    ("HT_DBH", "Height-DBH (HT_DBH)"),
    ("CR",     "Crown ratio (CR)"),
    ("HCB",    "Crown base height (HCB)"),
    ("Mort",   "Mortality (Mort)"),
    ("DG_Kue", "Diameter growth (DG Kuehne)"),
]

# Color palette: one per component
COLORS = {
    "HG":     "#2E5C8A",  # blue
    "HT_DBH": "#5DA1D8",  # light blue
    "CR":     "#4DAF7C",  # green
    "HCB":    "#A6CB66",  # light green
    "Mort":   "#C84B31",  # red
    "DG_Kue": "#D77949",  # orange
}

# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------

CELL_RE = re.compile(r"([-+]?[0-9.eE]+)\s*\[\s*([-+]?[0-9.eE]+)\s*,\s*([-+]?[0-9.eE]+)\s*\]")


def parse_cell(s):
    """Parse 'mean [low,high]' string into (mean, low, high)."""
    if pd.isna(s):
        return (np.nan, np.nan, np.nan)
    m = CELL_RE.search(str(s))
    if not m:
        return (np.nan, np.nan, np.nan)
    return float(m.group(1)), float(m.group(2)), float(m.group(3))


# ---------------------------------------------------------------------------
# Load
# ---------------------------------------------------------------------------

df = pd.read_csv(INPUT_PATH)
print(f"Loaded {INPUT_PATH.name}: {len(df)} modifier rows, {df.shape[1] - 1} components")

# Keep only modifiers we have labels for, in the documented order
df = df[df["alpha"].isin(ALPHA_LABELS.keys())].copy()
df["order"] = df["alpha"].map({k: i for i, k in enumerate(ALPHA_LABELS)})
df = df.sort_values("order").reset_index(drop=True)

# Build the long form for plotting
records = []
for _, row in df.iterrows():
    alpha = row["alpha"]
    label = ALPHA_LABELS[alpha]
    order = row["order"]
    for comp_key, comp_disp in COMPONENT_DISPLAY:
        mean, low, high = parse_cell(row[comp_key])
        records.append({
            "alpha": alpha,
            "label": label,
            "order": order,
            "comp_key": comp_key,
            "comp_disp": comp_disp,
            "mean": mean,
            "low": low,
            "high": high,
            "excludes_zero": (low > 0) or (high < 0) if not np.isnan(low) else False,
        })

long_df = pd.DataFrame(records)

# ---------------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------------

n_rows = len(ALPHA_LABELS)
n_comp = len(COMPONENT_DISPLAY)
band_height = 1.0
comp_dy = band_height * 0.7 / (n_comp - 1)
y_offsets = np.linspace(-band_height * 0.35, band_height * 0.35, n_comp)

fig, ax = plt.subplots(figsize=(9, 7.5))

for _, rec in long_df.iterrows():
    if np.isnan(rec["mean"]):
        continue
    color = COLORS[rec["comp_key"]]
    comp_idx = [k for k, _ in COMPONENT_DISPLAY].index(rec["comp_key"])
    y = (n_rows - 1 - rec["order"]) + y_offsets[comp_idx]

    ax.errorbar(
        rec["mean"], y,
        xerr=[[rec["mean"] - rec["low"]], [rec["high"] - rec["mean"]]],
        fmt="o",
        color=color,
        markersize=5,
        markerfacecolor=color if rec["excludes_zero"] else "white",
        markeredgecolor=color,
        markeredgewidth=1.0,
        elinewidth=1.2,
        capsize=0,
        alpha=0.9,
    )

ax.axvline(0, color="grey", linewidth=0.6, linestyle="--", zorder=-2)
ax.set_yticks(range(n_rows))
ax.set_yticklabels([ALPHA_LABELS[a] for a in reversed(list(ALPHA_LABELS.keys()))])
ax.set_xlabel("Posterior mean with 90% CI", fontsize=11)
ax.set_title(
    "Figure 3. CONUS Phase 4 modifier coefficients across six components (lambda = 10)",
    fontsize=12, fontweight="bold", loc="left", pad=10,
)

# Horizontal banding to separate modifier rows visually
for i in range(n_rows):
    if i % 2 == 0:
        ax.axhspan(i - 0.5, i + 0.5, color="#f5f5f5", zorder=-5)

ax.set_ylim(-0.6, n_rows - 0.4)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="x", linewidth=0.3, color="#dddddd", zorder=-3)
ax.set_axisbelow(True)
ax.tick_params(axis="y", length=0)

# Legend
legend_handles = [
    Line2D([0], [0], marker="o", color=color, markerfacecolor=color,
           markeredgecolor=color, markersize=6, linewidth=1.2,
           label=disp)
    for (key, disp), color in zip(COMPONENT_DISPLAY,
                                  [COLORS[k] for k, _ in COMPONENT_DISPLAY])
]
ax.legend(
    handles=legend_handles,
    loc="lower right",
    fontsize=8.5,
    frameon=True,
    framealpha=0.92,
    edgecolor="#cccccc",
    ncol=2,
    title="Component (filled = excludes 0)",
    title_fontsize=9,
)

fig.text(
    0.01, 0.005,
    "Filled markers = 90% CI excludes zero. Open markers = CI spans zero.",
    fontsize=8, color="#666666",
)

plt.tight_layout(rect=[0, 0.02, 1, 1])

png_path = OUT_DIR / "fig3_modifier_coefficients.png"
pdf_path = OUT_DIR / "fig3_modifier_coefficients.pdf"
fig.savefig(png_path, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(pdf_path, bbox_inches="tight")

# Brief excludes-zero count per component for the Discussion section
print("\nExcludes-zero count per component (out of 8 modifiers):")
for key, disp in COMPONENT_DISPLAY:
    n = long_df.loc[long_df["comp_key"] == key, "excludes_zero"].sum()
    n_total = long_df.loc[long_df["comp_key"] == key, "mean"].notna().sum()
    print(f"  {disp:34s} {n}/{n_total}")

print(f"\nWrote:\n  {png_path}\n  {pdf_path}")
