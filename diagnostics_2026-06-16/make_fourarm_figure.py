#!/usr/bin/env python3
"""Headline four-arm figure: within-framework |bias| reduction, engine vs projector (2026-06-18).
Communicates complementarity: the engine keyword calibration fixes size/density (QMD, TPH); the fvs-conus
equations fix level/scatter (BA, volume). Two panels, each its own framework, default vs adjusted |bias|."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

metrics = ["BA", "TPH", "QMD", "VOL"]
# Engine (FVS), median |bias| over 8 variants, OOS: default arm A -> calibrated arm B
eng_def = [11.3, 20.7, 15.7, 13.4]
eng_cal = [7.6, 18.9, 2.2, 8.3]
# Projector (fvs-conus), NE 21,811 undisturbed conditions: |bias| default -> fvs-conus equations (arm C)
prj_def = [12.3, 7.2, 5.3, 15.7]
prj_cal = [7.2, 6.2, 3.1, 9.2]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
x = np.arange(len(metrics)); w = 0.38
C_DEF, C_ADJ = "#9aa7b4", "#1f6f54"

def panel(ax, def_, adj_, title, adj_label):
    b1 = ax.bar(x - w/2, def_, w, label="default", color=C_DEF, edgecolor="white")
    b2 = ax.bar(x + w/2, adj_, w, label=adj_label, color=C_ADJ, edgecolor="white")
    for b in list(b1) + list(b2):
        ax.text(b.get_x()+b.get_width()/2, b.get_height()+0.3, "%.1f" % b.get_height(),
                ha="center", va="bottom", fontsize=8.5)
    ax.set_title(title, fontsize=12, fontweight="bold")
    ax.set_xticks(x); ax.set_xticklabels(metrics)
    ax.set_ylim(0, 24); ax.grid(axis="y", alpha=0.25)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, fontsize=9, loc="upper right")

panel(axes[0], eng_def, eng_cal, "FVS engine: keyword calibration\n(brms maxSDI + density recruitment + BAIMULT), 8 variants OOS", "calibrated")
panel(axes[1], prj_def, prj_cal, "fvs-conus projector: species-free equations\n(NE, 21,811 undisturbed conditions)", "fvs-conus")
axes[0].set_ylabel("median |bias|  (%)", fontsize=11)
fig.suptitle("Complementary gains: keyword calibration fixes QMD/TPH (size, density); fvs-conus equations fix BA/volume (level)",
             fontsize=11.5, y=1.02)
fig.text(0.5, -0.04, "Each arm shown as |bias| reduction within its own framework (engine over-predicts, projector under-predicts undisturbed), "
         "the framework-invariant comparison.", ha="center", fontsize=8.5, color="#555")
fig.tight_layout()
fig.savefig("fourarm_headline_20260618.png", dpi=300, bbox_inches="tight")
fig.savefig("fourarm_headline_20260618_thumb.png", dpi=70, bbox_inches="tight")
print("wrote fourarm_headline_20260618.png")
