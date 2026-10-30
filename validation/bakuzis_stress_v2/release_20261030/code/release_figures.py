#!/usr/bin/env python3
"""release_figures.py (2026-10-29)
Graphical abstract and release figures for the fvs-modern CONUS configuration catalog deposit.
Inputs: data/catalog_gates_by_arm.csv, data/stress_by_arm.csv, data/stress_by_variant.csv.
Usage: python3 release_figures.py --data data --out fig [--doi 10.5281/zenodo.NNNNNNN]
Headless (Agg), 300 dpi. Error bars: Wilson 95 percent intervals on exceedance shares where shown;
the stand bootstrap probabilities are themselves the uncertainty summary and carry no second interval.
"""
import argparse, os
os.environ.setdefault("MPLCONFIGDIR", "/tmp/mplcfg")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd, numpy as np

GREEN, GOLD, GREY, RED, BLUE = "#1a3d28", "#c5a55a", "#8a8f87", "#a23b2a", "#2f5f8a"
ARMS = ["default_native", "default_gompit", "calibrated_native", "calibrated_gompit"]
ARML = {"default_native": "default, native", "default_gompit": "default, GOMPIT",
        "calibrated_native": "calibrated, native", "calibrated_gompit": "calibrated, GOMPIT"}
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.facecolor": "white", "figure.facecolor": "white", "axes.grid": False})

def wilson(x, n, z=1.96):
    p = x / n
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h

def load(d):
    g = pd.read_csv(os.path.join(d, "catalog_gates_by_arm.csv"))
    g = g[g.arm.isin(ARMS)].copy()
    sa = pd.read_csv(os.path.join(d, "stress_by_arm.csv"))
    sv = pd.read_csv(os.path.join(d, "stress_by_variant.csv"))
    return g, sa, sv

def panel_variants(ax, g, sv):
    rec = g.groupby("group").admissible.any()
    sv = sv.set_index("group").loc[sorted(sv.group)]
    cols = [GREEN if rec[v] else GREY for v in sv.index]
    ax.bar(sv.index, sv.boot_any_admissible_prob, color=cols, width=0.7)
    ax.scatter(sv.index, sv.grid_any_share, color=GOLD, zorder=3, s=22, label="share of 144 rule cells")
    ax.set_ylim(0, 1.3); ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0]); ax.set_ylabel("Probability of an admissible engine arm")
    ax.set_title("A. Variant verdicts: 500 stand bootstraps (bars), 144 rule cells (points)", loc="left", fontsize=9.5, fontweight="bold")
    from matplotlib.patches import Patch
    h, l = ax.get_legend_handles_labels()
    ax.legend(h + [Patch(color=GREEN), Patch(color=GREY)], l + ["RECOMMENDED (7 variants)", "NOT RECOMMENDED (11)"], frameon=False, loc="upper center", fontsize=7.5, ncol=3)
    ax.tick_params(axis="x", labelsize=8)

def panel_gates(ax, g):
    for a, m in zip(ARMS, ["o", "s", "^", "D"]):
        s = g[g.arm == a]
        x = s.reineke_flag_rate.fillna(-0.02).values * 100
        y = s.share_ba999_50.values * 100
        lo, hi = wilson(s.x_ba999_50.values, s.n_traj.values)
        col = [GREEN if v else GREY for v in s.admissible]
        ax.errorbar(x, y, yerr=[y - lo * 100, hi * 100 - y], fmt="none", ecolor="#cccccc", elinewidth=0.6, capsize=0, zorder=1)
        ax.scatter(x, y, marker=m, c=col, s=26, edgecolor="white", linewidth=0.4, zorder=2, label=ARML[a])
    ax.axvline(10, color=RED, lw=0.8, ls="--"); ax.axhline(5, color=RED, lw=0.8, ls="--")
    ax.set_xlabel("Reineke flag rate among density-gated trajectories (%)")
    ax.set_ylabel("Stands above FIA BA p99.9 by year 50 (%, Wilson 95%)")
    ax.set_xlim(-4, 90); ax.set_ylim(-1, 40)
    ax.text(10.8, 38, "10% Reineke gate", color=RED, fontsize=7.5, va="top"); ax.text(60, 6, "5% magnitude gate", color=RED, fontsize=7.5)
    ax.text(-3.5, 36, "x < 0: not testable", fontsize=7, color=GREY, rotation=90, va="top")
    ax.set_title("B. The two gates that decide (72 arms)", loc="left", fontsize=9.5, fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5, loc="upper right", title="arm (green = admissible)", title_fontsize=7.5)

def panel_stability(ax, sa):
    sa = sa[sa.arm.isin(ARMS)]
    col = {"STABLE admissible": GREEN, "FRAGILE admissible": GOLD, "FRAGILE not admissible": BLUE, "STABLE not admissible": GREY}
    for k, c in col.items():
        s = sa[sa.verdict_stability == k]
        ax.scatter(s.grid_admissible_share, s.boot_admissible_prob, c=c, s=28, edgecolor="white", linewidth=0.4, label=f"{k} ({len(s)})", zorder=2)
    fr = sa[sa.verdict_stability == "FRAGILE admissible"].sort_values(["grid_admissible_share", "boot_admissible_prob"]).reset_index(drop=True)
    for i, r in fr.iterrows():
        dx, dy = (5, -3) if i % 2 == 0 else (-5, 3)
        ax.annotate(f"{r.group} {ARML[r.arm].replace('calibrated', 'cal').replace('default', 'def')}", (r.grid_admissible_share, r.boot_admissible_prob),
                    fontsize=5.5, xytext=(dx, dy), textcoords="offset points", ha="left" if dx > 0 else "right", color="#333333")
    ax.axhline(0.9, color=GREY, lw=0.6, ls=":"); ax.axhline(0.1, color=GREY, lw=0.6, ls=":"); ax.axvline(0.5, color=GREY, lw=0.6, ls=":")
    ax.set_xlabel("Share of 144 rule cells admitting the arm"); ax.set_ylabel("Bootstrap probability admissible (500 draws)")
    ax.set_xlim(-0.03, 1.08); ax.set_ylim(-0.03, 1.05)
    ax.set_title("C. Verdict stability (0 of 72 overturned)", loc="left", fontsize=9.5, fontweight="bold")
    ax.legend(frameon=False, fontsize=7, loc="lower right", title="fragile admissible arms labelled", title_fontsize=6.5)

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--data", default="data"); ap.add_argument("--out", default="fig"); ap.add_argument("--doi", default="")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    g, sa, sv = load(a.data)
    # graphical abstract
    fig = plt.figure(figsize=(8.2, 7.2), dpi=300)
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.42, wspace=0.28, left=0.07, right=0.98, top=0.86, bottom=0.08)
    ax1 = fig.add_subplot(gs[0, :]); ax2 = fig.add_subplot(gs[1, 0]); ax3 = fig.add_subplot(gs[1, 1])
    panel_variants(ax1, g, sv); panel_gates(ax2, g); panel_stability(ax3, sa)
    fig.text(0.07, 0.955, "A configuration catalog for fvs-modern across 18 FVS variants", fontsize=14, fontweight="bold", color=GREEN)
    fig.text(0.07, 0.925, "4,748 FIA stands x 4 survival and parameter arms x 100 years, gated on the FIA envelope and Reineke self-thinning.", fontsize=8.6, color="#333333")
    fig.text(0.07, 0.903, "7 of 18 variants recommended; 0 of 72 verdicts change across 144 rule cells and 500 stand bootstraps.", fontsize=8.6, color="#333333")
    foot = "Weiskittel, A.R. (2026). CRSF, University of Maine. Catalog version 1.2.0, October 30, 2026."
    if a.doi: foot += f"  https://doi.org/{a.doi}"
    fig.text(0.07, 0.015, foot, fontsize=7, color=GREY)
    fig.savefig(os.path.join(a.out, "00_fvs-modern_catalog_graphical_abstract.png"), dpi=300); plt.close(fig)
    # standalone figures for the report and deck
    for fn, f, size in [("fig_variant_stability.png", panel_variants, (8, 3.4)), ("fig_gate_space.png", panel_gates, (5.2, 4.6)), ("fig_arm_stability.png", panel_stability, (5.2, 4.6))]:
        fig, ax = plt.subplots(figsize=size, dpi=300)
        (f(ax, g, sv) if f is panel_variants else f(ax, g) if f is panel_gates else f(ax, sa))
        ax.set_title(ax.get_title()[3:], loc="left", fontsize=9.5, fontweight="bold")
        fig.tight_layout(); fig.savefig(os.path.join(a.out, fn), dpi=300); plt.close(fig)
    print("wrote", os.listdir(a.out))

if __name__ == "__main__":
    main()
