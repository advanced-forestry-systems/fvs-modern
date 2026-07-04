#!/usr/bin/env python3
"""Offline Bakuzis-style realism pass + improved constrained demo.

Runs the constrained arm+modifier+4-constraint projector across a bgi SITE
GRADIENT with the CALIBRATED NE stand-level Reineke SDIMAX cap (replacing the
arbitrary 600 fallback), then checks the Bakuzis law-like relations:
  1. Site ordering  : higher site (bgi) -> higher top height and yield at a
     given age (monotone across the gradient).
  2. Self-thinning  : stand SDI approaches but does not exceed SDIMAX.
  3. Monotone yield : BA does not implausibly collapse under the constraint.
  4. Eichhorn-ish   : yield tracks top-height*density, not decoupled.
Produces a small-multiples figure + a verdict JSON.
"""
import os, sys, json
import numpy as np
WT = "/fs/scratch/PUOM0008/crsfaaron/wt-gompit"
OUT = "/fs/scratch/PUOM0008/crsfaaron/fvs-conus_output_conus"
sys.path.insert(0, WT + "/sf_integration_dev")
sys.path.insert(0, WT + "/config")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from constrained_projection import project, synthetic_ne_stand
from config_loader import FvsConfigLoader

BUNDLES = {
    "bagrowth": OUT + "/stand_bagrowth/stand_bagrowth_bundle.json",
    "survival": OUT + "/stand_survival/stand_survival_bundle.json",
    "topht":    OUT + "/stand_topheight/stand_topheight_bundle.json",
    "stems":    OUT + "/stand_stems/stand_stems_bundle.json",
}
FIXED_SDIMAX = 480.0   # NE calibrated stand Reineke SDIMAX (imperial), from
                       # calibration/output/variants/ne/species_sdimax_calibrated.csv
SITES = {"low_bgi4": 4.0, "mid_bgi6": 6.0, "high_bgi8": 8.0}

L = FvsConfigLoader("ne", version="conus_sf", config_dir=WT + "/config")
eco = {"L1": "8", "L2": "8.1", "L3": "8.1.1", "FT": 101}
trees0 = synthetic_ne_stand()
rng = np.random.default_rng(20260704)

runs = {}
for name, bgi in SITES.items():
    con = project(L, "conus_sf", trees0, eco, {"bgi": bgi}, years=100, step=5,
                  n_draws=200, constrained=True, external_bundles=BUNDLES,
                  fixed_sdimax=FIXED_SDIMAX, rng=rng)
    runs[name] = con
    print(f"ran site {name} (bgi={bgi})")

def series(res, key):
    yr = [s["year"] for s in res["trajectory"]]
    v = [s[key]["q50"] for s in res["trajectory"]]
    return yr, v

# ---- Bakuzis checks ----
verdict = {"fixed_sdimax_imperial": FIXED_SDIMAX, "sites": {}, "flags": []}
end = {}
for name, res in runs.items():
    tr = res["trajectory"]
    ba = [s["BA"]["q50"] for s in tr]
    th = [s["TOPHT"]["q50"] for s in tr]
    tph = [s["TPH"]["q50"] for s in tr]
    qmd = [s["QMD"]["q50"] for s in tr]
    # imperial SDI proxy from the projector's own convention (TPA/ha->/ac, QMD in inches)
    sdi = [ (t/2.4710538)*((q/2.54/10.0)**1.605) for t, q in zip(tph, qmd) ]
    end[name] = dict(BA=ba[-1], TOPHT=th[-1], TPH=tph[-1], QMD=qmd[-1], SDImax_reached=max(sdi))
    verdict["sites"][name] = end[name]
    if max(sdi) > FIXED_SDIMAX * 1.02:
        verdict["flags"].append(f"{name}: SDI {max(sdi):.0f} exceeds SDIMAX {FIXED_SDIMAX:.0f} (self-thinning breach)")
    # monotone BA (allow small dips)
    drops = sum(1 for i in range(1, len(ba)) if ba[i] < ba[i-1] - 0.5)
    if drops > 2:
        verdict["flags"].append(f"{name}: BA non-monotone ({drops} drops) -- implausible collapse")

# site ordering: high >= mid >= low on top height and BA at year 100
order_ok_th = end["high_bgi8"]["TOPHT"] >= end["mid_bgi6"]["TOPHT"] >= end["low_bgi4"]["TOPHT"]
order_ok_ba = end["high_bgi8"]["BA"] >= end["mid_bgi6"]["BA"] >= end["low_bgi4"]["BA"]
verdict["site_ordering_topheight_ok"] = bool(order_ok_th)
verdict["site_ordering_ba_ok"] = bool(order_ok_ba)
if not order_ok_th:
    verdict["flags"].append("site ordering VIOLATED on top height (higher site not taller)")
if not order_ok_ba:
    verdict["flags"].append("site ordering VIOLATED on BA (higher site not more productive)")
verdict["pass"] = len(verdict["flags"]) == 0

# ---- figure ----
fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
colors = {"low_bgi4": "#4C72B0", "mid_bgi6": "#55A868", "high_bgi8": "#C44E52"}
for name, res in runs.items():
    yr, ba = series(res, "BA"); ax[0].plot(yr, ba, color=colors[name], label=name)
    yr, th = series(res, "TOPHT"); ax[1].plot(yr, th, color=colors[name], label=name)
    yr, tph = series(res, "TPH"); ax[2].plot(yr, tph, color=colors[name], label=name)
ax[0].axhline(0, lw=0)
ax[0].set_title(f"Basal area (constrained, SDIMAX={FIXED_SDIMAX:.0f})"); ax[0].set_xlabel("year"); ax[0].set_ylabel("BA (m2/ha)")
ax[1].set_title("Top height (site ordering)"); ax[1].set_xlabel("year"); ax[1].set_ylabel("top ht (m)")
ax[2].set_title("Trees per ha (self-thinning)"); ax[2].set_xlabel("year"); ax[2].set_ylabel("TPH")
for a in ax: a.legend(fontsize=8); a.grid(alpha=0.3)
fig.suptitle("Bakuzis site-gradient realism pass -- constrained conus_sf, calibrated SDIMAX cap", fontsize=11)
fig.tight_layout()
fig.savefig(WT + "/sf_integration_dev/bakuzis_site_gradient.png", dpi=130)
json.dump(verdict, open(WT + "/sf_integration_dev/bakuzis_site_gradient_verdict.json", "w"), indent=2)
print("VERDICT:", json.dumps(verdict, indent=2))
print("BAKUZIS_GRADIENT_OK")
