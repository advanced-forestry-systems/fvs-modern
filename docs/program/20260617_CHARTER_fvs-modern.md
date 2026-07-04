# fvs-modern — program charter and roadmap
2026-06-17. Track 1 of the two-track FVS modernization program. Companion: CHARTER_fvs-conus and
INTEGRATION_roadmap.

## One-sentence scope

Refine and modernize the existing USDA FVS codebase by replacing its 1970s-1990s parameters with a
Bayesian recalibration of all component equations and variant-specific parameters (notably maximum SDI),
fit to contemporary national FIA remeasurement data, with full uncertainty quantification, improving
predictions across species and region over CONUS while keeping the operational FVS engine and workflow
intact.

## Why it matters / positioning

fvs-modern is the deployable, agency-compatible track. It does not change FVS's equation FORMS or its
keyword/engine interface, so it drops into existing FVS workflows; it changes the NUMBERS (posterior
medians + intervals) and a small set of structural levers that the engine already exposes. It is the
fastest path to better, uncertainty-aware FVS predictions in operational use, and it is the validation
and benchmarking backbone for the more ambitious fvs-conus track.

## Components calibrated (per CALIBRATION.md)

Seven Bayesian component models per variant: diameter growth (Wykoff ln DDS), height-diameter
(Chapman-Richards), height increment (where applicable), mortality (annualized logistic), crown-ratio
change, stand density (SDIMAX), and self-thinning slope. Three runtime options: default, calibrated
(national FIA posteriors), custom (user data).

## What this session established (the adjustment-layer findings)

- DISTURBANCE-AWARE BENCHMARK (the methodological keystone): FIA-vs-FVS comparisons must be stratified by
  FIA COND treatment/disturbance. The naive pooled "+14-28% over-prediction" is mostly unsimulated
  harvest; on undisturbed plots the structured error is QMD over (+10%), TPH under (-14%), volume over
  (+15-24%), basal area mild (+7%).
- FOUR DEPLOYABLE LEVERS, by time scale: (1) max SDI from the brms site-specific (plot-level) model —
  the dominant long-term control (100-yr BA -5 to -37% where FVS defaults exceed the empirical limit);
  (2) recruitment/ingrowth injection at the per-variant FIA rate (sign-aware) — FVS adds no background
  ingrowth undisturbed; (3) signed BAIMULT for diameter growth / standing volume; (4) per-species,
  size-aware HT-DBH correction — worth 5-10% of volume in DBH-only inventories, hidden by unbiased top
  height. ACD/ADK run as NE sub-variants via customR; the layer applies unchanged (brms SDImax is keyed
  by plot).
- Full sign-aware stress test (15 variants, COND-undisturbed, one cycle, in-sample): median |bias|
  default->calibrated BA 9.9->7.9, TPH 19.2->9.1, QMD 11.1->6.4, volume 15.4->12.8.

## What's left (manuscript-blocking, from the red-team)

1. Out-of-sample validation: spatially-blocked held-out fold; relabel all numbers in-sample -> OOS.
   (Held-out validation executed this session for ne/sn/kt/pn — see held_out.csv.)
2. Removal-simulation converse test: simulate recorded harvest on harvested plots; show growth unbiased.
3. Volume in both height configs (measured vs imputed) + diameter/height decomposition; reconcile HT-DBH.
4. True maximum size-density boundary (stochastic frontier / quantile regression) replacing the p95;
   multi-cycle SDImax trajectory check.
5. Master results table + bootstrap CIs; reconcile variant counts; explain the Lake States SDImax
   anomaly; per-variant brms match rate (currently ~63% for NE, fallback otherwise).
6. Density-dependent recruitment form; multi-cycle (rotation) stress test of the combined levers.
7. brms SDImax model card (priors, predictors, R-hat/ESS, PPC); volume-definition spec + sensitivity;
   biomass on the FVS side via FFE.
8. CRNMULT keyword defect (crown multiplier inert) — isolated engine bug to fix.

## Deliverables / artifacts

- Repo: holoros/fvs-modern, branch conus-sf-integration-2026-05-21. CALIBRATION.md (pipeline +
  three options), DISTURBANCE_AWARE_VALIDATION.md (protocol + four levers), diagnostics_2026-06-16/
  (all validation scripts, per-variant config, brms SDImax, figures, calib_final.csv), manuscript/
  (draft + disturbance-aware section + red-team review + this charter).
- Manuscript pipeline #1: "National recalibration of FVS using Bayesian methods" (drafting).
- Existing Zenodo concept: fvs_perseus_conus (software). Deposit the calibration config + validation
  dataset as a new VERSION at submission.

## Recommended next session (fvs-modern)

Execute red-team items 1-3 (held-out fold extension to all variants; removal-simulation test; volume
height-config decomposition), then build the master results table with CIs and reframe the manuscript.
