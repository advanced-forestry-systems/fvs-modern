# Constrained projection results note (completes #89 stand layer)

Date: 2026-07-04. Branch: finish/stand-mort-survival (off feature/gompit-projection-wiring).

## What landed
All 6 shared modifier components (dg/hg/htdbh/hcb/cr/mortality) and all 4 García
state-space stand constraints (survival, BA-growth, top height, stems) are now landed
into the 25 CONUS variant configs via 62c. mort modifier + stand survival were the last
two (their 8h cloglog refits finished 2026-07-03 evening).

## BUG FOUND AND FIXED (load-bearing): stand-survival intercept refold
71_fit_stand_survival.R fit with an exposure offset log(YEARS)-3.9 but, unlike
76_fit_stand_stems.R, never refolded the intercept back by -3.9. The consumer
(stand_survival_eta in stand_constraint.py) reads the Intercept as the raw log annual
hazard, so the unrefolded value implied exp(-0.137)=0.87/yr baseline hazard => ~98.7%
stand mortality per 5-yr step. Corrected: Intercept -0.137 -> -4.037
(exp=1.77%/yr, matching the driver's own '~1.8%/yr' comment); baseline 5-yr stand
mortality now 8.4%. Fixed in the produced bundle (constant -3.9 shift on mean/q2.5/q97.5;
draws unchanged; original backed up as *.PREREFOLD.json) and in the 71 driver for
reproducibility.

## Coefficient sanity (all landed)
- mort modifier (cloglog log-hazard): management exp(-0.266)=0.77 (thinning lowers
  hazard), disturbance exp(+0.289)=1.33 (disturbance raises mortality). Sensible.
- stand survival: Intercept -4.04 (1.8%/yr), rd/ln_qmd/ba/bgi covariates finite.
- modifiers dg/hg/htdbh/hcb within 0.94-1.43; cr management x1.75 but sd 0.59 (CI spans 1).
- stand bagrowth/topht/stems as previously verified.

## Constrained vs unconstrained (conus_sf, synthetic NE stand, 100 yr, 400 draws)
Constrained tracks realistic stand targets; unconstrained tree-level equations drift to
degenerate values, and the four constraints rescue the trajectory. Year-100 medians:
- BA:    constrained 92.4 [79.3,101.5] m2/ha  vs  unconstrained 0.4 [0.1,1.6]
- TPH:   constrained 564 [533,569]           vs  unconstrained 6.7 [1.3,24.2]
- TOPHT: constrained 19.6 [19.0,25.3] m       vs  unconstrained 77.7 [63.1,92.6]
- QMD:   constrained 45.7 [42.4,49.2] cm      vs  unconstrained 27.2 [21.3,37.4]

## Honest caveats / follow-ups (not blocking the fix, but before manuscript use)
1. The harness renders only BA and TPH panels; top height and QMD are in the JSON but
   not plotted. Add the two panels for a full 4-panel figure.
2. Constrained top height is flat at ~20 m over 100 yr. Likely a synthetic-init artifact
   (init stand has QMD 16 cm but top height 20.7 m, internally slender). Re-run on a real
   FIA-derived NE init stand before drawing realism conclusions.
3. SDIMAX density ceiling used a fixed fallback of 600 (self-thinning posterior rds not
   passed). Wire --sdimax_samples so the ceiling carries the Bayesian posterior.
4. Constrained BA reaches ~92 m2/ha at year 100, high for NE; a Bakuzis realism pass
   across a site gradient (calibration/R/16_+17_) is the next gate.

Artifacts: constrained_vs_unconstrained_all4.png, constrained_projection_all4.json.
