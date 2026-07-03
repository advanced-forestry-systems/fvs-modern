# What's left, the three-way comparison (default / keyword-calibrated / fvs-conus), and where crown width fits
2026-06-17

## The three-way comparison (it's already wired in)

Three approaches now exist and can be compared:

- DEFAULT FVS (per-variant species defaults).
- KEYWORD-CALIBRATED FVS (this work): brms site-specific max SDI + sign-aware ingrowth injection +
  signed BAIMULT + per-species HT-DBH correction. Validated disturbance-aware (COND-undisturbed),
  one cycle, FVS engine. Result: QMD 11.1->6.4, BA 9.9->7.9, TPH 19.2->9.1, volume 15.4->12.8 (median
  |bias|, in-sample).
- FVS-CONUS species-free equations (the Bayesian trait-driven refit): banked per variant; benchmark
  already computed in fvs-conus/output/ (allvar_calibration.csv; comparisons_overstory/manuscript_tables/
  fia_benchmark_results.csv; NEonly multi-metric tables, 50k+ conditions).

The fvs-conus benchmark (all conditions, condition-level, NOT disturbance-stratified) shows the key,
honest result:

| metric | default FVS | fvs-conus equations |
|---|---|---|
| BA bias (overall) | +13.0% | +17.3% |
| merch volume (CFNET) bias | +22.6% | +29.2% |
| top height bias | -2.9% | -1.7% |
| RMSE (allvar_calibration) | native | LOWER in ~most variants (e.g. NE 34.6->30.5, NC 75.2->54.0, ACD 52.8->36.9) |

Interpretation: the fvs-conus species-free equations REDUCE scatter (RMSE) and species inconsistency and
slightly improve height, but they do NOT reduce stand-level BA/volume BIAS - they make it modestly worse
on the current benchmark. This is exactly what the disturbance-aware reframing predicts: that benchmark
is not disturbance-stratified, so it inherits harvest inflation, and the equations alone do not touch
stand density (SDImax) or recruitment, which is where the bias lives.

Conclusion: the two efforts are COMPLEMENTARY, not competing. fvs-conus gives species-free consistency,
lower scatter, and better height; the keyword calibration (especially brms SDImax + ingrowth) addresses
the stand-level density/bias the equations leave untouched. The right product is fvs-conus equations
RUNNING UNDER the disturbance-aware calibration layer, benchmarked on the disturbance-clean basis.

## Caveat on comparability (must fix before claiming a clean 3-way)

The three are currently on different bases: keyword-calibrated = COND-undisturbed, plot-level, one cycle,
FVS engine; fvs-conus = all conditions, condition-level, tree-list projector (no disturbance split). A
clean three-way requires running all three on ONE basis - the disturbance-clean (COND-undisturbed) plot
set, same metrics (BA/TPH/QMD/volume), same horizon. That harmonization is the single most valuable
remaining analysis and the headline figure for the manuscript.

## Where crown width / MCW fits (Marshall note, 2026-06-17)

David Marshall's note recovers Maximum Crown Width (MCW) from each variant's FVS CCF coefficients:
MCW = sqrt(CCF / 0.001803026); for the linear form, B0 = sqrt(R1/k), B1 = sqrt(R3/k), k=0.001803026;
power and (R1+R2+R3)*DBH forms handled too. This matters because CCF/crown width drives the competition
term feeding diameter growth and mortality - the crown component of the calibration (the CRNMULT keyword
was the one defective lever). Key points: Douglas-fir has 4+ different MCW curves across 17 variants, so
a CONUS-wide MCW needs one equation per species (245 FVS / 471 FIA species; 232 to map by genus or
softwood/hardwood); Russell & Weiskittel (2010) Acadian MCW equations exist and should be used for the
Acadian region. Action: implement MCW recovery for all variants (trivial from the CCF coefficients),
assess the cross-variant spread per species, select/fit a CONUS-consistent MCW per species, and feed it
into both the fvs-conus crown component and the competition term of the benchmark. This is a clean,
self-contained next component.

## What's left (consolidated, prioritized)

A. HARDEN THE CALIBRATION (red-team roadmap, manuscript-blocking):
   1. Held-out spatial validation fold; relabel all numbers in-sample -> out-of-sample.
   2. Removal-simulation converse test (simulate recorded harvest on harvested plots -> show growth
      unbiased): the decisive proof of the disturbance-artifact headline.
   3. Volume in BOTH height configurations (measured vs imputed) + diameter/height decomposition;
      reconcile the HT-DBH narrative.
   4. True max size-density boundary (frontier/quantile) replacing p95; multi-cycle SDImax check.
   5. Master results table + bootstrap CIs; reconcile variant counts; explain the Lake States SDImax
      anomaly; per-variant brms match rate.
   6. Density-dependent recruitment form; multi-cycle stress test.
   7. brms SDImax model card (priors, diagnostics); volume-definition spec + sensitivity; biomass via FFE.

B. THREE-WAY HARMONIZATION (this session's finding):
   8. Re-run the fvs-conus benchmark on the disturbance-clean basis; produce the unified default /
      keyword-calibrated / fvs-conus / (combined) comparison, all metrics, same plots. Headline figure.
   9. Test the COMBINED product: fvs-conus equations + brms SDImax + ingrowth on the disturbance-clean
      basis - the hypothesis that they are complementary.

C. CROWN WIDTH / MCW (new component, Marshall):
   10. Implement MCW recovery for all variants; assess per-species cross-variant spread; build a
       CONUS-consistent per-species MCW (with R&W Acadian); wire into the crown/competition term.

D. REMAINING ENGINE/COVERAGE ITEMS:
   11. ACD/ADK regional benchmark (NE+customR) - the NY run stalled on I/O; finish it.
   12. fvs2py in-process injection (Product B engine path) - still maintainer-level; only needed to run
       the fvs-conus equations INSIDE the FVS engine rather than the standalone projector.
   13. CRNMULT keyword defect fix (crown multiplier inert) - isolated bug.

E. PUBLICATION / DEPOSIT:
   14. Manuscript: integrate the disturbance-aware section + the three-way comparison; reframe to
       "disturbance-aware benchmark + prototype adjustment layer + species-free equation comparison."
   15. Zenodo: DEFER until A1-A3 and B8 are done; then deposit the calibration config + validation
       dataset as a NEW VERSION of the existing fvs_perseus_conus concept DOI (new_version.py).

## Bottom line

The fvs-conus equations are wired in and compared: they cut scatter and fix height/species consistency
but not stand bias, while the keyword calibration fixes stand bias on disturbance-clean plots - so the
finished product combines them. The two highest-value next steps are the held-out + removal-simulation
validation (A1-A2) and the disturbance-clean three-way harmonization (B8-B9). Crown width (C10) is a
clean parallel component. Zenodo waits until the validation is honest.
