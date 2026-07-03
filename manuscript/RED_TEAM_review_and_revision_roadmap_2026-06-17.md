# Red-team review of the FVS recalibration manuscript + revision roadmap
2026-06-17. Independent adversarial review (subagent, FEM/Forest Science reviewer stance) plus disposition
and the corrected full stress-test table. Verdict: MAJOR REVISION bordering on reject-and-reframe. The
core empirical contribution (disturbance-stratified benchmark) is real and publishable; the calibration
claims outrun the evidence and must be reframed and re-validated out-of-sample.

## Full stress test (corrected, sign-aware, all 15 variants, COND-undisturbed, one cycle)

This supersedes earlier partial tables. Default -> calibrated bias (%); inj = ingrowth injected (sign-aware).

| variant | n | inj | BA def>cal | TPH def>cal | QMD def>cal | VOL def>cal |
|---|---|---|---|---|---|---|
| ne | 55 | 1 | 9.9>9.7 | -10.9>-4.1 | 11.2>4.7 | 23.9>22.7 |
| acd | 70 | 1 | 5.4>3.8 | -21.4>-16.2 | 11.1>6.4 | 15.3>12.8 |
| sn | 51 | 1 | 16.4>16.1 | -19.4>-6.5 | 22.3>12.1 | 18.4>15.8 |
| ls | 32 | 0 | 19.8>14.7 | 1.2>-7.8 | 8.9>10.3 | 22.4>19.1 |
| cs | 19 | 0 | 22.3>20.0 | 8.3>5.8 | 7.6>7.6 | 31.2>28.8 |
| ie | 26 | 1 | 3.8>2.2 | -23.6>-12.8 | 11.7>3.5 | -8.2>-10.5 |
| kt | 30 | 1 | 7.4>1.5 | -47.0>-25.2 | 28.9>6.0 | 15.4>4.8 |
| ci | 23 | 1 | 1.7>-1.4 | -25.4>-11.8 | 6.3>-3.7 | -10.4>-14.1 |
| cr | 8 | 0 | 33.2>26.7 | 46.9>38.0 | -10.2>-10.3 | 38.5>32.6 |
| ut | 8 | 1 | 13.4>12.1 | -19.2>-5.6 | 16.3>6.9 | 37.0>32.7 |
| ca | 24 | 0 | 9.5>6.0 | 7.2>2.0 | 9.8>10.0 | 6.0>3.6 |
| nc | 30 | 1 | 7.5>10.0 | -14.2>2.3 | 16.2>7.2 | -3.7>-0.8 |
| ec | 44 | 1 | -12.4>-1.1 | -34.3>-9.1 | 12.2>-4.4 | -16.4>-6.1 |
| wc | 44 | 1 | 5.9>3.0 | -15.8>-16.2 | 4.7>-0.9 | 1.0>-0.6 |
| pn | 44 | 1 | 12.3>7.9 | -16.0>-19.0 | 7.9>2.6 | 15.4>12.8 |

Median |bias| default -> calibrated: BA 9.9 -> 7.9, TPH 19.2 -> 9.1, QMD 11.1 -> 6.4, volume 15.4 -> 12.8.
Sign-aware injection (vs the earlier uniform version) is what improves the TPH median (19.2 -> 9.1 rather
than 15.8 -> 11.8). NOTE small n for cr (8), ut (8), cs (19); cr/ca/cs/ls were not injected.

## The five most serious issues (reviewer) and our disposition

1. IN-SAMPLE / TRAIN-TEST LEAKAGE (disqualifying as written). Every lever (ingrowth rate, BAIMULT,
   HT-DBH ratios) is derived from FIA and scored on the same FIA plots. The "stress test" is a
   goodness-of-fit, not a validation; a multiplier tuned to the mean residual will halve it by
   construction. DISPOSITION: ACCEPT. Required fix before any "validated" claim: spatially-blocked
   held-out fold (by EVALID/state or hex), derive levers on the calibration fold, report bias reduction
   only out-of-sample. All current numbers must be relabeled "in-sample / apparent." HIGH PRIORITY.

2. APPLES-TO-ORANGES HEIGHT COMPARISON. The HT-DBH critique blanks input heights (so FVS imputes), but
   the volume/BA validation uses FIA-provided heights (curve bias absent by construction). The two
   documents then contradict: htdbh says "HT-DBH inflates volume +5-10%, essential"; the capstone says
   "height is accurate, HT-DBH is NOT the volume lever." DISPOSITION: ACCEPT - both are true in their own
   configuration and we conflated them. Fix: state height configuration for every run; run volume
   validation in BOTH (DBH-only and height-supplied); decompose volume bias into diameter vs
   height-imputation on the same plots. Reconcile to: "when heights are supplied FVS self-calibrates and
   the curve bias is absent (TopHt +/-3%); in DBH-only inventories the curve introduces a size-dependent
   bias that inflates volume +5-10%." HIGH PRIORITY.

3. INTER-DOCUMENT NUMERICAL INCONSISTENCIES. Reviewer caught: (a) TPH/volume figures garbled in the
   manuscript section; (b) BAIMULT quoted as both 0.70 (shown to work) and 0.90 (deployed in config);
   (c) NE single-stand seed result (QMD +12->+3) vs all-variant run (+11.2->+4.7) quoted interchangeably;
   (d) the LAKE STATES SDImax anomaly - default 100-yr SDI 321 is BELOW the FIA limit 337 yet shows
   -26% BA, contradicting the "largest where default exceeds limit" mechanism. DISPOSITION: ACCEPT all.
   Fix: single master results table (this file's stress table is the start); regenerate every quoted
   number; reconcile BAIMULT (the deployable default is 0.90 conservative; 0.70 is the stronger setting
   tested - state both and which produced which number); EXPLAIN LS (likely the brms site SDImax for LS
   plots sits well below the p95=337 used in the long-term run, plus n is modest and the density-dependent
   mortality scales with SDImax so a stand below the limit still thins harder under a lower cap - but the
   magnitude needs a direct check). HIGH PRIORITY.

4. 100-YR SDImax LEVERAGE IS NEVER VALIDATED (model-vs-model only); p95 of current stand SDI is not the
   self-thinning limit; dense-stand selection (sdi0 >= 0.45 x maxSDI) and tiny n (ca 24, kt 33) inflate
   the effect. DISPOSITION: ACCEPT. Fix: (a) fit a true maximum size-density boundary (stochastic
   frontier / quantile regression on the QMD-TPH frontier) instead of p95; (b) validate trajectories
   against multi-cycle FIA panels where they exist; (c) drop the causal "is why FVS over-predicts on
   multi-decade projections" - we have no multi-decade observation. Reframe to "lowers the 100-yr
   asymptote; whether that is more correct is untested." MEDIUM-HIGH.

5. COND PROXY IS COARSE AND THE CONVERSE TEST IS MISSING. TRTCD=10 alone undercounts harvest; DSTRBCD
   misses sub-threshold disturbance; interval-timing of the codes vs the growth period is unstated;
   "undisturbed FIA plots" != "stands FVS is built to project" (which include managed stands). And the
   paper never runs the decisive converse test: SIMULATE the recorded removals on harvested plots and
   show FVS growth is unbiased there too. DISPOSITION: ACCEPT. Fix: broaden disturbance definition + a
   sensitivity analysis; state interval-timing handling; and actually run the removal-simulation test -
   that is the real proof that harvest, not growth, drives the pooled bias. MEDIUM-HIGH.

## Medium issues (accepted, for revision)

- brms match ~63%; 37% use the variant-median fallback, so SDImax is NOT plot-level for over a third of
  data - report with/without fallback and per-variant match rate.
- Small and unevenly reported n (cr=8, ut=8, ci=23); no confidence intervals anywhere - add bootstrap CIs
  to every bias estimate; directional tallies (18/19, 15/19) hinge on variants with n in the teens.
- ec/wc/pn share one OR+WA plot pool - three "variants" are not independent; fix all "X of N variants"
  tallies and the variant count (22 vs 19 vs 15 must be reconciled and the attrition explained).
- Ingrowth as fixed %/decade of INITIAL TPA is density-independent and will over-recruit over a rotation,
  and may fight the lowered SDImax - needs a density-dependent recruitment form and a multi-cycle test.
- Volume definition matching (FVS MCuFt vs FIA VOLCFNET) drove an earlier +52% -> +24% correction, so it
  is highly sensitive - specify top dib, stump ht, min merch DBH, sound vs gross, and give a sensitivity.
- Biomass is asserted throughout but never computed on the FVS side (needs FFE) - drop or run FFE.
- The brms SDImax model is a black-box backbone - report priors, predictors, R-hat/ESS, PPCs.

## Overclaims to soften (reviewer rewordings adopted)

- "fully calibrated FVS" / "validated" -> "a four-lever adjustment evaluated in-sample over one
  remeasurement interval."
- "the apparent over-prediction is largely a disturbance artifact" -> keep but add "provisional; we did
  not directly simulate the recorded removals" (until issue 5 is done).
- "that long-term over-stocking is why FVS over-predicts on multi-decade projections" -> "substituting an
  empirical SDImax lowers 100-yr BA up to 37%; whether this is more accurate is untested."
- "ingrowth injection closes the gap" -> "reduces the gap in-sample."
- "HT-DBH curve is biased" (stated as fact) -> reconcile with the height-configuration caveat above.
- "roughly self-replace" -> "net density change is near zero or positive in most undisturbed variants."

## Revision roadmap (ordered)

1. Held-out spatial validation fold; relabel all current numbers in-sample; re-report out-of-sample.
2. Removal-simulation converse test on harvested plots (the decisive proof for the headline).
3. Volume in both height configurations + diameter/height decomposition; reconcile the HT-DBH narrative.
4. True max size-density boundary (frontier/quantile) replacing p95; multi-cycle SDImax trajectory check.
5. Master results table + bootstrap CIs; reconcile variant counts and the LS anomaly; per-variant brms
   match rate.
6. Density-dependent recruitment form; multi-cycle stress test of the combined levers.
7. brms model card (priors/diagnostics); volume-definition spec + sensitivity; biomass via FFE or drop.
8. Reframe title/abstract: "A disturbance-aware benchmark and a prototype adjustment layer for FVS"
   rather than "national recalibration, fully calibrated."

## Bottom line

The disturbance-stratification insight is the keeper and is genuinely useful to the FVS community. The
adjustment layer is a promising prototype, not a finished calibration. Do not submit (or deposit) under
the current "fully calibrated / validated" framing; execute items 1-3 (held-out fold, removal-simulation,
height decomposition) and the paper becomes solid.
