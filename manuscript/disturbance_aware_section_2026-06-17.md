# Disturbance-aware validation and a four-lever calibration of FVS (draft section for the national recalibration manuscript)
2026-06-17. Intended for insertion into manuscript/fvs_combined_draft.md (Results + Discussion).

## Methods addition: disturbance-stratified FIA benchmark

We evaluated default and calibrated FVS against FIA remeasurement plots, stratifying every plot by its
recorded management/disturbance history from the FIA COND table: harvested (any cutting treatment,
TRTCD = 10), naturally disturbed (any DSTRBCD > 0), or undisturbed (neither). Because the default FVS
projection is undisturbed by construction (it does not cut), comparisons were made on the undisturbed
stratum unless otherwise noted; pooled comparisons that ignore disturbance are reported only to expose
the artifact they create. Bias is reported as 100 x (predicted - observed) / observed for basal area
(BA), trees per hectare (TPH), quadratic mean diameter (QMD), and merchantable cubic volume (FVS MCuFt
vs FIA VOLCFNET, matched definitions). Maximum stand density index was supplied per plot from a Bayesian
(brms) site-specific SDImax model rather than the FVS species defaults.

## Result 1: the apparent FVS over-prediction is largely a disturbance artifact

Across 22 variants the pooled FVS basal-area bias was +14 to +28%, the value usually cited as "FVS
over-prediction." Stratifying by COND history dissolves most of it: harvested plots carry a +33 to +55%
apparent bias purely because the default run never removes the harvested trees, while on undisturbed
plots the median basal-area bias is small (a few percent in most regions). Calibrating growth downward to
erase the pooled bias would therefore be fitting a harvest signal and would push the undisturbed
projection into under-prediction. The first requirement for any FVS calibration study is a
disturbance-aware benchmark.

## Result 2: on undisturbed plots the residual error is structured, not a uniform over-prediction

Using the unbiased COND-undisturbed stratum (19 variants, n >= 15 each), three coherent signals emerge:
basal area median +7%, QMD median +10% (over-predicted in 18 of 19 variants), and TPH median -14%
(under-predicted in 15 of 19). FVS carries too few, too-large trees. Top height was within +/-3% in most
variants, which masks the underlying error.

## Result 3: the density deficit is missing recruitment

Decomposing the undisturbed TPH change shows real stands recruit 17 to 70% of initial density per decade
and roughly self-replace (ingrowth approximately balances mortality), whereas default FVS only loses
stems (net negative everywhere) because its establishment model is disturbance-triggered and adds no
background ingrowth in undisturbed projections. The per-variant ingrowth rate predicts the TPH bias
(Lake States, lowest ingrowth, is the one variant FVS tracks; Kootenai and Southern, highest, are the
worst). Injecting a recruitment cohort at the observed per-variant rate closes the gap (e.g. NE TPH
-11 -> -1, QMD +12 -> +3); the standard ESTAB keyword does not, because it requires a disturbance.

## Result 4: maximum SDI governs the long-term trajectory

Over a single remeasurement interval SDImax has little leverage, but over a 100-year projection it is the
dominant control. Replacing the FVS default species SDImax with the site-specific brms maximum reduces
100-year basal area by 5 to 37%, largest where the FVS default exceeds the FIA-observed self-thinning
limit (Kootenai default 100-yr SDI 448 vs observed 392, -37%; Inland California 505 vs 459, -35%;
PNW Coast 619 vs 473, -20%). Several western variants permit long-term over-stocking beyond the
empirical limit. The brms SDImax is a plot-level, site-driven quantity, not species-specific.

## Result 5: volume and biomass require an HT-DBH correction that top height hides

Merchantable volume is over-predicted +15 to +24% in the East/Lake States/PNW, driven by the QMD
over-prediction (volume scales with DBH^2). Height growth and top height are well calibrated (TopHt bias
+/-3%), but the HT-DBH curve, evaluated tree-by-tree against FIA measured heights, is biased in a
variant- and size-specific way: pooled it runs from -7% at 1-3 in to +9% at 19-40 in. Inland Empire and
Kootenai over-predict height across all sizes (+10 to +35%), Southern under-predicts large trees, and
PNW has a shape error (under small, over large). In DBH-only inventories, where FVS imputes height, this
inflates merchantable volume +5 to +10%. Eighty-two per-species correction ratios were derived; because
several variants show a shape error, refitting the curve coefficients is preferred over a uniform
multiplier. REGHMULT (height growth) is not the right lever for the static imputation bias.

## Result 6: the four-lever adjustment, evaluated in-sample over one interval

Combining the levers (site-specific brms SDImax; sign-aware recruitment injection at the per-variant FIA
rate; a signed diameter-growth multiplier; and an HT-DBH curve correction) and validating default vs
calibrated across all variants on the undisturbed stratum reduces the QMD bias (median |bias| ~11% ->
~4%) and the TPH bias wherever it was under-predicted, while leaving the already-good variants alone (the
recruitment injection is applied only where TPH is under-predicted). Across 15 variants the median
absolute bias falls from 11.1 to 6.4% for QMD, 9.9 to 7.9% for basal area, 19.2 to 9.1% for TPH (the
sign-aware injection is what produces the TPH gain), and 15.4 to 12.8% for merchantable volume (full
four-metric per-variant table: calib_final.csv). These are IN-SAMPLE figures: the adjustment levers are
derived from FIA and evaluated on the same plots over a single ~5-10 yr remeasurement interval, so they
are apparent (goodness-of-fit) reductions and an upper bound on out-of-sample performance.

## Discussion framing

The study reframes the long-standing claim that FVS over-predicts growth. On the stands FVS is built to
project, undisturbed and self-replacing, its growth is close; the headline over-prediction in naive FIA
comparisons is mostly disturbance the model does not simulate. The genuine, calibratable errors are
structural and decompose by scale: maximum SDI sets the multi-decade density ceiling, recruitment sets
decadal stem density, the diameter-growth multiplier sets standing size and therefore volume, and the
HT-DBH curve sets the height on which volume and biomass depend. All four are grounded in FIA and the
brms SDImax model, and all four deploy through existing FVS machinery (SDIMAX, treelist recruitment,
BAIMULT, and HT-DBH coefficients). The Acadian and Adirondacks sub-variants run under the Northeast
engine through the customR interface, and the calibration layer applies to them unchanged because the
brms SDImax is keyed by plot rather than species.

## Caveats

This is a prototype adjustment layer, not a finished, out-of-sample-validated calibration. Key
limitations, to be resolved before submission: (1) Levers are derived from and evaluated on the same FIA
plots; reported reductions are in-sample and require a spatially-blocked held-out fold for an honest
out-of-sample estimate. (2) The volume/BA validation uses FIA-supplied heights (where FVS self-calibrates
the HT-DBH curve), while the HT-DBH bias is shown only with heights blanked; the two configurations must
be reported together and the volume bias decomposed into diameter and height-imputation components. (3)
The 100-yr SDImax leverage is an FVS-vs-FVS simulation, not validated against multi-decade observations,
and the revised SDImax is a p95 of current stand SDI rather than a fitted maximum size-density boundary;
the causal "over-stocking is why FVS over-predicts long-term" is therefore softened to "lowers the 100-yr
asymptote." (4) The disturbance-artifact claim should be confirmed by the converse test - simulating the
recorded removals on harvested plots and showing FVS growth is then unbiased - which is not yet done. (5)
COND under-counts disturbance (TRTCD=10 only, DSTRBCD above threshold only); brms SDImax is plot-level for
~63% of plots and variant-median otherwise; ec/wc/pn share one OR+WA plot pool (not independent); several
variants have small n (cr, ut = 8); no confidence intervals are yet reported; biomass is not yet computed
on the FVS side (needs FFE); and the fixed-rate ingrowth form is density-independent and untested over a
rotation. See the red-team review and revision roadmap (20260617_RED_TEAM_review_and_revision_roadmap.md)
for the ordered fixes.
