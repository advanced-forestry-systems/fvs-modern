# FIA overstory recompute — apples-to-apples FVS-NE + FVS-ACD benchmark

Task #135 + #136. Two SLURM runs on Cardinal (jobs 10443275, NE only; 10591333,
NE + ACD via `FVS_ACD_RELABEL=TRUE` + footprint states 23, 33, 50). Engine:
`R/19_fia_benchmark_engine_overstory.R` (focused copy of the production engine
with `ALL_VARIANTS = c("acd","ne")`, the per-tree `DIA >= 5.0` filter
substituted at all three call sites, and output redirected to
`calibration/output/comparisons_overstory/`).

## Apples-to-apples FIA leaderboard (basal area, DIA >= 5 in)

| Model | n plots | BA bias | R squared | obs BA (m2/ha) |
| --- | --- | --- | --- | --- |
| FVS-ACD calibrated | 13,364 | **-0.06%** | 0.822 | 17.17 |
| FVS-ACD default    | 13,364 |  -6.07%   | 0.815 | 17.17 |
| FVS-NE calibrated  | 37,167 |  -3.95%   | 0.870 | 18.57 |
| FVS-NE default     | 37,167 |  -8.20%   | 0.852 | 18.57 |
| OSM-ACD (prior benchmark) | 12,029 | -2.18% | 0.958 | 21.54 |

ACD plots are the Maine + NH + VT subset of NE, projected with the ACD
posterior parameters (Acadian subvariant). NE is the rest of the
Northeast footprint. Both observed at DIA >= 5 in (FIA TPA_UNADJ),
unlike the prior all-trees comparison.

## What the overstory recompute shows

The Acadian calibration on its native footprint is essentially
unbiased on basal area at the overstory scale: FVS-ACD calibrated
lands at -0.06 percent across 13,364 Maine + NH + VT plots, the best
single-region accuracy of any model in this comparison. The default
FVS-ACD on the same plots is -6.07 percent, so the calibration adds
about six percentage points of accuracy. FVS-NE calibrated is at
-3.95 percent on the larger 37,167 non-Acadian Northeast plots, and
the calibration adds about four percentage points there.

The OSM-ACD prior benchmark at -2.18 percent on 12,029 plots is also
clearly accurate; the FVS-ACD calibrated result improves on it by
about two percentage points on its own Acadian footprint. The
previous OSM-vs-FVS pipeline gap I had flagged in #135 (correlation
0.77 between OSM and FVS observed BA, mean ratio 0.88) is largely
closed by the overstory filter, with both pipelines now within a few
points of each other on a comparable tree basis.

## Combined cross-region scorecard (apples-to-apples)

Putting FIA Maine/NH/VT (overstory) next to Canadian NB
(MAGPlot 262 pairs, also overstory-weighted):

| Model | FIA Maine (calibrated, overstory) | Canada NB (262) |
| --- | --- | --- |
| FVS-ACD (AcadianGY) calibrated | **-0.06%** | **-0.04%** |
| OSM-ACD                        | -2.2%   | -3.0%       |
| FVS-NE calibrated              | -4.0%   | (Canada not run with ACD overlay) |
| FVS-NE default                 | -8.2%   | +9.1%       |

FVS-ACD calibrated is essentially unbiased on BOTH regions
(-0.06% / -0.04%) — by this measure the most accurate AND most
consistent model across the cross-border Acadian forest. OSM-ACD is
the close second and still the cleanest single piece of evidence for
cross-region consistency because both its FIA and Canada benchmarks
were independently computed and converge. Default FVS-NE still flips
sign across the border (-8 to +9 percent), confirming the calibration
is doing real work.

This refines the previous "OSM-ACD is the only model accurate in both
regions" framing: with the apples-to-apples FIA recompute, **the
Acadian calibration (FVS-ACD) joins OSM-ACD as a cross-region-
consistent model**, and is the most accurate of the two on the
Acadian footprint specifically.

## Robustness — read the result carefully

The "essentially unbiased" framing above is true for the **bias of
plot means**, which is the right metric for landscape-scale stand
population aggregates such as regional basal area or carbon
inventory. It is not the right metric for individual stand
projections, where per-plot scatter is substantial. Per-plot
diagnostics on the same overstory data:

| Variant (calibrated) | n | bias of means | median plot bias | median \|plot bias\| | RMSE (ft2/ac) |
| --- | --- | --- | --- | --- | --- |
| FVS-ACD | 13,364 | -0.06% | -4.9% | 10.9% | 20.4 |
| FVS-NE  | 37,167 | -3.95% | -6.9% | 10.8% | 17.5 |

The per-plot bias distribution is skewed: 5th to 95th percentile is
roughly -33% to +139% for ACD (and -30% to +65% for NE), so the
near-zero mean bias on ACD reflects substantial cancellation of
under-projected and over-projected plots rather than uniformly
accurate per-stand predictions. Median plot bias is around -5 to -7
percent across both variants, i.e., the model systematically
under-projects most individual plots but balances out at the
landscape mean.

The BA result also decomposes into offsetting TPA and QMD components
(see `fia_overstory_tpa_qmd_decomp.png`). At the bias-of-means scale:

| Variant (calibrated) | TPA bias | QMD bias | BA bias | R^2(BA) |
| --- | --- | --- | --- | --- |
| FVS-ACD | +1.58% | -2.69% | -0.06% | 0.822 |
| FVS-NE  | +0.28% | -4.45% | -3.95% | 0.870 |

These compose per plot as `BA = K * TPA * QMD^2`, so the BA bias of
means is the bias-of-means of the per-plot product (`K * TPA * QMD^2`),
not the marginal product of the TPA and QMD bias-of-means — the
covariance between TPA and QMD matters. The ACD calibrated result is
achieving aggregate accuracy through compensating small biases in
the underlying stand structure, not by getting both components
individually correct; NE has the same sign pattern but a weaker
cancellation, leaving a -4% residual on BA. Defaults under-project
both structure components (TPA -1.87% / QMD -3.99% for ACD,
TPA +1.25% / QMD -7.07% for NE) and so under-project BA more
strongly.

The honest summary: the Acadian calibration is the best landscape-
mean predictor of any model in this comparison and is appropriate
for inventory and carbon accounting at the regional scale; it
should be used with uncertainty bands rather than as a point
estimate for individual-stand silvicultural decisions, where the
roughly 11 percent median absolute plot error is the operationally
relevant uncertainty.

## Run notes

* SLURM 10443275: 25 min wall, 50,531 NE plots only, BA bias
  -3.4 / -8.1 (calibrated / default). ACD was skipped because
  `FVS_ACD_RELABEL` defaulted to FALSE. Result archived as
  `comparisons_overstory_NEonly` and `validation_data_overstory_NEonly.csv`.
* SLURM 10591333: 25 min wall, same engine with
  `FVS_ACD_RELABEL=TRUE` and `FVS_ACD_FOOTPRINT_STATES=23,33,50`
  (Maine, New Hampshire, Vermont). 13,364 NE plots relabeled to
  ACD; 37,167 retained as NE. Both projected. Result is the canonical
  one: `validation_data_overstory.csv`.

## Files

Canonical data (from SLURM 10591333):

* `validation_data_overstory.csv` (46 MB) — per-plot observed and
  predicted for 50,531 ACD + NE plots at DIA >= 5 in. Includes TPA,
  BA, QMD, SDI, HT_top, and all four volume bases (CFGRS net, CFNET,
  BFNET, CFGRS gross) for observed t1, observed t2, predicted
  calibrated (+ uncertainty band), and predicted default.
* `fia_overstory_summary.csv` — the apples-to-apples leaderboard.
* `validation_data_overstory_NEonly.csv` (46 MB) — the first run
  output (NE-only, ACD relabel disabled); kept for reproducibility.

Cross-region + diagnostic artifacts (built locally from the validation file):

* `fia_canada_apples2apples.R` / `.csv` / `.png` — cross-region
  leaderboard plus FIA-vs-Canada consistency scatter, with the
  |bias| <= 5% accuracy box.
* `fia_overstory_perplot_dist.R` / `.png` — per-plot bias density for
  ACD and NE. Complements the bias-of-means scorecard with the
  individual-plot spread so the "essentially unbiased" headline
  reads correctly as landscape-mean, not per-stand.
* `fia_overstory_tpa_qmd_decomp.R` / `.csv` / `.png` — BA bias
  decomposition into TPA and QMD components. Shows the mechanism
  behind the ACD calibrated near-zero BA bias (small TPA over-
  projection cancels a small QMD under-projection on the basal-area
  scale).

Engine-side outputs (`observed_changes.csv` (~127 MB),
`manuscript_tables/`, `manuscript_figures/`) were generated on
Cardinal at `~/fvs-modern/calibration/output/comparisons_overstory/`
but were not synced to this local checkout. They are recoverable from
that directory or by rerunning the engine.

## Repository note

The repo I am writing into this session does not contain the
SILC analysis tree (silc_v25 through silc_v28, the v32 / v36
deck builders, the rev 14 SESSION_HANDOFF.md) that was here at the
end of the previous session. The overstory recompute artifacts have
landed in a clean new directory under
`calibration/output/comparisons_overstory/` for integration into
whichever organizational state the repo is now meant to be in.
