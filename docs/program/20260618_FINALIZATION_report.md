# Finalization of fvs-modern and fvs-conus: consolidated results

Date 2026-06-18. Branch holoros/fvs-modern conus-sf-integration-2026-05-21. Produced by an autonomous OODA
run that consolidates the committed result CSVs into one reproducible package (this report, the script
final_consolidate.R, and three headless figures). All metrics are out-of-sample where labeled, with
percentile-bootstrap 95 percent confidence intervals computed by resampling conditions.

## Summary

The program rests on a settled result and a validated adjustment layer. The widely cited FVS
over-prediction is mostly a disturbance-benchmark artifact, not a growth-equation bias, established three
independent ways. A keyword adjustment layer (brms site-specific max SDI, density-dependent recruitment,
per-species basal area increment) reduces size and density bias out-of-sample. The fvs-modern and
fvs-conus tracks are complementary in principle, but the per-species emulation of the fvs-conus growth
equations does not add a median gain over the density layer; a fair test of the true equations awaits the
in-engine injection (Route A), which is now diagnosed with a confirmed fix path.

## 1. The disturbance artifact (established, three ways)

| evidence | result |
|---|---|
| COND stratification, 22 variants | undisturbed median BA bias +1.8 percent; pooled +14; harvested +42 |
| removal-simulation converse test | harvested BA bias collapses to undisturbed when recorded harvest is simulated: ne +51.2 to +1.3, sn +110.1 to +0.1, pn +41.9 to +12.9 |
| fvs-conus projector, same conditions | same pattern; harvest inflation appears in the projector too |

Implication: validate on undisturbed plots or simulate the recorded removals before comparing.
Recalibrating growth equations downward does not help undisturbed prediction. Figure: disturbance_artifact_all25.png.

## 2. Four-arm comparison, FVS engine, out-of-sample (21 variants)

All four arms run in one framework so the projector-versus-engine sign difference does not confound the
comparison. Arm A default, B density layer (brms max SDI + density recruitment + BAIMULT), C per-species
BAIMULT emulating the fvs-conus growth equations, D combined.

| metric | A default | B density | C growth emul. | D combined | best |
|---|---|---|---|---|---|
| BA  | 6.1 | 4.9 | 6.7 | 5.6 | B |
| TPH | 16.3 | 9.3 | 17.2 | 9.3 | B and D tie |
| QMD | 12.5 | 5.9 | 15.1 | 5.7 | D (inside noise of B) |
| VOL | 12.0 | 10.7 | 13.0 | 11.9 | B |

Values are median absolute percent bias across variants. The density layer is the workhorse. The growth
emulation (C) does not help at the median and runs worse than default on BA, QMD, and VOL, so the combined
arm D tracks B and inherits little from C. The hypothesis that the growth signal stacks with the density
layer is not supported by these medians.

![Figure 1. Four-arm median absolute bias by metric, FVS engine, out-of-sample, 21 variants.](fig_fourarm_byarm.png)

Caveat: arm C is a per-species BAIMULT proxy for the fvs-conus equations, not the equations themselves; it
captures the diameter-growth level effect only. A fair arm C and D require the true equations in the
engine (Section 6).

## 3. Out-of-sample transfer of the size levers (held-out spatial fold)

Calibration derived on county-hashed fold A, applied unchanged to held-out fold B. QMD bias collapses
toward zero in every variant with adequate sample; CIs are 95 percent bootstrap.

| variant | QMD default | QMD calibrated | calibrated 95% CI |
|---|---|---|---|
| ne | +12.1 | +2.1 | [-2.3, +7.2] |
| sn | +11.3 | -2.1 | [-8.3, +4.4] |
| kt | +23.0 | -5.6 | [-13.5, +2.8] |
| pn | +17.5 | +1.7 | [-4.6, +8.2] |
| nc | +21.9 | -0.5 | [-14.0, +12.7] |
| ec | +23.4 | -5.8 | [-11.1, -0.3] |
| wc | +13.9 | -2.2 | [-8.0, +3.8] |
| cr | -1.0 | -15.4 | [-25.4, -4.2] |

cr (n small) over-corrects and is the documented exception; the cross-variant pattern is the robust
result. Density-dependent recruitment removed the earlier fixed-rate failure (the Southern out-of-sample
TPH over-correction of +16.1 became +1.4).

![Figure 2. QMD bias by variant, out-of-sample, default versus calibrated, with 95 percent bootstrap CIs.](fig_qmd_oos.png)

## 4. brms max SDI coverage per variant

Fraction of FIA plots in each variant's states carrying a brms site-specific max SDI, with the SDImax
distribution (English). State-level mapping conflates variants sharing states; treat as coverage, not a
per-variant rate.

| variant | match % | SDImax median [p10-p90] |
|---|---|---|
| acd | 73.4 | 386 [251-568] |
| ne | 63.2 | 370 [229-543] |
| sn | 41.1 | 337 [209-495] |
| ls | 28.4 | 322 [191-491] |
| ut | 18.5 | 265 [147-471] |

Lake States SDImax sits inside the cross-variant range, so the flagged LS anomaly is not a wild value; it
is flagged for the LS-specific calibration rerun. Conditions without a brms match fall back to FIA SICOND
or default variant SDIMAX.

## 5. Crown-width unification (Marshall)

MCW recovered from the CCF coefficients for the 15 western CCF variants plus AK; one CONUS-consistent
consensus curve MCW = B0 + B1 D was selected per species as the cross-variant median (44 species, 25
informed by more than one variant). Douglas-fir collapses from five disagreeing curves to one consensus
(4.71 + 1.499 D). Eastern variants carry no parabolic MCW (they zero the base coefficients); Russell and
Weiskittel 2010 is the recommended source there. Remaining step: emit the consensus coefficients into the
engine competition term and the fvs-conus crown component.

![Figure 3. Cross-variant maximum crown width spread at 20 in DBH by species; dot is the consensus mean.](fig_mcw_spread.png)

## 6. In-engine injection (Route A): diagnosed, fix path confirmed

The capability gap (in-process per-tree read and write) is closed: fvs2py now exposes get_tree_attr and
set_tree_attr, the variant library was rebuilt to clear the keyword-reader EOF, and in-process load plus
tree-attribute access are verified. The remaining segfault was root-caused by gdb to base/extree.f90 line
38, where the example-tree index array is unset because the in-process database read fails first
(SQLite "unrecognized token"). The confirmed fix is to load the stand through fvsAddTrees
(base/apisubs.f90 line 844) in memory, bypassing the in-process database read; this is the focused next
implementation step, with stand setup without a database as the one open unknown.

## 7. Finalization readiness

| track | state |
|---|---|
| Disturbance-aware benchmark | complete, three independent confirmations |
| Keyword adjustment layer (size/density) | validated out-of-sample with CIs |
| Density-dependent recruitment | prototype, transfers for most variants, PN excepted |
| Four-arm A/B/C/D | complete; growth emulation non-stacking documented |
| Crown-width consensus | complete for western variants; application step remains |
| In-engine true four-arm (Route A) | diagnosed, fvsAddTrees path confirmed, build pending |
| fvs-modern manuscript | disturbance-aware and four-arm sections integrated |
| fvs-conus manuscript v2 | pending (Section 3.7, Discussion, Figure 1) |
| Zenodo | deferred until validation is out-of-sample-honest, then new version of the concept DOI |

Reproduce: final_consolidate.R reads the committed CSVs and regenerates the tables and figures headlessly.

[DATA_STATE]: seven committed summary CSVs consolidated; 21-variant four-arm, 8-variant engine A/B with volume and CIs, NE projector within-framework deltas, 4-variant held-out, brms coverage for 15 variants, 44-species MCW consensus.
[OUTCOME_VERIFICATION]: median QMD |bias| 12.5 to 5.9 (B) out-of-sample; arm D 5.7 does not improve on B beyond noise; size levers transfer with bootstrap CIs in 7 of 8 variants (cr the small-n exception).
[IMPACT_UTILITY]: manuscript-preparer (fvs-modern report sections ready; fvs-conus v2 pending), data-curator and zenodo-deposit (defer DOI until OOS-honest; then new version of fvs_perseus_conus).
[NEXT_AUTONOMOUS_STEP]: build the fvsAddTrees in-memory stand loader to enable the true in-engine arms C and D, then re-run the four-arm and refresh this report.
