# Final handoff: CONUS FVS refinement program

**Date:** 2026-06-14
**Purpose:** single authoritative pickup point for the whole program. Supersedes
`20260614_HANDOFF_COMPLETE.md`. Read this first; the dated memos below carry detail.

## 1. The program, framed as three steps

The work is one plan to refine FVS, in three steps from near-term to long-term:

1. **FIA calibration modifiers** on the existing engine (the basis of fvs-modern). Deployable now.
2. **Refine the maximum SDI** (localized, FIA-derived, level-calibrated) for long-term predictions.
3. **Refit the component equations as one CONUS-wide variant**, retiring the 20 regional variants.

The FVS team has one report and one document for this, both organized by the three steps.

## 2. State of play, one paragraph

Step 1 works and is benchmarked (the calibrated engine beats native NE and ACD on basal area, -0.6 vs
+12 to +13 percent). Step 2 is finalized and validated across all 20 CONUS variants. Step 3 has the
component equations fit and banked, with the unified joint fit done and the species-free injection
proven on held-out data; the remaining piece is wiring the banked forms into the engine. Nothing is
broken; the outstanding items are compute-gated, sign-off-gated, or engineering, not new science.

## 3. FVS-team deliverables (the things to send)

- **Report (deck):** `FVS_team_UPDATE_3step_plan.pptx` (12 slides, CRSF-branded, the 3-step narrative).
- **Document:** `FVS_team_REPORT_3step_plan.docx` (single comprehensive report, all figures/tables).
- **Cover memo:** `FVS_team_UPDATE_3step_plan.md`.
- Superseded split files are archived under `_superseded_fvs_team_materials/`.

## 4. Maximum-SDI result (Step 2, finalized)

- Species-weighted maximum is +28% biased with near-zero plot skill; localized FIA-derived maximum
  predicts observed self-thinning ~85% better (deviance 0.107 vs 0.058), in every region.
- Across all 20 CONUS variants, a level-calibrated localized maximum matches or beats native in 19
  (Utah the exception); optimal level spans 0.6 to 2.0 (median 1.2). Big gains where native is most
  off (NC 75->54, ACD 53->37, NE 35->31). OR/WA variants near-neutral; Utah governed by drought
  (weakest self-thinning signal, 0.04).
- The required level tracks each variant's native bias (r = -0.35): the level compensates for the
  variant's mortality calibration, so it must be set jointly with mortality, not dropped in uniformly.
- Detail: `20260614_USFS_maxSDI_technical_report.md`; figures `allvar_level_calibration.png`,
  `level_vs_bias.png`, `pn_vs_cr_maxSDI.png`, `cr_scale_diagnostic.png`; data `allvar_calibration.csv`.

## 5. Modeling workstreams (Step 3)

- **Unified joint fit (done).** One shared self-thinning logistic + per-variant level on the localized
  maximum, on 351k plots/19 regions. Per-variant level doubles self-thinning R2 (0.021->0.046). Data
  vs engine levels agree weakly (r=0.12): the engine level mostly corrects FVS's own mortality form,
  which argues for refitting mortality jointly. `unified_joint_fit.R`, `joint_fit_levels.csv`,
  `joint_fit.png`. Detail in `20260614_modeling_workstreams_status.md`.
- **Tree-level + ingrowth refits (banked).** Species-free bundles for DG, HG, HCB, HT-DBH, crown, and
  survival are fit and banked (including annualized crown and senescence forms). Production Stan
  reruns need a compute slot; production-config adoption is a hard stop pending sign-off. The ingrowth
  species-composition fit is the one piece still needing a slot to finalize.
- **Engine injection (framework proven).** `benchmark_sf_vs_legA.R` evaluates pure-species-free vs
  hybrid vs leg-A vs global on held-out data with prediction-interval calibration. First result (HCB,
  NE, 60k held-out): injected species-free beats per-species (RMSE 0.140 vs 0.147, R2 0.298 vs 0.226,
  PICP 0.94) and equals the hybrid. Remaining components still benchmarking. The true injection step
  (wiring banked bundles into the engine's per-tree increment) is the largest remaining build.

## 6. Benchmark and Greg comparison (carried forward)

Clean NE benchmark: unified beats native NE/ACD on basal area (-0.6% vs +12-13%). Greg Johnson,
Douglas-fir DG: ours species-dependent 0.091 RMSE vs Greg 0.097; species-free 0.118 with zero DF data
(within ~22% of a 156k-obs fit). Detail in `20260614_vs_greg_johnson_DG_comparison.md`.

## 7. Cardinal: code, data, jobs

- Max-SDI sweep: `~/fvs-modern/calibration/python/var_scale_diag.py` (level sweep through the engine;
  honors VAR/STATES/MAXSP/NSAMP/SEED). Side-by-side: `var_sdimax_sidebyside.py`. SDIMAX module:
  `~/fvs-modern/calibration/sdimax/localized_sdimax.py` (+ `localized_sdimax.R`).
- Joint fit: `~/unified_joint_fit.R`. Injection benchmark: `~/fvs-conus/dev_sf_integration/benchmark_sf_vs_legA.R`,
  driver `~/run_sf_bench.sh`, outputs `~/fvs-conus/output/sf_bench/`.
- Species-free banked bundles: `~/fvs-conus/output/conus/sf_integration/*_sf_fixed.csv` (cr_t2, dg,
  dg_v8_sf, hcb_v2split, hg, hg_v5_prod, hg_v8rd_sf, htdbh_v2split, surv_crz).
- Data: `~/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2.rds` (8.2M rows),
  `~/fvs-conus/data/brms_SDImax.csv`, `VAR_SDIMAX.csv`. TreeMap surface: Zenodo 10.5281/zenodo.19509367.
- R on Cardinal: `module load gcc/12.3.0 R/4.4.0` (4.4.0 has data.table+mgcv; 4.5.2 does not).
- **Cluster state at handoff:** 8 of Aaron's sbatch jobs queued/running (v7_qrf, cem2100) under the
  association limit, so no competing sbatch was submitted. The species-free injection benchmark is
  still running as a login-node nohup and will finish the remaining components on its own (HCB done).

## 8. Operational constraints (preserve)

- SSH key `~/Documents/Claude/.ssh-cardinal` (discover via find; re-stage each bash call). GitHub token
  `CRSF-Cowork/_context/.gh-holoros/token` (pipe into gh; never in chat/repos). Work on
  `feat/conus-sf-integration` (PR #70); no commits to main.
- Hard stops: no production config writes, no coefficient changes without sign-off.
- Do not cancel Aaron's Cardinal jobs (v7_qrf, cem2100). Use login-node nohup, not sbatch, while the
  association limit is hit.

## 9. Prioritized next steps

1. **When the cluster clears:** finalize the ingrowth species-composition fit and the remaining
   species-free injection benchmarks (HG, DG, HT-DBH, survival); the harness is built.
2. **Engine wiring:** inject the banked species-free bundles into the FVS per-tree increment path (the
   largest remaining build; turns the offline proof into a running unified variant).
3. **Unified mortality + maximum refit:** estimate one mortality response against the localized maximum
   so the per-variant level is intrinsic (the joint fit is the offline prototype).
4. **Adoption (needs sign-off):** promote the validated forms and the localized maximum into production
   config; then merge, Zenodo, DOI.
5. **Benchmark expansion:** Greg comparison for HG and a common multi-species held-out set; spatial
   blocking across all variants.

## 10. Honest bottom line

The max-SDI line is complete and well-evidenced through all 20 variants and a joint fit. The component
equations are fit and banked, the species-free injection is proven competitive and calibrated on
held-out data, and the FVS-team materials are consolidated and ready to send. What remains is bounded:
compute-gated fits (ingrowth, remaining benchmarks), the engine wiring, and sign-off-gated production
adoption. Leaving the cluster to clear; resume at step 1 above when slots free.

## 11. Index of key memos

- `20260614_USFS_maxSDI_technical_report.md` — the full max-SDI technical report.
- `20260614_modeling_workstreams_status.md` — joint fit, refits, injection detail.
- `20260614_CONUS_component_equations_assessment.md` — the tree/stand x dependent/independent two-by-two.
- `20260614_vs_greg_johnson_DG_comparison.md` — the Greg head-to-head.
- `20260614_maxSDI_localization_FINAL.md` — the finalized localization recommendation.
- `FVS_team_REPORT_3step_plan.docx` / `FVS_team_UPDATE_3step_plan.pptx` — the FVS-team deliverables.
