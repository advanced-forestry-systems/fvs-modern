# Pre-registration: CONUS model run-off on identical FIA stands

Written October 30, 2026, before any full-run result. The pilot (section 9) was used only to check that each contestant runs, that units and schema agree, and to time the work. Sections 1 to 8 are locked. Any later change goes in section 11 with the date and the reason, before the affected results are looked at.

Root: `/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/runoff/` (written `runoff/` below). Catalog: `/fs/scratch/PUOM0008/crsfaaron/bakuzis_v2_20261027/catalog/` (written `catalog/`, read only).

## 1. Question

For each of the 18 FVS variants that have FIA fitting data (ACD, BM, CA, CI, CR, CS, EC, EM, IE, LS, NE, PN, SN, SO, TT, UT, WC, WS), which model configuration gives (a) admissible 100-year behavior and (b) the most accurate one-remeasurement projection, when every contestant projects the same FIA stands? The October 27 projector scores (`proj_rescore2/`, 1,500 stands per variant, seed 20261027) were drawn from a different stand sample than the engine catalog. This run-off removes that confound.

## 2. Contestants and exact configurations

### A. fvs-modern engine (four arms)

- Libraries: `/fs/scratch/PUOM0008/crsfaaron/wt-bakuzis-main/lib/FVS<var>.so` (NE library md5 prefix 27ed5279dae0), driven by fvs2py (`wt-bakuzis-main/deployment/fvs2py`) through `calibration/python/bakuzis_uncertainty_comparison.run_scenario`.
- Parameter sets: `default` (no calibration keywords) and `calibrated` (`config/calibrated/<var>.json` through `FvsConfigLoader`).
- Survival: `native` (FVS variant mortality, `FVS_GOMPIT` unset) and `gompit` (`FVS_GOMPIT=1`, `FVS_GOMPIT_COEF=wt-bakuzis-main/config/greg_mortality_coefficients.csv`).
- Arms: `default_native`, `default_gompit`, `calibrated_native`, `calibrated_gompit`.
- Long-term: not rerun. The trajectories are the catalog runs_v5 output (`catalog/runs_v5/`, concatenated in `catalog/catalog_long_v5.csv`, script `catalog/conus_catalog_runs.py`, array `catalog/catalog_runs_v6.slurm`).
- Near-term: new runs with `runoff/scripts/runoff_engine_runs.py --mode near`.

### B. fvs-conus projector (four arms)

- Script: `/fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/conus_eq_projector_v2.R`. This is the version that produced the Paper 2 Table 6 inputs rescored in `proj_rescore2/` (`fvs_stress/conus_eq_proj/conus_eq_array.sbatch` and the two GOMPIT array scripts call v2).
- Arms, named by the projector's CONFIG label:
  - `conus_b1`: `--mode=free`, the species-free b1 diameter growth, fvs-conus logit mortality.
  - `conus_b2`: `--mode=dependent`, the v8 species-aware b2 diameter growth, logit mortality.
  - `conus_b1_gompit` and `conus_b2_gompit`: the same two growth models with `--mort=gompit --cch_module=conus_eq_proj/cch_module.R`, using the coefficients in `/fs/scratch/PUOM0008/crsfaaron/conus_mort/full_out/greg_mortality_coefficients.csv`.
- Here "native" mortality means the fvs-conus logit model, not FVS variant mortality.
- Ingrowth is off: `--ingrowth=runoff/b_inputs/ingrowth_zero.rds`, the banked lookup with every rate set to 0 (`runoff/scripts/make_zero_ingrowth.R`). This matches the engines' no-regeneration protocol.
- The projector's SDImax self-thinning ramp (its FIX 3) stays on, because it is part of the model.
- Site covariates come from the nearest pairs plot by public PLOT coordinates, inside v2.
- Driver: `runoff/scripts/runoff_b_project.sh`. Post-processing: `runoff/scripts/runoff_b_post.R`.

### C. USDA FVS, unmodified (one arm, `usfs_native`)

- Source: USDAForestService/ForestVegetationSimulator at commit b183c10f (2025-11-06). This is the fork base recorded in `wt-bakuzis-main/.upstream-sync.json`. The volume/NVEL submodule (FMSC VolumeLibrary) is at its pinned commit 24705360.
- The clone is `runoff/upstream_build/usfs_fvs`. The build copy is `runoff/upstream_build/src_b183c10f`. No source file was edited.
- Build: the upstream `bin/makefile` (`make FVS<var>.so`). Its FFLAGS and CFLAGS were overridden to the flags fvs-modern uses (`wt-bakuzis-main/deployment/scripts/build_fvs_libraries.sh` lines 68 and 69), so compiler settings match contestant A. Upstream's default FFLAGS add `-ffpe-trap=invalid,zero,underflow,overflow,denormal`, which would abort runs on floating-point underflow.
- Libraries: `runoff/upstream_build/lib/FVS<var>.so` for 17 variants. Build logs are `runoff/upstream_build/build_<var>.log` and `runoff/logs/upstream_build_16012719.out`.
- ACD: upstream ships `acd/` source but no `bin/FVSacd_sourceList.txt`. A list derived the way fvs-modern derives its own fails to compile: `acd/badist.f` uses `BAU`, which upstream's commons do not declare (`runoff/logs/upstream_build_acd_16012803.out`). So C on the ACD stands runs on the upstream NE library (`--engine-variant ne`), the variant a USDA user in Maine would run.
- Parameters: default, with no calibration keywords. Mortality: FVS variant mortality. The fvs2py runner and keyword generators are the same as for A.

## 3. Stand sample

- The catalog sample, unchanged: `catalog/conus_catalog_runs.py::sample_plots` on `catalog/fia_variant_assignment.csv` with `--n-stands 300 --seed 20261028`. ACD uses NE-assigned plots in ME, NH and VT.
- The plot CN is recovered from the public PLOT table by (STATECD, COUNTYCD, PLOT, INVYR1), because the pairs PLT_CN is a double-cast denormal.
- Stand construction is the same for every contestant:
  - Condition trees only, live trees only.
  - Null TPA_UNADJ records are dropped.
  - TPA is divided by CONDPROP_UNADJ.
  - The traj id is C%04d by sample order.
- Contestant B receives exactly these trees (`runoff_engine_runs.py --b-inputs`) as original FIA SPCD, DBH in inches, height in feet, CR in percent and TREE_COUNT.
- The sample has 4,748 stands: ACD 300, BM 300, CA 142, CI 300, CR 300, CS 300, EC 300, EM 300, IE 300, LS 300, NE 300, PN 38, SN 300, SO 300, TT 142, UT 300, WC 226, WS 300 (`catalog/catalog_long_v5.csv`, default_native arm). The catalog v4 and v5 runs share this sample. v5 is used because all 4,748 stands pass the year-0 seed gate in every engine arm (`catalog/catalog_seedgate_v5.csv`), against 4,485 in v4.

## 4. Long-term protocol

- 100 years at a 5-year step (20 cycles). No regeneration keywords.
- `NOAUTOES` is appended in EM and IE (`FVS_NOAUTOES=1`) for engines A and C. The projector's ingrowth is off.
- No management.
- A stand whose last tree dies ends its trajectory (rows with TPH = 0 are dropped for every contestant). B's count of such stands is reported.
- Common volume for every contestant: VOL = 0.40 × BAPH × HT. HT is top height, defined by the FVS ATTOPHT rule: the TPA-weighted mean height of the 40 largest-DBH trees per acre. For the engines it comes from FVS ATTOPHT. For B it is computed from the projector tree lists. The engines' native total cubic volume is kept as VOL_NATIVE and reported but not scored.

## 5. Near-term protocol (Paper 2 convention)

- Each plot is projected in a single cycle whose length is its own remeasurement interval: k = round(YEARS) from the pairs file, or INVYR2 − INVYR1 when YEARS is missing. The engines use NUMCYCLE 1 and TIMEINT 0 k. The projector uses `--ncycles=1 --cyclelen=k`, one call per distinct k.
- The observations are in `runoff/obs/nearterm_obs_frame.csv` (script `runoff/scripts/runoff_nearterm_obs.R` v4, log `runoff/logs/obs_16012840.out`). They come from `/users/PUOM0008/crsfaaron/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds`, restricted to the catalog's (plot_key, INVYR1).
- Expansion is TPH_UNADJ1 / CONDPROP_UNADJ_cond1. DBH is in inches and HT in metres; both were checked against the file's stand columns in `runoff/logs/obsdiag_16012770.out`. Recomputed time-1 BA matches the file's BA1 within 2% on 91.5% of plots.
- Primary target: time-2 survivor metrics (trees live at time 1 and live at time 2). The pairs file has no ingrowth rows, and no contestant regenerates. Recomputed survivor BA2 is at or below the file's BA2 on 98.7% of plots, so the file's BA2 includes ingrowth.
- Sensitivity target: the file's BA2 and TPH2, which count all trees live at time 2.
- Exclusions, fixed now:
  - Plots with no time-2 record.
  - Plots where more than 5% of time-1 BA has unknown time-2 fate.
- Across the 17,769-plot frame, 14,931 plots are usable (obs_ok). 2,536 are excluded for unknown fate and 302 have no time 2.

## 6. Metrics

### Near-term

For each contestant arm and variant, on BA (m² ha⁻¹), TPH (ha⁻¹), QMD (cm) and VOLC (m³ ha⁻¹):

- bias = mean(pred − obs)
- percent bias = 100 × mean(pred − obs) / mean(obs)
- RMSE

Each gets a stand-bootstrap 95% percentile interval (2,000 resamples, seed 20261028) from `runoff/scripts/runoff_score_near.R`. The primary set is "common": stands that ran in every arm of the variant with an engine seed gate of OK. The "all" set is a sensitivity analysis. The year-0 check reports engine year-0 BA ÷ BALIVE and year-0 BA ÷ observed BA1.

### Long-term

All long-term metrics use the existing scorer and gates, unchanged:

- `bakuzis_stress_v2.R`: the A and C long table is built by `runoff/scripts/runoff_long_prep.R` (`--source=long`), and B is scored with `--source=projector` with HT_DOM_M added. The scorer reports Bakuzis relations T0 to T9 with evaluated, not-evaluated and flagged counts. Site and density relations are not evaluated because the sample has a single site and density level.
- `catalog/catalog_gates.R` with `catalog/fia_variant_thresholds_plot_acd.csv` reports, by arm and variant:
  - share of stands above FIA BA p99.9 at 50 and 100 years
  - share above SDI p99.9 at 50 and 100 years
  - Reineke flag rate, evaluated count and median slope
  - median peak BA against FIA BA p50
  - the `admissible` verdict

## 7. Decision rule for "best per variant"

Fixed before results (`runoff/scripts/runoff_decide.R`):

1. **Admissible.** The admissible set is the arms for which `catalog_gates.R` returns admissible = TRUE on the run-off stands. That requires all of:
   - share of max BA over 0 to 50 years above FIA BA p99.9 is under 5%, with a Wilson upper bound under 10%
   - the same for peak SDI against SDI p99.9
   - Reineke flag rate among RD-gated evaluated trajectories under 10%
   - median peak BA from 0 to 100 years at or above FIA BA p50
2. **Most accurate.** Among admissible arms, the best has the smallest absolute near-term BA percent bias on the common stand set.
3. **Tie.** An admissible arm ties the best if the paired stand-bootstrap 95% interval of (|pbias_arm| − |pbias_best|) contains 0 (2,000 resamples, seed 20261028). Among tied arms, the simplest wins, in this order:
   `usfs_native` < `default_native` < `default_gompit` < `calibrated_native` < `calibrated_gompit` < `conus_b1` < `conus_b1_gompit` < `conus_b2` < `conus_b2_gompit`.
   This is the catalog's order (default < calibrated < conus_b1 < conus_b2; native < gompit), with unmodified USDA FVS first.
4. **No admissible arm.** The variant is NOT RECOMMENDED for long-horizon use. The best near-term arm is named for single-cycle use only.

TPH, QMD and VOLC biases, RMSE and the sensitivity target are reported but do not enter the rule.

## 8. Fairness controls and known asymmetries

### Controls

- Same stands, built by one function.
- Same year-0 trees. In the pilot, the median ratio of year-0 BA to BALIVE was 1.001 for the engines and 1.00003 for B.
- Same thresholds file, scorer, gate script and decision script.
- Same unit conversions: TPA→TPH 2.4710538, ft²/ac→m²/ha 0.2295684, ft→m 0.3048. QMD is computed from BA and TPH for every contestant.
- Same top-height rule and the same form-factor volume.
- Regeneration is off everywhere.
- The engines use the same fvs2py runner, keyfile generators, species crosswalk and relabel rule for A and C.
- The same compiler flags for A and C.

### Known asymmetries (reported, not corrected)

- Missing tree heights are filled by FVS dubbing in A and C and by the fvs-conus height-diameter model in B.
- Species outside the variant's list are relabeled to 299, 998 or 999 in the engines. B keeps the original SPCD and borrows the stand's modal in-set species for growth.
- "Native" mortality differs by family: FVS variant mortality in A and C, fvs-conus logit in B.
- B caps diameter growth at 2.0 in per year and applies its SDImax ramp.
- FVS cycle lengths are whole years, so k is rounded for every contestant.
- C on ACD stands uses the upstream NE library, because the upstream ACD source does not build.

## 9. Pilot (10 NE and 10 PN stands, catalog traj C0001 to C0010)

These are the first 10 stands of the catalog sample. Near-term n is 8 for NE and 9 for PN, because the other stands have no usable time 2.

- Engine runs: `runoff/pilot/logs/eng_16012843_*.out`.
- Projector runs: `runoff/pilot/logs/b_16012845_*.out`.
- Scores:
  - `runoff/pilot/score/near/near_scores.csv` and `near_year0_check.csv`
  - `runoff/pilot/score/gates/catalog_gates_by_arm.csv`
  - `runoff/pilot/score/{eng,b}/pilot_*_scorecard_summary.csv`

Results:

- **Runs.** All stands ran in every contestant. There were no RUN_FAIL or SEED_FAIL statuses, and the seed gate was OK on 100% of stands.
- **Schema and units.** Every contestant's long_metric file has arm, group, traj, age, TPH, BAPH, QMD, HT, VOL and SDI in metric units. Year-0 medians in NE are BAPH 27.1 m²/ha for the engines and 27.0 for B, and TPH 1,575 against 1,576.
- **C against A default_native.** In a long-term screen of the first 10 stands of all 18 variants (`runoff/pilot/c_screen/c_vs_a_default_native_first10.csv`), C is bit-identical to A default_native over 100 years in 16 variants (share of identical BA rows 1.000, max |ΔBA| 0). It differs in SN (77.6% of rows identical, max |ΔBA| 2.53 m²/ha) and ACD (expected, since C runs the upstream NE library). Near-term results for C and A default_native are identical in NE and PN.
- **B with ingrowth off.** Some stands lose every tree by year 100: 1 of 10 in PN in each arm, and 1 of 10 in NE in each GOMPIT arm (`runoff/pilot/logs/rescore_16012884.out`).
- **Admissibility.** With 10 stands, no arm can be admissible: a Wilson upper bound under 10% needs about 35 or more stands with zero exceedances. The pilot gates are a schema check only.

## 10. Expected compute, from pilot timings

Pilot timings: `runoff/pilot/logs/time_*.txt` and `cscreen_16012906_*.out`.

| Step | Tasks | Pilot time | Full run estimate |
|---|---|---|---|
| A near | 36 (18 variants × 2 survival arms, both parameter sets per task) | 15 to 20 s per 10 stands | Under 1 h per task |
| C long | 18 | 15 to 56 s per 10 stands | Up to 2.5 h per task |
| C near | 18 | Similar to A near | Under 1 h per task |
| B | 144 (18 × 4 arms × long and near) | Long 8 to 10 s, near 23 to 32 s per 10 stands; max RSS 2.7 GB | Under 30 min each |

- The C long upper bound comes from the catalog: runs_v5 took 1.5 to 2.5 h per 300-stand task with two configs (`sacct` for job 15998817). The 10-stand C screen took 15 to 56 s per variant, including reading the state FIA tables.
- Requested wall times: 6 h for engines, 6 h for B, 4 h for scoring.
- Scoring: one job, under 4 h.
- Upper bound: about 250 core-hours. That is B at 144 tasks × 0.5 h × 2 cores = 144, C long at 18 × 2.5 h = 45, and A near plus C near at 54 × 1 h = 54. The pilot rates suggest the real figure is far smaller (tens of core-hours).

## 11. Change log

- 2026-10-30: pre-registration written. Observation script v4 (survivor target, unknown-fate exclusion) was fixed before any contestant prediction was compared with observations.

## 12. Amendment A1 (2026-10-30 06:39 EDT, before any full-run submission; decisions from the coordinator, delegated by Aaron)

1. C runs in full on all 18 variants. On ACD stands its arm is labelled `usfs_native(ne_lib)` (upstream NE library).
2. B runs both projector versions as separate arms, `conus_v2_*` (conus_eq_projector_v2.R) and `conus_v6_*` (conus_eq_projector_v6.R; called with --horizon/--cyclelen/--treelists, ingrowth off through the same zero lookup). Pilot (NE, PN, 10 stands, `runoff/pilot/logs/b26_*`, `v6cmp_16013107.out`): v6 runs and seeds bit-identical year-0 trees to v2, so both arms are kept. v6 behaviour flagged: its stand-level BA envelope pulls stands down sharply (pilot near-term mean BA NE 25.75 -> about 15.0 to 15.8 m2/ha over about 5.5 yr; PN 70.4 -> 16.9 to 21.1 over about 10 yr). B task list: 288 tasks (18 x 4 arms x long/near x v2/v6).
3. B ingrowth off (primary); no ingrowth-on arm.
4. Simplicity order: usfs_native = usfs_native(ne_lib) < default_native < default_gompit < calibrated_native < calibrated_gompit < conus_v2_b1 < conus_v2_b1_gompit < conus_v2_b2 < conus_v2_b2_gompit < conus_v6_b1 < conus_v6_b1_gompit < conus_v6_b2 < conus_v6_b2_gompit.
5. Near-term target: survivors primary; the file's BA2/TPH2 (ingrowth included) as sensitivity (metrics BA_total, TPH_total in near_scores.csv).
6. Long single cycles: no cap. FVS single cycles of 13, 14, 17 and 19 yr ran with years_run = k and plausible increments in A (default, calibrated) and C (`runoff/pilot/logs/longk_16013055.out`, `longk2_16013095.out`). Only 51 of 4,222 sample stands with usable time 2 have k > 12.
7. Calibrated-parameter path of A (gate check): `run_scenario` uses BOTH the keyword overlay from `config_loader.generate_keywords` (BAIMULT and MORTMULT from the precomputed calibration_multipliers, plus SDIMAX where present) AND the in-process fvs2py `apply_calibrated_config` / `apply_to_fvs` (FVS(config_version="calibrated") at stop point 7), which in this config sets only `spsdi` (13 of 18 variants) and nothing in ACD, BM, IE, PN, SO (`runoff/logs/calpath_16013046.out`, `runoff/pilot/logs/longk_16013055.out`). Per the coordinator's gate, the full run was NOT submitted pending a decision on this path. The same path produced the catalog long-term calibrated arms.

## 13. Amendment A2 (2026-10-30 07:20 EDT, before full-run submission; coordinator decision)

1. Calibrated-path gate resolved as Option 1: launch with the path as is. Rationale (coordinator, independent audit): the keyword overlay carries the full calibrated block (NE: 58 BAIMULT and 84 MORTMULT values different from 1 after the run); the in-process fvs2py overlay adds only SDImax (spsdi); this keeps the run-off consistent with the shipped catalog. Open audit item D2 (in-process spsdi zero fill for NA species) is checked separately in runoff/d2check/.
2. The conus_v6_* arms are excluded from the decision rule. They are run, scored and reported, labelled "diagnostic: stand BA cap suspected miscalibrated" in the decision table. The winner is chosen among usfs_native / usfs_native(ne_lib), the four A arms and the four conus_v2 arms (runoff/scripts/runoff_decide.R, md5 prefix e1f113960f0c).

## 14. Amendment A3 (2026-10-30 11:24 EDT, after the full run finished and before any rescoring; coordinator decision)

1. fvs-modern main d08b8b96 fixes the D2 zero fill (in-process spsdi set to 0 for NA species; material in CA, CR, CS, EM, LS, WC per runoff/d2check/d2_compare_by_variant.csv) and three other defects. Fixed worktree: /fs/scratch/PUOM0008/crsfaaron/wt-fix-20261030.
2. Primary analysis: A's calibrated near-term arms (calibrated_native, calibrated_gompit) are rerun with the fixed worktree (FVS_PROJECT_ROOT, FVS_LIB_DIR, FVS_CONFIG_DIR, PYTHONPATH -> wt-fix-20261030) into runoff/full/engines_fix/ (array job 16017891, 36 tasks) and replace the in-process-path calibrated near-term rows. A default near-term rows, C and B are unchanged.
3. A long-term arms come from catalog/catalog_long_v6.csv (calibrated arms regenerated with the fixed code by another agent; default arms reused from v5) instead of catalog_long_v5.csv.
4. The Option 1 results (in-process path; full/a_near calibrated rows and catalog_long_v5 calibrated arms) are kept and reported as a labelled sensitivity (arms suffixed _inproc), never entering the decision rule.
5. Decision rule, metrics, thresholds, simplicity order and the v6-diagnostic rule (A2) are unchanged. Rescoring output: runoff/full/score_fix/ (primary) and runoff/full/score_fix_sens/ (sensitivity).
