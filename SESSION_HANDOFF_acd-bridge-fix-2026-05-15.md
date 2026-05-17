# Session Handoff — acd-bridge-fix-2026-05-15

Date written: 2026-05-16
Author of session: Aaron Weiskittel (working with Claude assistance)
Branch: acd-bridge-fix-2026-05-15
Remote: github.com:holoros/fvs-modern.git
Latest commit at handoff: 44d1f0a (NE vs ACD runtime A/B) plus pending commits described below.

## TL;DR

ACD was a sub-variant of NE that behaved oddly when called. After
~15 commits this branch:

1. Fixes the ACD calibration post-pass (parse bug, NE param
   fallback, ingrowth_lookup cache).
2. Defaults NSVB volume/biomass/carbon equations to on across all
   seven Eastern variants (ACD, NE, CS, LS, SN, KT, EM) via the
   LFIANVB flag in grinit.f90.
3. Fixes seven F77 to F90 conversion bugs in the standalone FVS
   build path that were blocking native binary execution.
4. Delivers three runtime A/B comparisons proving the fork carries
   genuine variant divergence at three resolution levels: dispatch,
   estimator, and calibration.
5. Confirms with Bayesian posteriors and FIA benchmarks that NE
   currently outperforms ACD calibration on the Acadian footprint,
   with a clear fix path identified.

The fork is now in a state where calibrated ACD and calibrated NE
projections can be A/B tested at runtime against the same FIA input
deck, and the build pipeline is solid enough to support future
variant-level recalibration work without further F77/F90 hazard.

## Branch lineage

```
44d1f0a NE vs ACD runtime A/B: confirm subvariant divergence
cc00b92 NSVB vs CRM runtime A/B: definitive divergence demonstrated
26cc488 spctrn.f90: handle ACD variant in species fallback + guard JSPIN(0)
c75e576 build_fvs_executables.sh: shield variant PRGPRM/ESPARM from base/ shadows
6c029bd filopn.f90: open ISTDAT with <basename>.tre in non-interactive path
6a27c2b filopn.f90: init FVS unit numbers before OPEN — unblocks keyrdr.f90
3db614a fvs.f90: initialize FVS unit-number COMMON variables on entry
1cd784f errgro.f90: defensive JOSTND guard to avoid UNFORMATTED I/O abort
e0aeeec autopilot 2: extend NSVB to KT/EM, fix exe build, document divergence
3625a01 Extend NSVB default to NE / CS / LS / SN variants
254cce3 ACD: default FIA NSVB volume/biomass/carbon equations to ON
86a7110 ACD stand-level calibration factors + bias diagnostic suite
a8ee4e9 ACD post-pass: cache ingrowth_lookup so full-fidelity ingrowth is reusable
f5bb1bf ACD subvariant bridge: fix post-pass parse + relabel + NE param fallback
```

Plus pending commits at handoff time:
- sumout.f90: deactivate dead CS/LS/NE/SN CASE branches (fixes FVSne summary writer crash)
- calibrated_ne_vs_acd directory: SLURM script, CSV, README
- Updated ne_vs_acd README reflecting fixed FVSne output

## Key technical findings

### NSVB vs CRM (commit cc00b92)

Two FVSacd executables built from identical source except for the
LFIANVB default produce different .sum / .out / FVSOut.db. CFV
differs at cycle 0 (1530 NSVB vs 1526 CRM), BFV at cycle 0 is the
headline (292 NSVB vs 0 CRM because CRM applies stricter merch-spec
thresholds). TPA, BA, QMD remain bit-identical, confirming the
divergence is isolated to the volume estimator path. Pattern matches
Westfall et al. 2024 NSVB equations vs the Component Ratio Method.

### NE vs ACD runtime (commit 44d1f0a, then revised with sumout fix)

Stand S248112 UNTHINNED CONTROL at year 2090 (age 160):

| Metric | NE | ACD | Delta |
| --- | --- | --- | --- |
| TPA | 111 | 94 | NE +17% |
| BA  | 194 | 169 | NE +15% |
| QMD | 17.9 | 18.2 | ACD +1.7% |
| CFV | 7,638 | 6,727 | NE +13% |
| BFV | 43,258 | 38,540 | NE +12% |

ACD produces a lower-stocked, higher-QMD stand at maturity, the
regional spruce-fir signature expected from the Acadian variant. NE
keeps more trees; ACD has higher mortality but larger residuals.

### Calibrated NE vs Calibrated ACD on Acadian-footprint FIA

SLURM job 9718410 with FVS_ACD_RELABEL=TRUE. 30,146 ME/NH/VT
conditions (13,693 ACD + 16,453 NE post relabel).

| Variant | BA %RMSE | BA bias_pct |
| --- | --- | --- |
| ACD_acadian | 28.52 | +51.56 |
| NE_acadian  | 23.19 | +27.32 |

NE outperforms ACD on the Acadian footprint by 5.3 pp %RMSE and
~24 pp less positive bias. Driver: ACD diameter-growth HMC fit did
not converge (rhat 1.95, converged FALSE). Fix path: re-fit with
longer warmup, smaller step size, adaptive reparameterization.

All 2,261 common posterior variables between NE and ACD differ.
The posteriors are not aliases of each other; ACD just needs better
mixing.

## Files modified or added

### Source code (src-converted/)

- acd/grinit.f90, ne/grinit.f90, cs/grinit.f90, ls/grinit.f90,
  sn/grinit.f90, kt/grinit.f90, em/grinit.f90 — LFIANVB = .TRUE.
- base/errgro.f90 — defensive JOSTND guard
- base/fvs.f90 — unit-number init on entry
- base/filopn.f90 — unit-number init before OPEN, ISTDAT tree-file
  open in non-interactive path
- base/varver_stub.f90 — new stub for VARVER subroutine
- vls/spctrn.f90 — AC alias + JSPIN(0) guard
- vbase/sumout.f90 — deactivate dead CS/LS/NE/SN CASE branches (NEW)

### Build pipeline

- deployment/scripts/build_fvs_executables.sh — conditional
  econ_stubs inclusion, MAXSP shadow shield for variant builds

### Calibration pipeline

- calibration/R/19_fia_benchmark_engine.R — FVS_ACD_RELABEL env gate,
  ACD param fallback to NE, ingrowth_lookup cache
- calibration/R/27_acd_post_pass.R — marker-based slicing, fallback,
  ingrowth cache reading
- calibration/slurm/submit_benchmark_acd_relabel.sh — NEW

### Analysis artifacts (calibration/analysis/acd_stand_level_2026-05-16/)

- ne_vs_acd/ — runtime NE vs ACD outputs + README
- nsvb_runtime_ab/ — runtime NSVB vs CRM outputs + README
- calibrated_ne_vs_acd/ — calibrated NE vs ACD FIA benchmark CSV +
  SLURM script + README (NEW)
- ACD_bias_findings.md — diagnostic narrative
- stand_level_calibration.R — 81-row stand-level calibration factors
- acd_calibration_factors.csv + guidance

## How to resume

1. Pull the branch on Cardinal: `cd ~/fvs-modern-acdbridge && git
   pull origin acd-bridge-fix-2026-05-15`.
2. Rebuild executables: `bash deployment/scripts/build_fvs_executables.sh
   . lib acd` and `... lib ne`.
3. Smoke-test: run the four-corner integration test described in
   `calibration/analysis/acd_stand_level_2026-05-16/ne_vs_acd/README.md`.
4. Open PR or merge into main when ready.

## Open follow-ups (priority order)

1. Re-fit ACD diameter-growth HMC with longer warmup and smaller
   step size. Target: rhat < 1.05. The Stan model file is at
   calibration/stan/wykoff_dg/.
2. Widen ACD footprint to include New York Adirondack counties via
   county-level filter in 19_fia_benchmark_engine.R. FIPS county
   list documented in calibrated_ne_vs_acd/README.md.
3. Repeat calibrated A/B with the 81-row ACD post-pass calibration
   factors applied to see whether the stand-level corrections close
   the bias gap.
4. STOP 10 at the very end of standalone FVS runs — post-completion
   IC code, not a crash. Investigation: main.f90 SELECT CASE on
   fvsGetICCode. Low priority; doesn't affect output.
5. Consider merging acd-bridge-fix-2026-05-15 into main after the
   ACD HMC re-fit completes.

## Autopilot round 4 — 2026-05-16 (later)

Three follow-ups from this handoff's prioritized list landed in
commit 8510b6f:

1. **HMC sampler is env-tunable.** 02c_fit_dg_hmc_small.R now
   reads `FVS_HMC_WARMUP`, `FVS_HMC_SAMPLING`, `FVS_HMC_ADAPT_DELTA`,
   `FVS_HMC_TREEDEPTH`, `FVS_HMC_CHAINS`. The `project_root` resolver
   is now robust to non-interactive Rscript invocation (the prior
   sys.frame(1) failure mode is fixed).

2. **NY Adirondack footprint added.** 19_fia_benchmark_engine.R
   reads `FVS_ACD_NY_COUNTIES` (default: FIPS COUNTYCD 19, 31, 33,
   35, 41, 49, 89, 115 — the 8-county Adirondack Park footprint)
   and joins COUNTYCD into the plots lookup. When relabel is on, NY
   plots in those counties are now retagged as ACD alongside ME/NH/VT.

3. **ACD post-pass factors wired in.** Setting `FVS_ACD_POSTPASS=TRUE`
   in 19_fia_benchmark_engine.R applies population-level multipliers
   (BAPH 1.0168, QMD 1.0071, TOPHT 1.0117, TPA 1.0147) to calibrated
   ACD predictions before metrics. Multipliers individually
   overridable via FVS_ACD_POSTPASS_BAPH/QMD/TOPHT/TPA. Per-stratum
   factors (the 81-row table) are still a follow-up.

4. **Longer ACD HMC re-fit running.** SLURM job 9733825 launched at
   2026-05-16 20:14 UTC with warmup=2000, sampling=1500,
   adapt_delta=0.99, max_treedepth=12, max_obs=10000. Target:
   rhat < 1.05 for diameter_growth_posterior. Expected runtime
   ~4 to 6 hours. Result will land at
   calibration/output/variants/acd/diameter_growth_posterior.csv;
   prior posterior is snapshotted as
   diameter_growth_posterior.csv.refit_pre_*.

5. **STOP 10 at end of standalone runs is expected.** errgro.f90
   sets ICCODE=1 whenever a non-fatal note is logged (e.g.,
   "INPUT SPECIES CODE (WH) WAS SET TO (HI)" species substitutions
   in net01.key). main.f90 translates that to STOP 10 at exit. Not
   a bug. Closed.

## Next round priorities

1. Wait for HMC job 9733825 to finish; verify rhat < 1.05; if not,
   try adapt_delta=0.995 and warmup=3000.
2. Once ACD posterior is converged, rerun the calibrated A/B with
   FVS_ACD_RELABEL=TRUE and the new posteriors. Expect ACD %RMSE to
   drop closer to NE's 23.19.
3. Same run with `FVS_ACD_POSTPASS=TRUE` to see the combined effect
   of converged posterior + population multipliers.
4. With NY-county filter (`FVS_ACD_NY_COUNTIES` default) the next
   relabel run will pull in additional Adirondack plots —
   re-baseline the Acadian-footprint metrics.
5. Build out per-stratum FVS_ACD_POSTPASS variant that looks up the
   correct row in acd_calibration_factors.csv by FT_GROUP / SI_class
   / BA_t1_class / interval_class. Estimated +10 lines of dplyr.

## Autopilot round 5 — 2026-05-16/17

Round 4 closed several items but left an unconverged HMC and only
the population-level post-pass. Round 5 closes the loop.

1. **HMC re-fit budget right-sized.** Job 9733825 (warmup=2000,
   sampling=1500, adapt_delta=0.99, max_treedepth=12) hit the 6-hour
   walltime mid-sampling and timed out. The original ACD posterior
   was preserved (snapshot at
   `diameter_growth_posterior.csv.refit_pre_20260516_201255`).
   `refit_acd_dg.sh` updated to warmup=1500, sampling=1000,
   adapt_delta=0.98, max_treedepth=11, walltime=12h. Resubmitted as
   job **9812192** (RUNNING).

2. **Per-stratum post-pass implemented.** `19_fia_benchmark_engine.R`
   gained a stratified mode via `FVS_ACD_POSTPASS_MODE=stratified`.
   Each ACD-tagged row gets four candidate multipliers (one per
   stratum: FT_GROUP, BA_t1_class, SI_tercile, interval_years) and
   the geometric mean across the four is applied. The lookup reads
   `acd_calibration_factors.csv`. Population mode remains the
   default; stratified is opt-in.

3. **SLURM dependency chain for three back-to-back A/Bs.**
   `run_ab_after_hmc.sh` runs three benchmark passes once HMC
   completes successfully, each writing a tagged
   `fia_benchmark_pctrmse_*.csv`:

   - `refit_only`: converged ACD posterior, no post-pass
   - `refit_postpass_pop`: + population multipliers
   - `refit_postpass_strat_ny`: + stratified post-pass + NY Adirondack
     county filter (the full stack)

   Submitted as job **9812377** with `--dependency=afterok:9812192`.

## What lands automatically once HMC 9812192 returns

Job 9812377 fires three A/Bs back-to-back (~2 hours total) and
saves six tagged CSVs under
`calibration/output/comparisons/manuscript_tables/`:

```
fia_benchmark_pctrmse_refit_only.csv
fia_benchmark_pctrmse_refit_postpass_pop.csv
fia_benchmark_pctrmse_refit_postpass_strat_ny.csv
fia_benchmark_results_refit_only.csv
fia_benchmark_results_refit_postpass_pop.csv
fia_benchmark_results_refit_postpass_strat_ny.csv
```

The headline comparison is the ACD row %RMSE across these three
files vs the pre-refit ACD baseline (28.52%).

## Next round priorities

1. When 9812192 finishes, read the rhat status from
   `calibration/logs/acd_posterior_check_9812377.log` (written by
   the A/B chain's preflight check).
2. If rhat is still > 1.05, retry with adapt_delta=0.99 and
   warmup=2500 (24h walltime).
3. Compare the three A/B CSVs side-by-side and pick the configuration
   that yields the lowest ACD %RMSE.
4. Draft a manuscript-section update if the calibrated A/B story
   meaningfully improved (target: ACD %RMSE matches or beats NE's
   23.19).
5. Open a PR from `acd-bridge-fix-2026-05-15` to main if the
   converged-posterior + post-pass results are sufficient.

