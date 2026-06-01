# Session Handoff — 2026-05-31

Scope of this session: (1) landed the issue #54 calibration fix as a PR, (2) built
and staged a full-CONUS FIADB stress-test harness, (3) began incorporating Greg
Johnson's new CONUS mortality equation by re-fitting it on Cardinal, and (4)
stress-tested the fvs-conus variant projection outputs. Work spans three places:
the `fvs-modern` repo (GitHub holoros), the `~/fvs-conus` project on OSC Cardinal,
and Cardinal scratch.

---

## 1. Issue #54 — calibration multipliers (DONE, in review)

**Problem found:** the calibrated configs shipped only 3 of 5 per-species
multiplier components (dropped height-diameter and crown-ratio), so the new
validator/CI guardrail failed all 25 configs and the calibrated keyword
propagation was incomplete (the FVS-ACD "BAMAX only" symptom in #54).

**Fix shipped as PR #62** (https://github.com/holoros/fvs-modern/pull/62):
- Backfilled `htdbh_multiplier` + `cr_multiplier` into all 25 calibrated configs.
- Gated `dds`/`htg` (DG/HI) to the authoritative `equation_availability_full.csv`
  (DG adopted in 7 variants: ACD, CA, CS, KT, LS, NC, ON; HI in 6: BC, CI, EM,
  IE, KT, WS), so they are 1.0 where not adopted.
- `config/validate_calibrated.py`: 25 configs, **0 hard errors** (6 residual
  warnings = HI-adopted variants with sparse fits, expected).
- Branch `fix/issue-54-calibration-multipliers` (commit 6c31835), based on
  `feature/silc-v10-mortcal-yr100` for a clean 1-commit diff. Old local branch
  `fix/issue-54-keyword-multipliers` (superseded 3-component version) deleted.

**Open item:** PR base is the silc branch, not `main` (cherry-pick onto main
conflicts — main has diverged ~14 commits). Retarget to main after silc lands, or
merge there first. #54 will need a manual close (won't auto-close from a non-main
base). PR description saved at `PR_DESCRIPTION_issue54_5comp.md`.

---

## 2. Full-CONUS FIADB stress-test harness (STAGED, not submitted)

The engine is robust (an earlier benchmark projected 111,777 conditions with 0
failures), but the run harness re-scanned the 791 MB stand-init on every task and
only logged failures. New hardened kit (in repo at `calibration/stress/`, staged
on Cardinal at `/fs/scratch/PUOM0008/crsfaaron/fvs_stress/`):

- `prep_split_standinit.py` — split `ENTIRE_FVS_STANDINIT_PLOT.csv`
  (1,874,031 stands) by variant once.
- `build_manifest.py` — emit the SLURM array manifest (~381 tasks at batch 5000).
- `run_stress_task.py` — one (variant, batch) task with an explicit per-stand
  failure ledger.
- `submit_conus_stress.slurm` — array job (scratch-first).
- `summarize_stress.py` — aggregate ledgers into pass/fail + failure rate.
- `STRESS_TEST_PLAN.md` — per-variant counts and the exact run book.

~19 variants present (SN largest at 524k). Validated end to end on a sample (prep,
manifest, one live task with ledger). **Submit when ready** (held back because the
cluster is busy with the `asym_agb_analysis` jobs):
`sbatch --array=0-380%16 /fs/scratch/PUOM0008/crsfaaron/fvs_stress/submit_conus_stress.slurm`

---

## 3. CONUS mortality — Greg's new gompit (IN PROGRESS)

**Greg's model** (Johnson/Marshall/Weiskittel, 2026-05-26), per species:
`P_surv = 1 - exp(-exp(b0 + b1*(cr+0.01)^b2 + b3*cch^b4))`, cr = crown ratio,
cch = crown closure at tree tip. Supersedes the older DBH/CR/BAL gompit. The PDF
gives the form + parameter distributions but not the coefficient table.

**Decisions taken:** re-fit his form on Cardinal; swap mortality only; regional
pilot comparison first.

**Key finds:**
- `cch` already exists in `~/fvs-conus/data/conus_remeasurement_pairs_metric_cond_v2.rds`
  (CCH1/CCH2 on the 8.22M-pair national panel; cch range [0, 8.09], matching the
  PDF). No derivation needed.
- A naive per-species NLL fit is ill-identified: as b2 -> 0 the cr term flattens
  and b0/b1 blow up (one species diverged to b0 ~ -3.7e15). Fixed with a
  **profiled fit**: solve linear (b0,b1,b3) via a cloglog GLM with log(years)
  offset, search the exponents (b2,b4), bound b2 to Greg's interior support,
  reject separated GLM solutions, and fall back to base-rate-only (flagged).
- Validation (600k subsample, 34 species): 34/34 improved over base rate, no
  divergence, sane coefficients (b0 ~ -4..-10, b1 ~ 0.5..4.9, b2 in [-0.6,-0.2],
  b4 ~ 0..1.7), correct biological signs.

**Code (repo):** `calibration/R/35_fit_greg_mortality_cch.R` (reference, with
synthetic-recovery test) and `35b_fit_greg_mortality_conus_profiled.R`
(production, profiled + hardened). Staged on Cardinal at
`/fs/scratch/PUOM0008/crsfaaron/conus_mort/`.

**Full fit:** runs as a SLURM job over all ~133 species (>=5000 obs) of the 7.8M
filtered panel. NOTE the memory lesson: load the 172-col rds then slim to 5 cols;
use a keyed `setkey(d, SPCD)` lookup (NOT `split()` — it OOMs) and ~90G/4 cores.
Latest job: **11118369**. Output: `conus_mort/full_out/greg_mortality_coefficients.csv`.

**Next steps (not yet done):**
1. Land coefficients into `categories_conus.mortality` per variant with
   `conus_mort/land_mortality_coeffs.py` (maps SPCD -> FVS slot via FIAJSP).
2. **Engine gap:** the FVS Fortran mortality routine cannot evaluate gompit(cr,cch)
   natively. Either implement it in the `cn` variant Fortran, or apply the new
   survival in the Python projection path. This is the main remaining build.
3. Regional pilot (NE/LS/CS): project current vs new mortality, compare AGB/BA.

---

## 4. fvs-conus variant stress / realism test (DONE — a check, not a full run)

Ran `~/fvs-conus/R/18_projection_realism_check.R` over the existing calibrated
CONUS projections: **24/25 variants realistic**, only **TT flagged** ("QMD did
not grow"). BA growth mean 0.29 ft2/ac/yr (range -0.22..1.27), final BA
49.2 (TT) .. 123.6 (CS) ft2/ac. Output:
`~/fvs-conus/output/comparisons/manuscript_tables/projection_realism_detail.csv`.
Follow-up: investigate TT's stagnant QMD. A broader stress run (full perseus
projection across CONUS) can reuse the harness in section 2.

---

## Cardinal reference

- SSH: key in `CRSF-Cowork/_context/.ssh-cardinal/`; `module load gcc/12.3.0 R/4.4.0`;
  Rscript at `/apps/spack/0.21/.../r/gcc/12.3.0/4.4.0-bgpim4f/bin/Rscript`.
- Mortality work: `/fs/scratch/PUOM0008/crsfaaron/conus_mort/`
- Stress harness: `/fs/scratch/PUOM0008/crsfaaron/fvs_stress/`
- Other active jobs this session were the unrelated `asym_agb_analysis` arrays.

## Other fvs-modern refinements noted (not yet actioned)

Junk `$HOME_OUT` / `$O` / `$OUTDIR` dirs in `~/fvs-modern` (unquoted shell var
bug); the loose benchmark R validation script produces 0 validation pairs ->
NaN -> ggplot crash; `conus_stress.sh` writes to HOME instead of scratch.

---

## UPDATE (2026-05-31, later) — mortality fit complete + landed

**Full re-fit done** (after the prep->slim->fit two-step fixed the OOM): 133
species, **131 improved over base rate, 2 base-rate fallbacks** (degenerate
species correctly caught). Coefficients sane: b0 in [-10.6, -2.8], b2 mostly at
the -0.2 bound, b4 in [0.01, 3]. Total NLL 3.20M vs base 3.37M (~5% reduction).
Coefficients: `/fs/scratch/PUOM0008/crsfaaron/conus_mort/full_out/greg_mortality_coefficients.csv`.

**Model comparison (7.6M obs, `35c_compare_greg_mortality.R`):**
- AUC **0.740 (new) vs 0.673 (base rate)**; log-loss 0.420 vs 0.443; new model
  well calibrated (pred survival 0.8198 vs observed 0.8193).
- Captures the cch signal: observed survival falls 0.830 -> 0.777 across cch
  quartiles; the new model tracks it (0.834 -> 0.779) while the base rate goes
  the WRONG way (0.794 -> 0.835). This is the crown-closure-at-tip mortality
  signal from Greg's paper, which the base rate misses entirely.

**Landed** into `categories_conus.mortality` for NE (75/108 species), LS (52/68),
CS (65/96) via `calibration/python/land_mortality_coeffs.py`. Unmapped (rare,
< 5000 obs) species fall back to existing per-variant mortality by design. These
3 configs are modified locally (uncommitted).

**Engine wiring (started):** `calibration/python/greg_mortality.py` — the
`GregMortality` applicator that loads the re-fit coefficients and computes
per-tree annual hazard / period survival, with `apply_to_treelist()` to scale
TPA. Species without a fitted row return None so the caller keeps native
mortality. Validated (exact formula match; survival falls monotonically with
crown closure; unfit species preserved). **Remaining piece:** compute `cch`
(crown closure at tree tip) per tree each cycle inside perseus — FVS does not
expose cch in its treelist output, so it must be recomputed from the cycle
treelist (crown-width profile). Once that lands, the projection-level (AGB/BA)
old-vs-new comparison can run via `process_plot` with a mortality override.

**Committed:** PR #64 (`feature/conus-mortality-gompit`) — the fitters,
comparison, landing script, stress harness, applicator, and this handoff. The
landed NE/LS/CS configs were kept out of the PR pending a surgical-write fix to
`land_mortality_coeffs.py` (it currently reformats the whole JSON).

### cch provenance — the projection-comparison blocker (investigated 2026-05-31)

The projection-level (AGB/BA) old-vs-new comparison needs `cch` recomputed each
cycle. I traced its definition and it is **not reproducible from the fvs-conus
repo**:
- `30_build_conus_dataset.R` reads `CCH` directly from the change/tree zips
  (`read_tree_from_zip`, kept-cols list incl. "CCH", line ~887) — i.e. cch is
  precomputed by an UPSTREAM pipeline that generated those zips. No crown-profile
  formula exists in the fvs-conus R scripts.
- Helper scripts disagree: `41_compile_new_states.R` uses `CCH = pmin(CCFL,100)`
  (a rough proxy for new states), while `30d_fix_units.R` treats CCH on the same
  scale as HT. The v2 panel's CCH1 came from the 30_build (zip) path.
- Empirical check on 300k panel rows: CCH1 in [0, 7.1], median 0.19; cor(CCH1,
  CCFL1)=0.73, cor(CCH1, HT1)=-0.46, cor(CCH1, DBH1)=-0.42. So cch is a genuine
  crown-closure-above-the-tip metric (understory trees high, dominants low) on an
  unknown scale — NOT a simple CCFL/100 transform.

**To unblock the projection comparison, pick one:**
1. Get the cch formula from the upstream change-zip generator (whoever built
   them) — fastest and exact.
2. Fit a proxy: compute crown-competition-in-taller-trees from raw FIA tree lists
   (FVS species MCW equations) and regress against the stored CCH1 to recover the
   scaling. Self-contained but a modeling task, and the proxy will carry error
   into the cch^b4 term.
Do NOT run the projection comparison with a guessed cch — the gompit is sensitive
to it and the result would be silently wrong.
