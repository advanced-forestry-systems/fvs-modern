# ACD bridge follow-up: calibrated-row diagnostic + samples-link helper

Follow-up to merged PR #20. Targets the one open item: ACD's own row
is still missing from the calibrated FIA benchmark even though all 12
standard variants land cleanly (PR #20 shipped a 6.3% OVERALL BA RMSE
reduction across 12 variants + OVERALL).

## What this PR adds

1. **ACD NA diagnostic** (`calibration/R/19_fia_benchmark_engine.R`,
   gated by `FVS_ACD_DEBUG=TRUE`). Before the validation_data NA-filter,
   reports how many ACD rows carry NA in each prediction column
   (BA_pred_calib, BA_pred_default, TPA_*, QMD_*) plus a sample row.
   This pinpoints exactly which column drops ACD from the calibrated
   A/B, so the fix can target the real source instead of guessing.

2. **`calibration/slurm/link_variant_samples.sh`** — one-shot helper
   that symlinks the heavy `_samples.rds` Bayesian posteriors from the
   canonical `fvs-conus/output/variants/` into a fresh worktree's
   `calibration/output/variants/`. Without these files, every variant
   loads `params$dg = NULL` and the calibrated A/B produces all-NA
   predictions. This was the root cause of the rounds 10-15 debugging
   detour in the original effort. Tested: links 204 sample files
   across 25 variants. ACD borrows NE's diameter_growth + height_diameter
   samples (the ACD-uses-NE fallback handles the rest).

3. `.gitignore` entries so the absolute-path sample symlinks never get
   committed.

## Diagnostic run

SLURM job 10123603 running with `FVS_ACD_DEBUG=TRUE`,
`FVS_ACD_RELABEL=TRUE`, CI brackets + bootstrap disabled for speed.
Output will show the ACD NA breakdown; the actual one-line fix lands
as a follow-up commit on this branch once the column is identified.

## Why a separate PR

PR #20 is already merged and stands on its own (12 variants, 6.3%
improvement, all build/runtime fixes). This PR is the incremental
"ACD calibrated row" closeout plus the worktree-robustness helper that
prevents the missing-samples regression from recurring.

## Branch

- Branch: `acd-bridge-followup-2026-05-20` (off current main)
- Base: `main`
