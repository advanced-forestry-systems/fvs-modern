# fvs2py engine-interface fixes: smoke and script checks (PR 143, October 30, 2026)

Evidence for the summary and shared-engine notes in `docs/user-guide/outputs.md` and
`docs/user-guide/running-projections.md`. All runs were on OSC Cardinal (glibc 2.34) with the
variant libraries of the PR 143 worktree. No plot identifiers, PLT_CN values or coordinates are
in these files.

## Smoke comparison (summary width, first commit of PR 143)

`smoke_summary.py` projects 5 real FIA stands per variant (NE and ACD; scenario index 13 of
`bakuzis_uncertainty_comparison.enumerate_scenarios`, stands from `generate_real_stand` with seed
42 + scenario number) under the default and calibrated MAP configurations through
`bakuzis_uncertainty_comparison.run_scenario`. It writes the summary table of every run.

| File | Code |
|---|---|
| `smoke_before.csv` | main before PR 143 (20 summary columns) |
| `smoke_after.csv` | PR 143 first commit (22 summary columns) |
| `smoke_after_review.csv` | PR 143 with the review fixes (byte-identical to `smoke_after.csv`) |

`smoke_comparison.csv` (written by `review2_checks.slurm`, which calls `compare_pairs.py` on the three CSVs in this folder; the labels in the file are the ones in that script): each file has 420 rows (2 variants x 2
configurations x 5 stands x 21 summary rows per run). Before versus first commit, the 23 shared columns (20
summary columns plus variant, replicate, config) are identical; only `scuft` and `rscuft` are new.
In `smoke_after.csv`, `scuft` runs from 249 to 8,450 ft3/ac and `rscuft` is 0 (no harvest).
First commit versus review fix: all 25 columns identical.

## Script checks for the review fixes

The review fix closes a collected `FVS` object's library handle (weakref finalizer). Scripts that
build one `FVS` object per stand and never close it now load a fresh engine for each stand.
`review_checks.slurm` ran each script below on the PR 143 first commit ("old") and on the review
fix ("new") and compared the outputs. `review_script_checks.csv` holds the result: all six
comparisons are identical, with maximum absolute difference 0.

- `smoke_summary.py` (calls `bakuzis_uncertainty_comparison.run_scenario`), 420 rows.
- `smoke_both.py` (the same stands through `run_scenario` of both
  `bakuzis_uncertainty_comparison.py` and `bakuzis_100yr_comparison.py`), 840 rows.
- `validation/bakuzis_stress_v2/engine_bakuzis_v2.py` on FVSne.so, 756 rows (36 of 36 scenarios
  pass its seeding and site gates in both runs).
- `validation/bakuzis_stress_v2/ensemble_regen4/confirm/confirm_real.py d12` (worktree and output
  paths substituted): `d12_runs.csv`, 26 rows, and `d1_echo_decode.csv`, 456 rows.

The file names in `review_script_checks.csv` refer to the job's scratch outputs on Cardinal, which
are not copied here.

## Second review round (concurrency fix)

`_check_and_register` now loads the library, sets `owner._lib` and attaches the finalizer under
the registry lock, and undoes the registration if loading fails. Before, the library was loaded
after the lock was released, so a second constructor in another thread could see a registered
owner whose `_lib` was still None and treat it as not live. `review2_checks.slurm` ran on the PR
head (Cardinal, glibc 2.34) and wrote:

- `test_engine_fixes_run.txt`: the 12 tests of `test_engine_fixes.py`, all passed (including the
  new `test_failed_load_does_not_leave_a_registration`).
- `concurrency_20_repeats.txt`: `test_concurrent_construction_registers_one_owner` in 20 separate
  pytest runs, 20 of 20 passed.
- `smoke_comparison.csv`, regenerated; it is byte-identical to the earlier copy.

`review_checks.slurm` is the first round's job as run; its last call wrote a scratch-only
comparison with other labels, which is not part of this folder.
