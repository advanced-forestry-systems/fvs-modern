# Full-CONUS FVS stress test — plan and run book

Goal: exercise the modernized FVS engine over the **entire FIADB CONUS plot set**
(1,874,031 stands in `ENTIRE_FVS_STANDINIT_PLOT.csv`), projecting each stand 100
years under both the default and the calibrated configurations, and produce an
explicit pass/fail ledger so any crash, build failure, or projection error is
captured with its `STAND_CN`, variant, and config.

## Scale

| Variant | Stands | Batches (5k) | | Variant | Stands | Batches (5k) |
|---|---|---|---|---|---|---|
| SN | 523,936 | 105 | | WS | 13,977 | 3 |
| LS | 345,167 | 70  | | EC | 12,100 | 3 |
| CR | 315,529 | 64  | | CA | 11,357 | 3 |
| CS | 231,414 | 47  | | TT | 8,765  | 2 |
| NE | 139,189 | 28  | | WC | 8,629  | 2 |
| EM | 99,222  | 20  | | PN | 8,401  | 2 |
| UT | 66,364  | 14  | | BM | 7,925  | 2 |
| AK | 19,832  | 4   | | NC | 5,131  | 2 |
| IE | 15,059  | 4   | | (blank VARIANT: 12,848 — excluded) |
| SO | 14,822  | 3   | | | | |
| CI | 14,364  | 3   | | **Total run** | **~1.86M** | **~381 tasks** |

~19 variants are present in the CONUS stand-init. Each task projects up to 5,000
stands × 2 configs × 20 five-year cycles.

## Why a dedicated harness (not conus_100yr_projection.py directly)

`conus_100yr_projection.py` re-scans the full 791 MB stand-init on *every* task to
filter one variant's batch — ~381 full-file scans and heavy shared-FS contention.
This kit pre-splits the stand-init by variant once, so each task reads only its
(small) variant file, and adds an explicit failure ledger (the original only logs
failures, it does not collect them).

## Run book

```bash
S=/fs/scratch/PUOM0008/crsfaaron/fvs_stress         # staged here
FIA=/fs/scratch/PUOM0008/crsfaaron/FIA
module load gcc/12.3.0 python/3.12

# Stage 0 — one-time split by variant (~1-2 min I/O; run via salloc or a tiny sbatch,
# not a bare login-node shell if you want to be strict about OSC limits):
python3 $S/prep_split_standinit.py \
    --standinit $FIA/ENTIRE_FVS_STANDINIT_PLOT.csv \
    --out-dir   $S/standinit_by_variant --batch-size 5000

# Stage 1 — build the array manifest (prints the exact --array range, ~0-380):
python3 $S/build_manifest.py \
    --counts   $S/standinit_by_variant/counts.tsv \
    --manifest $S/manifest.tsv --batch-size 5000

# Stage 2 — submit the array, throttled to backfill politely around other jobs:
sbatch --array=0-380%16 $S/submit_conus_stress.slurm

# Stage 3 — after completion, aggregate the ledger:
python3 $S/summarize_stress.py --output-dir $S/out
#   -> stress_summary.csv, stress_failures.csv, and the headline failure rate
```

`%16` keeps it within ~16 concurrent tasks so it does not starve the active
`asym_agb_analysis` jobs. Raise the throttle once the cluster frees up.

## Calibration status (calibrated arm)

Cardinal's `~/fvs-modern` already carries the corrected 5-component calibrated
configs for all 25 standard variants (verified: `config/validate_calibrated.py`
reports 0 hard errors on the standard set; DG/HI gated to the availability table,
same as PR #62). The calibrated arm of this run is therefore meaningful as-is.
Only the DEV-only `ne.sf_preview.json` still fails validation (missing htdbh/cr) —
it is not a CONUS standard variant and is not exercised by this run.

## Validation already performed (staging, not the full run)

- All four Python stages `py_compile` clean locally and on Cardinal.
- `prep_split_standinit.py` + `build_manifest.py` verified on a 2,000-row sample.
- `run_stress_task.py` executed one task end-to-end on Cardinal: engine loads
  (NSBE 465 species, FIA tree CSV read), ledger + output CSV written. The sample
  batch (leading CS/Illinois stands) matched 0 trees — a real no-tree-condition
  property of those plots, which the ledger records; the engine itself is proven
  by the earlier benchmark that projected 111,777 conditions with 0 failures.

## Cost / runtime note

381 tasks; per-task compute is small (a 5k-stand batch is ~minutes), but stand↔tree
joins and per-state TREE CSV reads dominate. At `%16` throttle expect roughly a day
of wall time. Outputs and intermediates stay on scratch per OSC storage policy;
only the summary CSVs need promoting to home.
