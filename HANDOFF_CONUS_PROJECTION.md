# Handoff — CONUS-wide FVS projection + Greg mortality (2026-06-01)

Start-here doc for a fresh session. Everything below is committed to GitHub PR #64
(branch `feature/conus-mortality-gompit`, repo holoros/fvs-modern) unless noted.

## TL;DR — where we are

1. **Greg's CONUS mortality (gompit on crown ratio + crown closure at tree tip):**
   fully re-fit, validated, applicator + projection comparison built. DONE.
2. **CONUS-wide FVS projection by variant from FIADB:** the harness + engine work,
   but the run was broken by a **tree-data mismatch**. Root cause found and fixed;
   the missing **FVS TreeInit tables are re-downloading now**; a corrected runner
   is built and a proof job is running. RESUME HERE.

## The CONUS projection blocker + fix (the important part)

- The standinit `ENTIRE_FVS_STANDINIT_PLOT.csv` keys stands by the **15-digit FVS
  `STAND_CN`** (e.g. 750085992290487), NOT the FIA `PLT_CN`.
- The raw DataMart `<ST>_TREE.csv` uses a short `PLT_CN` (e.g. 11839) and does not
  contain the standinit stands -> the old harness matched ~0.7% of stands.
- The correct tree lists are the **FVS-native `<ST>_FVS_TREEINIT_PLOT.csv`**,
  bundled inside each state's DataMart `<ST>_CSV.zip`. MT and TX both join the
  standinit **3000/3000** sampled. This is the fix.

## Cardinal locations (user crsfaaron, account PUOM0008)

Connect: SSH key at `CRSF-Cowork/_context/.ssh-cardinal/id_ed25519_cardinal`
(`ssh -F ~/.ssh/config cardinal`). R: `module load gcc/12.3.0 R/4.4.0`; Rscript:
`/apps/spack/0.21/cardinal/linux-rhel9-sapphirerapids/r/gcc/12.3.0/4.4.0-bgpim4f/bin/Rscript`.

```
/fs/scratch/PUOM0008/crsfaaron/
  FIA/ENTIRE_FVS_STANDINIT_PLOT.csv          # CONUS standinit, 1,874,031 stands (STAND_CN-keyed)
  FIA_fresh/
    dl_fvs.sh                                # resumable downloader (RE-RUN to resume)
    dl_fvs.log                               # progress log
    treeinit/<ST>_FVS_TREEINIT_PLOT.csv      # THE FIX: FVS-native tree lists per state
  fvs_stress/
    standinit_by_variant/  (counts.tsv)      # standinit pre-split by variant (19 variants)
    manifest.tsv                             # 381 array tasks: idx<TAB>variant<TAB>batch<TAB>size
    run_conus_task_fvstreeinit.py            # CORRECTED runner (uses FVS TreeInit)
    submit_conus_fvs.slurm                   # array submit for the corrected runner
    run_stress_task.py / submit_conus_stress.slurm   # OLD runner (wrong tree source - do not use)
    out_fvs/                                 # output + ledgers from corrected runner
  conus_mort/                                # all Greg-mortality work (below)
~/fvs-modern/lib/FVS*.so                     # 25 compiled FVS variant libraries
~/fvs-conus/data/raw_fia/                    # has MT_FVS_TREEINIT_PLOT.csv, TX_CSV.zip
```

## RESUME STEPS for the CONUS projection

1. **Check the download finished:**
   `ls /fs/scratch/PUOM0008/crsfaaron/FIA_fresh/treeinit/*_FVS_TREEINIT_PLOT.csv | wc -l`
   (target ~49). If interrupted, re-run `bash /fs/scratch/PUOM0008/crsfaaron/FIA_fresh/dl_fvs.sh`
   (it skips completed states). URL base: `https://apps.fs.usda.gov/fia/datamart/CSV/<ST>_CSV.zip`.
2. **Check the proof job** (SN batch 0, corrected runner): job `11204332`,
   ledger `fvs_stress/out_fvs/ledger_sn_b0.json`. Confirm `n_stands_with_trees`
   is now a large fraction of 5000 (the fix worked) and `n_output_rows` > 0.
3. **Launch the full CONUS array** (throttle to backfill around other jobs):
   `sbatch --array=0-380%16 /fs/scratch/PUOM0008/crsfaaron/fvs_stress/submit_conus_fvs.slurm`
4. **Aggregate** when done: each task writes `out_fvs/conus_<variant>_b<batch>.csv`
   (STAND_CN, STATE, YEAR, PROJ_YEAR, VARIANT, CONFIG=default|calibrated, AGB_TONS_AC)
   and `ledger_*.json`. Sum/compare default vs calibrated AGB by variant.

Note: the corrected runner loads each state's TreeInit file once and groups by
STAND_CN; large southern states (TX/SN region) make the first task slow. If memory
is tight, raise `--mem` in `submit_conus_fvs.slurm` (currently 12G).

## Greg mortality work (DONE) — locations

- Coefficients: `/fs/scratch/PUOM0008/crsfaaron/conus_mort/full_out/greg_mortality_coefficients.csv`
  (133 species, 131 improved over base rate). Slim panel: `conus_mort/mort_slim.rds`.
- Code (repo `calibration/`): `R/35b_fit_greg_mortality_conus_profiled.R` (fitter),
  `R/35c_compare_greg_mortality.R` (model comparison: AUC 0.74 vs 0.67 base),
  `R/35d_validate_cch.R` (cch port validated, Spearman 0.93),
  `python/cch_organon.py` (CCH from CAL_CCH.for), `python/greg_mortality.py`
  (gompit applicator), `python/project_mortality.py` (projection glue),
  `python/land_mortality_coeffs.py` (lands into categories_conus.mortality),
  `python/project_compare.py` (100yr mortality comparison: Greg retains +41.5%
  BA at yr100 vs base rate; trajectory in `calibration/output/`).
- Open item: make `land_mortality_coeffs.py` write surgically (it reformats the
  whole JSON) before landing into tracked configs.

## Other open items
- PR #62 (issue #54 calibration fix) is open; decide merge/retarget vs `main`.
- Full mortality detail in `SESSION_HANDOFF_2026-05-31.md` (UPDATES 1-7).
