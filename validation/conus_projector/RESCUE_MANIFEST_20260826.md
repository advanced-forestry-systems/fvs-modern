# FVS Rescue Manifest — 2026-08-26

Mechanical rescue of scratch-only, version-control-absent FVS work to durable Cardinal storage
(`~/fvs-modern-rescue-20260826/`). All originals left untouched on scratch. Nothing committed
to any git branch. This is a copy operation only.

## 1. conus_eq_projector_v4.R (three scratch copies plus two PREV snapshots)

| Rescued file | Source path | Source mtime | md5sum |
|---|---|---|---|
| conus_eq_projector_v4/v4_from_fvs_stress_20260708.R | /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/conus_eq_projector_v4.R | 2026-07-08 17:16:55 | 3d191c0362ddbef4839ff856e493eef6 |
| conus_eq_projector_v4/v4_from_fvs_stress_v4run_20260804.R | /fs/scratch/PUOM0008/crsfaaron/fvs_stress_v4run/conus_eq_projector_v4.R | 2026-08-04 16:41:51 | 4f4faabeb9deba4638f5476cabbf3a34 |
| conus_eq_projector_v4/v4_PREV_20260728.R | /fs/scratch/PUOM0008/crsfaaron/fvs_stress_v4run/conus_eq_projector_v4.R_PREV_20260728 | 2026-07-25 19:07:44 | 3d026fb37f987e6a85efa6b803d853a0 |
| conus_eq_projector_v4/v4_PREV_20260804.R | /fs/scratch/PUOM0008/crsfaaron/fvs_stress_v4run/conus_eq_projector_v4.R_PREV_20260804 | 2026-07-28 16:58:49 | dc3f5090f179efb3b69a00eba9385948 |
| conus_eq_projector_v4/v4_m5sizing_from_m5fix_20260804.R | /fs/scratch/PUOM0008/crsfaaron/m5_fix_20260804/conus_eq_projector_v4_m5sizing.R | 2026-08-04 22:37:57 | 0805311d8097aea0f9aeef7cb02bb602 |

All five md5sums verified byte-identical between source and rescued copy. See README.md in
`conus_eq_projector_v4/` for the diff summary across the three primary (non-PREV) copies.

## 2. conus_eq_projector_greg_par.R

| Rescued file | Source path | Source mtime | md5sum |
|---|---|---|---|
| conus_eq_projector_greg_par/greg_par_from_fvs_stress_20260708.R | /fs/scratch/PUOM0008/crsfaaron/fvs_stress/conus_eq_proj/conus_eq_projector_greg_par.R | 2026-07-08 17:16:54 | e7080d315e872d0ed83829820c409ffa |

Verified. Location anchor `conus_eq_projector_greg.R` (already durably captured elsewhere)
confirmed present in the same source directory; `conus_eq_projector_greg_par.R` was the only
uncaptured sibling matching the greg_par pattern (`greg_par_tail.R` is a distinct, unrelated file
and was not rescued as it was outside scope).

## 3. wt-gompit git worktree, uncommitted changes

Worktree: `/fs/scratch/PUOM0008/crsfaaron/wt-gompit`, branch `analysis/greg-review`, HEAD
`6cb6ba9664587ee56624db32c423d4a16677d6d`
("Final stress test artifacts: offline projector (3 defects: 2 guards + realized top-height
non-monotonicity), engine robust (90pct die-off not in fixed build; residual composition
first-cycle bump), config matrix 125/125 pass").

`git status --short` enumerated 130 entries, all prefix `??` (untracked/new). Zero tracked-modified
files, so `git diff` for tracked files produced an empty patch (0 bytes), saved for completeness at
`wt-gompit-uncommitted/wt-gompit-tracked-changes.patch`.

All 130 untracked files copied preserving relative path structure under
`wt-gompit-uncommitted/`:
- `config/calibrated/` — 125 files, 25 variant codes (acd, ak, bc, bm, ca, ci, cr, cs, ec, em, ie,
  kt, ls, nc, ne, oc, on, op, pn, sn, so, tt, ut, wc, ws) x 5 timestamped
  `.json.pre_stand_YYYYMMDD_HHMMSS` snapshots each.
- `sf_integration_dev/` — 5 files: `constrained_all4_conus.png`, `constrained_all4_organon.png`,
  `constrained_projection_all4_conus.json`, `constrained_projection_all4_organon.json`,
  `constrained_projection_result.json`.

Full per-file source path, md5sum, and mtime verification: all 130 files checked individually
against the live worktree copy, 0 mismatches (see `wt-gompit-uncommitted/wtgompit_file_manifest.tsv`
for the complete per-file table with md5sums).

Checked the three JSON output artifacts in `sf_integration_dev/` for coordinate-like fields
(`lat`, `lon`, `latitude`, `longitude`, `coord`, case-insensitive) — none found. No FIA true plot
coordinates present in any rescued file.

## Verification summary

- 6 files (v4 x5, greg_par x1): md5 source-vs-rescued match, 6/6.
- 130 files (wt-gompit untracked): md5 source-vs-rescued match, 130/130, 0 mismatches.
- 1 patch file (wt-gompit tracked changes): empty, as expected (0 tracked-modified files).
- Total files in rescue bundle: 137 (6 + 130 + 1 patch).
- Originals: untouched on scratch (copy only, no scratch deletions performed).
- Git: nothing committed to any branch in fvs-modern or the wt-gompit worktree.
