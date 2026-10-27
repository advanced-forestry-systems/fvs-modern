# Autopilot run book, 27 October 2026 (prepared while the computer link was down)

Everything in this folder was built in the cloud container after the pop-os bridge dropped at about 11:30 local, which took firebreather and Cardinal offline with it. Nothing here has been executed on a cluster yet. The steps below are in the order they should run once the link is back; each one is reversible and nothing deletes or overwrites a prior version.

## Files

| File | What it is | Where it goes |
|---|---|---|
| `bakuzis_stress_v2.R` | the v2 scorer (already on PR 112) | unchanged |
| `engine_bakuzis_v2.py` | engine harness with SITECODE in place of the non-keyword SETSITE, a year-0 seeding gate (20 percent), a site-index readback gate, and A/B over two survival libraries | `fvs-modern/validation/bakuzis_stress_v2/` |
| `fia_ba_sdi_percentiles.R` | observed BA and SDI p99 and p99.9 by variant from the CONUS pairs RDS, for the magnitude test | same |
| `2026-10-27_bakuzis-reruns_LAUNCHER.sh` | three opt-in phases on Cardinal (engine A/B, FIA ensemble seeding audit, projector v4 rescore) | same |
| `20261027_fvs-modern_fvs-conus_bakuzis_assessment.md` | the memo, voice-passed (one British spelling and six aphoristic colons recast) | replaces the PR 112 copy |
| `2026-10-27_paper1-combined-draft_HOLD-INSERTS.md` | four additive inserts for `manuscript/fvs_combined_draft.md` | `fvs-modern/manuscript/` after the draft is backed up as `_PREV_20261027` |
| `Paper2_CONUS-FVS_benchmarking_2026-10-27_DRAFT.md`, `CHANGES_2026-10-27_Paper2.md` | Paper 2 with the hold extended in the header, Section 3.4 and limitations; nothing removed | `MAINE/.../FVS_CONUS/` beside the 4 August version |
| `Paper2_FEM_cover_letter_2026-10-27_DRAFT.md` | the R26 to R40 cover letter rewrite on the repaired-build narrative (gross basis, mortality, Greg 19 percent) | same folder; replaces the stale 17 July letter once Aaron approves |

## Steps when the bridge returns

1. Commit this folder to `CRSF-Cowork/outputs/2026-10/bakuzis_stress_v2/autopilot/` and the two Paper 2 files to `FVS_CONUS/` (device_commit_files; fresh filenames, never overwrite the 4 August files).
2. `fire_upload` the folder into `jobs/bakuzis_fvsconus_20261027/autopilot/`, check md5 on arrival.
3. Push to GitHub from firebreather with `gh_api_push.py --repo advanced-forestry-systems/fvs-modern --base validation/bakuzis-stress-v2-20261027 --branch validation/bakuzis-stress-v2-20261027-reruns`, adding `engine_bakuzis_v2.py`, `fia_ba_sdi_percentiles.R`, the launcher, the hold-inserts note, and the memo with `--allow-overwrite` on the memo path. Dry run first. Token `~/.config/gh-tokens/advanced-forestry-systems` (the holoros token refused blob writes on 27 October; check its Contents permission on GitHub).
4. On Cardinal, in the wt-engine worktree: copy the folder, then `bash 2026-10-27_bakuzis-reruns_LAUNCHER.sh` (plan only), then `PHASE1=1`, read `seed_site_gate_*.csv` before reading any scorecard. A SITE_FAIL on every row means SITECODE did not take either and the NE site mechanism needs a different keyword; a SEED_FAIL means the engine is not seeding what the generator wrote.
5. `PHASE2=1` prints year-0 BA against its label per cell for the banked ensemble. The expected finding is a ratio near 0.1; then read `calibration/python/inventory_keyfile.py` for the per-record TPA written to TREEINIT and fix the expansion factor before any re-run of `bakuzis_uncertainty_comparison.py`.
6. `PHASE3=1` needs the projector v4 `*_metrics.csv` files under `$SCR/conus_eq_proj`; if array 13309592 has read out, point `SCR` at its output directory.
7. Apply the Paper 1 inserts to the combined draft on the firebreather clone (back up first), and re-render both manuscripts to PDF for the visual pass before anything is sent.

## Open items not advanced here

The Paper 2 Zenodo deposit and evidence map (R26 to R40 remainder) need the repaired-build artifacts, which are on Cardinal. The holoros token permission needs the GitHub settings page. The Table 6 scorer either turns up in a Cardinal scratch directory (`/fs/scratch/PUOM0008/crsfaaron/fvs_stress/`) or is retired; search `fvs_stress` for a script that computes "inversions" before writing a replacement.
