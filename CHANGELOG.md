# Changelog

All notable changes to fvs-modern are documented in this file. The format is
based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this
project adheres to calendar-based versioning (YYYY.MM).

## [Unreleased]

### Added
- `calibration/R/63_build_variant_ecodiv_coverage.R`: one-time generator
  for `calibration/data/variant_ecodiv_coverage.csv`. Reads the area
  weighted variant x ecodivision overlay
  (`fvs_variant_ecodivision_map.csv`), filters water bodies, normalizes
  weights to sum to 1.0 per variant, and adds state-based proxies for
  the six variants without polygon geometry (ACD->NE, BC->PN, KT->EM,
  OC->CA, ON->LS, OP->PN). Verified equivalent to a pandas reference
  implementation; 96 rows covering all 25 variants.
- `calibration/data/variant_ecodiv_coverage.csv`: the produced table.
  Consumed by `62_conus_to_variant_json.R` to write coverage-weighted
  ecodivision random effect blocks. Removes the warned-against
  "all ecodivisions equally weighted" fallback path.
- `categories_conus` block in all 25 `config/calibrated/{variant}.json`
  files. Each block carries: fixed-effects summary (mean, sd, q025,
  median, q975, rhat) per model parameter; species random intercepts
  filtered to the variant's `FIAJSP` set with `species_missing` recording
  SPCDs not present in the CONUS posterior; ecodivision intercept block
  (empty for cspi_traits1 fits by model design); modifier block at the
  selected lambda; relative `draws_csv` path to the full posterior CSV.
  Pipeline version `fvs-conus phase 4`, integration date 2026-05-11,
  selected lambda 10.

### Changed
- `calibration/R/61_extract_conus_summaries.R`: broadened the
  fixed-effect name pattern to cover the unified `cspi_traits1`
  parameterization. Adds `h0..h5` (HCB), `m0..m6` (mortality),
  `r0..r6` (crown recession), `a_<name>$` for HT-DBH coefficients
  (`a_bal`, `a_ba`, `a_cspi`, `a_bard`, `a_blrd`), bare `gamma`, and
  `trait_effect` to the allow list. Removed `trait_effect` from the
  auxiliary block list (it is a real fixed effect in cspi_traits1
  models). Extended `load_meta_lookups()` to read `meta$prep_meta$sp`
  in addition to the older `meta$prep$sp_levels` and top-level
  `meta$species` conventions, so species intercept CSVs now write a
  real `SPCD` column rather than just `idx`.
- `calibration/R/62_conus_to_variant_json.R`: missing
  `_ecodiv_intercepts.csv` is now treated as "model has no ecodivision
  RE" rather than fatal. Matches the cspi_traits1 model design where
  hierarchical structure uses `L1`/`L2`/`L3` plot grouping instead of
  Bailey ecodivision.

### Build system
- `deployment/scripts/build_fvs_libraries.sh` and
  `deployment/scripts/build_fvs_executables.sh`: formally retired
  ifort support. Both scripts now force `FC=gfortran` and emit a
  warning if an environment exports anything else. The rationale,
  documented in the script header, is the `!DEC$ ATTRIBUTES ALIAS`
  decorator in `src-converted/base/fvs.f90` and similar locations:
  ifort honors the uppercase ALIAS name as the exported symbol, which
  does not match the lowercase + underscore mangling that downstream
  Fortran callers and Python ctypes wrappers expect.

## [2026.05.2] — 2026-04-26

### Fixed
- `deployment/scripts/build_fvs_libraries.sh`: corrected gfortran
  include path order. Variant-specific include directories
  (`-I$var`, `-I$var/common`) now precede `-Ibase` so that headers
  like `PRGPRM.f90` resolve to the variant's own copy with the
  correct `MAXSP` parameter. Without this, variant-shared modules
  such as `vwc/morts.f90` (used by both PN and WC, sized to
  `MAXSP=39` from `pn/common/PRGPRM.f90`) silently picked up
  `base/PRGPRM.f90`'s `MAXSP=23` and failed to compile their
  39-element DATA arrays. Foundational fix that benefits every
  variant build.
- Repaired seven F77-to-F90 conversion bugs in `src-converted/econ/`:
  `ecvol.f90`, `ecinit.f90`, `echarv.f90`, `ecstatus.f90`,
  `eccalc.f90`, `ecsetp.f90`, `ecin.f90`. The conversion script had
  mangled `C`-prefixed Fortran 77 comment blocks into trailing-`&`
  continuations of the preceding statement, leaving the `econ`
  extension uncompilable. After repair, the FVS-PN, FVS-SN, and
  FVS-IE shared libraries build cleanly and load via Python
  ctypes; twelve previously-undefined symbols (`morcon_`, `ecvol_`,
  `ecinit_`, `echarv_`, `ecstatus_`, `setpretendstatus_`,
  `eccalc_`, `ecsetp_`, `eckey_` plus three downstream) are now
  resolved.

### Added
- `scripts/fix_f77_comment_bug.py`: heuristic auto-repair tool for
  the comment-conversion pattern. Walks any Fortran file and
  rewrites trailing-`&` plus stripped-comment-fragment pairs into
  proper `!` comments. Used to repair `eccalc.f90` and `ecsetp.f90`
  during this release.
- `calibration/python/marshall_to_fia_csv.py`: adapter that
  converts Marshall-format curated FIA CSVs (in
  `calibration/data/marshall_csvs/`) to standard FIA DataMart
  schema readable by `fia_stand_generator.py`. Used to convert
  per-state inventory data for OR, WA, AL, FL, GA, MS, SC, TN, ID,
  MT during this release. Ready to drive PN/SN/IE FIA Bakuzis
  evaluation once the runtime keyword-reader regression is closed.
- `calibration/python/PN_SN_LIBRARY_DIAGNOSIS.md`: living document
  capturing the FVS-PN, FVS-SN, FVS-IE library state. Fully traces
  the conversion-bug pattern and lists every file repaired.
- `calibration/slurm/submit_bakuzis_fia_pnsnie.sh` (Cardinal-side):
  SLURM array template for the future PN, SN, IE Bakuzis run.

### Known
- FVS-PN, FVS-SN, FVS-IE libraries now load cleanly but the
  Fortran runtime hits an EOF in `base/keyrdr.f90` line 47 when
  reading any keyword file (both INVENTORY-mode and canonical
  upstream USDA test files). This is a deferred runtime-layer
  issue from the F77-to-F90 conversion that requires comparison
  against a USDA reference binary or deeper Fortran I/O debugging
  with verbose tracing. Tracked in `PN_SN_LIBRARY_DIAGNOSIS.md`.

## [2026.05.1] — 2026-04-26

### Added
- FIA-derived stand generator (`calibration/python/fia_stand_generator.py`)
  for Bakuzis matrix scenarios. Stratified random sample of FIA condition
  pairs by site index band and basal area band, returning up to five real
  plots per (variant, species_group, site_class, density_class) cell.
  Variant-to-state coverage for NE (CT MA ME NH NY RI VT) and ACD
  (ME NH VT). Inventory year forced to a canonical 2000 so all sampled
  plots project from the same calendar baseline.
- INVENTORY-mode keyfile generator (`calibration/python/inventory_keyfile.py`)
  for non-eastern variants whose DBS readers do not handle Maine-style
  state-county fields. Generates self-contained keyfiles with embedded
  TREEDATA records in fixed-width format.
- `--use-fia` and `--fia-n-plots` flags wired into
  `calibration/python/bakuzis_uncertainty_comparison.py`. The runner
  iterates over `(stand_id, stand_df, tree_df, replicate_idx)` tuples
  and adds a `replicate` column to ensemble output.
- Variant-specific species groups (`SPECIES_GROUPS_PN`,
  `SPECIES_GROUPS_SN`, `SPECIES_GROUPS_IE`) and per-variant location
  defaults in the Bakuzis runner.
- `calibration/python/run_bakuzis_post_array.sh` post-array wrapper for
  Cardinal that runs the aggregator, prints law compliance, and creates
  a tarball for fast scp pull.
- `calibration/python/marshall_to_fia_csv.py` adapter that converts
  Marshall-format FIA CSVs to standard FIA DataMart format. Ready for
  PN/SN/IE Bakuzis extension once the FVS-PN/SN library load issues
  are resolved.
- `calibration/python/PN_SN_LIBRARY_DIAGNOSIS.md` documenting the
  undefined `morcon_` symbol in FVS-PN/SN shared libraries (PN missing
  morts.f90 source entirely; SN missing it from the link line).
- Per-variant figure rendering in
  `calibration/python/bakuzis_uncertainty_figures.py`. NE and ACD
  trajectory, divergence, and band-growth panels are now produced
  automatically for every variant present in the summary.

### Changed
- `calibration/python/bakuzis_uncertainty_aggregate.py` extended for
  the FIA replicate dimension. Default and calibrated MAP rows are
  aggregated across replicates so each cell yields a single median
  plus an SD across plots. Posterior rows are pooled across the cross
  product of draws and replicates so the credible band reflects both
  parameter and stand-to-stand uncertainty.
- The aggregator now computes horizon from each replicate's own start
  year and quantizes to the nearest 5-year cycle. This aligns plots
  with different FIA INVYRs cleanly on a 0 to 100 year axis and
  removes the sawtooth pattern that would otherwise appear when
  plots with different inventory years are mixed by calendar year.
- Manuscript section 4.6 rewritten with FIA-derived numbers (33 of 36
  cells populated per variant, 3 unfilled where FIA bins lacked
  plots). ACD posterior achieves 50 percent Eichhorn compliance
  versus 20 percent point estimates, corroborating the synthetic-stand
  finding on real stands. NE Sukachev 0 percent now confirmed as a
  variant property, not synthetic artifact.
- Manuscript section 3.3 methods describes the FIA stratified sampling
  path (state coverage by variant, fixed inv_year baseline, five
  replicates per cell).
- Manuscript section 5.4 limitations updated to reflect 33-of-36 cell
  coverage per variant and the open path to extend to PN/SN/IE
  pending per-state CSV availability and shared library fixes.

### Fixed
- `docs/getting_started.md` Step 7: corrected the calibrated config
  filename from `config/calibrated/ne_calibrated.json` (does not exist)
  to `config/calibrated/ne.json` and `config/calibrated/ne_draws.json`,
  and added a working Python example showing how to load the calibrated
  overlay via `FVS(..., config_version="calibrated")`.
- `README.md` Bayesian calibration section: reconciled the component
  count (six → seven) to match `CALIBRATION.md` and `CLAUDE.md`. The
  seven components are diameter growth, height to diameter allometry,
  height increment (six variants only), mortality, crown ratio change,
  stand density, and the self thinning slope.
- `CALIBRATION.md`: replaced hardcoded OSC username (`crsfaaron`) and
  project code (`PUOM0008`) in the Cardinal example with `<OSC_USER>`
  and `<OSC_PROJECT>` placeholders. The scripts already read paths from
  environment variables; this change makes the user facing example
  portable to other OSC users and other SLURM resources.
- `config/config_loader.py` `FvsConfigLoader.compare()`: previously
  crashed with `ValueError: could not convert string to float: 'NA'`
  on most real calibrated variants because the helper tried to coerce
  species code tables and sentinel strings to `float64`. Now skips non
  numeric lists and mismatched length arrays with a conservative
  fallback. Confirmed working on NE, AK, IE, CS, LS, CA, and SN.
- `config/config_loader.py` `calibration_metadata` property and
  `summary_table()`: posterior_to_json writes the metadata block under
  `calibration` (with key `date`), but the helpers looked for
  `calibration_metadata` and `calibration_date`, so neither the date
  nor the component list rendered. Both now accept either key
  convention.

### Added
- `CALIBRATION.md`: worked `compare_configs("ne")` sample output so
  readers can see the shape of the comparison table without running it
  first, and a note clarifying that `fvs.uncertainty_summary` is a
  pandas DataFrame with a two level MultiIndex on the columns.

## [2026.04.7] — 2026-04-21

### Fixed
- `CITATION.cff`: bumped `version` from `2026.04.5` to `2026.04.7` and
  `date-released` from `2026-04-13` to `2026-04-21`. The v2026.04.6 tag
  shipped with stale citation metadata; anyone using GitHub's "Cite this
  repository" widget, Zotero, or cffconvert would have seen the older
  version string. Zenodo prefers `.zenodo.json` over `CITATION.cff` on
  ingestion so DOI minting is unaffected by the previous discrepancy.
- `RELEASE_READINESS_2026-04-21.md`: the "Latest published tag is
  v2026.04.5" line updated to reflect the v2026.04.6 push, and the
  Zenodo status line replaced with a pointer to the outstanding
  operational action on the Zenodo webhook.

### Changed
- `.gitignore` extended to cover keyword-file test outputs with `.sng`
  and `.sum` extensions; PNAS v6 manuscript response scratch at the
  repo root (`package.json`, `package-lock.json`, `node_modules/`,
  `apply_v6_patches.gs`, `build_letters.js`, `pnas_v6_figures/`); and
  dated internal status docs under `docs/` (`perseus_projection_status_*.md`
  and `docs/*_status_*.md`). None of these belong in the public source
  tree; they were working artifacts from parallel manuscript and
  calibration workstreams.

### Known operational actions
- Zenodo GitHub webhook for `holoros/fvs-modern` has never been enabled.
  A 2026-04-21 audit against the Zenodo public API confirmed zero
  fvs-modern records despite seven tags (v2026.04.0 through v2026.04.6).
  Enable at `https://zenodo.org/account/settings/github/` by toggling
  the repo ON while logged in as @holoros. After enabling, this tag
  (v2026.04.7) will be the first to mint a DOI. Past tags can be
  backfilled by manually uploading the source tarball for each prior
  tag through the Zenodo UI if historical version DOIs are desired.

## [2026.04.6] — 2026-04-21

### Added
- Standalone executable builder (`deployment/scripts/build_fvs_executables.sh`)
  for subprocess based projection pipelines. Supports all 25 variants with
  fallback linking against the shared library for unresolved symbols.
- F77 continuation fixer utility (`scripts/fix_f77_continuations.py`) for
  incremental cleanup of column 6 continuation markers remaining in .f90 files.
  Handles pure F77, half converted, and double operator patterns.
- CI expanded from 2 to 5 jobs: Fortran syntax check, NE shared library build
  with API symbol verification, executable build with smoke test (informational),
  NVEL coefficient audit, and script linting (shellcheck + ruff).
- NVEL upstream audit script (`scripts/audit_nvel_upstream.py`) comparing volume
  coefficients against the Open-FVS submodule to detect drift.
- Getting Started walkthrough (`docs/getting_started.md`) covering library build,
  ctypes verification, and keyword file simulation.
- Comprehensive FVS ecosystem landscape scan in README Related Projects section.
- `.github/CODEOWNERS` so GitHub review requests are auto-routed to the
  maintainer for Fortran source, calibration, deployment, CI, and release
  metadata paths.
- `docs/release_v2026_04_6_notes.md` release notes for the v2026.04.6 tag.

### Changed
- Executable builder includes 11 additional include directories (vvolume, vie,
  strp, clim, vdbs, vdbsqlite, econ, wpbr, wsbwe, pg, canada) for broader
  variant coverage beyond NE.
- Nine NVEL volume library files converted from fixed form to free form .f90,
  eliminating compilation stubs (volinit2, volinitnvb, fmatv, fmppget, etc.).
- Five previously stubbed subroutines (SUMOUT, OPADD, OPCSET, OPGET3, FILOPN)
  restored to full implementations.
- ACD variant .so now builds successfully on Cardinal (NVEL include path fix).
  2026-04-21 verification on Cardinal login01 under gcc/12.3.0 produced
  `FVSacd.so` (sha256 `357ac26b51a5dc18804d7f764c43b609bbac7b1f46bc3991cf94dadddf9105af`)
  that loads via `ctypes.CDLL` and exposes all four public API entry points
  (`fvssetcmdline_`, `fvssummary_`, `fvsdimsizes_`, `fvstreeattr_`). ACD is
  promoted out of advisory status and ships as a fully supported variant
  in the REST surface and shared-library matrix.
- `KNOWN_ISSUES.md`: iet03 reclassified from "segfault" to "resolved (crash
  fixed; baseline refresh pending)" after the SUMOUT/OPADD/OPCSET/OPGET3/FILOPN
  restorations. The test now exits cleanly with `STOP 10` on 2026-04-21; only
  the numeric summary diverges from the stale 2025-04-25 baseline.
- README variant counts reconciled: 23 US + 2 Canadian regional variants,
  all 25 fully supported in the REST surface description after the ACD
  Cardinal advisory closed.
- CLAUDE.md build command example corrected to `build_fvs_libraries.sh
  src-converted ./lib ne ak ie` (removed duplicated script path).

### Fixed
- FILOPN double open: JOTREE unit was opened twice when running via keyword file
  path, causing gfortran runtime errors on some platforms.
- NVEL include path ordering in build scripts so volume coefficients resolve
  correctly for all variants.
- CI symbol verification updated to check actual API entry points
  (fvssetcmdline_, fvssummary_) rather than the PROGRAM symbol.
- iet03 regression: previously segfaulted during FFE snag and fuels
  initialization, now exits cleanly (`STOP 10`). The fix rides on the
  SUMOUT/OPADD/OPCSET/OPGET3/FILOPN restorations landing before the FFE
  initializer fires. Tracked by GitHub issues #3 and #5; baseline refresh
  still pending before we claim byte-exact parity.
- Scrubbed stale Cowork session paths (`/sessions/kind-upbeat-darwin/...`)
  from eight committed docs and conversion manifests so the public checkout
  no longer leaks internal sandbox locations.

### Removed
- Stale CI workflow snapshot and internal process documents (email drafts,
  Cardinal handoff notes, release runbook) from docs/.
- `NSBE_DEPLOYMENT.md` and `NSBE_IMPLEMENTATION_SUMMARY.md` moved out of the
  public tree into `docs/internal/` (now gitignored). The underlying NSBE
  volume and biomass coefficient data remains in `data/NSBE/`.

### Fixed (continued)
- BC Canadian variant: declared `JDST(3)` in the BC-specific
  `src-converted/canada/bc/ESHAP.f90` common block and extended the
  `/ESHAP/` common list so `esinit.f90`, `varget.f90`, and `varput.f90`
  all resolve the symbol without relying on implicit typing. All three
  files now compile cleanly under `-fimplicit-none` and the rebuilt
  `FVSbc.so` loads via ctypes with the full API surface.

### Known issues
- `iet03` numeric baseline still references the 2025-04-25 snapshot; the
  simulation itself now completes cleanly but the summary diverges as
  expected from the fire and SDI plumbing changes (for example, 2003 TPA
  is 441 in the current run vs. 266 in the baseline). Baseline refresh
  is deferred pending reviewer sign-off on which run is the reference;
  tracked in GitHub issues #3 and #5.
- `src-converted/fire/base/common/FMCOM.f90` is an include-only common
  block fragment that references `MAXSP`, `MAXTRE`, and `MXFLCL` as
  array dimensions. It relies on callers including `PRGPRM.f90` and
  `FMPARM.f90` first. All production callers do, so the shared-library
  build is green; a future pass to make the file self-sufficient under
  `-fimplicit-none` remains tracked as IMPLICIT NONE hardening.

### Planned
- Refresh the iet03 summary baseline against the 2026-04-21 build and close
  issues #3 and #5.
- Make `FMCOM.f90` self-sufficient so it compiles cleanly when pulled
  in through `build_fvs_executables.sh` with `-fimplicit-none`.
- Add CONUS unified variant based on ORGANON growth forms
  (`conus-variant` branch with Greg Johnson).
- Incrementally apply F77 continuation fixes to expand exe build coverage
  beyond the current 4 variants (NE, ACD, CS, LS) that build cleanly.
- Ingrowth submodel and crown ratio revision for the next calendar tag.

## [2026.04.5] — 2026-04-13

### Fixed
- Re-tag only. The v2026.04.4 webhook delivery to Zenodo failed with a
  `ConnectTimeout` from Zenodo's worker to `api.github.com` (server-
  side transient). Zenodo's webhook receiver is idempotent on
  `delivery_id` and returned `The release has already been received`
  on retry, so a new tag is required to give Zenodo a fresh delivery
  to process. No source, calibration, deployment, or metadata changes.

## [2026.04.4] — 2026-04-13

### Added
- `.zenodo.json` with native Zenodo metadata (upload_type, license,
  creators with bare ORCID, keywords, related_identifiers). Zenodo uses
  `.zenodo.json` in preference to CITATION.cff on GitHub release
  ingestion. After v2026.04.1, .2, and .3 all returned "Citation
  metadata load failed" with distinct error_ids despite cffconvert
  validating the CITATION.cff against CFF schema 1.2.0, switching to the
  native format bypasses whatever parser quirk was hitting the CFF path.

## [2026.04.3] — 2026-04-13

### Fixed
- CITATION.cff: removed `license-url` field and the machine-readable
  `references:` block. Zenodo's GitHub release ingestion continued to
  return "Citation metadata load failed" on v2026.04.1 and v2026.04.2
  despite the file validating cleanly under cffconvert 2.0 against CFF
  schema 1.2.0. Empirically, Zenodo's parser accepts only a minimal
  subset of CFF fields. The upstream Weiskittel et al. 2011 citation is
  preserved in README.md and the v3 manuscript.

## [2026.04.2] — 2026-04-13

### Fixed
- `src-converted/base/CONTRL.f90`: stripped CRLF line terminators and
  converted fixed-form continuation (`     &` in column 6) to free-form
  trailing `&`. The file is consumed by 40+ base modules via `INCLUDE`,
  so the defect cascaded to `volkey.f90` and `evldx.f90` under the
  `gfortran -ffree-form` CI sanity check.
- `src-converted/base/errgro.f90`: renamed legacy `INCLUDE 'PRGPRM.F77'`,
  `'CONTRL.F77'`, and `'PLOT.F77'` references to the actual `.f90`
  include targets in the tree, and applied the same CRLF and
  continuation fixes as CONTRL.
- CI `syntax-check` job in `.github/workflows/ci.yml` now passes on all
  six sanity-check base files (main, myopen, volkey, algevl, evldx,
  errgro). Preexisting CI red since the PR #118 upstream sync is
  cleared.

## [2026.04.1] — 2026-04-13

### Fixed
- CITATION.cff `license` field switched from a list (`MIT`, `CC0-1.0`) to a
  single SPDX identifier (`MIT`) plus `license-url`. Zenodo's CITATION.cff
  ingestion rejected the list form with "Citation metadata load failed" on
  the v2026.04.0 release, preventing DOI minting. Dual-license intent is
  preserved in the adjacent comment block and the LICENSE file.

## [2026.04.0] — 2026-04-13

### Added
- Bayesian calibration pipeline producing posterior distributions for seven
  component models across 25 regional variants.
- UncertaintyEngine with 500 posterior draws per component preserving
  within-draw parameter correlations.
- Python ctypes wrapper (`fvs2py`) and FastAPI REST interface (`microfvs`).
- Turnkey deployment for Docker, AWS, Fedora/RHEL, macOS, and WSL2.
- Upstream sync tooling (`check_upstream.sh`, `sync_upstream.sh`) with
  per-commit traceability.
- Code of Conduct, changelog, and funding metadata.

### Changed
- Synchronized 74 Fortran files from upstream commit 66a5e14e (PR #118),
  adding SDI reporting updates, modern gfortran FORMAT descriptor fixes,
  Reineke diameter plumbing in CONTRL globals, and NSVB decay class = 3.
- Dual-licensed release: CC0-1.0 for original USDA FVS Fortran source,
  MIT for new Python/R/deployment code. CITATION.cff reflects the SPDX list.

### Fixed
- FORMAT descriptor regressions (`A8'`, `A10'`, `I3'`) that broke under
  gfortran 15.2; patches upstreamed in PR #118 so local sed hooks are
  retired in favor of a no-op stub.
- Untracked file triage: 58 stale top-level variant directories, PERSEUS
  and Bakuzis output, and ephemeral slide exports are now gitignored.

### Known Issues
- Standalone `iet03` keyword file segfaults (1 of 67 regression tests).
  Root cause under investigation; does not affect library load or rFVS
  `.Fortran()` path.

## [2025.12] — initial public draft

### Added
- Free-form Fortran 90 conversion of 2,247 source files across base plus
  25 regional variants and optional extension modules (fire, estb, volume).
- Regression harness comparing library loads and keyword-file runs against
  upstream baselines.
