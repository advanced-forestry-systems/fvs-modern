# Release and Zenodo deposit runbook (fvs-modern, deployment)
2026-06-18. The gated sequence to take fvs-modern from the current integration branch to a deployable,
archived release. Each step has a gate; do not advance past a failed gate.

## Gate 0: full build + regression (in progress)

Job 11758564 builds all variants (build_fvs_libraries.sh src-converted ./lib) and runs the regression
suite (run_regression_tests.sh ./lib src-converted/tests). The scheduled task harvest-build-regress
writes docs/program/20260618_build_regression_report.md with the .so count and the pass rate against the
66/67 baseline. PROCEED only if the build produces the expected libraries and the regression pass rate is
at or above the 66/67 baseline (the one known exception, iet03, is documented). If the pass rate drops,
fix before releasing.

## Step 1: land the branch on main

    cd ~/fvs-modern
    git checkout main && git pull --ff-only
    git merge --no-ff conus-sf-integration-2026-05-21 -m "Merge conus-sf-integration: disturbance-aware benchmark, four-arm, forest-type+ecoregion components, MCW consensus"
    # resolve the known working-tree cleanup (the SILC-output deletions) before or during the merge

## Step 2: update CITATION.cff and tag

    # bump version 2026.05.3 -> 2026.06.1 and date-released to the release date in CITATION.cff
    git add CITATION.cff && git commit -m "Release v2026.06.1: update CITATION version and date"
    git tag -a v2026.06.1 -m "v2026.06.1: disturbance-aware benchmark, OOS-validated adjustment layer, forest-type+ecoregion fvs-conus components, MCW consensus, full all-variant build + regression"
    git push origin main --tags

## Step 3: build the software artifact

    # a clean source archive of the tagged release (the deployable code, not the multi-GB data)
    cd ~/fvs-modern
    git archive --format=tar.gz --prefix=fvs-modern-2026.06.1/ -o ~/zenodo_staging/fvs-modern/fvs-modern-2026.06.1.tar.gz v2026.06.1
    # files_to_upload.txt: the tarball, README.md, CITATION.cff, the build_regression report

## Step 4: Zenodo new version of the concept DOI

The fvs-modern software already has concept DOI 10.5281/zenodo.19802673 (current version 2026.05.3). This
release is a NEW VERSION of that concept, not a new deposit.

    cd ~/zenodo_staging/fvs-modern
    module load python/3.11
    python new_version.py \
        --token-file ~/.zenodo_token \
        --parent-doi 10.5281/zenodo.19802673 \
        --metadata zenodo_metadata_fvsmodern_v2026.06.json \
        --files-list files_to_upload.txt \
        --publish

A sandbox dry run first (--sandbox, no --publish) is recommended since the metadata schema changed
(added forest type and ecoregion keywords, updated description). Confirm a 200 before the production run.

## Step 5: backfill the version DOI

The minted version DOI lands in CITATION.cff (identifiers and version), the README badge if a
version-specific badge is wanted (the concept badge already resolves to latest), and both manuscripts'
data-and-code-availability sections. Commit the backfill to main.

## Prerequisites and flags

- ~/.zenodo_token must exist on Cardinal (mode 600). If absent, issue a key at
  https://zenodo.org/account/settings/applications/tokens and save it with printf to that path.
- ORCID: CITATION.cff records 0000-0002-9249-1686 for this software; the lab default in other contexts is
  0000-0003-2534-4478. Confirm which ORCID to use before publishing (the metadata draft uses the
  CITATION.cff value to stay consistent with the existing record).
- The deposit is the deployable SOURCE release. The large fvs-conus data and posteriors are a separate
  data deposit (new version of the fvs_perseus_conus concept DOI), made when the validation is
  out-of-sample-honest, not bundled into the software release.

## Status

Gate 0 running (job 11758564); harvest scheduled. Steps 1 to 5 are staged and ready: the deposit metadata
draft (zenodo_metadata_fvsmodern_v2026.06.json) and this runbook are committed. Execute Step 1 onward once
the regression report confirms a clean build.
