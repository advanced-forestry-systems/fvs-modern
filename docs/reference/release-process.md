# Release process

The repository has no release script or release workflow. This page documents what the repository shows about versioning, the changelog, citation metadata and Zenodo, and the steps that past releases followed according to their notes.

## Versioning

Versions are calendar based. `CHANGELOG.md` states that the project follows calendar-based versioning (YYYY.MM), and `CLAUDE.md` states that releases are tagged `vYYYY.MM.N`. The changelog entries use the three-part form, from `[2026.04.0]` through `[2026.05.2]`, with an earlier `[2025.12]` initial public draft. The release notes in `docs/release_v2026_04_6_notes.md` and `docs/release_v2026_04_7_notes.md` name their tags `v2026.04.6` and `v2026.04.7`.

!!! warning "Metadata out of step at this commit"
    `CITATION.cff` carries `version: "2026.05.3"` and `date-released: "2026-04-26"`, while `CHANGELOG.md` has no `2026.05.3` entry; its newest released entry is `[2026.05.2]`, dated 2026-04-26, followed by `[Unreleased]`. The `2026.05.x` entries also carry April dates. Reconcile `CITATION.cff`, `CHANGELOG.md` and the tag name before cutting the next release.

## CHANGELOG

`CHANGELOG.md` follows Keep a Changelog 1.1.0. New work goes under `## [Unreleased]` with `### Added`, `### Changed`, `### Fixed` and, in some entries, `### Known` or `### Known Issues` subsections. At release time the unreleased block becomes `## [<version>]` with the release date.

## Citation and Zenodo metadata

Two files describe the software for citation:

- `CITATION.cff` (CFF 1.2.0) with title, author, ORCID, `repository-code`, `license: MIT`, `version`, `date-released` and a DOI identifier.
- `.zenodo.json` with title, description, `upload_type: software`, `license: MIT`, `access_right: open`, creators and keywords.

The changelog records why both exist. Zenodo's GitHub ingestion rejected `CITATION.cff` with "Citation metadata load failed" on v2026.04.0 to v2026.04.3: a list-valued `license` was replaced by a single SPDX identifier in 2026.04.1, and `license-url` and the `references:` block were removed in 2026.04.3. `.zenodo.json` was added in 2026.04.4 because Zenodo uses it in preference to `CITATION.cff`. Keep `CITATION.cff` to that minimal field set; comments in the file record these constraints.

The DOI in `CITATION.cff` and the README badge is `10.5281/zenodo.19802673`, which `CITATION.cff` describes as the Zenodo archive for v2026.05.3. The README also cites a separate data record for the calibrated parameter posteriors and validation outputs, `10.5281/zenodo.20760580`, and the configuration catalog work is archived at `10.5281/zenodo.23777723`.

!!! note "Zenodo webhook"
    `docs/release_v2026_04_7_notes.md` records that, as of 2026-04-21, the Zenodo GitHub integration had never been enabled for the repository and no tag had been ingested, and gives the enablement steps. The repository does not record whether the later DOI was minted through the webhook or uploaded manually.

## Cutting a release

The following steps are assembled from the release notes and changelog entries; no script automates them.

1. Make sure `main` is green in CI (`ci.yml`, `cross-platform.yml`) and that `python3 config/validate_calibrated.py` passes.
2. Move the `[Unreleased]` entries in `CHANGELOG.md` to a new `## [YYYY.MM.N]` heading with the date.
3. Update `version` and `date-released` in `CITATION.cff` to match. The v2026.04.7 release exists only to correct a tag that shipped with stale values in this file.
4. Optionally add release notes under `docs/release_<tag>_notes.md`, following `docs/release_v2026_04_7_notes.md`.
5. Tag `vYYYY.MM.N` and publish a GitHub release. A published release triggers `docker-publish.yml`, which pushes images tagged with the semver tag, and, once the Zenodo integration is enabled, the release notes say Zenodo reads `.zenodo.json` and mints a version DOI under a stable concept DOI.

If Zenodo misses a delivery, the changelog entry for 2026.04.5 records that its webhook receiver treats a repeated delivery as already received, so the fix used then was a new tag rather than a retry.

!!! warning "Unverified"
    Whether releases are cut from `main` directly or from a release branch, and who holds tagging rights, is not recorded in the repository. The repository URL in `CITATION.cff`, the README badges and the Docker image names is `holoros/fvs-modern`; whether the canonical location has moved to `advanced-forestry-systems/fvs-modern` is not recorded in the snapshot used for this page.
