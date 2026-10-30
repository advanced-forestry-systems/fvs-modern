# Contributing

The authoritative guide is `CONTRIBUTING.md` at the repository root, with `CODE_OF_CONDUCT.md` and `SECURITY.md` beside it. This page summarizes it and adds the conventions visible elsewhere in the repository.

## Workflow

Fork the repository, create a feature branch from `main` (the guide's example is `git checkout -b feature/my-improvement`), make and test the change, push to your fork and open a pull request against `main`. `ci.yml` runs on every pull request to `main`, and `cross-platform.yml` runs on every pull request.

Branches named `upstream-sync/YYYYMMDD` are opened by the `upstream-sync.yml` workflow and labelled `upstream-sync`; they are bot pull requests carrying `docs/upstream_report.md`, not human work.

## Pull requests

`CONTRIBUTING.md` asks that a pull request address one concern, describe what changed and why, and reference related issues. The pull request template (`.github/pull_request_template.md`) has sections for Summary, Changes and Variants tested, and a checklist:

- Fortran changes pass `gfortran -fsyntax-only -ffree-form -std=legacy`.
- At least one variant builds.
- No credentials, tokens or large data files are included.
- Documentation is updated where applicable.

If you change `deployment/scripts/build_fvs_libraries.sh`, build at least one variant. A quick single-variant build from the guide:

```bash
bash deployment/scripts/build_fvs_libraries.sh src-converted /tmp/test-build ie
```

Build variants serially in one source tree; see [Building the engine](building-the-engine.md#build-serially).

## Code style

From `CONTRIBUTING.md`:

- Fortran: free form `.f90`, 132-column width, lowercase keywords, `IMPLICIT NONE` in new routines. When converting computed GOTOs, use `SELECT CASE` where branches are self-contained and keep labelled GOTOs for fall-through.
- Shell: `set -euo pipefail` and `#!/usr/bin/env bash`.
- R: tidyverse style, `library()` calls at the top.

`CLAUDE.md` adds: do not modernize the vendored SQLite amalgamation (`dbsqlite/`), the vendored NVEL library (`volume/NVEL/`) or `src-converted/archive/`; new Python code targets 3.9 or later with type hints and `pathlib.Path`; new code reads `FVS_PROJECT_ROOT`, `FVS_LIB_DIR` and `FVS_FIA_DATA_DIR` from the environment instead of hardcoding HPC paths, which stay confined to `calibration/osc/` and `calibration/slurm/`; and large binaries (documents, slide decks, data dumps, large images) belong in a release asset or Zenodo deposit, not in version history.

## Upstream awareness

The fork tracks the two USFS repositories. Before editing a file in `src-converted/`, check how often it changes upstream (`docs/update_lifecycle.md`). `CONTRIBUTING.md` notes that `base/` and `vbase/` are relatively stable while variant directories change more often upstream, so edits there are more likely to conflict at the next sync.

## Reporting issues

Use GitHub Issues with the templates in `.github/ISSUE_TEMPLATE/` (`bug_report.md`, `feature_request.md`). Include the affected variants, steps to reproduce, expected and actual behavior, and your platform and compiler version.
