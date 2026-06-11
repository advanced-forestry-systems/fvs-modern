# History Rewrite & Large-File Plan

## Problem

The repository's git history carries large binary artifacts that were committed,
revised, and sometimes deleted over time. They stay in history forever and every
clone pays for them.

Measured on a fresh `--mirror` clone of `origin` (2026-06-11):

| State | Mirror size |
|-------|-------------|
| Current history (origin) | **210 MB** |
| After stripping derivable artifact paths | **17 MB** |

That is a **~92% reduction**. The heaviest offenders in history are repeated
revisions of `manuscript/*.docx` (a single SI file appears at ~28 MB across
several revisions), slide decks under `calibration/slides/*.pptx`, generated
figures under `calibration/output/**/*.{pdf,png}`, large CSV dumps, and the
`calibration/figshare/*.zip` bundle.

## What stays vs. goes

Keep in history (real source needed to build):
- `src-converted/dbsqlite/sqlite3.c` / `.h` (vendored SQLite, ~7.8 MB)
- `src-converted/volume/NVEL/` source (vendored volume library)

Strip from history (derivable, archivable, or distributable elsewhere):
- `manuscript/` (drafts, SI, figures) -> move final versions to a release/Zenodo
- `calibration/slides/` (`.pptx`)
- `calibration/output/` (generated figures, benchmark reports)
- `calibration/figshare/` (zip bundles)
- All `*.docx`, `*.pptx`, `*.xlsx`, `*.pdf` repo-wide

Optional further trim (vendored Windows build artifacts, ~3 MB):
- `src-converted/volume/NVEL/volCStest/` (`bin_/`, `publish*.zip`)

## Procedure (run when ready -- this rewrites history)

This is irreversible and changes every commit hash. It requires a force-push and
forces all collaborators to re-clone. Coordinate first.

1. Announce a freeze. Make sure open PRs are merged or noted; branches will need
   to be re-based or re-created from the rewritten history.

2. Back up the current remote:
   ```bash
   git clone --mirror https://github.com/holoros/fvs-modern.git fvs-backup.git
   tar czf fvs-backup-$(date +%Y%m%d).tar.gz fvs-backup.git   # store off-repo
   ```

3. Before stripping, preserve anything worth keeping as a release asset or
   Zenodo deposit (final manuscript PDF, key figures, the figshare bundle).

4. Rewrite on a fresh mirror:
   ```bash
   pip install git-filter-repo
   git clone --mirror https://github.com/holoros/fvs-modern.git fvs-rewrite.git
   cd fvs-rewrite.git
   git filter-repo --force \
     --path manuscript/ \
     --path calibration/slides/ \
     --path calibration/output/ \
     --path calibration/figshare/ \
     --path-glob '*.docx' --path-glob '*.pptx' --path-glob '*.xlsx' \
     --path-glob '*.pdf' \
     --invert-paths
   git reflog expire --expire=now --all && git gc --prune=now --aggressive
   du -sh .            # expect ~17 MB
   ```

5. Inspect the result (log, key files present, build still works from a normal
   clone of the rewritten mirror) before pushing.

6. Force-push all refs and tags:
   ```bash
   git push --force --mirror https://github.com/holoros/fvs-modern.git
   ```

7. Have every collaborator re-clone. Old local clones will diverge and must be
   discarded, not merged.

## Prevent recurrence

- Tighten `.gitignore` so generated figures/outputs and office docs never get
  staged (`*.docx`, `*.pptx`, `*.xlsx`, `perseus_*.{png,pdf,xlsx}`, `*.mod`,
  build scratch dirs).
- Adopt Git LFS for any binary that genuinely must be versioned, or attach it to
  a tagged release / Zenodo deposit instead.
- Treat `manuscript/` and `calibration/output/` as build products, not source.
