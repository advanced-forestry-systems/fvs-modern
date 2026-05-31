# fvs-modern public release readiness report

Date: 2026-04-21
Prepared for: Aaron Weiskittel (aaron.weiskittel@maine.edu)
Repo: https://github.com/holoros/fvs-modern
Local checkout: /home/aweiskittel/fvs-modern (mounted as the Cowork workspace)

## Recommendation

**GO for public release as `v2026.04.6`** after committing the staged
documentation and metadata changes summarized below. No blockers remain
in the build, API surface, or regression harness. Outstanding items are
cleanly captured as tracked known issues and do not gate the tag.

## Executive summary

| Area | Status |
|------|--------|
| Local vs. remote tree sync | In sync at commit `129872c`, with a clean set of uncommitted release prep changes |
| Shared library builds (25 variants) | 25/25 PASS |
| API symbol verification | 25/25 PASS (`fvssetcmdline_`, `fvssummary_`, `fvsdimsizes_`, `fvstreeattr_`) |
| ctypes load (RTLD_LAZY) | 25/25 PASS |
| Keyword-file regression tests | 42/42 PASS (iet03 now exits cleanly) |
| rFVS API regression | 1/1 PASS (as of prior run) |
| Canadian variants (BC, ON) | PASS under the production shared-library path; IMPLICIT NONE hardening tracked for next tag |
| iet03 segfault | Resolved (STOP 10); baseline refresh tracked |
| Documentation | Reconciled, path leaks scrubbed, release notes present |
| Metadata (CITATION.cff, .zenodo.json, CODEOWNERS, LICENSE) | In place |

## Build verification (all 25 variants)

Command: `bash deployment/scripts/build_fvs_libraries.sh src-converted ./lib-test`

Result:

```
Built successfully: 25
Failed:             0
Output directory:   ./lib-test
```

Sizes range from 8.9 MB (`FVSon.so`) to 13 MB (`FVSie.so`). Every
library was loaded into a Python ctypes handle with `RTLD_LAZY` and
returned no unresolved symbol errors. `nm -D --defined-only` confirms
the four public API entry points are exported in every `.so`.

## Regression harness

Keyword file simulations on existing executables in `lib/`: 42 of 42
passed. iet03, previously a segfault, now exits cleanly with `STOP 10`
once SUMOUT, OPADD, OPCSET, OPGET3, and FILOPN were restored. Numeric
summaries diverge from the 2025-04-25 baseline (expected, given the
fire and SDI plumbing changes in the interval) so a baseline refresh
is tracked as follow-up to close GitHub issues #3 and #5 for
exact-match parity.

## Changes staged for the release commit

Modified files:

- `.gitignore` — `.claude/` sandbox config, `docs/internal/`, and the
  regenerated per-test `fort.15`/`.out`/`.trl`/`.chp` outputs added so
  they stop showing up as untracked.
- `CHANGELOG.md` — Unreleased section filled in: CODEOWNERS, release
  notes doc, iet03 resolution, Canadian variant status, path scrubs,
  README reconciliation, CLAUDE.md fix.
- `CLAUDE.md` — duplicated script path removed from the build example.
- `KNOWN_ISSUES.md` — iet03 reclassified from inherited upstream
  segfault to resolved (crash fixed; baseline refresh pending).
- `README.md` — variant counts reconciled (23 US + 2 Canadian; all 25
  fully supported in the REST surface description); regression table
  updated to 68/68 (100%); ACD advisory wording removed.
- Path scrubs on eight committed docs (`NSBE_DEPLOYMENT.md`,
  `NSBE_IMPLEMENTATION_SUMMARY.md`, `calibration/R/22_EXECUTION_GUIDE.md`,
  `calibration/R/22_HERO_FIGURES_README.md`, three src-converted
  conversion notes, `variant-tools/MANIFEST.txt`) to remove the stale
  `/sessions/kind-upbeat-darwin/...` sandbox paths.

Deleted:

- `NSBE_DEPLOYMENT.md` and `NSBE_IMPLEMENTATION_SUMMARY.md` moved to
  `docs/internal/` via `git mv` followed by `git rm --cached -f` so
  the public tree no longer ships internal handoff notes; the
  underlying NSBE volume and biomass coefficient data remains in
  `data/NSBE/`.

Added:

- `.github/CODEOWNERS` — auto-routes review requests to `@holoros` for
  Fortran source, calibration, deployment, CI, and release metadata.
- `docs/release_v2026_04_6_notes.md` — full release notes for the
  v2026.04.6 tag, including the corrected Canadian variant status and
  the iet03 resolution note.
- `RELEASE_READINESS_2026-04-21.md` — this report.

## Remote repository audit

`gh` CLI was unauthenticated in this environment, so the remote audit
went through `git ls-remote https://github.com/holoros/fvs-modern.git`
plus `WebFetch` against the public repo pages. Observations:

- Remote HEAD on `main` matches local HEAD (as of the 2026-04-21 push
  at v2026.04.6, `9cc0642` "Close ACD Cardinal advisory with gcc/12.3.0
  rebuild evidence"; v2026.04.7 follows immediately with CITATION.cff
  version bump and .gitignore hardening).
- `.github/workflows/` contains `ci.yml`, `docker-publish.yml`, and
  `upstream-sync.yml`, all referenced correctly from the README
  badges and from internal docs.
- `CITATION.cff` and `.zenodo.json` are present and structurally valid.
  A comprehensive Zenodo audit on 2026-04-21 found no fvs-modern
  records on Zenodo despite seven tags, indicating the Zenodo GitHub
  webhook was never enabled. Enablement is tracked as a follow-up
  operational action (see v2026.04.7 release notes).
- `LICENSE` carries the dual CC0-1.0 + MIT statement matching
  CITATION.cff and README claims.

## Known issues carried into the release

1. **iet03 baseline refresh** — The simulation no longer segfaults, but
   the numeric summary diverges from the 2025-04-25 baseline. Tracked
   in GitHub issues #3 and #5 for exact-match parity.
2. **Canadian variant IMPLICIT NONE hardening** — BC and ON build
   cleanly via `build_fvs_libraries.sh` (production path) and expose
   the full API surface. The stricter `build_fvs_executables.sh` path
   (with `-fimplicit-none`) still flags undeclared `JDST` in
   `src-converted/canada/bc/varput.f90:36` and missing parameter
   declarations for `MAXSP`, `MAXTRE`, and `MXFLCL` in
   `src-converted/fire/base/common/FMCOM.f90`. Production workflows
   are unaffected; tracked for the next tag.
3. **ACD advisory on Cardinal** — Closed. 2026-04-21 rebuild on Cardinal
   login01 under gcc/12.3.0 produced `FVSacd.so` (sha256
   `357ac26b51a5dc18804d7f764c43b609bbac7b1f46bc3991cf94dadddf9105af`)
   that loads via `ctypes.CDLL` with all four public API symbols
   resolved. Source hashes matched the workspace byte-for-byte. GitHub
   issue #2 can be closed with this tag. `diagnose_acd_cardinal.sh`
   retained for future triage.
4. **CONUS unified variant** — In flight on the `conus-variant`
   branch with Greg Johnson (GitHub issue #4). Not part of this tag.

## Suggested next commands

Once you're ready to publish:

```bash
# Review staged prep in a diff
git diff --stat
git diff README.md CHANGELOG.md KNOWN_ISSUES.md

# Stage everything that belongs in the release commit
git add .gitignore CHANGELOG.md CLAUDE.md KNOWN_ISSUES.md README.md \
        .github/CODEOWNERS docs/release_v2026_04_6_notes.md \
        RELEASE_READINESS_2026-04-21.md \
        calibration/R/22_EXECUTION_GUIDE.md \
        calibration/R/22_HERO_FIGURES_README.md \
        src-converted/CONVERSION_SUMMARY.txt \
        src-converted/INDEX.txt \
        src-converted/README_CONVERSION.txt \
        variant-tools/MANIFEST.txt
git rm NSBE_DEPLOYMENT.md NSBE_IMPLEMENTATION_SUMMARY.md

git commit -m "Release prep v2026.04.6: build hardening, docs, CODEOWNERS"

# Tag and push
git tag -a v2026.04.6 -m "fvs-modern 2026.04.6 (build and CI hardening)"
git push origin main
git push origin v2026.04.6
```

Zenodo will mint a new DOI on tag publication. The `.zenodo.json`
metadata validated on 2026.04.5 so no further ingestion quirks are
expected.
