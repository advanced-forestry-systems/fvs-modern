# fvs-modern release readiness — 19 November 2026

This record distinguishes a working engine build from a validated growth-and-yield release. It describes `main` at `d0ff5bc35bc4ae6fe8204d6c4125c8c62f7b66ab` and the open research branches on this date. Update the commit and evidence before using it for a later release.

## What is on `main`

The platform build and load checks pass on Linux, macOS, and Windows after PRs 147–149. Input and output guards from PR 149 are present. These checks establish that the engine can build and reject specified invalid inputs. They do not establish biological forecast accuracy. The last published GitHub release is `v2026.05.4`; subsequent `main` work has not been tagged as a new release.

## Open changes and their gates

| Change | Current result | Gate before adoption |
|---|---|---|
| [PR 143](https://github.com/advanced-forestry-systems/fvs-modern/pull/143), Python engine interface | Fixed shared-library owner tracking, added a regression test, and passed 13 tests using a matching FVSne library on Cardinal. The updated branch passes Linux, macOS, and Windows CI. | Independent code review and merge; retest the merged commit. |
| [PR 144](https://github.com/advanced-forestry-systems/fvs-modern/pull/144), FIA plot-size expansion | Corrects the condition expansion divisor by tree plot size; the reported reruns compare seeded basal area with FIA `BALIVE`. Updated onto current `main`; current platform results are on the linked PR. | Review the expansion rule and rerun provenance, then merge before the catalog. |
| [PR 145](https://github.com/advanced-forestry-systems/fvs-modern/pull/145), catalog v8 | Stacked on PR 144. The primary recommendation rule is unchanged in the submitted comparison, but alternative-rule and stress labels change. Updated onto the new PR 144 head; current platform results are on the linked PR. | Merge only after PR 144, review the v8 tables and release package against the rerun. |
| [PR 146](https://github.com/advanced-forestry-systems/fvs-modern/pull/146), Greg 2026 refit | Updated onto current `main`; Linux, macOS, and Windows checks pass. The 36-scenario Maine stress result reports northern-hardwood and oak-pine dominant heights about 25–35% above native FVSne. Its site-ordering flags reflect the absence of a site term. | Reconcile height response and intended site scope against held-out measurements before selecting this arm as a default. Review the coefficient provenance and climate definitions. |

PR 92 is a large, conflicting integration branch. The smaller research PRs 100 and 102–110 overlap parts of that proposal. Decide their disposition by comparing each feature and its tests with current `main`; do not treat an open PR or a catalog label as evidence of release readiness.

## Scientific release gate

1. Freeze the exact engine route, coefficient files and hashes, climate inputs, species mapping, initial-condition expansion, random seed, code commit, and geographic scope. Keep experimental Greg and CONUS arms explicitly labeled until validated.
2. Use independent FIA remeasurement plots and available long-term field plots. Split by plot or site before fitting or tuning. Report observed-versus-predicted basal area, trees per hectare, quadratic mean diameter, dominant height, and volume at actual remeasurement intervals, with bias, MAE/RMSE, interval coverage, and plot-cluster bootstrap intervals. Set equivalence margins before examining the test set.
3. Audit sensitivity by species, forest type, site, state, starting density, and age. Explain failed strata and long-horizon behavior rather than collapsing them into an overall score.
4. Run a shared seeded cohort through the intended fvs-modern configuration and the standalone [FVS-CONUS projector](https://github.com/holoros/fvs-conus-pipeline/blob/main/RELEASE_READINESS_2026-11-19.md). Compare tree and stand trajectories at years 0, 5, 50, and 100; distinguish defined model differences from implementation drift.
5. Regenerate affected manuscript figures, tables, and release artifacts from the frozen route. Record the source/data manifest and review decision before a new tag or DOI deposit.

The existing Bakuzis stress tests are useful diagnostics, but they are not an independent remeasurement validation or a substitute for the gates above.
