# Implementation notes, CONUS configuration catalog (October 28 to 29, 2026)

## Deviations from the pre-registered plan

1. Reineke gate not applied below 10 evaluated trajectories (NA passes). Pre-registered as "undefined"; the code lets an undefined rate through, which admits TT default GOMPIT on 5 evaluated trajectories. Kept as written and stated in the README and the memo rather than changed after the fact.
2. Stress-test alternative Reineke gate (Wilson upper bound under 25 percent, at least 30 evaluated) added after the interim scoring (01:45 EDT, 150 to 200 stands per variant) and before the final scoring. Primary verdicts unchanged; both carried.
3. `INSUFFICIENT SAMPLE` status for a variant whose best arm has fewer than 100 stands; `RECOMMENDED (low sample)` when the recommended arm has fewer than 100. Applies to AK (97, projector only). PN (36) is NOT RECOMMENDED on the gates themselves, not on sample size.
4. Trajectory ids collided across variants in the engine runs (C0001 in every variant); fixed at concatenation by prefixing the variant. Found when NE and PN returned no per-trajectory Reineke rows in the test scoring.
5. ACD thresholds set equal to NE (no ACD plots in the pairs file), stated as an assumption.
6. Runs re-submitted twice: once for the pairs file's double-cast PLT_CN (recovered from the PLOT table by state, county, plot and inventory year), once for NaN STDAGE crashing the stand builder. The second array (15968171) was cut at 3 h wall time with no checkpoints; the third (15969413) wrote checkpoints every 25 stands under a 10 h limit.
7. Western variants ran at about 8 s per projection on shared Cardinal nodes against 0.5 s for NE; the cause was not investigated.

## Discovered unknowns

- Crosswalk gaps in fvs2py's species mapping: 165 of 300 SN stands, 76 BM and 73 CS did not run (SPCD 403, 409, 682, 15, 475 among others). A per-variant crosswalk audit is the first thing to do before the southern catalog row is trusted.
- Year-0 basal area above BALIVE by a median factor of 1.25 (CS) and 1.42 (SN) on the seed-gate failures in those two variants only; every other variant's failures are under-seeded (0.75 to 0.79). Not resolved.
- Automatic establishment in EM and IE makes the within-stand Reineke slope uninformative there; NOAUTOES or an ingrowth-aware slope is needed before those rows mean anything.
- Dominant height is non-monotone in most trajectories (T1 HT_dom flags 4,071 of 4,399 in the engine arms), which is expected when tall trees die and is not gated.
- Per-variant near-term accuracy exists only for the Northeast cohort.

## Failed repairs

- Parallel library builds on firebreather with four private source copies still produced a 548-object NE library that segfaulted; rebuilt serially (555 objects). Cause not isolated beyond the builder's shared `/tmp/stub_*.err` and in-place header shadowing.

## Deferred

- Posterior draws for the catalog arms. Regeneration treatment. Per-variant near-term scorecard. Projector arms rerun on the common stand sample.
