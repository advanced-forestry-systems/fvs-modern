# Implementation notes, CONUS configuration catalog (October 28 to 29, 2026)

## Deviations from the pre-registered plan

1. Reineke gate not applied below 10 evaluated trajectories (NA passes). Pre-registered as "undefined"; the code lets an undefined rate through, which admits TT default GOMPIT on 5 evaluated trajectories. Kept as written and stated in the README and the memo rather than changed after the fact.
2. Stress-test alternative Reineke gate (Wilson upper bound under 25 percent, at least 30 evaluated) added after the interim scoring (01:45 EDT, 150 to 200 stands per variant) and before the final scoring. Primary verdicts unchanged; both carried.
3. `INSUFFICIENT SAMPLE` status for a variant whose best arm has fewer than 100 stands; `RECOMMENDED (low sample)` when the recommended arm has fewer than 100. PN (36 stands) fails the BA gate only through the Wilson bound (one stand over on the GOMPIT arms), which is a sample-size failure in effect; AK and NC have no FIA thresholds and are NOT EVALUATED.
4. Trajectory ids collided across variants in the engine runs (C0001 in every variant); fixed at concatenation by prefixing the variant. Found when NE and PN returned no per-trajectory Reineke rows in the test scoring.
5. ACD thresholds set equal to NE (no ACD plots in the pairs file), stated as an assumption.
6. Runs re-submitted twice: once for the pairs file's double-cast PLT_CN (recovered from the PLOT table by state, county, plot and inventory year), once for NaN STDAGE crashing the stand builder. The second array (15968171) was cut at 3 h wall time with no checkpoints; the third (15969413) wrote checkpoints every 25 stands under a 10 h limit.
7. Western variants ran at about 8 s per projection on shared Cardinal nodes against 0.5 s for NE; the cause was not investigated.

8. Output red team (October 29, after the first catalog was pushed as PR 120) confirmed four defects, all repaired and the catalog regenerated: (a) the under-stocking gate on the 100 year peak was vacuous (peak includes year 0; every arm passed by construction), now the median year-50 basal area against the FIA p50; (b) an untestable Reineke gate (fewer than 10 evaluated) passed as NA, now NOT TESTED and inadmissible under both rules (this removes TT default GOMPIT, which had been RECOMMENDED on 5 evaluated trajectories); (c) seed-gate failures were scored, now excluded, and the Reineke counts use the same stand set; (d) calibrated arms identical to default (ACD, IE: the calibrated config fails to load because per-species BAMAX is not a stand-level keyword; SO: the config carries nothing for the variant) were counted as evidence, now flagged `calibrated_distinct = FALSE` and excluded. The selector now carries the Reineke-not-testable, not-distinct and rule-sensitive notes and keeps the NE near-term preference when the caller constrains the arm.
9. The exact 0.75 x BALIVE under-seeding on failing stands in 15 variants was the condition proportion: tree expansion is per plot acre and BALIVE per condition acre. Trees are now expanded by CONDPROP_UNADJ. The 1.25 to 1.64 x over-seeding in SN and CS was the runner's all-plot fallback when the condition's trees were not found; removed.
10. Species crosswalk: unmapped species are relabelled to the variant's other-softwood or other-hardwood code (FIA 299 or 998, or 999 where that is the only catch-all) instead of aborting the stand; counts and the basal-area share relabelled are in the seed gate file.

## Discovered unknowns

- Species crosswalk gaps (SN 165 of 300 stands, BM 76, CS 73 in the first run) are handled by relabelling; which species carry the relabelled basal area is reported, not audited.
- Automatic establishment in EM and IE makes the within-stand Reineke slope uninformative there; NOAUTOES or an ingrowth-aware slope is needed before those rows mean anything.
- Dominant height is non-monotone in most trajectories (T1 HT_dom flags 4,071 of 4,399 in the engine arms), which is expected when tall trees die and is not gated.
- Per-variant near-term accuracy exists only for the Northeast cohort.

## Failed repairs

- Parallel library builds on firebreather with four private source copies still produced a 548-object NE library that segfaulted; rebuilt serially (555 objects). Cause not isolated beyond the builder's shared `/tmp/stub_*.err` and in-place header shadowing.

## Deferred

- Posterior draws for the catalog arms. Regeneration treatment. Per-variant near-term scorecard. Projector arms rerun on the common stand sample.
