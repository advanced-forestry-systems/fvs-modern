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

11. ACD and IE calibrated arms: `config_loader.generate_keywords` raised on the variants' own per-species BAMAX array (blkdat, not a calibrated quantity) and aborted the whole block, so the calibrated arms ran as default. The loader now skips BAMAX with a warning (`patch_config_loader.py`, applied to `config/config_loader.py`); ACD now emits 68 calibrated keyword lines and IE 49, and both variants were rerun (array 15984133).
12. EM and IE rerun with NOAUTOES (the runner appends it to the keyword block when `FVS_NOAUTOES=1`). Stems now decline over 50 years; the Reineke failures remain.
13. Over-seeding diagnostic (`diag_overseed.py`, `diag_overseed.csv`): for the failing SN, CS and NE stands the condition's tree sum per condition acre reproduces BALIVE (first case 90.4 against 90.3 ft2 ac-1), so the 1.3 to 1.8 x year-0 excess is inside the engine's year-0 accounting for those stands and is unrelated to relabelled species, sapling share or condition proportion. Next test: run one failing stand with a single cycle and compare the engine tree list at year 0 against the input records.
14. Shipped under the primary rule (Aaron, October 29, 14:36 EDT). Stress test after shipping (`stress_catalog.R`, job 15990793): 144-cell rule grid and 500-draw stand bootstrap, zero disagreements with the shipped verdicts, 8 of 23 admissible arms and 8 of 49 not-admissible arms labelled FRAGILE; the Reineke tolerance is the only axis that changes recommended variants. No verdict changed.
15. Over-seeding resolved (October 29, evening). FIA live records with a null TPA_UNADJ (P2A_GRM_FLG = Y) were filled with 6.0 TPA by `fia_stand_generator.build_tree_init`; the engine reproduces its input exactly (`diag_y0/`, jobs 15998331 and 15998396). Generator and runner now drop them and the seed gate carries `n_null_tpa` and a 1 ft² ac⁻¹ absolute floor (PR 129). Catalog regenerated on all 18 variants (array 15998817, finalize 16002593, stress 16003647, scorecards 16003648): 4,748 of 4,748 stands pass the year-0 gate, 3,617 records dropped from 371 stands, zero changes in the 72 arm verdicts, the recommendation table and the stress labels. Deviation 13's next test is closed.
16. P2A_GRM_FLG (October 30). The field is FIA's periodic to annual growth, removal and mortality flag (FEMC metadata for the FIA TREE table, https://www.uvm.edu/femc/data/archive/project/federal-forest-inventory-analysis-data-for/dataset/tree-level-data-gathered-through-forest/metadata). No FIADB guide reachable from here states that TPA_UNADJ is null for these records (the v4.0 and v6.0 Phase 2 guides fetched truncate before the TREE table), so the null expansion factor is reported as an observed property of the data (3,617 live records in 371 sampled eastern stands, `catalog_seedgate_v5_noid.csv`), not as a documented rule.

## Discovered unknowns

- Species crosswalk gaps (SN 165 of 300 stands, BM 76, CS 73 in the first run) are handled by relabelling; which species carry the relabelled basal area is reported, not audited.
- Dominant height is non-monotone in most trajectories (T1 HT_dom flags 4,071 of 4,399 in the engine arms), which is expected when tall trees die and is not gated.
- Per-variant near-term accuracy exists only for the Northeast cohort.

## Failed repairs

- Parallel library builds on firebreather with four private source copies still produced a 548-object NE library that segfaulted; rebuilt serially (555 objects). Cause not isolated beyond the builder's shared `/tmp/stub_*.err` and in-place header shadowing.

## Deferred

- Posterior draws for the catalog arms. Regeneration treatment. Per-variant near-term scorecard. Projector arms rerun on the common stand sample.
