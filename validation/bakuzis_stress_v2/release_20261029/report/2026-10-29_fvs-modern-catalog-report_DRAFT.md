# A configuration catalog for fvs-modern across the FVS variants of the conterminous United States: Bakuzis admissibility gates, stress test, and selector

Aaron R. Weiskittel, Center for Research on Sustainable Forests, University of Maine. CRSF Technical Report, October 29, 2026.

## Summary

fvs-modern can run every Forest Vegetation Simulator (FVS) variant with either the variant's native mortality or GOMPIT crown closure survival, and with either the fvs-modern default parameters or a Bayesian calibrated maximum a posteriori (MAP) set. A user needs to know which of those four combinations stays biologically plausible over the horizons FVS is used for. This report describes a configuration catalog that answers that question per variant from evidence, the pre-registered rule behind it, the red team repairs that changed it, the stress test run after it shipped, and the selector that reads it. For each of 18 variants, up to 300 FIA plots were projected 100 years under all four arms and gated on the FIA fitting envelope (the plot level basal area p99.9 and the stand density index p99.9), on Reineke self-thinning, and on a stocking lower bound. Seven variants have at least one admissible fvs-modern arm (ACD, BM, CS, LS, NE, SN, WS) and eleven do not (CA, CI, CR, EC, EM, IE, PN, SO, TT, UT, WC), nearly all on the self-thinning axis alone. A stress test of 144 alternative rules and 500 stand bootstraps changed none of the 72 arm verdicts; eight admissible and eight not-admissible arms are fragile, and the catalog carries that label. The catalog, the gate table, the scripts, and the per variant Bakuzis scorecards are on the fvs-modern main branch (commit 1a520f5e) and in the Zenodo record that accompanies this report.

## Why a catalog

The Bakuzis matrix (Leary 1997) tests whether a growth model's long projections obey the law-like relations among stand properties: height and age ordering by site, Reineke's stems against mean diameter, the Sukachev effect, Eichhorn's rule, and yield against density. The earlier fvs-modern work applied that matrix to a synthetic 36 cell design for the Northeast and Acadian variants and to the Northeast Paper 2 cohort, and the fvs-conus projector was scored the same way. Those tests settled the two northeastern variants and said nothing about the other twenty. At the same time, the fvs-modern survival switch and the calibrated parameter set had never been tested together outside the Northeast, and a code defect found on October 28 (a six argument GOMPSURV interface updated in one variant's mortality routine and not in the other 21) meant that every GOMPIT result produced on a variant other than NE, CS and LS between July 25 and October 28 came from an undefined read. With that fixed and all 22 variant libraries rebuilt, the question became whether a user could pick a survival model and a parameter set for their variant on evidence rather than by habit. The catalog is that evidence.

## Stand sample and projection arms

The sample frame is the plot level roll up of the CONUS remeasurement pairs file used to fit the fvs-conus projector, which assigns each plot to an FVS variant. Up to 300 plots per variant were drawn with a fixed seed (17 variants carry plots; PN has 38, CA and TT 142, WC 226, the rest 300). The Acadian variant has no assignment in that file and uses the Northeast plots in Maine, New Hampshire and Vermont, which is an assumption the catalog states. Five variants with no plots in the frame (AK, KT, NC, OC and OP) are reported as not evaluated rather than filled from another source. Each plot's first measurement condition was written to FVS stand and tree initialization tables with the public FIA coordinate, the condition's site species and the variant's Forest Service region. Trees are those of the sampled condition expanded to the condition acre by CONDPROP_UNADJ, which is the definition of the FIA live basal area the gates compare against, and species the compiled variant does not carry are relabelled to its other softwood or other hardwood code as FVS itself does on input. A year zero gate records whether each stand entered the engine within 20 percent of its FIA live basal area; only stands that pass enter the gates. That gate passes 95 to 100 percent of stands in 15 variants, 89 percent in LS, 78 in CS and 54 in SN, where the failing stands enter above their FIA basal area at a median of 1.46 times with the input verified, an engine side year zero effect that remains open. The catalog rests on 4,485 seed passing stands.

Every stand was projected 100 years at a 5 year step under four arms: default and calibrated MAP parameters, each with native and GOMPIT survival, with no regeneration keywords and no posterior draws. EM and IE run with NOAUTOES because their establishment models otherwise add stems in every projection; no other variant establishes regeneration. The runs were Cardinal array jobs 15976824 and 15984133 after library build 15967530 with the GOMPSURV fix. The fvs-conus projector arms (b1 and b2 by native and GOMPIT) are listed in the gate table for comparison from an October 27 rescore of 1,500 stands per variant, with the caveat that they were neither run on the common sample nor through this engine.

## The pre-registered rule

The rule was fixed before any catalog run was scored and red teamed blind before the runs started; the seven confirmed design defects (one common stand sample for every engine arm, a 50 year primary horizon with Wilson bounds, an under stocking lower bound, one shared variant assignment, and confirmation that the projector's GOMPIT is R side and unaffected by the Fortran defect) were folded in. An arm is admissible for a variant when all of the following hold on seed passing stands: fewer than 5 percent of stands exceed the FIA basal area p99.9 by year 50 and the Wilson 95 percent upper bound of that share is under 10 percent; the same for the SDI p99.9; the Reineke flag rate among trajectories evaluated above a relative density of 0.55 is under 10 percent on at least 10 evaluated trajectories, with fewer evaluated counted as not tested and inadmissible; the median year 50 basal area is at or above the FIA p50; and a calibrated arm's trajectories differ from the default arm's. Among admissible arms the simplest wins (default before calibrated, native before GOMPIT), except where per variant near term evidence separates the survival models, which at present is only the Northeast cohort of Paper 2, so NE prefers GOMPIT among its admissible arms. A variant with no admissible arm is not recommended for long horizon use and the least bad arm by 50 year basal area exceedance is named for information.

The thresholds are plot level percentiles of the fitting data by variant. They are age censored population extremes, so a variant whose inventory is young or managed has a lower ceiling than its unmanaged stands could reach; the 50 year gate limits how far that bias reaches, and where a variant has fewer than 1,000 fitting plots the p99.9 is the sample maximum, which the gate table discloses. An alternative Reineke gate, the Wilson upper bound of the flag rate under 25 percent on at least 30 evaluated trajectories, was added after the interim scoring showed the western pattern and before the final scoring; both verdicts are carried and the catalog flags the eight variants where they differ. The catalog ships under the primary rule, a decision of October 29, 2026, because the alternative admits western arms by declining to test their self-thinning.

A first catalog was pushed on the morning of October 29 and red teamed blind. Four defects that each made the gates easier than the documentation claimed were confirmed and the catalog was regenerated rather than patched: the under stocking bound used the 100 year peak, which includes year zero and could not fail; a Reineke gate with fewer than 10 evaluated trajectories counted as a pass; stands that failed the year zero gate were scored; and the calibrated arms for ACD and IE were identical to default because the keyword loader raised on those variants' per species BAMAX array and aborted the calibrated block. The loader now skips that array with a warning, both variants were rerun, and a distinctness check guards every variant. The same review found the cause of the two seeding failure modes carried from earlier work: partial conditions expanded per plot acre rather than per condition acre (the exact 0.75 under seeding), and an all plot fallback in the runner when a condition's trees were not found (the southern over seeding). Both are gone.

## Results

Table 1 gives the verdict by variant and Figure 1 places all 72 engine arms in the space of the two gates that decide. Under the primary rule seven variants have admissible fvs-modern arms: ACD, BM, CS, NE and SN on all four, LS on the two calibrated arms, and WS on calibrated native alone. Eleven do not. In CA, CI, CR, EC, IE, SO, TT and UT every arm passes the magnitude and stocking gates and fails only on the Reineke axis; in EM and WC the arms that pass magnitude fail on Reineke; PN fails the basal area gate on 38 stands. The magnitude test is passed almost everywhere: the 50 year share above the basal area p99.9 is at or below 1.7 percent in 14 of 18 variants, the SDI ceiling is never exceeded by more than 4.7 percent of stands, and the year 50 lower bound passes in every arm. The exceptions are the native arms in EM (22 percent default, 10 percent calibrated, with 8 percent on default GOMPIT), LS (14 and 10 percent on the default arms), WC (11 and 9 percent native) and PN (34 and 32 percent native).

Two patterns sit behind the Reineke failures. In CA, TT, UT, WC, EM and CR the default native arm does not self-thin: median within stand slopes of log stems on log quadratic mean diameter are -0.83, -0.84, -1.12, -1.11, -1.04 and -1.46 against Reineke's -1.605, with flag rates of 27 to 84 percent, and calibrated GOMPIT moves the medians to -1.59 to -1.90 while the flag rates stay at 17 to 41 percent because the slopes are dispersed across stands. In EM and IE the regeneration confound is gone with NOAUTOES, stems now decline, and the flag rates remain at 13 to 77 percent in EM and 32 to 56 percent in IE, so the self-thinning failure stands on its own. Across the West only 18 to 202 trajectories per variant reach the relative density gate on the best populated arm, and in SO and TT the flag rate rests on 17 to 18 and 7 to 31 trajectories. The projector arms are a different population: their median 50 year basal area exceedance across variants is 19 percent for b1 native, 6 percent for b2 native and 2 to 3 percent under GOMPIT, and none is admissible under the primary rule.

**Table 1.** Catalog verdict by FVS variant under the primary rule (`config/catalog/catalog_recommendation_by_variant.csv`). Stands are seed passing stands on the engine arms. The bootstrap probability is the share of 500 stand resamples in which the variant keeps at least one admissible engine arm; the grid share is the share of 144 alternative rule cells that admit one (`stress_by_variant.csv`).

| Variant | Stands | Recommended arm | Status | Admissible engine arms | Bootstrap probability | Grid share | Rule sensitive |
|---|---|---|---|---|---|---|---|
| ACD | 299 | default, native | Recommended | all four | 1.00 | 0.50 | no |
| BM | 298 | default, native | Recommended | all four | 1.00 | 0.75 | yes |
| CA | 142 | calibrated, GOMPITᵃ | Not recommended | none | 0.00 | 0.00 | yes |
| CI | 297 | calibrated, nativeᵃ | Not recommended | none | 0.49 | 0.67 | yes |
| CR | 300 | calibrated, GOMPITᵃ | Not recommended | none | 0.35 | 0.52 | yes |
| CS | 233 | default, native | Recommended | all four | 1.00 | 1.00 | no |
| EC | 300 | calibrated, GOMPITᵃ | Not recommended | none | 0.16 | 0.33 | yes |
| EM | 300 | calibrated, GOMPITᵃ | Not recommended | none | 0.00 | 0.03 | no |
| IE | 298 | calibrated, GOMPITᵃ | Not recommended | none | 0.00 | 0.00 | no |
| LS | 266 | calibrated, native | Recommended | calibrated native, calibrated GOMPIT | 1.00 | 0.54 | no |
| NE | 285 | default, GOMPIT | Recommended | all four | 0.99 | 0.58 | no |
| PN | 38 | calibrated, GOMPITᵃ | Not recommended | none | 0.00 | 0.00 | no |
| SN | 161 | default, native | Recommended | all four | 1.00 | 1.00 | no |
| SO | 300 | calibrated, GOMPITᵃ | Not recommended | none | 0.42 | 0.33 | yes |
| TT | 142 | projector b2, GOMPITᵃ | Not recommended | none | 0.24 | 0.17 | yes |
| UT | 300 | calibrated, GOMPITᵃ | Not recommended | none | 0.01 | 0.00 | no |
| WC | 226 | calibrated, GOMPITᵃ | Not recommended | none | 0.03 | 0.13 | yes |
| WS | 300 | calibrated, native | Recommended | calibrated native | 0.90 | 0.50 | yes |

ᵃ Least bad arm by 50 year basal area exceedance, named for information; no arm is admissible. AK and NC have no fitting plots and are NOT EVALUATED.

![Figure 1](../fig/fig_gate_space.png)

**Figure 1.** The 72 engine arms in the space of the two gates that decide the catalog: the Reineke flag rate among trajectories above the relative density gate against the share of stands above the FIA basal area p99.9 by year 50, with Wilson 95 percent intervals. Green symbols are admissible under the primary rule. Arms at a negative flag rate have fewer than 10 evaluated trajectories and are not testable.

## Stress test of the shipped catalog

After the catalog shipped, the verdicts were stress tested against the analyst choices in the rule that could have gone another way, from the trajectory level data and without new engine runs (`stress_catalog.R`, Cardinal job 15990793). The first part is a rule grid of 144 cells: horizon 50 or 100 years, exceedance tolerance 3, 5 or 10 percent, Wilson upper bound 10 or 15 percent, basal area ceiling at the p99.9 or the p99, Reineke flag tolerance 10, 15 or 25 percent, and a minimum of 10 or 30 evaluated trajectories. The second is a stand bootstrap, 500 resamples with replacement per variant under the shipped rule. Recomputing the shipped rule inside the stress script reproduced all 72 catalog verdicts, which is also an independent check of the gate script.

The bootstrap splits the 23 admissible arms into 15 stable (admissible in at least 90 percent of resamples and at least half the grid) and 8 fragile, and the 49 not admissible arms into 41 stable and 8 fragile (Table 2, Figure 2). At the variant level, which is what the selector acts on, the seven recommended variants keep at least one admissible engine arm in 100 percent of resamples for ACD, BM, CS, LS and SN, 99 percent for NE and 90 percent for WS. The fragile not admissible arms are all among the thirteen that a 25 percent Reineke tolerance admits, and the bootstrap puts seven of the eight below even odds (CI default native sits at 0.51), so the primary rule's refusal to recommend them is not an artifact of one draw of stands. CA, IE, PN and UT are admitted by no cell of the grid. The grid shows where the rule does its work: the Reineke tolerance is the only axis that changes recommended variants (25 percent admits 13 further arms in CI, CR, EC, SO, TT, WC and WS); the p99 ceiling admits only five arms anywhere, so the shipped p99.9 is the permissive end of the magnitude axis; the 100 year horizon removes one arm (LS calibrated native) and no variant; the 3 percent tolerance removes four arms and no variant; the minimum of 30 evaluated trajectories changes nothing. In NE the stable arm is default native (98 percent) while the shipped preference for GOMPIT rests on the Paper 2 near term cohort, so the catalog's NE choice trades long horizon stability for near term evidence and the selector's note says so.

**Table 2.** The sixteen arms whose verdict is fragile under the stand bootstrap or the rule grid (`stress_by_arm.csv`). Bootstrap probability is the share of 500 stand resamples in which the arm is admissible under the shipped rule; grid share is the share of 144 rule cells that admit it.

| Variant | Arm | Shipped verdict | Bootstrap probability | Grid share | Stability |
|---|---|---|---|---|---|
| ACD | calibrated, GOMPIT | admissible | 0.58 | 0.33 | fragile |
| ACD | default, GOMPIT | admissible | 0.58 | 0.33 | fragile |
| CS | default, GOMPIT | admissible | 0.77 | 0.50 | fragile |
| LS | calibrated, native | admissible | 0.84 | 0.25 | fragile |
| NE | calibrated, GOMPIT | admissible | 0.57 | 0.58 | fragile |
| NE | calibrated, native | admissible | 0.74 | 0.54 | fragile |
| NE | default, GOMPIT | admissible | 0.84 | 0.33 | fragile |
| WS | calibrated, native | admissible | 0.90 | 0.50 | fragile |
| CI | default, native | not admissible | 0.51 | 0.67 | fragile |
| CR | calibrated, native | not admissible | 0.24 | 0.47 | fragile |
| EC | default, native | not admissible | 0.16 | 0.33 | fragile |
| SO | calibrated, GOMPIT | not admissible | 0.44 | 0.33 | fragile |
| SO | calibrated, native | not admissible | 0.19 | 0.17 | fragile |
| SO | default, GOMPIT | not admissible | 0.46 | 0.33 | fragile |
| SO | default, native | not admissible | 0.24 | 0.17 | fragile |
| TT | calibrated, GOMPIT | not admissible | 0.25 | 0.17 | fragile |

![Figure 2](../fig/fig_arm_stability.png)

**Figure 2.** Stability of every shipped verdict: the share of 144 rule cells admitting each arm against the bootstrap probability that it is admissible under the shipped rule. No verdict is overturned; the eight fragile admissible arms are labelled.

## The selector

The catalog lives in `config/catalog/` as two tables and is read by `fvs2py.catalog.select_config(variant, survival=None, parameter_set=None)`. The call returns the arm, the configuration version to pass to the engine, the environment variables that switch GOMPIT on, the admissible set, and a note when the caller forces an arm the catalog did not admit, when the Reineke axis was not testable, when a calibrated arm is not distinct from default, or when the verdict is rule sensitive. When a variant has no admissible arm the selector returns the least bad engine arm with the not recommended status attached, falling back from a projector arm where that is the least bad overall. AK and NC return not evaluated. The stress labels sit beside the catalog in `validation/bakuzis_stress_v2/catalog/stress/` and a user choosing a fragile admissible arm should read them.

## Limitations and open items

Per variant near term accuracy is the gap. The only near term scorecard is the Northeast cohort of Paper 2, so outside NE the tie break among admissible arms is simplicity alone and the catalog says so in its `nearterm_evidence` column. The calibrated MAP is a point estimate here; the Paper 1 ensemble shows the posterior band for NE and ACD only. None of the arms establishes regeneration, so a 100 year plateau or decline in stems is partly a design choice. The common sample is unstratified, so the site ordering, Eichhorn, Sukachev and yield density relations are confounded by stand selection and are reported in the per variant scorecards but not gated; the synthetic 36 cell matrix remains the clean test of those relations. The 263 stands that fail the year zero gate, 256 of them eastern and nearly all entering the engine above their FIA basal area with verified input, are excluded and unexplained. The thresholds are the sample maximum in small variants. The projector b1 GOMPIT arm has no NE or SN rows in the October 27 rescore. These are stated in `implementation-notes.md` with the thirteen numbered deviations from the pre-registered plan.

## Data and code availability

The catalog tables, gate table, stress tables, per variant scorecards, scripts, run index and the running memo (Sections 1 to 16) are on the fvs-modern main branch of the advanced-forestry-systems organization at commit 1a520f5e (pull requests 118 to 126), under `config/catalog/` and `validation/bakuzis_stress_v2/`. The Zenodo record that accompanies this report (concept DOI 10.5281/zenodo.23777723, version 1.0.0 at 10.5281/zenodo.23777724) carries the same tables and scripts with their checksums, this report, and the presentation. FIA plot identifiers and coordinates do not appear in any deposited file; the stand sample is described by variant and count only.

## References

Crookston, N.L., and G.E. Dixon. 2005. The forest vegetation simulator: a review of its structure, content, and applications. Computers and Electronics in Agriculture 49: 60 to 80.

Leary, R.A. 1997. Testing models of unthinned red pine plantation dynamics using a modified Bakuzis matrix of stand properties. Ecological Modelling 98: 35 to 46.

Reineke, L.H. 1933. Perfecting a stand-density index for even-aged forests. Journal of Agricultural Research 46: 627 to 638.
