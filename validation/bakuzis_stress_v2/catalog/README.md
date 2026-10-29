# CONUS configuration catalog for fvs-modern (October 28, 2026)

This folder holds the evidence and the rule behind `config/catalog/`, which `fvs2py.catalog.select_config(variant)` reads. The catalog answers one question per FVS variant: which survival model (the variant's native mortality or GOMPIT crown-closure survival) and which parameter set (fvs-modern default, fvs-modern Bayesian calibrated MAP, or for comparison the CONUS trait-mediated projector arms b1 and b2) stays inside the observed envelope of the FIA fitting data over a 50 and a 100 year projection, and which of the admissible arms is the simplest.

## Common stand sample

For every variant, up to 300 FIA plots were drawn (seed 20261028) from the fitting data's own variant assignment, the plot-level roll-up of `conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds` (`fia_variant_assignment.csv`, 17 variants with plots; PN 38, CA 142, TT 142, WC 226, all others 300). The Acadian variant uses the Northeast-assigned plots in Maine, New Hampshire and Vermont, an assumption, since the pairs file carries no ACD assignment. Variants with no plots in the pairs file (KT, OC, OP, NC for the engine runs) are reported as not evaluated rather than filled from another source. Each plot's first-measurement condition was written to FVS_StandInit and FVS_TreeInit through the repaired `fia_stand_generator.py` (PR 117: `basal_area_factor -1`, `inv_plot_size 1`, `brk_dbh 999`, so TPA_UNADJ passes through), with the public FIA plot coordinate, the condition's site species (SISP) where present, and the variant's USFS region. The same stands were projected 100 years at a 5 year step under four engine arms, default and calibrated MAP by native and GOMPIT survival, with no regeneration keywords and no posterior draws. A year-0 gate records whether each stand entered the engine within 20 percent of its FIA live basal area.

The CONUS projector arms (b1 and b2 by native and GOMPIT) come from the October 27 rescore (`cardinal_rescore/`), 1,500 stands per variant drawn from the projector's own stand set. They are listed in the same gate table for comparison, with the caveat that they were not run on the common sample and not through this engine, so differences between them and the engine arms confound engine with parameter set.

## Thresholds

Observed envelopes are plot-level percentiles of the fitting data by variant (`fia_variant_thresholds_plot.csv`): basal area p50, p90, p99 and p99.9 in m² ha⁻¹ and the SDI p99.9 in trees ha⁻¹ (the file's per-acre SDI converted). They are age-censored population extremes, so a variant whose inventory is young or managed (SN, CS, PN plantations) has a lower ceiling than its unmanaged stands could reach; the 50 year gate limits how far that bias reaches, and the 100 year shares are reported beside it.

## Pre-registered rule (fixed before any catalog run was scored; the plan red team's repairs are folded in)

An arm is admissible for a variant when all four hold: fewer than 5 percent of stands exceed the BA p99.9 by year 50 and the Wilson 95 percent upper bound of that share is under 10 percent; the same for the SDI p99.9; the Reineke flag rate among trajectories evaluated above the relative density gate (SDI/SDImax 0.55) is under 10 percent, not gated when fewer than 10 trajectories are evaluated (which is how TT default GOMPIT is admitted on 5 evaluated trajectories); and the median peak basal area over 100 years is at or above the FIA p50, which catches an arm that passes the ceilings by killing the stand. Among admissible arms the simplest wins (default before calibrated before b1 before b2; native before GOMPIT), except where per-variant near-term evidence separates the survival models, which at present is only the Northeast cohort of Paper 2 (basal area bias +10.78 percent native against +4.62 percent GOMPIT), so NE prefers GOMPIT among its admissible arms. A variant with no admissible arm is NOT RECOMMENDED for long-horizon use and the least-bad arm by 50 year BA exceedance is named for information.

## Stress-test alternative and rule sensitivity

The protocol's alternative Reineke gate, added after the interim scoring showed the western pattern and before the final scoring, uses the Wilson 95 percent upper bound of the flag rate under 25 percent with at least 30 evaluated trajectories and does not gate below that count. `catalog_gates_by_arm.csv` carries both verdicts (`admissible`, `admissible_alt`) and `catalog_recommendation_by_variant.csv` flags `rule_sensitive` where the admissible sets differ (eight variants). The primary rule is what `recommended` and `status` report; the alternative admits western arms by declining to test their self-thinning, which is why it is shown and not shipped as the default.

## Results in one paragraph (October 29)

Eighteen variants ran (4,299 stands). Under the primary rule eight have admissible fvs-modern arms (ACD, BM, CS, NE, SN on all four; LS on the two calibrated arms; TT default GOMPIT; WS calibrated native) and ten do not (CA, CI, CR, EC, EM, IE, PN, SO, UT, WC), nearly all on the Reineke gate alone; the magnitude gates pass in 13 of 18 variants at or below 1.3 percent of stands above the BA p99.9 by year 50, and no arm under-stocks. Under the alternative, CA, CI, CR, EC, SO and WC gain admissible arms and EM, IE, PN and UT do not. No CONUS projector arm is admissible under the primary rule. Section 15 of the memo has the numbers.

## Files

- `conus_catalog_runs.py`, `catalog_runs.slurm`: the engine runs (Cardinal array 15967542, after library build 15967530 with the GOMPSURV fix).
- `catalog_gates.R`, `finalize_catalog.slurm`: scoring with `bakuzis_stress_v2.R --source=long` and the gate table (job 15967981).
- `fia_variant_assignment.csv`, `fia_variant_thresholds_plot.csv`: the stand sample frame and the thresholds.
- `gates/catalog_gates_by_arm.csv`, `gates/catalog_recommendation_by_variant.csv`: the catalog itself (copied to `config/catalog/`).
- `catalog_seedgate_all.csv`, `score_catalog/`: the year-0 gate and the full v2 scorecards.

## What the catalog does not cover

Per-variant near-term accuracy. The only near-term scorecard is the Northeast cohort in Paper 2, so outside NE the tie-break is simplicity alone and the catalog says so in `nearterm_evidence`. Posterior uncertainty: the calibrated MAP is a point estimate here; the Paper 1 ensemble shows the 100 draw band for NE and ACD. Regeneration: none of the arms establishes regeneration, so a 100 year plateau or decline in stems is partly a design choice. Species groups: the common sample is unstratified, so the Eichhorn, site ordering and Sukachev relations are not part of the catalog gates.
