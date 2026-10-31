# D5 close-out, October 30, 2026

Two open items from the Bakuzis stress test v2 (Paper 1, Section 4.6) and the
CONUS catalog run-off.

## A. Narrow ACD within-stand band (`a_acd_band/`)

In the Paper 1 ensemble 81 of 180 ACD stands gave the same year 100 basal area
under all 100 posterior draws. This is a property of the ACD variant, not a
keyword failure.

* The draws are read. Each draw keyfile carries 74 to 91 BAIMULT and MORTMULT
  records, they are echoed in the FVS output (`a_keyword_echo.csv`), and tree
  lists differ between draws from cycle 1 (`a_treelists_cycles.csv`).
* Year 100 stem counts differ among draws in 80 of the 81 stands, and basal
  area diverges at some year in 80 of the 81 (`yr100_draw_spread_by_stand.csv`).
* `src-converted/acd/morts.f90` applies the mortality multiplier only to the
  background rate (line 522, `IF(RIP .EQ. RN) X=1.0`), which governs a stand
  below 55 percent of its maximum SDI. Under GOMPIT the multiplier is dropped
  (lines 524 to 529;
  `vls/morts.f90` lines 599 to 608).
* The ACD posterior has no SDIMAX term, so every draw approaches the same
  density ceiling. In a three stand probe (`a_runs_summary.csv`) BAIMULT 1.5
  changed year 100 basal area by at most 1 ft2/ac, MORTMULT 2 changed no
  stand's year 100 basal area, and an all-species SDIMAX of 150 or 600 moved
  all three stands from 178 to 197 ft2/ac under the MAP parameters to 84 to
  107 or 219 to 278 ft2/ac.

The rerun ensemble (below) still has 81 of 180 ACD stands and 14 of 180 NE
stands with a zero year 100 band.

## B. Plot-size expansion (`b_plotsize/`)

The stand generator and the catalog runner expanded every tree by
1/CONDPROP_UNADJ. FIA live basal area (BALIVE) expands microplot, subplot and
macroplot trees by the condition's share of that plot size (MICRPROP_UNADJ,
SUBPPROP_UNADJ, MACRPROP_UNADJ). `condition_expansion_divisor()` in
`calibration/python/fia_stand_generator.py` now does this, and the catalog
runner imports it. A tree is assigned to a plot size only when its TPA_UNADJ
is within 1 percent of 74.965, 6.018 or 0.9997; a tree that fits none, or
whose proportion is missing, is divided by CONDPROP_UNADJ (then 1) and
counted (`return_counts=True`; `condition_trees` logs a warning and the
catalog runner prints a line).

`b_expansion.py` calls the shipped function (review rerun, October 30, 2026,
after PR 144 review). Ratio of seeded BA to BALIVE on the multi-condition
stands: 125 of 125 ensemble and 971 of 972 catalog stands within 1e-4 of 1,
and the remaining catalog stand at 1.00012 (BALIVE 0.55 ft2/ac, so a
rounding residual; the CONDPROP expansion gives the same value). Fallback
counts (`n_no_plot_size`, `n_no_proportion` in the check files):

* Ensemble: 949 trees in 19 of 360 stands fit no plot size (2 of them
  multi-condition); 0 trees lack a proportion.
* Catalog: 4,714 trees in 231 of 4,748 stands fit no plot size (1
  multi-condition); 2 trees in 2 stands lack a proportion (single-condition).
* In all three multi-condition stands with fallback trees CONDPROP_UNADJ
  equals SUBPPROP_UNADJ and MICRPROP_UNADJ, and no stand's seeded basal area
  differs from the first D5 rule (TPA_UNADJ above 50 micro, below 2 macro):
  `shipped_vs_threshold_pct` is 0 on every row, so the reruns below stand.

| Sample | Off BALIVE by >2% before | Stands moved >2% | Largest move |
|---|---|---|---|
| Ensemble (360 stands) | 14 | 14 (3.9%) | 12.8% |
| Catalog (4,748 stands) | 263 of 972 multi-condition | 266 (5.6%) | 25.5% |

Both exceed the 1 percent rerun trigger, so both were rerun.

## Ensemble rerun (`ensemble_rerun/`)

`ensemble_full5.slurm` then `finish_ensemble5.slurm` (seed gate, aggregation,
`bakuzis_stress_v2.R`, `ens5_numbers.R`). All 360 stands now enter inside their
selection bin. `ens5_numbers_old_vs_new.csv` lists every Paper 1 number, old
and new; `ens4_vs_ens5_check.txt` confirms the laws, exceedances and scorer
verdicts. The rerun's law compliance table and species-group
scorecard summary are byte-identical to the tracked
`ensemble_regen4/agg/bakuzis_laws_compliance.csv` and
`ensemble_regen4/score_ens_full4b/ens_full4_20261030b_scorecard_summary.csv`,
so they are not repeated here.

## Catalog rerun (`catalog_rerun/`)

`catalog_runs_v8.slurm` then `finalize_catalog_v8.slurm` (concat, scorer,
gates, `compare_gates_v6_v8.R`). No primary-rule
verdict changed (`catalog_v6_v8_verdict_changes.csv`, which compares the
primary-rule admissibility, recommended arm and admissible arms, has no rows),
and the recommended arm, admissible arms and status of every variant are
unchanged (7 RECOMMENDED, 11 NOT RECOMMENDED, 2 NOT EVALUATED). The full
comparison, which adds the alternative rule and the stress test (stress v8,
run for PR 145; `validation/bakuzis_stress_v2/catalog/catalog_v6_v8_verdict_changes_full.csv`
in PR 145), finds five changed fields: under the alternative rule WC
default GOMPIT and WC calibrated GOMPIT are no longer admissible, so WC has no
alternative-rule admissible arm and is no longer rule sensitive, and WS
calibrated native moved from FRAGILE to STABLE admissible. Stress v8 has 17
stable admissible, 6 fragile admissible, 41 stable not admissible and 8
fragile not admissible arms (v6: 16, 7, 41, 8); in both, every primary-rule
admissible arm is admissible in at least 56 percent of 500 stand bootstraps and
every other arm in fewer than 50 percent. The largest change in an admissible
arm's median year 50 basal area is 0.69 m2/ha. Seed-gate rows within 2 percent
of BALIVE rose from 90.73 to 95.66 percent; all 18,992 rows pass in both.

## Identifiers

No file here contains FIA plot coordinates, PLT_CN or CN values. Stands are
indexed by sample keys only. FIA data stayed on OSC Cardinal.

`MD5SUMS` covers every file in this directory.
