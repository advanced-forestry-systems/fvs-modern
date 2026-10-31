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
  (lines 524 to 527; `vls/morts.f90` lines 600 to 604).
* The ACD posterior has no SDIMAX term, so every draw approaches the same
  density ceiling. In a three stand probe (`a_runs_summary.csv`) BAIMULT 1.5
  changed year 100 basal area by at most 1 ft2/ac, MORTMULT 2 changed nothing,
  and an all-species SDIMAX of 150 or 600 moved one stand from 178 to 107 or
  278 ft2/ac.

The rerun ensemble (below) still has 81 of 180 ACD stands and 14 of 180 NE
stands with a zero year 100 band.

## B. Plot-size expansion (`b_plotsize/`)

The stand generator and the catalog runner expanded every tree by
1/CONDPROP_UNADJ. FIA live basal area (BALIVE) expands microplot, subplot and
macroplot trees by the condition's share of that plot size (MICRPROP_UNADJ for
TPA_UNADJ > 50, MACRPROP_UNADJ for TPA_UNADJ < 2, SUBPPROP_UNADJ otherwise).
`condition_expansion_divisor()` in `calibration/python/fia_stand_generator.py`
now does this, falling back to CONDPROP_UNADJ and then 1, and the catalog
runner imports it. With the fix the ratio to BALIVE is 1.0000 on all 125
ensemble and all 972 catalog multi-condition stands (`b_expansion_report.txt`).

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
gates, `compare_gates_v6_v8.R`). No verdict field changed
(`catalog_v6_v8_verdict_changes.csv` has no rows), the recommended arm,
admissible arms and status of every variant are unchanged (7 RECOMMENDED,
11 NOT RECOMMENDED, 2 NOT EVALUATED), and the largest change in an admissible
arm's median year 50 basal area is 0.69 m2/ha. Seed-gate rows within 2 percent
of BALIVE rose from 90.73 to 95.66 percent; all 18,992 rows pass in both.

## Identifiers

No file here contains FIA plot coordinates, PLT_CN or CN values. Stands are
indexed by sample keys only. FIA data stayed on OSC Cardinal.

`MD5SUMS` covers every file in this directory.
