# Validation

Long-horizon validation in fvs-modern asks whether 100-year projections stay inside the range of stand structure that FIA data support and obey law-like stand relations. The current pipeline is `validation/bakuzis_stress_v2/`. Its output is the configuration catalog in `config/catalog/`, which records, per variant, which combination of parameter set and survival model passed pre-registered gates. This page describes the scorer, the catalog pipeline, the stress test and the scorecards, and how to regenerate `config/catalog/`.

!!! warning "Unverified"
    The `validation/` tree was not part of the source snapshot used to write this page. File names, the run design and the decision rules below come from the merged validation work (PRs 118 to 137) and are cross-checked against `config/catalog/README.md`, `deployment/fvs2py/fvs2py/catalog.py` and the published tables. Command-line options of the validation scripts are not documented here because they could not be read.

## Pipeline files

All paths are under `validation/bakuzis_stress_v2/`.

| File | Role |
|---|---|
| `bakuzis_stress_v2.R` | The versioned Bakuzis scorer. Scores each trajectory against the stand relations listed below. |
| `catalog/conus_catalog_runs.py` | Engine runs on Cardinal: up to 300 FIA plots per variant, 4 arms, 100 years on a 5-year step, with `NOAUTOES` for EM and IE. |
| `catalog/catalog_gates.R` | Applies the pre-registered decision rule and writes the gate and recommendation tables. |
| `catalog/stress_catalog.R` | Stress test of the decision rule: a 144-cell grid of rule settings and a 500-stand bootstrap. |
| `catalog/bakuzis_scorecards.R` | Per-variant Bakuzis scorecards. |

The four engine arms are the combinations of `default` or `calibrated` parameters with `native` or `gompit` survival. The gate table also carries rows for `conus_b1` and `conus_b2` arms, which `fvs2py.catalog` describes as projector arms from the fvs-conus pipeline, not run through this engine and listed for comparison only. Their `n_traj` values differ from the engine arms (for example 1,500 against 300 for NE).

The current catalog tables are released as version 1.2.0 of the catalog record, packaged in `validation/bakuzis_stress_v2/release_20261030/` (PR 140; version 1.1.0 stays in `release_20261029/`). Cite the record by its concept DOI, [10.5281/zenodo.23777723](https://doi.org/10.5281/zenodo.23777723). The version 1.2.0 deposit had not been published on Zenodo when this page was written, so the concept DOI resolves to the latest published version until it is.

A separate, pre-registered comparison of fvs-modern, the fvs-conus projector and unmodified USDA FVS on the same stands is described on [Model run-off (October 2026)](runoff.md).

## Stand sample and seeding

Stands come from FIA through `calibration/python/fia_stand_generator.py`, which converts FIA condition and tree records into FVS `standinit` and `treeinit` rows. Since PR 129, `build_tree_init` drops live tree records with a null `TPA_UNADJ` (records with `P2A_GRM_FLG = Y`, carried for growth, removal and mortality accounting only) and logs how many it dropped. The generator previously filled them with 6.0 trees per acre, and its comment records that this over-seeded 250 of 256 eastern seed-gate failures in the October 2026 catalog run. The catalog runner records the number of null-TPA records dropped (`n_null_tpa`) in its seed gate. `config/catalog/README.md` states that after this fix every sampled stand passes the year-0 gate and that no verdict changed from the previous catalog version.

## Scorer relations

The scorer reports, for each relation, the number of trajectories on which the relation was evaluated alongside the number flagged. The relations, as described in the combined manuscript draft (`manuscript/fvs_combined_draft.md`, section 3.3), are:

- the initial step shock;
- monotone height and monotone quadratic mean diameter;
- site ordering of height and of volume;
- the Reineke self-thinning slope, on trajectories that reach a relative density of 0.55, with an acceptance window of -2.2 to -1.2, and the pooled tau = 0.95 frontier;
- maximum basal area against the 99th and 99.9th percentiles of the FIA fitting data, and peak stand density index against its 99.9th percentile;
- Eichhorn's rule;
- Sukachev's per-tree volume ordering;
- yield against initial density;
- a basal area plateau;
- step inversions of basal area and of volume.

!!! warning "Unverified"
    The scorer labels these relations T0 to T9. The mapping from each label to a relation, and the exact test statistic and threshold for each, are defined in `bakuzis_stress_v2.R`, which was not available in the snapshot. Read the scorer source before citing a relation by its T label.

## Decision rule

The catalog ships under the pre-registered primary rule (decision of October 29, 2026, per `config/catalog/README.md`). An arm is admissible for a variant when all of the following hold:

1. Fewer than 5 percent of stands exceed the FIA basal area 99.9th percentile by year 50, and fewer than 5 percent exceed the FIA SDI 99.9th percentile by year 50, each with a Wilson 95 percent upper bound below 10 percent.
2. The Reineke flag rate is below 10 percent, evaluated on at least 10 trajectories that reach the relative density gate.
3. The median year-50 basal area is at least the FIA 50th percentile.
4. The calibrated arm is distinct from the default arm, that is, the calibrated configuration actually loaded.

The alternative rule replaces condition 2 with a Wilson upper bound on the Reineke flag rate below 25 percent on at least 30 evaluated trajectories. The recommendation table reports the alternative in `admissible_arms_alt` and sets `rule_sensitive` to `TRUE` where the two rules give different admissible sets.

Both rules can be reproduced from `config/catalog/catalog_gates_by_arm.csv` alone; the column expressions are on the [Configuration files](configuration-files.md#catalog_gates_by_armcsv) page. Variants without FIA thresholds (AK and NC at this commit) are reported as `NOT EVALUATED`.

## Regenerating config/catalog end to end

The catalog tables must be regenerated only through the pipeline; `config/catalog/README.md` says not to edit them by hand. The sequence is:

1. Build the engine libraries serially (see [Building the engine](building-the-engine.md)), and make sure `config/calibrated/<variant>.json` and `config/greg_mortality_coefficients.csv` are the versions you intend to evaluate. The gompit arms run with `FVS_GOMPIT=1` and `FVS_GOMPIT_COEF` pointing at the coefficient file.
2. Run `catalog/conus_catalog_runs.py` to project the sampled FIA stands for every variant and arm.
3. Score the trajectories with `bakuzis_stress_v2.R`.
4. Apply the gates with `catalog/catalog_gates.R`. The current tables (fifth version, October 30, 2026, PR 137; catalog release 1.2.0) were produced under `gates_v6/`; that version reran only the calibrated arms, with the overlay fix from PR 135, and carried the default arms over from the fourth version.
5. Run `catalog/stress_catalog.R` and `catalog/bakuzis_scorecards.R` for the rule stress test and the scorecards.
6. Copy `catalog_recommendation_by_variant.csv` and `catalog_gates_by_arm.csv` into `config/catalog/` and update `config/catalog/README.md` with the date, version and what changed.
7. Check the result through the runtime reader:

```bash
python3 deployment/fvs2py/fvs2py/catalog.py ne sn acd
```

The module's `__main__` block prints, per variant, the selected arm, `config_version`, the environment variables, whether the arm is admissible, the status and any notes.

!!! warning "Unverified"
    Step order 3 to 5 and the location of intermediate outputs follow the pipeline description, not a run of the scripts. Check the README in `validation/bakuzis_stress_v2/catalog/` for the exact invocations, working directories and HPC job templates.

## Consuming the catalog

```python
from fvs2py.catalog import select_config
choice = select_config("ne")
choice = select_config("sn", survival="gompit")
```

`select_config` reads both tables from `config/catalog/` (override with `catalog_dir` or the `FVS_CATALOG_DIR` environment variable). It returns a `CatalogChoice` with the arm, `config_version` (`calibrated` or `default`), the environment variables to set before the first run (`FVS_GOMPIT=1` and `FVS_GOMPIT_COEF=<config dir>/greg_mortality_coefficients.csv` for gompit arms), `admissible`, `status` and a `note`. It raises `KeyError` for a variant that is absent or not evaluated. When the catalog names a `conus_*` arm as least bad, it returns the engine arm with the lowest `share_ba999_50` (ties broken on `share_sdi_50`) and says so in the note; it never returns a projector arm as runnable.

## Other validation material

`calibration/stress/` holds a separate full-CONUS stress battery (plan in `calibration/stress/STRESS_TEST_PLAN.md`) that projects FIADB stands 100 years under default and calibrated configurations to produce a pass or fail ledger of crashes and projection errors. `calibration/python/bakuzis_100yr_comparison.py` is the earlier 36-scenario Bakuzis runner (4 species groups by 3 site classes by 3 densities). `evaluation/postfreeze_20260805/` holds the keyword-layout, input-guard and cycle-ceiling proofs referenced in `config_loader.py` and `tools/`.
