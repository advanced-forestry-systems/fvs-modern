# Configuration files

Parameter sets live as JSON under `config/`. The engine does not read them. `config/config_loader.py` reads them and either writes FVS keywords or sets species attributes on an fvs2py instance. This page documents the schema as the loader reads it, the difference between default and calibrated files, keyword generation, and the configuration catalog tables.

## Files

| File | Written by | Read by |
|---|---|---|
| `config/<variant>.json` | `modernization/extract_parameters.py` | `FvsConfigLoader` (version `default`), `calibration/R/06_posterior_to_json.R` (as template and prior source) |
| `config/calibrated/<variant>.json` | `calibration/R/06_posterior_to_json.R`; CONUS blocks by `62_conus_to_variant_json.R` and `62b_speciesfree_to_variant_json.R` | `FvsConfigLoader` (all versions except `default` and `custom`), `config/validate_calibrated.py` |
| `config/calibrated/<variant>_draws.json` | `06_posterior_to_json.R` | `config/uncertainty.py` (`UncertaintyEngine`) |
| `config/catalog/catalog_recommendation_by_variant.csv`, `config/catalog/catalog_gates_by_arm.csv` | the validation catalog pipeline | `fvs2py.catalog` |
| `config/greg_mortality_coefficients.csv` | Output name of `calibration/R/35_fit_greg_mortality_cch.R` | the engine, through `FVS_GOMPIT_COEF` |

`config/` currently holds 26 default files: `acd ak bc bm ca ci cn cr cs ec em ie kt ls nc ne oc on op pn sn so tt ut wc ws`.

## Schema as the loader reads it

A default file has five top-level keys:

```json
{
  "variant": "ne",
  "variant_name": "Northeastern US",
  "maxsp": 108,
  "files_parsed": ["..."],
  "categories": {
    "bark_ratio": {"BKRAT": [...]},
    "growth": {"B1": [...], "B2": [...]},
    "crown": {"BCR1": [...], "BCR2": [...], "BCR3": [...], "BCR4": [...]},
    "site_index": {"SDICON": [...]},
    "species_definitions": {"JTYPE": [...], "JSP": [...], "FIAJSP": [...], "PLNJSP": [...], "NSP": [...]},
    "other": {"...": [...]}
  }
}
```

`categories` maps a category name to a mapping of Fortran array names to values. Species-indexed arrays have `maxsp` entries in FVS species order. The extractor assigns arrays to categories by name pattern; its full category set is `bark_ratio`, `growth`, `crown`, `height_diameter`, `height_growth`, `site_index`, `volume`, `mortality`, `regeneration`, `species_definitions` and `other`, and empty categories are dropped. NE has six.

The loader uses these parts of the schema:

- `maxsp` (via `FvsConfigLoader.maxsp`).
- `categories.<category>.<name>` (via `get_param` and `get_species_params`, which raise `KeyError` listing what is available).
- `categories.species_definitions.FIAJSP`, the FIA species code per FVS slot, used to map SPCD-keyed calibration results onto FVS slots. `config/validate_calibrated.py` fails a calibrated config that lacks it.
- The stand density parameter, found by `_find_sdi_param` as the first of `SDICON`, `R5SDI`, `R4SDI`, `FMSDI`, `SDIDEF` present in `categories.site_index` or `categories.other`.
- The maximum basal area parameter, found by `_find_bamax_param` as the first of `BAMAXA`, `BAMAX1`, `BAMAX` in the same two categories.

## Default and calibrated

`FvsConfigLoader(variant, version=...)` accepts nine versions:

```python
VALID_VERSIONS = ("default", "calibrated", "conus", "hybrid", "custom",
                  "conus_sf", "conus_hybrid", "conus_organon", "conus_greg")
```

`default` reads `config/<variant>.json`. `custom` reads the file passed as `custom_config` and raises if none is given. Every other version reads `config/calibrated/<variant>.json`; the version selects which top-level block the access methods use (`categories`, `categories_conus` or `categories_conus_sf`). `config_dir` overrides the directory.

A calibrated file is the default file loaded as a template with fitted values written over it, so `categories` keeps the default schema. `06_posterior_to_json.R` updates arrays in categories including `height_diameter`, `mortality`, `crown`, `crown_change`, `self_thinning`, `species_definitions` and `other`, writes the stand density array it finds (`SDICON`, `R5SDI`, `R4SDI` or `FMSDI`) and the BAMAX array (`BAMAXA` or `BAMAX1`) when a stand density fit exists, and adds these top-level blocks:

| Block | Content |
|---|---|
| `calibration` | `date`, `method`, `fia_data_source`, `r_version`, `prior_type`, `components_updated`, `components_skipped`, `notes`. The loader also accepts the older key `calibration_metadata`. |
| `calibration_multipliers` | Five per-species arrays of length `maxsp` (`htdbh_multiplier`, `mort_multiplier`, `cr_multiplier`, `dds_multiplier`, `htg_multiplier`) and a `provenance` record, computed by `calibration/R/multipliers.R`. |
| `categories_conus` | CONUS Phase 4 per-species (Leg A) blocks per component plus `metadata`. Schema in `CONUS_BLOCK_SCHEMA.md`. |
| `categories_conus_sf` | Species-free (Leg B) blocks: fixed effects, trait coefficients and random-effect tables. |

An optional top-level flag `_emit_sdimax` (default true) turns SDIMAX keyword emission off for a variant.

`calibration_multipliers` definitions, from `multipliers.R`:

| Array | Component | Factor |
|---|---|---|
| `htdbh_multiplier` | Height-diameter | asymptote ratio `(a_pooled + r_SPCD__a[s]) / a_pooled` |
| `mort_multiplier` | Mortality | `(1 - plogis(b0 + r_SPCD[s])) / (1 - plogis(b0))` |
| `cr_multiplier` | Crown ratio | `exp(r_SPCD_grouped[s])` |
| `dds_multiplier` | Diameter growth | `exp(b0[s] - mean b0)` |
| `htg_multiplier` | Height increment | `exp(b0[s] - mean b0)` |

Height-diameter, mortality and crown ratio are SPCD-keyed and map to FVS slots exactly through `FIAJSP`. Diameter growth and height increment use a dense Stan species index and map exactly only when the fit wrote a `<component>_species_index.csv` crosswalk. Each component is gated by `calibration/data/equation_availability_full.csv`; where a component is not adopted for a variant its array is all 1.0.

!!! warning "Known issue: unity DG and HTG arrays"
    `KNOWN_ISSUES.md` records that `dds_multiplier` and `htg_multiplier` are exactly 1.0 for all 108 slots in `config/calibrated/ne.json` and `config/calibrated_sdifix/ne.json`, although the fits contain species-level terms, because the species index crosswalk was never written. A unity array is accepted by `_array_or_none` and therefore suppresses the legacy fallback computation.

`config/validate_calibrated.py`, run in CI, checks every `config/calibrated/*.json` except `*_draws.json`: `maxsp` present, `FIAJSP` present, all five arrays present with length `maxsp`, every value finite and positive, and, where the availability table marks an SPCD-keyed component as adopted, that its array is not all 1.0 (a hard error; for DG and HI the same condition is a warning). It skips the availability checks with a warning when `equation_availability_full.csv` is absent.

## Keyword generation

`FvsConfigLoader.generate_keywords(include_comments=True)` returns a keyword block built from the active config. `generate_calibration_keyfile(variant, output_path=None)` wraps it with version `calibrated` and optionally writes it to a file you can pull in with `ADDFILE`. `config/export_keywords.py` is the command-line front end:

```bash
python export_keywords.py --variant ne --output ne_calibrated.key
python export_keywords.py --variant ne --keyfile input.key --output output.key
python export_keywords.py --variant ne --keyfile input.key \
    --n-draws 50 --seed 42 --output-dir ensemble/
python export_keywords.py --all-variants --output-dir calibrated_keywords/
python export_keywords.py --list
```

With `--keyfile`, the block is inserted before the first line that starts with `PROCESS`, or appended if there is none.

The block is assembled in this order. All records use a 10-column keyword field followed by 10-column numeric fields.

1. **SDIMAX**, one record per species with a positive value: `f"{'SDIMAX':<10}{i + 1:10d}{val:10.1f}"`. Skipped when `_emit_sdimax` is false.
2. **BAMAX**, a stand-level keyword with no species field: `f"{'BAMAX':<10}{value:10.1f}"`. If the config holds a per-species BAMAX vector with more than one distinct positive value, the formatter raises `ValueError`; `generate_keywords` catches it, issues a warning (`BAMAX skipped for variant ...`) and continues with the rest of the block. Before this change the exception aborted the whole calibrated keyword block for ACD and IE.
3. **MORTMULT**, from `calibration_multipliers.mort_multiplier` when present, else from the legacy computation (odds ratio of mortality intercepts, clipped to 0.1 to 10). Only species with a multiplier more than 0.01 from 1 are written. MORTMULT is a scheduled activity, so field 1 is the date (left blank) and the record is `keyword, blank date, species, multiplier, 0.0, 999.0`.
4. **BAIMULT**, from `sqrt(clip(dds_multiplier, 0.01, 100))` when present, else from the legacy computation (square root of `exp(cal_B0 - def_B0)`, clipped). Same 0.01 threshold and scheduled layout (blank date, species, multiplier).
5. **HTGMULT**, from `htg_multiplier` when present, else from the legacy ratio of `HGLD` or `HG_B0`. Same layout as BAIMULT.

When a non-empty BAIMULT or HTGMULT vector yields no records because every value is within 0.01 of 1, the loader writes a note to stderr (and into the block when comments are on) saying the family is a no-op for the variant.

!!! note "Column layout"
    The docstrings in `config_loader.py` record that earlier writers padded keywords to 16 columns, which shifted values into the wrong fields and caused FVS to ignore or misread them without a visible error. The current writers use 10 columns. `KNOWN_ISSUES.md` still lists BAMAX, BAIMULT and HTGMULT as using the 16-column form; that entry is stale relative to the code.

`apply_to_fvs(fvs_instance)` is the fvs2py path. It sets the species attributes `spsdi`, `baimult`, `mortmult` and `htgmult` through `set_species_attr`, padding or trimming each multiplier array to the instance's `maxspecies`. For `spsdi` it reads the current FVS values and overwrites only species with a finite, positive calibrated value; since PR 135 (October 30, 2026) species with a missing value keep their FVS value instead of being written as 0. It uses the legacy multiplier computations, not `calibration_multipliers`. The formatter methods (`_format_sdimax_keywords`, `_format_bamax_keywords`, `_format_mortmult_keywords`, `_format_baimult_keywords`, `_format_htgmult_keywords`) are static methods so that `UncertaintyEngine.generate_keywords_for_draw` writes draw keywords with the same layout.

`FvsConfigLoader.compare(variant)` and `compare_configs(variant)` report default against calibrated values per numeric parameter, skipping arrays of different lengths.

## Catalog tables

`config/catalog/` holds the configuration catalog read by `fvs2py.catalog.select_config(variant)`. The tables are generated by the validation pipeline and copied here; the catalog README says not to edit them by hand. See [Validation](validation.md) for regeneration.

Arms are named `<parameter set>_<survival>`. Parameter sets are `default`, `calibrated`, `conus_b1` and `conus_b2`; survival is `native` (the variant's own mortality) or `gompit` (GOMPIT survival with `FVS_GOMPIT=1` and `FVS_GOMPIT_COEF` set). The `conus_*` arms are projector arms from the separate fvs-conus pipeline, not run through this engine and listed for comparison only. In the gates table the bare `conus_b1` and `conus_b2` labels carry no survival suffix.

### catalog_recommendation_by_variant.csv

One row per variant (`group`, upper case). 20 rows at commit fb4cf10.

| Column | Content |
|---|---|
| `group` | Variant code. |
| `recommended` | The recommended arm, or the least-bad arm when the variant is not recommended. Empty when not evaluated. |
| `status` | `RECOMMENDED`, `NOT RECOMMENDED for long horizons (least-bad arm named)` or `NOT EVALUATED (no FIA thresholds for this variant)`. |
| `admissible_arms` | Arms that pass the primary rule, separated by `; `. |
| `nearterm_evidence` | Per-variant near-term evidence; `none per variant` except NE. When the caller constrains the survival model or parameter set, `select_config` prefers a gompit arm if this field mentions the NE cohort. |
| `rule_sensitive` | `TRUE` when the admissible set differs under the alternative rule. |
| `admissible_arms_alt` | Arms that pass the alternative rule. |

At commit fb4cf10, 7 variants are `RECOMMENDED` (ACD, BM, CS, LS, NE, SN, WS), 11 are `NOT RECOMMENDED for long horizons` and 2 are `NOT EVALUATED` (AK, NC).

### catalog_gates_by_arm.csv

One row per arm and variant, 146 rows at commit fb4cf10, with these columns:

```text
arm,group,n_traj,fia_n_plots,ba_p999,ba_p99,ba_p50,sdi_p999,x_ba999_50,x_ba999_100,
x_ba99_50,x_ba99_100,x_sdi_50,x_sdi_100,med_peak_ba100,med_peak_sdi100,med_ba0,med_ba50,
med_stem_ret50,share_peak_below_p50,share_ba999_50,share_ba999_100,share_ba99_50,
share_ba99_100,share_sdi_50,share_sdi_100,wil_ba999_50,wil_sdi_50,reineke_eval,
reineke_flag,reineke_na,reineke_slope_med,reineke_flag_rate,reineke_eval_share,seed_n,
seed_ok_share,run_fail_share,calibrated_distinct,reineke_unknown,admissible,
reineke_flag_wil,admissible_alt,why_not,threshold_n,simplicity
```

The columns that carry the decision rules can be checked against the table itself. On every row where all inputs are present, `admissible` equals the conjunction

```text
share_ba999_50 < 0.05 and share_sdi_50 < 0.05
and wil_ba999_50 < 0.10 and wil_sdi_50 < 0.10
and reineke_flag_rate < 0.10 and reineke_eval >= 10
and med_ba50 >= ba_p50 and calibrated_distinct == TRUE
```

and `admissible_alt` equals the same conjunction with the two Reineke terms replaced by `reineke_flag_wil < 0.25 and reineke_eval >= 30`. This matches the pre-registered primary rule and the alternative rule described on the [Validation](validation.md) page. The other identifying columns are `arm`, `group`, `n_traj` (trajectories scored), `fia_n_plots` and `threshold_n`. `fvs2py.catalog` reads `why_not` (reason text), `reineke_unknown` (Reineke not testable), `calibrated_distinct` (`FALSE` means the calibrated config did not load and the arm equals default) and, when it must replace a projector arm with an engine arm, `share_ba999_50` and `share_sdi_50`.

!!! warning "Unverified"
    Definitions of the remaining diagnostic columns (`x_*`, `med_peak_*`, `med_ba0`, `med_stem_ret50`, `share_peak_below_p50`, the `*_100` and `ba99` shares, `reineke_flag`, `reineke_na`, `reineke_slope_med`, `reineke_eval_share`, `seed_n`, `seed_ok_share`, `run_fail_share`, `simplicity`) are written in `validation/bakuzis_stress_v2/catalog/`, which was not part of the source snapshot used for this page. Read them from that README rather than inferring them from column names.
