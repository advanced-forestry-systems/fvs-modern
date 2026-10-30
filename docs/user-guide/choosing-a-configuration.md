# Choosing a configuration

Every projection runs one parameter set with one survival model. fvs-modern offers two parameter sets and two survival models for each variant it has calibrated, and a configuration catalog that says which combination has passed long-horizon realism gates for that variant. This page explains the four combinations, how to read the catalog, and how to select a configuration with `fvs2py.catalog.select_config`.

## Parameter sets

**Default** (`config_version=None` or `"default"`) runs the coefficients compiled into the variant library. No overlay is applied.

**Calibrated MAP** (`config_version="calibrated"`) overlays the Bayesian calibration stored in `config/calibrated/<variant>.json`, which `CALIBRATION.md` describes as posterior medians from hierarchical models fitted to national FIA remeasurements. The catalog code calls this arm "fvs-modern Bayesian calibrated MAP". There are two ways the overlay reaches the engine, and they do not apply the same quantities:

| Path | Used by | What is changed |
|:---|:---|:---|
| In process (`FvsConfigLoader.apply_to_fvs`) | `fvs2py.FVS` with `config_version="calibrated"` or `"custom"` | Species attributes `spsdi`, and `baimult`, `mortmult`, `htgmult` where they can be derived from differences between calibrated and default coefficient arrays in `categories` |
| Keyword (`FvsConfigLoader.generate_keywords`) | microfvs, standalone executables, `config.config_loader.generate_calibration_keyfile` | `SDIMAX` records per species, a stand-level `BAMAX` record, and `MORTMULT`, `BAIMULT`, `HTGMULT` records from the `calibration_multipliers` block |

For NE the in-process path can set only `spsdi`, because the default NE configuration has no growth `B0` array and no mortality category to difference against; see [Running projections](running-projections.md#what-the-in-process-overlay-changes). The keyword path skips any multiplier within 0.01 of 1 and writes a message to standard error when a whole family is unity (`... This family is a NO-OP for this variant; the arm runs on the uncalibrated growth model.`).

!!! warning "Diameter growth is uncalibrated in most variants"
    `KNOWN_ISSUES.md` records that only 7 of 25 variants (`acd`, `ca`, `cs`, `kt`, `ls`, `nc`, `on`) carry non-unity per-species diameter growth multipliers. The other 18, including `ne`, run default diameter growth under the calibrated set because of a species-index mapping failure in the serializer. In `config/calibrated_sdifix/ne.json`, the NE configuration available in the snapshot used for this page, `dds_multiplier` and `htg_multiplier` are 1.0 for all 108 species while `mort_multiplier` ranges from 0.24 to 6.69. The same file records that the Acadian variant does not respond to FVS keyword multipliers (`CALIBRATION.md`, "Variant caveat: FVS-ACD"). Treat "calibrated" as a calibration of height-diameter, mortality, crown ratio and SDI maximum unless your variant is one of the seven.

!!! warning "Unverified"
    The calibrated configuration files `config/calibrated/<variant>.json` are tracked in the repository (25 files at commit fb4cf10) but were too large for the text snapshot used to write this page, so their contents were not checked directly. The posterior draw files `config/calibrated/<variant>_draws.json` are not tracked in the repository at that commit.

**Custom** (`config_version="custom", custom_config="my_ne.json"`) applies a JSON file you supply in the same schema as `config/<variant>.json`, through the same two paths. It is not part of the catalog.

## Survival models

**Native** survival is the variant's own mortality routine.

**GOMPIT** replaces native mortality, tree by tree and cycle by cycle inside the growth loop, with a per-species annual hazard, H = exp(b0 + b1(CR + 0.01)^b2 + b3·CCH^b4), where CR is crown ratio and CCH is crown closure at the tree tip, and period survival exp(-H × cycle length) (`src-converted/base/gompmort.f90`). Species without a fitted row in the coefficient file keep native background mortality. GOMPIT is switched on by environment variables, not by `config_version`:

```bash
export FVS_GOMPIT=1
export FVS_GOMPIT_COEF="$PWD/config/greg_mortality_coefficients.csv"
```

It can also be switched on with the `GOMPMORT` keyword, but the coefficient path still has to come from `FVS_GOMPIT_COEF`. When it is active the main output listing contains a line such as `GOMPMORT enabled: 92 fitted species read, 63 matched to this variant.` (the NE figures from the documentation build). If you see `GOMPMORT: ... NOT enabled.` instead, the coefficient path is blank or unreadable and the run used native mortality. The engine reads these variables once per process, so set them before the first `run()`.

`KNOWN_ISSUES.md` records two related limitations for Northeast long runs on the frozen v1 model: GOMPIT is density aware through crown closure but has no productivity term, and the Greg arms, in which the GOMPIT override drives mortality, lose 79 to 84 percent of stems over a long horizon with almost no difference across site classes.

## The configuration catalog

The catalog is two CSV files in `config/catalog/`, regenerated on October 30, 2026 (the fifth version, merged in PR 137) by the benchmark pipeline under `validation/bakuzis_stress_v2/catalog/` and not to be edited by hand (`config/catalog/README.md`):

- `catalog_recommendation_by_variant.csv`: one row per variant with `group`, `recommended`, `status`, `admissible_arms`, `nearterm_evidence`, `rule_sensitive` and `admissible_arms_alt`;
- `catalog_gates_by_arm.csv`: one row per variant and arm with the gate statistics (`n_traj`, `fia_n_plots`, FIA percentile thresholds, exceedance shares, Reineke diagnostics, `admissible`, `admissible_alt`, `why_not` and others).

Arms are named `<parameter set>_<survival>`: `default_native`, `default_gompit`, `calibrated_native`, `calibrated_gompit`. The files also carry `conus_b1`, `conus_b2`, `conus_b1_gompit` and `conus_b2_gompit` rows for the CONUS trait-mediated projector, which is not run through this engine and is listed for comparison only.

The manuscript draft in the repository describes the gates. Up to 300 FIA plots per variant (4,748 stands across 18 variants) were projected 100 years under each arm. An arm is admissible when fewer than 5 percent of its stands exceed the FIA 99.9th percentiles of basal area and stand density index by year 50, it self-thins inside the Reineke window on at least 90 percent of density-gated trajectories, and its median year 50 basal area is at or above the FIA median.

### Reading status

| `status` | Meaning | What `select_config` does |
|:---|:---|:---|
| `RECOMMENDED` | At least one engine arm passed every gate. `recommended` names the arm chosen from the admissible set. | Returns that arm with `admissible=True`. |
| `NOT RECOMMENDED for long horizons (least-bad arm named)` | No engine arm passed. `recommended` names the least-bad arm and `admissible_arms` is empty. | Returns the named arm with `admissible=False` and the reasons in `note`. If the named arm is a CONUS projector arm, returns the engine arm with the lowest year 50 basal area exceedance share, then SDI exceedance share, and says so in `note`. |
| `NOT EVALUATED (no FIA thresholds for this variant)` | The variant was not gated. `recommended` is blank. | Raises `KeyError`. |

A variant missing from the file altogether also raises `KeyError` (`variant KT is not in the catalog (no fitting plots assigned or not run)`).

"NOT RECOMMENDED for long horizons" is a statement about 100 year behaviour on the gates above. The `why_not` column gives the reason per arm; in the current files it is `Reineke flags` for every engine arm of CA, CI, CR, EC, IE, SO and UT and for three of the four TT engine arms (the fourth, `default_gompit`, reads `Reineke not testable (7 evaluated)`), and it adds `BA>p99.9 by 50 yr` for some arms of EM, LS, PN and WC.

### Stability labels

Two fields describe how stable a verdict is.

`rule_sensitive` is TRUE when the admissible set under the alternative Reineke gate (`admissible_arms_alt`) differs from the set under the primary rule. The catalog ships under the primary rule (`config/catalog/README.md`); the alternative is reported so that you can see which verdicts depend on that choice. For example, WS is RECOMMENDED with `calibrated_native` only, but `default_native` is also admissible under the alternative rule, and CI has no admissible arm under the primary rule but `default_native` under the alternative.

`nearterm_evidence` records whether any per-variant near-term evidence informed the choice. It is `Paper 2 NE cohort (stocking)` for NE and `none per variant` for every other variant. When a caller constrains the survival model or parameter set, `select_config` prefers the GOMPIT arm among admissible candidates only when this field mentions the NE cohort.

Per arm, `reineke_unknown` is TRUE when too few trajectories reached the relative density gate for the self-thinning test (`select_config` then adds `Reineke self-thinning not testable (N trajectories reached the relative-density gate).` to `note`), and `calibrated_distinct` is FALSE when the calibrated configuration failed to load for that variant so the calibrated arm equals the default. In the current files `reineke_unknown` is TRUE only for TT `default_gompit` (7 trajectories evaluated) and `calibrated_distinct` is TRUE everywhere.

## Current catalog by variant

The table reproduces `config/catalog/catalog_recommendation_by_variant.csv` as of commit fb4cf10 (catalog fifth version). The last column is what `select_config(variant)` returned for each variant when run against those files.

!!! note "Changes in the fifth version (October 30, 2026)"
    The calibrated arms of all 18 variants were rerun after the in-process overlay fix in PR 135 (see [Running projections](running-projections.md#what-the-in-process-overlay-changes)); the default arms were carried over unchanged. According to `config/catalog/README.md`, the least-bad arm named for EM is now `conus_b2_gompit` (was `calibrated_gompit`) and for WC `default_gompit` (was `calibrated_gompit`), both still NOT RECOMMENDED for long horizons; WC `calibrated_gompit` becomes admissible under the alternative rule; and no status and no primary-rule admissibility changed.

| Variant | `status` | `recommended` | `admissible_arms` | `rule_sensitive` | `admissible_arms_alt` | `select_config` returns |
|:---|:---|:---|:---|:---|:---|:---|
| ACD | RECOMMENDED | `default_native` | `default_native` `default_gompit` `calibrated_native` `calibrated_gompit` | FALSE | `default_native` `default_gompit` `calibrated_native` `calibrated_gompit` | `default_native` (`config_version="default"`, `admissible=True`) |
| BM | RECOMMENDED | `default_native` | `default_native` `default_gompit` `calibrated_native` `calibrated_gompit` | TRUE | `default_native` `default_gompit` `calibrated_native` `calibrated_gompit` `conus_b1_gompit` `conus_b2_gompit` | `default_native` (`config_version="default"`, `admissible=True`) |
| CA | NOT RECOMMENDED for long horizons (least-bad arm named) | `calibrated_gompit` | none | TRUE | `conus_b1_gompit` | `calibrated_gompit` (`config_version="calibrated"`, `admissible=False`) |
| CI | NOT RECOMMENDED for long horizons (least-bad arm named) | `calibrated_native` | none | TRUE | `default_native` | `calibrated_native` (`config_version="calibrated"`, `admissible=False`) |
| CR | NOT RECOMMENDED for long horizons (least-bad arm named) | `calibrated_gompit` | none | TRUE | `calibrated_native` `conus_b1_gompit` | `calibrated_gompit` (`config_version="calibrated"`, `admissible=False`) |
| CS | RECOMMENDED | `default_native` | `default_native` `default_gompit` `calibrated_native` `calibrated_gompit` | FALSE | `default_native` `default_gompit` `calibrated_native` `calibrated_gompit` | `default_native` (`config_version="default"`, `admissible=True`) |
| EC | NOT RECOMMENDED for long horizons (least-bad arm named) | `calibrated_gompit` | none | TRUE | `conus_b1_gompit` | `calibrated_gompit` (`config_version="calibrated"`, `admissible=False`) |
| EM | NOT RECOMMENDED for long horizons (least-bad arm named) | `conus_b2_gompit` | none | FALSE | none | `default_gompit` (`config_version="default"`, `admissible=False`) |
| IE | NOT RECOMMENDED for long horizons (least-bad arm named) | `calibrated_gompit` | none | FALSE | none | `calibrated_gompit` (`config_version="calibrated"`, `admissible=False`) |
| LS | RECOMMENDED | `calibrated_native` | `calibrated_native` `calibrated_gompit` | FALSE | `calibrated_native` `calibrated_gompit` | `calibrated_native` (`config_version="calibrated"`, `admissible=True`) |
| NE | RECOMMENDED | `default_gompit` | `default_native` `default_gompit` `calibrated_native` `calibrated_gompit` | FALSE | `default_native` `default_gompit` `calibrated_native` `calibrated_gompit` | `default_gompit` (`config_version="default"`, `admissible=True`) |
| PN | NOT RECOMMENDED for long horizons (least-bad arm named) | `calibrated_gompit` | none | FALSE | none | `calibrated_gompit` (`config_version="calibrated"`, `admissible=False`) |
| SN | RECOMMENDED | `default_native` | `default_native` `default_gompit` `calibrated_native` `calibrated_gompit` | FALSE | `default_native` `default_gompit` `calibrated_native` `calibrated_gompit` | `default_native` (`config_version="default"`, `admissible=True`) |
| SO | NOT RECOMMENDED for long horizons (least-bad arm named) | `calibrated_gompit` | none | TRUE | `conus_b1_gompit` | `calibrated_gompit` (`config_version="calibrated"`, `admissible=False`) |
| TT | NOT RECOMMENDED for long horizons (least-bad arm named) | `conus_b2_gompit` | none | FALSE | none | `default_native` (`config_version="default"`, `admissible=False`) |
| UT | NOT RECOMMENDED for long horizons (least-bad arm named) | `calibrated_gompit` | none | FALSE | none | `calibrated_gompit` (`config_version="calibrated"`, `admissible=False`) |
| WC | NOT RECOMMENDED for long horizons (least-bad arm named) | `default_gompit` | none | TRUE | `default_gompit` `calibrated_gompit` | `default_gompit` (`config_version="default"`, `admissible=False`) |
| WS | RECOMMENDED | `calibrated_native` | `calibrated_native` | TRUE | `default_native` `calibrated_native` | `calibrated_native` (`config_version="calibrated"`, `admissible=True`) |
| AK | NOT EVALUATED (no FIA thresholds for this variant) | (blank) | none | FALSE | none | raises `KeyError` |
| NC | NOT EVALUATED (no FIA thresholds for this variant) | (blank) | none | FALSE | none | raises `KeyError` |

For EM and TT the catalog's least-bad arm is the CONUS projector arm `conus_b2_gompit`, which cannot run in this engine, so `select_config` falls back to the engine arm with the lowest year 50 basal area exceedance share. For EM that is `default_gompit` (an exceedance share of 0.077 across 300 stands, against 0.097 for `calibrated_gompit` and 0.223 for both native arms). All four TT engine arms have the same year 50 basal area and SDI exceedance shares (one stand in 142), so the tie-break falls to file order and `select_config("tt")` returns `default_native`, the TT engine arm with the highest Reineke flag rate (0.84). Treat TT projections beyond the near term with particular caution.

## Using select_config

```python
from fvs2py.catalog import select_config, catalog_table

select_config(variant, survival=None, parameter_set=None,
              catalog_dir=None, config_dir=None)
```

| Argument | Meaning |
|:---|:---|
| `variant` | Variant code, any case. |
| `survival` | `"native"` or `"gompit"` to force the survival model. The result is the first admissible arm, in the order of `admissible_arms`, that matches; when none matches, the arm is built from the forced component and the other component of the recommended arm, with `admissible=False`. |
| `parameter_set` | `"default"` or `"calibrated"` to force the parameter set, by the same rule. When both constraints are given, both must match. |
| `catalog_dir` | Directory holding the two CSV files. Default: `FVS_CATALOG_DIR`, then `config/catalog` in the checkout. |
| `config_dir` | Directory used to build the `FVS_GOMPIT_COEF` path. Default: `FVS_CONFIG_DIR`, then `config/` in the checkout. |

It returns a `CatalogChoice` dataclass:

| Field | Content |
|:---|:---|
| `variant` | Upper-case variant code |
| `arm` | Arm name, for example `default_gompit` |
| `parameter_set` | `default` or `calibrated` |
| `survival` | `native` or `gompit` |
| `config_version` | `"calibrated"` or `"default"`, ready to pass to `FVS` |
| `env` | `{}` for native survival; for GOMPIT, `{"FVS_GOMPIT": "1", "FVS_GOMPIT_COEF": "<config_dir>/greg_mortality_coefficients.csv"}` |
| `admissible` | Whether the returned arm passed every gate for this variant |
| `status` | The variant's `status` string from the catalog |
| `admissible_arms` | List of admissible arms under the primary rule |
| `nearterm_evidence` | As in the catalog |
| `evidence` | The arm's full row from `catalog_gates_by_arm.csv` as a dictionary of strings |
| `note` | Empty for an admissible arm without caveats; otherwise the reasons, joined into one string |

`select_config` never substitutes an admissible arm silently for one the catalog did not admit: it returns the requested or least-bad arm with `admissible=False` and leaves the decision to you. It raises `ValueError` if the resulting arm is a CONUS projector arm, which can happen when you force a survival or parameter set on a variant whose only candidates are projector arms.

`catalog_table(catalog_dir=None)` returns every row of `catalog_gates_by_arm.csv` as a list of dictionaries. Running the module prints one line per variant code given on the command line:

```bash
PYTHONPATH=deployment/fvs2py python3 -m fvs2py.catalog ne ls ws
```

### Examples

```python
import os
from fvs2py import FVS
from fvs2py.catalog import select_config

choice = select_config("ne")
print(choice.arm, choice.config_version, choice.admissible)
# default_gompit default True

if not choice.admissible:
    print("Not admissible:", choice.note)

os.environ.update(choice.env)          # before the first run() in this process
fvs = FVS("lib/FVSne.so", config_version=choice.config_version)
```

Forcing a component:

```python
select_config("ne", survival="native").arm        # 'default_native'
select_config("ne", parameter_set="calibrated").arm  # 'calibrated_gompit'
select_config("sn", survival="gompit").arm        # 'default_gompit'

c = select_config("ca", survival="native")
c.arm, c.admissible                               # ('calibrated_native', False)
```

For CA the note reads `calibrated_native is not an admissible arm for CA: Reineke flags; Verdict is rule-sensitive (admissible set differs under the alternative Reineke gate). NOT RECOMMENDED for long horizons (least-bad arm named)`.

!!! note "GOMPIT is process-wide"
    `os.environ.update(choice.env)` affects every library loaded later in the same process, and the engine caches the GOMPIT decision after the first run. To compare native and GOMPIT arms, run them in separate processes.

## Choosing in practice

For a variant marked RECOMMENDED, use `select_config(variant)` and record the returned `arm` with your results. For a variant marked NOT RECOMMENDED, the returned arm is the least-bad engine configuration, not a validated one; the catalog says nothing about short-horizon accuracy, so use it for 100 year projections only with the `note` reported alongside. AK and NC are listed as NOT EVALUATED, and BC, KT, OC, ON and OP are absent from the file; for all seven `select_config` raises `KeyError`, and the choice between default and calibrated parameters rests on the calibration caveats above. See [Variants](variants.md) for the coverage of each code.
