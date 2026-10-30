# `categories_conus` JSON schema

Documented schema for the new per-variant config block carrying CONUS Phase 4 calibration parameters. Sits alongside the existing `categories` block in `config/calibrated/{variant}.json`. Existing readers continue to use `categories`; new readers can opt into `categories_conus` via the runtime switch in `config_loader.py` (see `parameter_source_switch` section below).

## Top-level shape

```json
{
  "variant": "ne",
  "variant_name": "Northeast",
  "maxsp": 39,
  "files_parsed": {...},
  "categories": {... existing per-variant fits, untouched ...},
  "calibration": {... existing metadata, untouched ...},

  "categories_conus": {
    "metadata": {
      "integration_date": "2026-05-09",
      "pipeline_version": "fvs-conus phase 4",
      "selected_lambda": 10,
      "components_present": ["diameter_growth", "height_growth"]
    },

    "diameter_growth":   { ... see component block below ... },
    "height_growth":     { ... },
    "height_diameter":   { ... },
    "height_crown_base": { ... },
    "mortality":         { ... },
    "crown_recession":   { ... }
  }
}
```

`categories_conus` is sparse: only components that have been integrated land here. The runtime fallback rule (in `config_loader.py`) reads from `categories` for any component missing in `categories_conus`.

## Per-component block

Each component shares the same shape (the inner content varies by model parameterization but the keys are stable):

```json
"diameter_growth": {
  "model": "dg_kuehne_cspi_traits1",
  "modifier_lambda": 10,
  "fixed_effects": {
    "param":  ["mu_b0", "b1", "b2", "b3", "b4", "b5", "b6", "b7", "b8", "b9", "b10", "K1", "K2", "sigma_sp", "sigma_eco", "sigma"],
    "mean":   [-1.36,    0.228, -3.7e-5, 0.869, 0.077, ...],
    "sd":     [0.041,    0.012,  4.2e-6, 0.045, 0.022, ...],
    "q025":   [...],
    "median": [...],
    "q975":   [...],
    "rhat":   [1.001, 1.000, 1.000, 1.001, ...]
  },
  "species_intercepts": {
    "SPCD": [12, 16, 18, 36, 95, 96, 97, 105, 261, 263, 318, 371, 375, 391, 400, 402, 462, 471, 491, 531, 541, 543, 552, 555, 591, 601, 602, 611, 612, 621, 651, 654, 661, 693, 701, 711, 720, 741, 762, 771, 802, 812, 823, 826, 832, 833, 837, 901, 922, 951, 972, 975],
    "mean": [-0.024, 0.073, -0.012, ...],
    "sd":   [0.018, 0.022, 0.015, ...]
  },
  "species_missing": [],

  "ecodiv_intercepts": {
    "ecodiv": ["210", "211", "212", "M211", "M212", "221", "222", "M221", ...],
    "mean":   [-0.015, -0.038,  0.012, -0.005, ...],
    "sd":     [0.022,  0.018,   0.025,  0.020, ...]
  },
  "ecodiv_weights": {
    "M211": 0.42,
    "M212": 0.31,
    "211":  0.18,
    "212":  0.09
  },

  "modifier": {
    "coef": ["beta_csi_dia", "beta_traits_wd", "beta_csi_x_wd", ...],
    "mean": [0.043, -0.012, 0.008, ...],
    "sd":   [0.011, 0.007,  0.004, ...]
  },

  "draws_csv": "calibration/output/conus/dg_kuehne_cspi_traits1_fixed_draws.csv",
  "notes": "Posterior summary; full draws referenced via draws_csv."
}
```

## Field reference

| Field | Type | Purpose |
|---|---|---|
| `model` | string | identifier matching the file stem in `calibration/output/conus/` |
| `modifier_lambda` | int | the regularization choice from the 5/10/20 sweep that was picked for this integration |
| `fixed_effects.param` | string[] | parameter names from the Stan model (e.g. b0..b10, K1, K2, sigma_*) |
| `fixed_effects.mean` etc. | float[] | posterior summary statistics; one entry per param |
| `species_intercepts.SPCD` | int[] | filtered to the variant's `categories.species_definitions.FIAJSP` set |
| `species_missing` | int[] | SPCDs in the variant's species set that have no CONUS posterior estimate (need a fallback rule at runtime) |
| `ecodiv_intercepts.ecodiv` | string[] | full set of ecodivisions from the CONUS fit, NOT pre-filtered |
| `ecodiv_weights` | object | per-variant ecodivision coverage weights (sum to 1.0); runtime computes a coverage-weighted RE |
| `modifier.coef` | string[] | modifier model coefficient names (the climate / trait interaction terms) |
| `draws_csv` | string | relative path to the full posterior draws CSV; runtime samplers read this on demand |

## Component-specific param names

The `fixed_effects.param` list varies per component. Reference values (subject to refinement when each component's final model is selected):

- **diameter_growth (dg_kuehne)**: `mu_b0, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10, K1, K2, sigma_sp, sigma_eco, sigma`
- **height_growth (hg_organon_fixedK)**: `mu_b0, b1, b2, b3, b4, b5, K1, K2, sigma_sp, sigma_eco, sigma`
- **height_diameter (htdbh_wykoff_lognormal)**: `mu_b1, b2, b3, sigma_sp, sigma_eco, sigma`
- **height_crown_base (hcb_organon)**: `mu_b0, b1, b2, b3, b4, sigma_sp, sigma_eco, sigma`
- **mortality (mort_logit_simple)**: `mu_alpha, beta_dbh, beta_ba, beta_csi, sigma_sp, sigma_eco`
- **crown_recession (cr_recession)**: `mu_b0, b1, b2, sigma_sp, sigma_eco, sigma`

## Why the schema is shaped this way

**Posterior summaries land in JSON, full draws stay as CSV.** Embedding 4,000 draws × 16 parameters × 25 variants × 6 components in JSON would inflate each variant config from ~50 KB to several MB and slow runtime parsing. The summary statistics (mean / sd / quantiles) are sufficient for point estimates and Wald-style intervals; the runtime samplers in `config/uncertainty.py` load the draws CSV when full posterior propagation is needed.

**Species intercepts pre-filtered to the variant's species set.** Each variant's `FIAJSP` list defines its in-scope species. Storing the full 100+ CONUS species per variant would waste space and confuse readers; storing only the variant's species keeps each block self-contained.

**Ecodivision intercepts kept whole, weights stored alongside.** A variant may span multiple ecodivisions in non-uniform proportions. Pre-collapsing the ecodivisions into a single weighted RE here would lose information; the runtime sampler does the weighting using `ecodiv_weights` so different downstream tasks (e.g. spatial visualization vs point projection) can use the data differently.

**Modifier kept at one chosen lambda.** The lambda 5/10/20 sweep is a development artifact; production needs one number. The integration script picks the one specified by `--modifier-lambda` (default 10). The choice should be driven by held-out FIA RMSE per component — TODO: build a `63_select_modifier_lambda.R` that picks per component instead of global.

## `parameter_source_switch` extension for `config_loader.py`

The runtime layer needs to know which block to read from. Proposed change to `config/config_loader.py`:

```python
class ParameterSource(Enum):
    DEFAULT = "default"        # built-in FVS defaults
    CALIBRATED = "calibrated"  # per-variant Bayesian fits (existing categories.*)
    CONUS = "conus"            # CONUS Phase 4 fits (new categories_conus.*)
    HYBRID = "hybrid"          # CONUS where available, fall back to calibrated
    CUSTOM = "custom"          # user override


def load_parameters(variant: str, source: ParameterSource = ParameterSource.CALIBRATED, components: Optional[List[str]] = None) -> Dict:
    """
    Load parameters for `variant` using the requested `source` strategy.

    For `CONUS` source: reads `categories_conus.{component}` blocks. Raises if
    the requested component is not present (no automatic fallback).

    For `HYBRID` source: reads `categories_conus.{component}` if present, else
    falls back to `categories.{component}`. Components missing from BOTH raise.
    """
    cfg = _load_variant_config(variant)
    components = components or _all_components(cfg)
    out = {}

    for comp in components:
        if source == ParameterSource.CONUS:
            block = cfg.get("categories_conus", {}).get(comp)
            if block is None:
                raise KeyError(f"CONUS block missing for component {comp} in variant {variant}")
            out[comp] = _build_runtime_block_from_conus(block)
        elif source == ParameterSource.HYBRID:
            block = cfg.get("categories_conus", {}).get(comp)
            if block is not None:
                out[comp] = _build_runtime_block_from_conus(block)
            else:
                out[comp] = cfg["categories"][comp]
        elif source == ParameterSource.CALIBRATED:
            out[comp] = cfg["categories"][comp]
        # ... DEFAULT and CUSTOM unchanged ...

    return out


def _build_runtime_block_from_conus(block: Dict) -> Dict:
    """
    Convert the JSON-side categories_conus block into the in-memory format
    the FVS Fortran shim expects. Posterior means become point estimates;
    the species_intercepts SPCD list becomes the per-species coefficient
    vector indexed by FIAJSP order; ecodiv_weights are applied to compute
    the variant's effective ecodivision RE.
    """
    fe = block["fixed_effects"]
    fixed = dict(zip(fe["param"], fe["mean"]))
    species_re = dict(zip(block["species_intercepts"]["SPCD"], block["species_intercepts"]["mean"]))
    weights = block.get("ecodiv_weights", {})
    ecodiv_means = dict(zip(block["ecodiv_intercepts"]["ecodiv"], block["ecodiv_intercepts"]["mean"]))
    if weights:
        eco_re = sum(weights[e] * ecodiv_means[e] for e in weights)
    else:
        eco_re = 0.0
    return {
        "fixed": fixed,
        "species_re": species_re,
        "eco_re": eco_re,
        "modifier": dict(zip(block["modifier"]["coef"], block["modifier"]["mean"])) if block.get("modifier") else {},
        "_draws_csv": block.get("draws_csv"),
    }
```

## `uncertainty.py` extension

The existing `UncertaintyEngine` reads 500 draws per component from the per-variant `*_draws.json` files. The CONUS variant draws csv is shaped differently (one CSV per CONUS fit, draws are global not variant-specific), so the engine needs a parallel code path:

```python
class CONUSUncertaintyEngine:
    """
    Posterior draw engine for CONUS Phase 4 fits. Reads the global draws CSV
    once, applies the variant-specific species RE and ecodiv weighting at
    sample time.
    """
    def __init__(self, variant: str, component: str, n_draws: int = 500):
        cfg = _load_variant_config(variant)
        block = cfg["categories_conus"][component]
        self.draws = pd.read_csv(_resolve_path(block["draws_csv"]))
        self.species_re = dict(zip(block["species_intercepts"]["SPCD"],
                                   block["species_intercepts"]["mean"]))
        self.eco_re = self._weighted_eco_re(block)
        self.draw_indices = np.random.choice(len(self.draws), n_draws, replace=False)

    def get_draw(self, i: int, spcd: int) -> Dict:
        row = self.draws.iloc[self.draw_indices[i]]
        return {
            "fixed": row.to_dict(),
            "species_re": self.species_re.get(spcd, 0.0),
            "eco_re": self.eco_re,
        }
```

## Migration order

Once HG (the closest-to-done component) lands:

1. `Rscript 61_extract_conus_summaries.R --component hg --model hg_organon_fixedK_cspi_traits1`
2. `Rscript 62_conus_to_variant_json.R --variant ne --component hg --modifier-lambda 10 --dry-run` — verify the block shape
3. Drop `--dry-run` and write the NE config
4. Add the `ParameterSource.CONUS` enum value to `config_loader.py` and the `_build_runtime_block_from_conus` function
5. Add a unit test that loads NE with `source=CONUS, components=["height_growth"]` and asserts the parameters are non-empty
6. Run a one-variant FVS smoke test with the CONUS HG parameters to verify the runtime path works
7. Repeat for the other 24 variants once NE is validated
8. Repeat the whole sequence for the next component

The schema lock-in risk lives in steps 4 and 5: any changes to the block shape after this point require both the R writer (62) and Python reader (config_loader) to be updated together.
