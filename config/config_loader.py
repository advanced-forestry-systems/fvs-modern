"""
FVS Parameter Configuration Loader

Provides runtime selection among three parameter sets:

  - "default":    Original FVS parameters extracted from Fortran source code
  - "calibrated": Bayesian posterior estimates fit to national FIA data
  - "custom":     User supplied JSON calibrated to independent data
                  (cooperative plot networks, regional inventories, etc.)

Works with both fvs2py (shared library API) and microfvs (keyfile/subprocess).

Usage patterns:

  1. Python API (fvs2py):
       # National FIA calibration
       fvs = FVS(lib_path, config_version="calibrated")

       # Custom calibration from your own plot network
       fvs = FVS(lib_path, config_version="custom",
                  config_dir="/path/to/my_calibration")

  2. Keyfile injection (microfvs or standalone FVS):
       loader = FvsConfigLoader("ne", version="calibrated")
       keywords = loader.generate_keywords()
       # Append keywords to .key file before running FVS

  3. Custom calibration from independent data:
       loader = FvsConfigLoader("ne", version="custom",
                                custom_config="/path/to/my_ne.json")
       loader.apply_to_fvs(fvs)  # or loader.generate_keywords()

  4. Comparison mode:
       diff = FvsConfigLoader.compare("ne")
       # Returns parameter by parameter differences with credible intervals
"""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Any, Optional

import numpy as np

logger = logging.getLogger(__name__)


class FvsConfigLoader:
    """Loads and applies FVS variant parameters from JSON configuration files.

    Supports three config versions:
      - "default":    Original FVS parameters extracted from Fortran source
      - "calibrated": Bayesian posterior estimates from national FIA data
      - "custom":     User supplied JSON calibrated to independent data
    """

    # Auto detect project root from this file's location
    _CONFIG_DIR = Path(__file__).parent

    # "default":     Original FVS parameters extracted from Fortran source code
    # "calibrated":  Per-variant Bayesian posterior estimates (categories.*)
    # "conus":       CONUS Phase 4 fits (categories_conus.*), strict — raises if a
    #                requested component is not present in categories_conus
    # "hybrid":      categories_conus.* where available, falls back to
    #                categories.* per component. Recommended when CONUS
    #                integration is partial.
    # "custom":      User supplied JSON
    # "conus_sf":    Species-free (Leg B) coefficients from categories_conus_sf.
    #                Every species effect is computed from traits at runtime
    #                (standardized traits times gamma) plus ecoregion / forest
    #                type random effects. Generalizes to species with no
    #                per-species fit.
    # "conus_hybrid": Per species, prefer the Leg A per-species block where the
    #                species has a reliable per-species fit (hybrid source map
    #                marks it leg_a), else fall back to the species-free (Leg B)
    #                trait effect. Recommended default once both legs are landed.
    VALID_VERSIONS = ("default", "calibrated", "conus", "hybrid", "custom",
                      "conus_sf", "conus_hybrid", "conus_greg")

    def __init__(
        self,
        variant: str,
        version: str = "default",
        config_dir: Optional[str | Path] = None,
        custom_config: Optional[str | Path] = None,
    ):
        """Initialize config loader for a specific variant.

        Args:
            variant: FVS variant code (e.g., 'ne', 'ca', 'sn')
            version: 'default', 'calibrated', or 'custom'
            config_dir: Override path to config directory (where default
                and calibrated JSONs live)
            custom_config: Path to a user supplied JSON config file.
                Required when version='custom'. The JSON must follow the
                same schema as config/{variant}.json. This allows users
                to calibrate FVS to their own plot network data (e.g.,
                cooperative inventories, long term silvicultural trials,
                or regional permanent sample plots).
        """
        self.variant = variant.lower()
        self.version = version.lower()

        if self.version not in self.VALID_VERSIONS:
            raise ValueError(
                f"version must be one of {self.VALID_VERSIONS}, got '{self.version}'"
            )

        if self.version == "custom" and custom_config is None:
            raise ValueError(
                "custom_config path is required when version='custom'. "
                "Provide the path to your calibrated JSON file."
            )

        self._custom_config_path = Path(custom_config) if custom_config else None

        if config_dir is not None:
            self._config_dir = Path(config_dir)
        else:
            self._config_dir = self._CONFIG_DIR

        self._config: dict[str, Any] | None = None
        self._default_config: dict[str, Any] | None = None

    @property
    def config_path(self) -> Path:
        """Path to the active config file."""
        if self.version == "custom" and self._custom_config_path is not None:
            return self._custom_config_path
        if self.version in ("calibrated", "conus", "hybrid",
                            "conus_sf", "conus_hybrid", "conus_greg"):
            # All live in the same variant JSON; the version selects which
            # top-level block to read (categories / categories_conus /
            # categories_conus_sf).
            return self._config_dir / "calibrated" / f"{self.variant}.json"
        return self._config_dir / f"{self.variant}.json"

    @property
    def default_path(self) -> Path:
        """Path to the default (uncalibrated) config file."""
        return self._config_dir / f"{self.variant}.json"

    @property
    def calibrated_path(self) -> Path:
        """Path to the calibrated config file."""
        return self._config_dir / "calibrated" / f"{self.variant}.json"

    @property
    def has_calibrated(self) -> bool:
        """Whether a calibrated config exists for this variant."""
        return self.calibrated_path.exists()

    @property
    def config(self) -> dict[str, Any]:
        """Loaded configuration dictionary."""
        if self._config is None:
            self._config = self._load_json(self.config_path)
        return self._config

    @property
    def default_config(self) -> dict[str, Any]:
        """Default configuration (always available)."""
        if self._default_config is None:
            self._default_config = self._load_json(self.default_path)
        return self._default_config

    def _load_json(self, path: Path) -> dict[str, Any]:
        """Load a JSON config file."""
        if not path.exists():
            raise FileNotFoundError(
                f"Config file not found: {path}\n"
                f"Available configs: {list(self._config_dir.glob('*.json'))}"
            )
        with open(path) as f:
            return json.load(f)

    # =========================================================================
    # Parameter Access
    # =========================================================================

    def get_param(self, category: str, name: str) -> list | float:
        """Get a parameter value from the active config.

        Args:
            category: Parameter category (e.g., 'growth', 'mortality', 'bark_ratio')
            name: Parameter name (e.g., 'B1', 'BKRAT', 'SDICON')

        Returns:
            Parameter value (list for species indexed, scalar otherwise)
        """
        cats = self.config.get("categories", {})
        if category not in cats:
            raise KeyError(f"Category '{category}' not in config. Available: {list(cats.keys())}")
        if name not in cats[category]:
            raise KeyError(f"Parameter '{name}' not in category '{category}'. Available: {list(cats[category].keys())}")
        return cats[category][name]

    def get_species_params(self, category: str, name: str) -> np.ndarray:
        """Get a species indexed parameter as a numpy array.

        Args:
            category: Parameter category
            name: Parameter name

        Returns:
            numpy array of length maxsp
        """
        vals = self.get_param(category, name)
        if isinstance(vals, list):
            return np.array(vals, dtype=np.float64)
        return np.array([vals], dtype=np.float64)

    @property
    def maxsp(self) -> int:
        """Number of species in this variant."""
        return self.config.get("maxsp", 0)

    @property
    def calibration_metadata(self) -> dict | None:
        """Calibration metadata (only present in calibrated configs).

        The posterior_to_json pipeline writes this under the top level
        key `calibration`. Older exports used `calibration_metadata`. We
        accept either for backward compatibility.
        """
        return self.config.get("calibration") or self.config.get("calibration_metadata")

    # =========================================================================
    # CONUS Phase 4 access (categories_conus block)
    # =========================================================================

    @property
    def has_conus_block(self) -> bool:
        """Whether the variant config carries a categories_conus block."""
        return "categories_conus" in self.config

    def conus_components_present(self) -> list[str]:
        """Component names available under categories_conus."""
        block = self.config.get("categories_conus", {})
        return [k for k in block.keys() if k != "metadata"]

    def get_conus_block(self, component: str) -> dict:
        """Return the categories_conus block for `component`.

        Behavior depends on self.version:
          - "conus":  reads categories_conus.{component}; raises KeyError
                       if the component is not present (no fallback)
          - "hybrid": reads categories_conus.{component} if present,
                       otherwise falls back to categories.{component}
          - other:    raises (this access path is CONUS-specific)
        """
        if self.version not in ("conus", "hybrid"):
            raise ValueError(
                f"get_conus_block requires version in ('conus','hybrid'); "
                f"got '{self.version}'"
            )
        cu = self.config.get("categories_conus", {})
        if component in cu:
            return cu[component]
        if self.version == "hybrid":
            cats = self.config.get("categories", {})
            if component not in cats:
                raise KeyError(
                    f"Component '{component}' not in categories_conus or categories"
                )
            # Legacy block shape; signal source so callers know it's not CONUS
            return {"_source": "legacy_categories", "data": cats[component]}
        raise KeyError(
            f"CONUS block for '{component}' not present in {self.variant}; "
            f"available: {self.conus_components_present()}"
        )

    def get_conus_runtime_block(self, component: str) -> dict:
        """Decompose a categories_conus block into the runtime form.

        Returns a dict with:
          fixed:      {param_name: posterior_mean}
          species_re: {SPCD: posterior_mean}
          eco_re:     scalar (variant-specific weighted ecodivision RE)
          modifier:   {coef_name: posterior_mean}
          _draws_csv: path to the full posterior draws CSV (for uncertainty)

        For legacy fallback blocks (when version='hybrid' and component
        wasn't integrated yet), returns the legacy categories.{component}
        dict directly with _source='legacy_categories' so downstream code
        can branch on the source.
        """
        block = self.get_conus_block(component)
        if block.get("_source") == "legacy_categories":
            return block

        # Guard: R/jsonlite writes empty named lists as JSON [] which
        # deserialize to Python lists; coerce those back to {} so the
        # .get() accessors below are safe for components with empty
        # modifier / ecodiv_intercepts blocks (all cspi_traits1 fits).
        def _asdict(x):
            return x if isinstance(x, dict) else {}
        fe = _asdict(block.get("fixed_effects", {}))
        params = fe.get("param", [])
        means = fe.get("mean", [])
        fixed = dict(zip(params, means))

        si = _asdict(block.get("species_intercepts", {}))
        spcds = si.get("SPCD") or si.get("idx") or []
        sp_means = si.get("mean", [])
        species_re = {int(s): float(m) for s, m in zip(spcds, sp_means)}

        ei = _asdict(block.get("ecodiv_intercepts", {}))
        eco_codes = ei.get("ecodiv") or ei.get("idx") or []
        eco_means = ei.get("mean", [])
        eco_lookup = dict(zip([str(c) for c in eco_codes], eco_means))

        weights = _asdict(block.get("ecodiv_weights", {}))
        if weights:
            eco_re = sum(
                float(w) * float(eco_lookup.get(str(eco), 0.0))
                for eco, w in weights.items()
            )
        else:
            # Unweighted fallback: simple mean across all ecodivisions
            eco_re = float(np.mean(eco_means)) if eco_means else 0.0

        mod = _asdict(block.get("modifier", {}))
        mod_coefs = mod.get("coef", [])
        mod_means = mod.get("mean", [])
        modifier = dict(zip(mod_coefs, mod_means))

        return {
            "_source": "categories_conus",
            "model": block.get("model"),
            "modifier_lambda": block.get("modifier_lambda"),
            "fixed": fixed,
            "species_re": species_re,
            "eco_re": eco_re,
            "modifier": modifier,
            "_draws_csv": block.get("draws_csv"),
            "_species_missing": block.get("species_missing", []),
        }

    def conus_summary(self) -> dict:
        """One-line summary of CONUS integration state for this variant.

        Useful for logging and for the comparison report. Returns:
          {
            "has_block": bool,
            "components": [str, ...],
            "metadata": {... categories_conus.metadata ...} or None
          }
        """
        block = self.config.get("categories_conus", {})
        return {
            "has_block": bool(block),
            "components": self.conus_components_present(),
            "metadata": block.get("metadata"),
        }

    # =========================================================================
    # Species-free (Leg B) access (categories_conus_sf block)
    # =========================================================================

    @property
    def has_conus_sf_block(self) -> bool:
        """Whether the variant config carries a categories_conus_sf block."""
        return "categories_conus_sf" in self.config

    def conus_sf_components_present(self) -> list[str]:
        """Component names available under categories_conus_sf."""
        block = self.config.get("categories_conus_sf", {})
        return [k for k in block.keys() if k != "metadata"]

    # -------------------------------------------------------------------------
    # Greg arm (categories_conus_greg) + keyword-selectable site driver
    # -------------------------------------------------------------------------
    def has_conus_greg_block(self) -> bool:
        """Whether the variant config carries a categories_conus_greg block."""
        return "categories_conus_greg" in self.config

    def get_conus_greg_block(self, component: str) -> dict:
        """Raw categories_conus_greg.components.{component} block, or KeyError."""
        greg = self.config.get("categories_conus_greg", {})
        comps = greg.get("components", {}) if isinstance(greg, dict) else {}
        if component not in comps:
            raise KeyError(
                f"Greg-arm block for '{component}' not present in {self.variant}; "
                f"available: {list(comps.keys())}"
            )
        return comps[component]

    def get_greg_driver_coefficients(self, component, driver=None):
        """Resolve the site-driver coefficient file for a Greg-arm component.

        The site driver is keyword-selectable (DGDRIVER for diameter_growth,
        MORTDRIVER for survival). If `driver` is None, use the component's
        default `site_driver`. Returns the coefficient filename (relative to
        the config dir). Raises if the driver is not an offered option.
        A/B evidence (2026-07-03): DG default 'cspi' (driver is a minor lever),
        survival default 'bgi' (driver matters, ~5.7% log-loss gain).
        """
        b = self.get_conus_greg_block(component)
        drv = driver or b.get("site_driver")
        opts = b.get("site_driver_options", []) or []
        cmap = b.get("coefficients_by_driver", {}) or {}
        if opts and drv not in opts:
            raise ValueError(
                f"driver '{drv}' not offered for {component}; options: {opts}"
            )
        if drv not in cmap:
            raise KeyError(
                f"no coefficient set for driver '{drv}' in {component}; "
                f"have: {list(cmap.keys())}"
            )
        return cmap[drv]

    # =========================================================================
    # Per-variant DG native-hook policy (DGDRIVER)
    # =========================================================================
    # Source of truth: config/dg_hook_policy.json. Decided from an in-engine
    # per-variant native-bias screen. Enable the Greg diameter-growth hook (with
    # the REFIT coefficient set, DGDRIVER code 2) only where native FVS DG
    # under-predicts; keep native FVS DG (no DGDRIVER keyword) elsewhere. The
    # keyword is parsed engine-side (initre.f90 option 149 -> IDGDRV; gregdghg.f90
    # maps code 2 -> config/greg_dg_coefficients_refit.csv). This is opt-in and
    # backward compatible: nothing here changes runtime DG behavior unless a run
    # explicitly emits the keyword via generate_keywords(apply_dg_hook_policy=True)
    # or the dgdriver_keyword() helper.

    _DG_HOOK_POLICY_FILE = "dg_hook_policy.json"
    _dg_hook_policy_cache: dict | None = None

    @classmethod
    def load_dg_hook_policy_table(cls, config_dir=None) -> dict:
        """Load and cache the full per-variant DG hook policy table.

        Returns the parsed dg_hook_policy.json (keys: '_meta', 'variants').
        """
        cdir = Path(config_dir) if config_dir is not None else cls._CONFIG_DIR
        path = cdir / cls._DG_HOOK_POLICY_FILE
        if cls._dg_hook_policy_cache is None or config_dir is not None:
            with open(path) as f:
                table = json.load(f)
            if config_dir is None:
                cls._dg_hook_policy_cache = table
            return table
        return cls._dg_hook_policy_cache

    def dg_hook_policy(self, variant: str | None = None) -> dict:
        """Return the DG hook policy dict for a variant.

        Args:
            variant: FVS variant code (e.g. 'ne'). Defaults to this loader's
                own variant.

        Returns:
            Dict with keys: hook (bool), dgdriver_code (2 for hook / None for
            native), classification ('HOOK'/'NATIVE'/'DEFER'), basis (str),
            native_bias (value or None), note (str).

        Raises:
            KeyError if the variant is not in the policy table.
        """
        v = (variant or self.variant).lower()
        table = self.load_dg_hook_policy_table(self._config_dir)
        variants = table.get("variants", {})
        if v not in variants:
            raise KeyError(
                f"variant '{v}' not in DG hook policy table; "
                f"have: {sorted(variants.keys())}"
            )
        return variants[v]

    def dgdriver_keyword(self, variant: str | None = None) -> str | None:
        """Return the FVS DGDRIVER keyword line for a variant per policy.

        HOOK variants -> the keyword string that enables the refit Greg DG hook
        (DGDRIVER code 2). NATIVE and DEFER variants -> None (emit nothing;
        native FVS DG is used).

        The keyword is fixed-format: keyword left-justified to 10 columns, the
        numeric field in the next 10-column field. A keyword-file builder
        appends the returned string (when not None) to the .key file.
        """
        pol = self.dg_hook_policy(variant)
        if not pol.get("hook"):
            return None
        code = pol.get("dgdriver_code")
        if code is None:
            return None
        # Fixed 10-col keyword + 10-col integer field (matches FVS keyword parser,
        # e.g. initre.f90 option 149 reads ARRAY(1) from the first field).
        return "%-10s%10d" % ("DGDRIVER", int(code))

    # =========================================================================
    # Per-variant mortality tier policy (MORTDRIVER)
    # =========================================================================
    # Source of truth: config/mort_hook_policy.json. Decided from a per-variant
    # 4-arm A/B on held-out FIA remeasurement pairs (native FVS, crown-only
    # gompit, size-extended gompit, size+BGI gompit; log-loss/AUC/bias). Unlike
    # DGDRIVER, this is NOT a swappable site-driver family: the merged GOMPSURV
    # form is an additive TIER structure (crown-only -> +size -> +size+BGI).
    # size+BGI (MORTDRIVER code 2) wins on both log-loss and AUC in all six
    # confirmed variants and is the recommended opt-in upgrade; crown-only
    # (code 0 / no keyword) remains the engine's byte-identical default. The
    # keyword is parsed engine-side (initre.f90 option 150 -> IMORTDRV) but is
    # a reproducibility-logging device only for now (GOMPLOAD logs IMORTDRV);
    # actual tier selection happens by which coefficient file is passed to
    # FVS_GOMPIT_COEF. This is opt-in and backward compatible: nothing here
    # changes runtime mortality behavior unless a run explicitly emits the
    # keyword via the mortdriver_keyword() helper (or an equivalent
    # generate_keywords()-style policy application).

    _MORT_HOOK_POLICY_FILE = "mort_hook_policy.json"
    _mort_hook_policy_cache: dict | None = None

    @classmethod
    def load_mort_hook_policy_table(cls, config_dir=None) -> dict:
        """Load and cache the full per-variant mortality tier policy table.

        Returns the parsed mort_hook_policy.json (keys: '_meta', 'variants').
        """
        cdir = Path(config_dir) if config_dir is not None else cls._CONFIG_DIR
        path = cdir / cls._MORT_HOOK_POLICY_FILE
        if cls._mort_hook_policy_cache is None or config_dir is not None:
            with open(path) as f:
                table = json.load(f)
            if config_dir is None:
                cls._mort_hook_policy_cache = table
            return table
        return cls._mort_hook_policy_cache

    def mort_hook_policy(self, variant: str | None = None) -> dict:
        """Return the mortality tier policy dict for a variant.

        Args:
            variant: FVS variant code (e.g. 'ne'). Defaults to this loader's
                own variant.

        Returns:
            Dict with keys: recommended_tier (str), mortdriver_code (int),
            log_loss/auc/bias (float, for the recommended tier), n_trees
            (int), confidence ('high'/'low'), note (str), and a nested
            'tiers' dict with native_fvs/crown_only/size/size_bgi arms.

        Raises:
            KeyError if the variant is not in the six-variant confirmed table
            (other CONUS variants are unscreened for mortality tiers).
        """
        v = (variant or self.variant).lower()
        table = self.load_mort_hook_policy_table(self._config_dir)
        variants = table.get("variants", {})
        if v not in variants:
            raise KeyError(
                f"variant '{v}' not in mortality tier policy table (unscreened); "
                f"confirmed variants: {sorted(variants.keys())}"
            )
        return variants[v]

    def mortdriver_keyword(self, variant: str | None = None) -> str | None:
        """Return the FVS MORTDRVR keyword line for a variant per policy.

        Variants with a recommended tier whose code > 0 (i.e. size or
        size+BGI) get the keyword string that records that tier choice.
        A recommended code of 0 (crown-only) returns None: crown-only is the
        engine default, so no keyword needs to be emitted.

        The keyword is fixed-format: keyword left-justified to 10 columns
        (MORTDRVR, the FVS 8-char-limited form of MORTDRIVER), the numeric
        field in the next 10-column field. A keyword-file builder appends the
        returned string (when not None) to the .key file.
        """
        pol = self.mort_hook_policy(variant)
        code = pol.get("mortdriver_code")
        if code is None or int(code) <= 0:
            return None
        # Fixed 10-col keyword + 10-col integer field (matches FVS keyword parser,
        # e.g. initre.f90 option 150 reads ARRAY(1) from the first field).
        return "%-10s%10d" % ("MORTDRVR", int(code))

    def get_conus_sf_block(self, component: str) -> dict:
        """Raw categories_conus_sf.{component} block, or raise KeyError."""
        sf = self.config.get("categories_conus_sf", {})
        if component not in sf:
            raise KeyError(
                f"Species-free block for '{component}' not present in "
                f"{self.variant}; available: {self.conus_sf_components_present()}"
            )
        return sf[component]

    def get_conus_sf_runtime_block(self, component: str) -> dict:
        """Decompose a categories_conus_sf block into runtime lookups.

        Returns a dict with:
          intercept     scalar global intercept (a0 / b0 / h0)
          fixed         {param: mean} for all global fixed effects (incl sigmas)
          covariate     {param: mean} fixed effects excluding intercept + scale
          gamma         {trait_col: gamma_mean}
          scale         {trait_col: (mean, sd)}  standardization constants
          trait_effect  {SPCD: precomputed trait effect}
          raw_traits    {SPCD: {trait_col: raw value}}
          re_L1/re_L2/re_L3/re_FT  {str(level): RE mean}
          source_map    {SPCD: "leg_a" | "leg_b"}
        """
        b = self.get_conus_sf_block(component)

        # Guard: R/jsonlite writes empty named lists as JSON [] which
        # deserialize to Python lists; coerce to {} so the .get()
        # accessors below are safe for SF blocks with empty
        # species / hybrid_source_map / trait_gamma (e.g. pure
        # species-free components with no Leg A overrides).
        def _asdict_sf(x):
            return x if isinstance(x, dict) else {}
        fe = _asdict_sf(b.get("fixed_effects", {}))
        fixed = dict(zip(fe.get("param", []), fe.get("mean", [])))
        # intercept = a0 / b0 / h0 if present, else first non-sigma param
        intercept_name = next(
            (p for p in ("a0", "b0", "h0") if p in fixed),
            next((p for p in fixed if not p.startswith(("sigma", "phi"))), None),
        )
        intercept = float(fixed.get(intercept_name, 0.0)) if intercept_name else 0.0
        covariate = {
            p: float(m) for p, m in fixed.items()
            if p != intercept_name and not p.startswith(("sigma", "phi", "lp__"))
        }

        tg = _asdict_sf(b.get("trait_gamma", {}))
        cols = tg.get("trait_col", [])
        gamma = dict(zip(cols, tg.get("gamma_mean", [])))
        scale = {c: (float(m), float(s)) for c, m, s in
                 zip(cols, tg.get("scale_mean", []), tg.get("scale_sd", []))}

        def _f(v):
            try:
                return float(v)
            except (TypeError, ValueError):
                return float("nan")  # bundle stores missing raw traits as "NA"

        sp = _asdict_sf(b.get("species", {}))
        spcds = [int(x) for x in sp.get("SPCD", [])]
        te_vals = sp.get("trait_effect_mean", [])
        trait_effect = {s: _f(v) for s, v in zip(spcds, te_vals)}
        raw_traits: dict[int, dict[str, float]] = {s: {} for s in spcds}
        std_traits: dict[int, dict[str, float]] = {s: {} for s in spcds}
        for c in cols:
            rc = sp.get(f"raw_{c}")
            sc = sp.get(f"std_{c}")
            if rc is not None:
                for s, v in zip(spcds, rc):
                    raw_traits[s][c] = _f(v)
            if sc is not None:
                for s, v in zip(spcds, sc):
                    std_traits[s][c] = _f(v)

        def re_lookup(tag):
            t = b.get(tag)
            if not t:
                return {}
            return {str(l): float(m) for l, m in zip(t.get("level", []),
                                                     t.get("mean", []))}

        hm = _asdict_sf(b.get("hybrid_source_map", {}))
        source_map = {int(s): str(src) for s, src in
                      zip(hm.get("SPCD", []), hm.get("source", []))}

        return {
            "_source": "categories_conus_sf",
            "model": b.get("model"),
            "intercept_name": intercept_name,
            "intercept": intercept,
            "fixed": {k: float(v) for k, v in fixed.items()},
            "covariate": covariate,
            "gamma": gamma,
            "scale": scale,
            "trait_effect": trait_effect,
            "raw_traits": raw_traits,
            "std_traits": std_traits,
            "re_L1": re_lookup("re_L1"),
            "re_L2": re_lookup("re_L2"),
            "re_L3": re_lookup("re_L3"),
            "re_FT": re_lookup("re_FT"),
            "source_map": source_map,
        }

    @staticmethod
    def sf_trait_effect(rt: dict, spcd: int,
                        raw_traits: Optional[dict] = None) -> float:
        """Species effect from traits: standardized traits times gamma.

        Uses the precomputed value when the species is in the bundle. For a
        species outside the bundle, supply its raw trait values and they are
        standardized with the stored constants. Returns 0.0 if neither is
        available (engine should then fall back).
        """
        spcd = int(spcd)
        if spcd in rt["trait_effect"]:
            return rt["trait_effect"][spcd]
        traits = (raw_traits or {}).get(spcd) or rt["raw_traits"].get(spcd)
        if not traits:
            return 0.0
        eta = 0.0
        for c, g in rt["gamma"].items():
            if c in traits and c in rt["scale"]:
                val = traits[c]
                if val != val:        # NaN (missing trait) -> standardized 0
                    continue
                m, s = rt["scale"][c]
                if s and s != 0:
                    eta += ((val - m) / s) * g
        return eta

    def resolve_species_source(self, component: str, spcd: int) -> str:
        """For conus_hybrid: 'leg_a' if the species has a reliable per-species
        fit, else 'leg_b' (species-free trait fallback)."""
        try:
            rt = self.get_conus_sf_runtime_block(component)
        except KeyError:
            return "leg_b"
        return rt["source_map"].get(int(spcd), "leg_b")

    def sf_linear_predictor(self, component: str, spcd: int,
                            eco_codes: dict, covariates: dict,
                            runtime: Optional[dict] = None) -> float:
        """Runtime species-free linear predictor (the Leg B evaluator).

        eta = intercept + trait_effect(spcd)
            + z_L1[L1] + z_L2[L2] + z_L3[L3] + z_FT[FT]
            + sum_p covariate[p] * covariates[p]

        Args:
          component: category key (height_growth, diameter_growth,
                     height_crown_base, crown_recession, height_diameter,
                     mortality), matching the Leg A get_conus_block API.
          spcd: FIA species code.
          eco_codes: {"L1":code, "L2":code, "L3":code, "FT":code} (codes as in
                     the bundle RE tables; ints or strings both accepted).
          covariates: {param_name: value} for the covariate fixed effects
                      (e.g., {"a1": ln_dbh, "a2": ln_ht, ...}). Missing params
                      contribute zero.
        Note: for components with extra structure beyond the standard B1
        (HG v5 BGI site slopes), this evaluates the standard part; the BGI
        site-slope refinement is a documented extension.
        """
        rt = runtime or self.get_conus_sf_runtime_block(component)
        eta = rt["intercept"] + self.sf_trait_effect(rt, spcd)
        for tag, key in (("re_L1", "L1"), ("re_L2", "L2"),
                         ("re_L3", "L3"), ("re_FT", "FT")):
            code = eco_codes.get(key)
            if code is not None:
                eta += rt[tag].get(str(code), 0.0)
        for p, coef in rt["covariate"].items():
            eta += coef * float(covariates.get(p, 0.0))
        return eta

    # =========================================================================
    # fvs2py Integration
    # =========================================================================

    def apply_to_fvs(self, fvs_instance) -> dict[str, bool]:
        """Apply calibrated parameters to a running fvs2py FVS instance.

        Uses the existing set_species_attr() API to modify species level
        attributes. This should be called after load_keyfile() but before
        run(), ideally at stop point 7 (after input read, before imputation).

        Args:
            fvs_instance: An initialized fvs2py.FVS object

        Returns:
            Dictionary of which attributes were successfully applied
        """
        applied = {}
        cats = self.config.get("categories", {})

        # SDI max: maps to 'spsdi' species attribute
        sdi_param = self._find_sdi_param(cats)
        if sdi_param is not None:
            try:
                # Replace 'NA' strings with 0 before conversion
                cleaned = [0.0 if isinstance(v, str) else v for v in sdi_param]
                arr = np.array(cleaned, dtype=np.float64)
                # Pad or trim to match variant's maxspecies
                dims = fvs_instance.dims
                maxsp = dims.get("maxspecies", len(arr))
                if len(arr) < maxsp:
                    arr = np.pad(arr, (0, maxsp - len(arr)), constant_values=0)
                elif len(arr) > maxsp:
                    arr = arr[:maxsp]
                fvs_instance.set_species_attr("spsdi", arr)
                applied["spsdi"] = True
                logger.info(f"Applied calibrated SDI max for {self.variant}")
            except Exception as e:
                logger.warning(f"Could not apply SDI max: {e}")
                applied["spsdi"] = False

        # Basal area increment multiplier from calibration
        # The calibrated config stores growth multipliers relative to defaults
        growth_mult = self._compute_growth_multipliers(cats)
        if growth_mult is not None:
            try:
                dims = fvs_instance.dims
                maxsp = dims.get("maxspecies", len(growth_mult))
                arr = self._pad_array(growth_mult, maxsp)
                fvs_instance.set_species_attr("baimult", arr)
                applied["baimult"] = True
                logger.info(f"Applied calibrated BA increment multipliers for {self.variant}")
            except Exception as e:
                logger.warning(f"Could not apply BA increment multipliers: {e}")
                applied["baimult"] = False

        # Mortality multiplier
        mort_mult = self._compute_mortality_multipliers(cats)
        if mort_mult is not None:
            try:
                dims = fvs_instance.dims
                maxsp = dims.get("maxspecies", len(mort_mult))
                arr = self._pad_array(mort_mult, maxsp)
                fvs_instance.set_species_attr("mortmult", arr)
                applied["mortmult"] = True
                logger.info(f"Applied calibrated mortality multipliers for {self.variant}")
            except Exception as e:
                logger.warning(f"Could not apply mortality multipliers: {e}")
                applied["mortmult"] = False

        # Height growth multiplier
        hg_mult = self._compute_height_multipliers(cats)
        if hg_mult is not None:
            try:
                dims = fvs_instance.dims
                maxsp = dims.get("maxspecies", len(hg_mult))
                arr = self._pad_array(hg_mult, maxsp)
                fvs_instance.set_species_attr("htgmult", arr)
                applied["htgmult"] = True
                logger.info(f"Applied calibrated height growth multipliers for {self.variant}")
            except Exception as e:
                logger.warning(f"Could not apply height growth multipliers: {e}")
                applied["htgmult"] = False

        return applied

    # =========================================================================
    # Keyfile Keyword Generation (for microfvs or standalone FVS)
    # =========================================================================

    def generate_keywords(self, include_comments: bool = True,
                          apply_dg_hook_policy: bool = False) -> str:
        """Generate FVS keyword block from the calibrated configuration.

        These keywords can be appended to any FVS keyfile to apply the
        calibrated parameters. Works with all FVS invocation methods.

        Args:
            include_comments: Whether to include explanatory comments
            apply_dg_hook_policy: Opt-in. When True, emit the per-variant
                DGDRIVER keyword from config/dg_hook_policy.json (DGDRIVER 2
                for HOOK variants; nothing for NATIVE/DEFER). Default False
                keeps existing behavior and the v2026.07-calibrated tag
                unaffected.

        Returns:
            String of FVS keywords ready to insert into a .key file
        """
        lines = []
        cats = self.config.get("categories", {})
        meta = self.config.get("calibration_metadata", {})

        if include_comments:
            lines.append(
                f"!! Bayesian calibrated parameters for variant {self.variant.upper()}"
            )
            if meta:
                lines.append(f"!! Calibration date: {meta.get('calibration_date', 'unknown')}")
                lines.append(f"!! Components: {meta.get('components_updated', [])}")
            lines.append("!!")

        # SDIMAX keyword: FORMAT FIXED 2026-06-29 (was DISABLED 2026-06-16, WO-1).
        # Root cause (initre.f90 option 89): the field ORDER is correct (species,
        # value); the real bug was the keyword written as "SDIMAX"+10 spaces, a
        # 16-char prefix that pushed species/value out of their fixed 10-col fields
        # so FVS misread them (garbage MAX SDI -> over-thinning). _format_sdimax_
        # keywords now left-justifies the keyword to 10 cols (matches the tested
        # %-10s%10d%10.1f in sdimax_binding_test.py). Localized max-SDI surfaces
        # exist (brms_SDImax_site_specific.csv, alphaearth maxsdi_*_lookup.csv);
        # NA dropouts handled by make_sdifix. Note sdimax is often non-binding
        # (growth-engine-dominated, per the 2026-06-10 audit). Re-enable per variant
        # via config "_emit_sdimax": true once the NE binding test confirms binding.
        emit_sdimax = bool(self.config.get("_emit_sdimax", True))
        sdi_values = self._find_sdi_param(cats) if emit_sdimax else None
        if sdi_values is not None:
            lines.append(self._format_sdimax_keywords(sdi_values, include_comments))

        # BAMAX keyword (for variants that use it)
        bamax_values = self._find_bamax_param(cats)
        if bamax_values is not None:
            lines.append(self._format_bamax_keywords(bamax_values, include_comments))

        # Per-species multipliers. Prefer the precomputed calibration_multipliers
        # block written by the R serializer (calibration/R/06_posterior_to_json.R
        # via multipliers.R); fall back to the legacy on-the-fly computation from
        # raw coefficients for older configs that predate that block. This is the
        # fix for the "at most one keyword block per variant" gap (issue #54): the
        # legacy path returned None for mortality/growth/height because the writer
        # and emitter used mismatched coefficient schemas.
        precomputed = self.config.get("calibration_multipliers") or {}

        # MORTMULT keyword: mortality rate multipliers per species
        mort_mult = self._array_or_none(precomputed.get("mort_multiplier"))
        if mort_mult is None:
            mort_mult = self._compute_mortality_multipliers(cats)
        if mort_mult is not None:
            lines.append(self._format_mortmult_keywords(mort_mult, include_comments))

        # BAIMULT keyword: diameter growth multipliers per species. The precomputed
        # value is on the DDS scale (exp(delta b0)); convert to the diameter-growth
        # scale via sqrt, matching the legacy _compute_growth_multipliers convention.
        dds = self._array_or_none(precomputed.get("dds_multiplier"))
        if dds is not None:
            growth_mult = np.sqrt(np.clip(dds, 0.01, 100.0))
        else:
            growth_mult = self._compute_growth_multipliers(cats)
        if growth_mult is not None:
            lines.append(self._format_baimult_keywords(growth_mult, include_comments))

        # HTGMULT keyword: height growth multipliers per species
        hg_mult = self._array_or_none(precomputed.get("htg_multiplier"))
        if hg_mult is None:
            hg_mult = self._compute_height_multipliers(cats)
        if hg_mult is not None:
            lines.append(self._format_htgmult_keywords(hg_mult, include_comments))

        # Per-variant DG native-hook policy (opt-in, backward compatible).
        # HOOK variants get 'DGDRIVER 2' (refit Greg DG coefficients);
        # NATIVE/DEFER variants emit nothing and use native FVS DG.
        if apply_dg_hook_policy:
            dg_kw = self.dgdriver_keyword()
            if dg_kw is not None:
                if include_comments:
                    pol = self.dg_hook_policy()
                    lines.append(
                        f"!! DG hook policy: {pol.get('basis','')} "
                        f"(class {pol.get('classification','')})"
                    )
                lines.append(dg_kw)

        return "\n".join(lines)

    @staticmethod
    def _array_or_none(seq) -> "np.ndarray | None":
        """Coerce a JSON list of per-species multipliers to a float ndarray.

        Returns None when the input is missing, non-numeric, empty, or all-NaN so
        the caller can fall back to the legacy computation path. NaNs map to 1.0
        (no-op multiplier).
        """
        if seq is None:
            return None
        try:
            arr = np.asarray(seq, dtype=np.float64)
        except (TypeError, ValueError):
            return None
        if arr.size == 0 or bool(np.all(np.isnan(arr))):
            return None
        return np.where(np.isnan(arr), 1.0, arr)

    def _format_sdimax_keywords(self, values: list, comments: bool) -> str:
        """Format SDIMAX keyword block."""
        lines = []
        if comments:
            lines.append("!! Species specific SDI maximums")
        for i, val in enumerate(values):
            if isinstance(val, str) or val is None:
                continue
            if val > 0:
                # SDIMAX keyword: species_index  sdi_value
                lines.append(f"{'SDIMAX':<10}{i + 1:10d}{val:10.1f}")  # fixed: 10-col keyword (was SDIMAX+10 spaces, misaligned fields)
        return "\n".join(lines)

    def _format_bamax_keywords(self, values: list, comments: bool) -> str:
        """Format BAMAX keyword block."""
        lines = []
        if comments:
            lines.append("!! Maximum basal area")
        # BAMAX uses a single value or per species
        if isinstance(values, (int, float)):
            lines.append(f"BAMAX           {values:10.1f}")
        else:
            for i, val in enumerate(values):
                if isinstance(val, str) or val is None:
                    continue
                if val > 0:
                    lines.append(f"BAMAX           {i + 1:10d}{val:10.1f}")
        return "\n".join(lines)

    def _format_mortmult_keywords(self, multipliers: np.ndarray, comments: bool) -> str:
        """Format MORTMULT keyword block."""
        lines = []
        if comments:
            lines.append("!! Mortality rate multipliers (calibrated / default)")
        for i, mult in enumerate(multipliers):
            if abs(mult - 1.0) > 0.01:  # Only include if meaningfully different from 1.0
                # MORTMULT fields: species  proportion  lower_dbh  upper_dbh
                lines.append(f"MORTMULT        {i + 1:10d}{mult:10.4f}       0.0     999.0")
        return "\n".join(lines)

    def _format_baimult_keywords(self, multipliers: np.ndarray, comments: bool) -> str:
        """Format BAIMULT (growth multiplier) keyword block."""
        lines = []
        if comments:
            lines.append("!! Diameter growth multipliers (calibrated / default)")
        for i, mult in enumerate(multipliers):
            if abs(mult - 1.0) > 0.01:
                # BAIMULT via READCORD or GROWTH multiplier approach
                # Using species level growth multiplier
                lines.append(f"BAIMULT         {i + 1:10d}{mult:10.4f}")
        return "\n".join(lines)

    def _format_htgmult_keywords(self, multipliers: np.ndarray, comments: bool) -> str:
        """Format height growth multiplier keyword block."""
        lines = []
        if comments:
            lines.append("!! Height growth multipliers (calibrated / default)")
        for i, mult in enumerate(multipliers):
            if abs(mult - 1.0) > 0.01:
                lines.append(f"HTGMULT         {i + 1:10d}{mult:10.4f}")
        return "\n".join(lines)

    # =========================================================================
    # Comparison / Diagnostics
    # =========================================================================

    @classmethod
    def compare(
        cls,
        variant: str,
        config_dir: Optional[str | Path] = None,
    ) -> dict[str, Any]:
        """Compare default and calibrated parameters for a variant.

        Returns:
            Dictionary with per parameter comparisons including:
              - default_value, calibrated_value
              - percent_change
              - credible_interval (if available)
        """
        default = cls(variant, "default", config_dir)
        calibrated = cls(variant, "calibrated", config_dir)

        comparisons = {}
        default_cats = default.config.get("categories", {})
        calibrated_cats = calibrated.config.get("categories", {})

        for cat_name, cat_params in default_cats.items():
            if cat_name not in calibrated_cats:
                continue

            for param_name, default_val in cat_params.items():
                if param_name not in calibrated_cats[cat_name]:
                    continue

                cal_val = calibrated_cats[cat_name][param_name]

                if isinstance(default_val, list) and isinstance(cal_val, list):
                    # Skip non numeric arrays (species code tables, flags, etc.)
                    def _all_numeric(seq):
                        return all(isinstance(e, (int, float)) and not isinstance(e, bool) for e in seq)
                    if not (_all_numeric(default_val) and _all_numeric(cal_val)):
                        continue
                    # Skip arrays of different lengths (variant specific schema
                    # differences; comparing element wise is not meaningful)
                    if len(default_val) != len(cal_val):
                        continue
                    if len(default_val) == 0:
                        continue
                    d_arr = np.array(default_val, dtype=np.float64)
                    c_arr = np.array(cal_val, dtype=np.float64)

                    # Compute pct change (avoid div by zero)
                    with np.errstate(divide="ignore", invalid="ignore"):
                        pct_change = np.where(
                            d_arr != 0,
                            100 * (c_arr - d_arr) / d_arr,
                            np.where(c_arr != 0, np.inf, 0),
                        )

                    comparisons[f"{cat_name}/{param_name}"] = {
                        "category": cat_name,
                        "parameter": param_name,
                        "n_values": len(default_val),
                        "default_mean": float(np.mean(d_arr[d_arr != 0])) if np.any(d_arr != 0) else 0,
                        "calibrated_mean": float(np.mean(c_arr[c_arr != 0])) if np.any(c_arr != 0) else 0,
                        "mean_pct_change": float(np.nanmean(pct_change[np.isfinite(pct_change)])),
                        "max_abs_pct_change": float(np.nanmax(np.abs(pct_change[np.isfinite(pct_change)]))) if np.any(np.isfinite(pct_change)) else 0,
                        "n_changed": int(np.sum(np.abs(pct_change[np.isfinite(pct_change)]) > 1)),
                    }
                elif isinstance(default_val, (int, float)) and isinstance(cal_val, (int, float)):
                    pct = (100 * (cal_val - default_val) / default_val) if default_val != 0 else 0
                    comparisons[f"{cat_name}/{param_name}"] = {
                        "category": cat_name,
                        "parameter": param_name,
                        "n_values": 1,
                        "default_value": default_val,
                        "calibrated_value": cal_val,
                        "pct_change": pct,
                    }

        # Add calibration metadata if present
        meta = calibrated.calibration_metadata
        if meta:
            comparisons["_metadata"] = meta

        return comparisons

    @classmethod
    def summary_table(
        cls,
        variant: str,
        config_dir: Optional[str | Path] = None,
    ) -> str:
        """Generate a human readable comparison summary.

        Returns:
            Formatted string table of parameter changes
        """
        diff = cls.compare(variant, config_dir)
        meta = diff.pop("_metadata", {})

        lines = [
            f"Parameter Comparison: {variant.upper()} (Default vs. Calibrated)",
            "=" * 70,
        ]

        if meta:
            # Accept both key conventions: newer posterior_to_json output
            # writes `date`; earlier exports used `calibration_date`.
            cal_date = meta.get("date") or meta.get("calibration_date") or "?"
            components = meta.get("components_updated") or []
            lines.append(f"Calibration date: {cal_date}")
            if components:
                lines.append(f"Components: {', '.join(components)}")
            lines.append("")

        lines.append(f"{'Parameter':<35} {'Default':>10} {'Calibrated':>10} {'Change':>10}")
        lines.append("-" * 70)

        for key, info in sorted(diff.items()):
            if info.get("n_values", 0) == 1:
                d = info.get("default_value", 0)
                c = info.get("calibrated_value", 0)
                pct = info.get("pct_change", 0)
                lines.append(f"{key:<35} {d:>10.2f} {c:>10.2f} {pct:>+9.1f}%")
            else:
                d = info.get("default_mean", 0)
                c = info.get("calibrated_mean", 0)
                pct = info.get("mean_pct_change", 0)
                n = info.get("n_changed", 0)
                lines.append(
                    f"{key:<35} {d:>10.2f} {c:>10.2f} {pct:>+9.1f}% ({n} spp changed)"
                )

        return "\n".join(lines)

    # =========================================================================
    # Internal Helpers
    # =========================================================================

    def _find_sdi_param(self, cats: dict) -> list | None:
        """Find the SDI maximum parameter (varies by variant)."""
        site_cats = cats.get("site_index", {})
        other_cats = cats.get("other", {})
        all_cats = {**site_cats, **other_cats}

        for name in ("SDICON", "R5SDI", "R4SDI", "FMSDI", "SDIDEF"):
            if name in all_cats:
                return all_cats[name]
        return None

    def _find_bamax_param(self, cats: dict) -> list | None:
        """Find the BAMAX parameter (varies by variant)."""
        site_cats = cats.get("site_index", {})
        other_cats = cats.get("other", {})
        all_cats = {**site_cats, **other_cats}

        for name in ("BAMAXA", "BAMAX1", "BAMAX"):
            if name in all_cats:
                return all_cats[name]
        return None

    def _compute_growth_multipliers(self, cats: dict) -> np.ndarray | None:
        """Compute diameter growth multipliers as ratio of calibrated to default.

        Returns multipliers where 1.0 = no change.
        """
        if self.version == "default":
            return None

        try:
            default_cats = self.default_config.get("categories", {})
            growth_cal = cats.get("growth", {})
            growth_def = default_cats.get("growth", {})

            # Use B1 coefficients as primary indicator of growth rate change
            if "B1" in growth_cal and "B1" in growth_def:
                cal = np.array(growth_cal["B1"], dtype=np.float64)
                default = np.array(growth_def["B1"], dtype=np.float64)

                # For Wykoff model: ln(DDS) = B0 + B1*ln(DBH) + ...
                # Multiplier on DDS scale = exp(cal_B0 - default_B0)
                # Simplified: use ratio of intercepts as overall growth multiplier
                if "B0" in growth_cal and "B0" in growth_def:
                    cal_b0 = np.array(growth_cal["B0"], dtype=np.float64)
                    def_b0 = np.array(growth_def["B0"], dtype=np.float64)
                    # exp(delta_B0) gives the multiplicative change on DDS scale
                    with np.errstate(invalid="ignore"):
                        multipliers = np.where(
                            def_b0 != 0,
                            np.exp(cal_b0 - def_b0),
                            1.0,
                        )
                    # Convert DDS multiplier to diameter multiplier (sqrt)
                    multipliers = np.sqrt(np.clip(multipliers, 0.1, 10.0))
                    return multipliers

        except Exception as e:
            logger.warning(f"Could not compute growth multipliers: {e}")

        return None

    def _compute_mortality_multipliers(self, cats: dict) -> np.ndarray | None:
        """Compute mortality multipliers from calibrated vs default coefficients."""
        if self.version == "default":
            return None

        try:
            default_cats = self.default_config.get("categories", {})
            mort_cal = cats.get("mortality", {})
            mort_def = default_cats.get("mortality", {})

            # Use intercept (B0) change to estimate overall mortality rate shift
            for b0_name in ("MORT_B0", "B0", "MRT_B0"):
                if b0_name in mort_cal and b0_name in mort_def:
                    cal_b0 = np.array(mort_cal[b0_name], dtype=np.float64)
                    def_b0 = np.array(mort_def[b0_name], dtype=np.float64)

                    # Logistic model: logit(p) = B0 + ...
                    # Odds ratio = exp(cal_B0 - def_B0)
                    with np.errstate(invalid="ignore"):
                        odds_ratio = np.where(
                            def_b0 != 0,
                            np.exp(cal_b0 - def_b0),
                            1.0,
                        )
                    return np.clip(odds_ratio, 0.1, 10.0)

        except Exception as e:
            logger.warning(f"Could not compute mortality multipliers: {e}")

        return None

    def _compute_height_multipliers(self, cats: dict) -> np.ndarray | None:
        """Compute height growth multipliers from calibrated config."""
        if self.version == "default":
            return None

        try:
            default_cats = self.default_config.get("categories", {})

            # Check for height growth parameters
            for cat_name in ("height_growth", "growth"):
                hg_cal = cats.get(cat_name, {})
                hg_def = default_cats.get(cat_name, {})

                for b0_name in ("HGLD", "HG_B0"):
                    if b0_name in hg_cal and b0_name in hg_def:
                        cal = np.array(hg_cal[b0_name], dtype=np.float64)
                        default = np.array(hg_def[b0_name], dtype=np.float64)

                        with np.errstate(divide="ignore", invalid="ignore"):
                            multipliers = np.where(
                                default != 0,
                                cal / default,
                                1.0,
                            )
                        return np.clip(multipliers, 0.1, 10.0)

        except Exception as e:
            logger.warning(f"Could not compute height growth multipliers: {e}")

        return None

    def _pad_array(self, arr: np.ndarray, target_len: int) -> np.ndarray:
        """Pad or trim an array to target length."""
        if len(arr) < target_len:
            return np.pad(arr, (0, target_len - len(arr)), constant_values=1.0)
        elif len(arr) > target_len:
            return arr[:target_len]
        return arr


# =============================================================================
# Convenience Functions
# =============================================================================

def load_config(variant: str, version: str = "default", **kwargs) -> FvsConfigLoader:
    """Convenience function to create a config loader.

    Args:
        variant: FVS variant code
        version: 'default' or 'calibrated'

    Returns:
        Initialized FvsConfigLoader
    """
    return FvsConfigLoader(variant, version, **kwargs)


def generate_calibration_keyfile(
    variant: str,
    output_path: Optional[str | Path] = None,
    config_dir: Optional[str | Path] = None,
) -> str:
    """Generate a standalone FVS keyword file with calibrated parameters.

    This file can be included in any FVS simulation via the ADDFILE keyword
    or by appending its contents to an existing keyfile.

    Args:
        variant: FVS variant code
        output_path: Where to save the keyword file (optional)
        config_dir: Override config directory path

    Returns:
        String contents of the keyword file
    """
    loader = FvsConfigLoader(variant, "calibrated", config_dir)
    keywords = loader.generate_keywords(include_comments=True)

    if output_path is not None:
        path = Path(output_path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with open(path, "w") as f:
            f.write(keywords)
        logger.info(f"Calibration keywords written to {path}")

    return keywords


def compare_configs(variant: str, **kwargs) -> str:
    """Print a comparison of default vs calibrated parameters.

    Args:
        variant: FVS variant code

    Returns:
        Formatted comparison table
    """
    return FvsConfigLoader.summary_table(variant, **kwargs)
