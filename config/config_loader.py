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
    VALID_VERSIONS = ("default", "calibrated", "conus", "hybrid", "custom")

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
        if self.version in ("calibrated", "conus", "hybrid"):
            # All three live in the same variant JSON; the version selects
            # which top-level block to read (categories vs categories_conus).
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

        fe = block.get("fixed_effects", {})
        params = fe.get("param", [])
        means = fe.get("mean", [])
        fixed = dict(zip(params, means))

        si = block.get("species_intercepts", {})
        spcds = si.get("SPCD") or si.get("idx") or []
        sp_means = si.get("mean", [])
        species_re = {int(s): float(m) for s, m in zip(spcds, sp_means)}

        ei = block.get("ecodiv_intercepts", {})
        eco_codes = ei.get("ecodiv") or ei.get("idx") or []
        eco_means = ei.get("mean", [])
        eco_lookup = dict(zip([str(c) for c in eco_codes], eco_means))

        weights = block.get("ecodiv_weights", {})
        if weights:
            eco_re = sum(
                float(w) * float(eco_lookup.get(str(eco), 0.0))
                for eco, w in weights.items()
            )
        else:
            # Unweighted fallback: simple mean across all ecodivisions
            eco_re = float(np.mean(eco_means)) if eco_means else 0.0

        mod = block.get("modifier", {})
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

    def generate_keywords(self, include_comments: bool = True) -> str:
        """Generate FVS keyword block from the calibrated configuration.

        These keywords can be appended to any FVS keyfile to apply the
        calibrated parameters. Works with all FVS invocation methods.

        Args:
            include_comments: Whether to include explanatory comments

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

        # SDIMAX keyword: sets maximum SDI per species
        sdi_values = self._find_sdi_param(cats)
        if sdi_values is not None:
            lines.append(self._format_sdimax_keywords(sdi_values, include_comments))

        # BAMAX keyword (for variants that use it)
        bamax_values = self._find_bamax_param(cats)
        if bamax_values is not None:
            lines.append(self._format_bamax_keywords(bamax_values, include_comments))

        # MORTMULT keyword: mortality rate multipliers per species
        mort_mult = self._compute_mortality_multipliers(cats)
        if mort_mult is not None:
            lines.append(self._format_mortmult_keywords(mort_mult, include_comments))

        # BAIMULT keyword: diameter growth multipliers per species
        growth_mult = self._compute_growth_multipliers(cats)
        if growth_mult is not None:
            lines.append(self._format_baimult_keywords(growth_mult, include_comments))

        # HTGMULT keyword: height growth multipliers per species
        hg_mult = self._compute_height_multipliers(cats)
        if hg_mult is not None:
            lines.append(self._format_htgmult_keywords(hg_mult, include_comments))

        return "\n".join(lines)

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
                lines.append(f"SDIMAX          {i + 1:10d}{val:10.1f}")
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
