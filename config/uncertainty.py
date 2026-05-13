"""
FVS Uncertainty Propagation Engine

Provides Monte Carlo uncertainty estimation for FVS projections by sampling
from Bayesian posterior distributions of calibrated parameters. Each draw
from the posterior represents one plausible parameterization of the growth,
mortality, height, crown ratio, and stand density models. Running FVS
repeatedly with different draws produces an ensemble of projections whose
spread quantifies parameter uncertainty.

The engine preserves the joint posterior structure: within a single draw,
all component model parameters come from the same MCMC iteration, so
correlations between (for example) growth and mortality coefficients are
maintained.

Usage patterns:

  1. Full Monte Carlo ensemble via fvs2py:
       fvs = FVS(lib_path, config_version="calibrated",
                  uncertainty=True, n_draws=100)
       fvs.load_keyfile("stand.key")
       fvs.run()
       ensemble = fvs.uncertainty_results  # list of DataFrames

  2. Full Monte Carlo ensemble via microfvs:
       result = run_fvs(stand_init, tree_init,
                        config_version="calibrated",
                        uncertainty=True, n_draws=100)

  3. Manual draw selection:
       engine = UncertaintyEngine("ne", config_dir="/path/to/config")
       params = engine.sample_draw(42)  # deterministic draw index
       engine.apply_draw_to_fvs(fvs, params)

  4. Summarize ensemble:
       summary = engine.summarize_ensemble(results_list,
                                           quantiles=[0.025, 0.5, 0.975])
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any, Optional

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


class UncertaintyEngine:
    """Samples parameter vectors from Bayesian posterior draws for FVS.

    The engine loads a draws JSON file produced by the calibration pipeline
    (06_posterior_to_json.R) and provides methods for:
      - Drawing a complete parameter set from the posterior
      - Applying that parameter set to a running FVS instance
      - Generating FVS keywords for a specific draw (keyfile injection)
      - Summarizing an ensemble of FVS runs into credible intervals

    The draws file contains N posterior samples (typically 500) for each
    component model. All components share the same draw index, preserving
    the joint posterior correlation structure.
    """

    def __init__(
        self,
        variant: str,
        config_dir: Optional[str | Path] = None,
        draws_path: Optional[str | Path] = None,
        seed: Optional[int] = None,
    ):
        """Initialize the uncertainty engine.

        Args:
            variant: FVS variant code (e.g., 'ne', 'ca', 'sn')
            config_dir: Path to config directory containing calibrated/ subdir
            draws_path: Direct path to a draws JSON file (overrides config_dir)
            seed: Random seed for reproducible draw sequences
        """
        self.variant = variant.lower()
        self._rng = np.random.default_rng(seed)

        if draws_path is not None:
            self._draws_path = Path(draws_path)
        elif config_dir is not None:
            self._draws_path = (
                Path(config_dir) / "calibrated" / f"{self.variant}_draws.json"
            )
        else:
            self._draws_path = (
                Path(__file__).parent / "calibrated" / f"{self.variant}_draws.json"
            )

        self._draws_data: dict[str, Any] | None = None
        self._n_draws: int = 0
        self._components: list[str] = []

    @property
    def draws_available(self) -> bool:
        """Whether posterior draws are available for this variant."""
        return self._draws_path.exists()

    @property
    def n_draws(self) -> int:
        """Number of posterior draws available."""
        if self._draws_data is None:
            self._load_draws()
        return self._n_draws

    @property
    def components(self) -> list[str]:
        """Component models with available draws."""
        if self._draws_data is None:
            self._load_draws()
        return self._components

    def _load_draws(self) -> None:
        """Load posterior draws from JSON."""
        if not self._draws_path.exists():
            raise FileNotFoundError(
                f"Posterior draws not found: {self._draws_path}\n"
                f"Run the calibration pipeline with draw export enabled, "
                f"or check that 06_posterior_to_json.R completed successfully."
            )

        logger.info(f"Loading posterior draws from {self._draws_path}")
        with open(self._draws_path) as f:
            self._draws_data = json.load(f)

        self._n_draws = self._draws_data.get("n_draws", 0)
        self._components = self._draws_data.get("components", [])

        # Validate that draw arrays are consistent lengths
        draws = self._draws_data.get("draws", {})
        for comp, comp_draws in draws.items():
            actual_n = len(comp_draws)
            if actual_n != self._n_draws:
                logger.warning(
                    f"Component '{comp}' has {actual_n} draws, expected {self._n_draws}. "
                    f"Using min({actual_n}, {self._n_draws})."
                )
                self._n_draws = min(self._n_draws, actual_n)

        logger.info(
            f"Loaded {self._n_draws} draws for {len(self._components)} components: "
            f"{', '.join(self._components)}"
        )

    def sample_draw_index(self) -> int:
        """Sample a random draw index from the posterior.

        Returns:
            Integer index into the draws array
        """
        if self._draws_data is None:
            self._load_draws()
        return int(self._rng.integers(0, self._n_draws))

    def get_draw(self, draw_idx: int) -> dict[str, dict[str, float]]:
        """Get a complete parameter set for a specific draw index.

        Returns a nested dictionary: component -> parameter_name -> value.
        All components use the same draw_idx, preserving the joint
        posterior correlation structure.

        Args:
            draw_idx: Index into the posterior draws (0 to n_draws - 1)

        Returns:
            Dictionary mapping component names to parameter dictionaries
        """
        if self._draws_data is None:
            self._load_draws()

        if draw_idx < 0 or draw_idx >= self._n_draws:
            raise IndexError(
                f"draw_idx {draw_idx} out of range [0, {self._n_draws})"
            )

        draws = self._draws_data["draws"]
        result = {}

        for comp_name, comp_draws in draws.items():
            if draw_idx < len(comp_draws):
                # Each draw is a dict of param_name -> value
                draw = comp_draws[draw_idx]
                # Flatten any nested lists (JSON serialization artifact)
                flat = {}
                for k, v in draw.items():
                    if isinstance(v, list) and len(v) == 1:
                        flat[k] = v[0]
                    else:
                        flat[k] = v
                result[comp_name] = flat

        return result

    def sample_draw(self, draw_idx: Optional[int] = None) -> dict[str, dict[str, float]]:
        """Sample a parameter set from the posterior.

        If draw_idx is None, a random index is chosen. Otherwise uses the
        specified index for reproducibility.

        Args:
            draw_idx: Optional specific draw index

        Returns:
            Dictionary mapping component names to parameter dictionaries
        """
        if draw_idx is None:
            draw_idx = self.sample_draw_index()
        return self.get_draw(draw_idx)

    # =========================================================================
    # Apply Draws to FVS
    # =========================================================================

    def _build_fia_to_fvs_map(
        self, default_config: dict[str, Any]
    ) -> dict[int, int]:
        """Build mapping from FIA species codes to FVS species indices (1 based).

        The posterior draws use FIA species codes in random effect names
        (e.g. r_SPCD[12,Intercept]), but FVS keywords need species by
        sequential FVS index. The FIAJSP array in the config provides
        the mapping: index i (0 based) holds the FIA code for FVS species i+1.

        Args:
            default_config: Config dictionary with categories.species_definitions

        Returns:
            Dict mapping FIA code -> FVS species index (1 based)
        """
        sp_defs = default_config.get("categories", {}).get("species_definitions", {})
        fiajsp = sp_defs.get("FIAJSP", [])
        fia_to_fvs = {}
        for i, code_str in enumerate(fiajsp):
            try:
                fia_to_fvs[int(code_str)] = i + 1
            except (ValueError, TypeError):
                continue
        return fia_to_fvs

    def _reconstruct_hierarchical_intercepts(
        self,
        params: dict[str, float],
        fixed_name: str,
        random_prefix: str,
        n_species: int,
        fia_to_fvs: dict[int, int],
    ) -> np.ndarray | None:
        """Reconstruct species level intercepts from hierarchical brms draws.

        brms models parameterize as: species_intercept = b_Intercept + r_SPCD[code,Intercept]
        Stan hierarchical models use: species_b0 = mu_b0 + sigma_b0 * z_b0[i]

        This method detects the parameterization style and reconstructs an
        array of species level values indexed by FVS species number.

        Args:
            params: Parameter dictionary from a single draw
            fixed_name: Name of fixed effect (e.g. "b_Intercept", "mu_b0")
            random_prefix: Prefix for random effects (e.g. "r_SPCD", "z_b0")
            n_species: Number of FVS species slots
            fia_to_fvs: FIA code -> FVS index mapping

        Returns:
            Array of length n_species with species level values, or None
        """
        import re

        result = np.full(n_species, np.nan, dtype=np.float64)

        # Style 1: brms with r_SPCD[code,Intercept] random effects
        if random_prefix == "r_SPCD":
            fixed = params.get(fixed_name)
            if fixed is None:
                return None

            r_pattern = re.compile(r"r_SPCD\[(\d+),Intercept\]")
            found_any = False
            for key, val in params.items():
                m = r_pattern.match(key)
                if m:
                    fia_code = int(m.group(1))
                    fvs_idx = fia_to_fvs.get(fia_code)
                    if fvs_idx is not None and fvs_idx <= n_species:
                        result[fvs_idx - 1] = fixed + val
                        found_any = True

            # Species without random effects get the fixed intercept
            if found_any:
                result[np.isnan(result)] = fixed
                return result
            return None

        # Style 2: Stan non centered parameterization z_b0[i]
        if random_prefix == "z_b0":
            mu = params.get("mu_b0")
            sigma = params.get("sigma_b0")
            if mu is None or sigma is None:
                return None

            z_pattern = re.compile(r"z_b0\[(\d+)\]")
            z_vals = {}
            for key, val in params.items():
                m = z_pattern.match(key)
                if m:
                    z_vals[int(m.group(1))] = val

            if not z_vals:
                return None

            # z_b0 indices are sequential (1..N_species in the Stan model)
            # Map them to FVS indices via FIAJSP ordering
            # The Stan model species ordering matches FIAJSP
            for stan_idx, z_val in z_vals.items():
                if stan_idx <= n_species:
                    result[stan_idx - 1] = mu + sigma * z_val

            result[np.isnan(result)] = mu
            return result

        return None

    def _reconstruct_brms_intercepts_double(
        self,
        params: dict[str, float],
        fixed_name: str,
        random_tag: str,
        n_species: int,
        fia_to_fvs: dict[int, int],
    ) -> np.ndarray | None:
        """Reconstruct species intercepts from brms with double underscore notation.

        Height diameter models use r_SPCD__a[code,Intercept] style.

        Args:
            params: Parameter dictionary from a single draw
            fixed_name: Name of fixed effect (e.g. "b_a_Intercept")
            random_tag: Tag between underscores (e.g. "a" for r_SPCD__a)
            n_species: Number of FVS species slots
            fia_to_fvs: FIA code -> FVS index mapping

        Returns:
            Array of length n_species with species level values, or None
        """
        import re

        fixed = params.get(fixed_name)
        if fixed is None:
            return None

        result = np.full(n_species, np.nan, dtype=np.float64)
        r_pattern = re.compile(
            rf"r_SPCD__{re.escape(random_tag)}\[(\d+),Intercept\]"
        )
        found_any = False
        for key, val in params.items():
            m = r_pattern.match(key)
            if m:
                fia_code = int(m.group(1))
                fvs_idx = fia_to_fvs.get(fia_code)
                if fvs_idx is not None and fvs_idx <= n_species:
                    result[fvs_idx - 1] = fixed + val
                    found_any = True

        if found_any:
            result[np.isnan(result)] = fixed
            return result
        return None

    def compute_multipliers_for_draw(
        self,
        draw: dict[str, dict[str, float]],
        default_config: dict[str, Any],
    ) -> dict[str, np.ndarray]:
        """Compute FVS multipliers from a posterior draw.

        Translates raw posterior parameter values into the multiplicative
        adjustments that FVS uses internally (GROWMULT, MORTMULT, HTGMULT,
        SDIMAX). This mirrors the logic in FvsConfigLoader but operates on
        a single draw rather than the posterior median.

        Handles two posterior parameterization styles:
          - Stan non centered: mu_b0 + sigma_b0 * z_b0[i] (diameter growth)
          - brms random effects: b_Intercept + r_SPCD[code,Intercept] (mortality,
            stand density) or b_a_Intercept + r_SPCD__a[code,Intercept] (height)

        FIA species codes in the random effects are mapped to FVS species indices
        using the FIAJSP array from the variant config.

        Args:
            draw: Parameter draw from get_draw() or sample_draw()
            default_config: The default (uncalibrated) config dictionary

        Returns:
            Dictionary with keys: growth_mult, mort_mult, hg_mult, sdi_values
        """
        result = {}
        default_cats = default_config.get("categories", {})
        maxsp = default_config.get("maxsp", 0)
        if maxsp == 0:
            sp_defs = default_cats.get("species_definitions", {})
            nsp = sp_defs.get("NSP", [])
            maxsp = int(nsp[0]) if nsp else 108

        fia_to_fvs = self._build_fia_to_fvs_map(default_config)

        # ---- Growth multiplier from diameter growth intercept shift ----
        if "diameter_growth" in draw:
            dg_params = draw["diameter_growth"]
            default_growth = default_cats.get("growth", {})

            # Try hierarchical Stan parameterization first (z_b0)
            cal_b0 = self._reconstruct_hierarchical_intercepts(
                dg_params, "mu_b0", "z_b0", maxsp, fia_to_fvs
            )
            # Fall back to simple indexed form
            if cal_b0 is None:
                cal_b0 = self._extract_species_vector(dg_params, "b0")

            # Get default intercepts: try B0, then B1, then WEIBB1
            def_b0_arr = None
            for bname in ("B0", "B1", "WEIBB1"):
                if bname in default_growth:
                    def_b0_arr = default_growth[bname]
                    break

            if cal_b0 is not None and def_b0_arr is not None:
                def_b0 = np.array(def_b0_arr, dtype=np.float64)
                n = min(len(cal_b0), len(def_b0))
                cal_b0 = cal_b0[:n]
                def_b0 = def_b0[:n]

                with np.errstate(invalid="ignore"):
                    multipliers = np.where(
                        def_b0 != 0,
                        np.sqrt(np.clip(np.exp(cal_b0 - def_b0), 0.1, 10.0)),
                        1.0,
                    )
                result["growth_mult"] = multipliers

        # ---- Mortality multiplier from mortality intercept shift ----
        if "mortality" in draw:
            mort_params = draw["mortality"]
            default_mort = default_cats.get("mortality", {})

            # Try brms hierarchical (b_Intercept + r_SPCD[code,Intercept])
            cal_b0 = self._reconstruct_hierarchical_intercepts(
                mort_params, "b_Intercept", "r_SPCD", maxsp, fia_to_fvs
            )
            # Fall back to simple indexed
            if cal_b0 is None:
                cal_b0 = self._extract_species_vector(mort_params, "b_Intercept")
            if cal_b0 is None:
                cal_b0 = self._extract_species_vector(mort_params, "b0")

            def_b0_name = None
            for name in ("MORT_B0", "B0", "MRT_B0"):
                if name in default_mort:
                    def_b0_name = name
                    break

            if cal_b0 is not None and def_b0_name is not None:
                def_b0 = np.array(default_mort[def_b0_name], dtype=np.float64)
                n = min(len(cal_b0), len(def_b0))

                with np.errstate(invalid="ignore"):
                    odds_ratio = np.where(
                        def_b0[:n] != 0,
                        np.exp(cal_b0[:n] - def_b0[:n]),
                        1.0,
                    )
                result["mort_mult"] = np.clip(odds_ratio, 0.1, 10.0)

        # ---- Height growth multiplier ----
        # Draws may use "height_diameter" or "height_increment" key
        hg_key = None
        for k in ("height_increment", "height_diameter"):
            if k in draw:
                hg_key = k
                break

        if hg_key is not None:
            hg_params = draw[hg_key]
            for cat_name in ("height_diameter", "height_growth", "growth"):
                hg_def = default_cats.get(cat_name, {})
                for b0_name in ("HD_A", "HGLD", "HG_B0"):
                    if b0_name in hg_def:
                        # Try brms double underscore notation (r_SPCD__a)
                        cal_vals = self._reconstruct_brms_intercepts_double(
                            hg_params, "b_a_Intercept", "a", maxsp, fia_to_fvs
                        )
                        if cal_vals is None:
                            cal_vals = self._extract_species_vector(hg_params, "b1")

                        if cal_vals is not None:
                            def_vals = np.array(hg_def[b0_name], dtype=np.float64)
                            n = min(len(cal_vals), len(def_vals))
                            with np.errstate(divide="ignore", invalid="ignore"):
                                mult = np.where(
                                    def_vals[:n] != 0,
                                    cal_vals[:n] / def_vals[:n],
                                    1.0,
                                )
                            result["hg_mult"] = np.clip(mult, 0.1, 10.0)
                            break
                if "hg_mult" in result:
                    break

        # ---- SDI max values from stand density ----
        if "stand_density" in draw:
            sdi_params = draw["stand_density"]

            # Try brms hierarchical: b_Intercept + r_SPCD[code,Intercept]
            # The stand density intercept is ln(SDImax), so exponentiate
            sdi_intercepts = self._reconstruct_hierarchical_intercepts(
                sdi_params, "b_Intercept", "r_SPCD", maxsp, fia_to_fvs
            )
            if sdi_intercepts is not None:
                # Intercepts are on log scale; convert to SDI values
                sdi_vals = np.exp(sdi_intercepts)
                # Clip to reasonable range (100 to 2000)
                result["sdi_values"] = np.clip(sdi_vals, 100, 2000)
            else:
                # Fall back to direct sdimax values
                sdi_vals = self._extract_species_vector(sdi_params, "sdimax")
                if sdi_vals is not None:
                    result["sdi_values"] = sdi_vals

        return result

    def apply_draw_to_fvs(
        self,
        fvs_instance,
        draw: dict[str, dict[str, float]],
        default_config: dict[str, Any],
    ) -> dict[str, bool]:
        """Apply a posterior draw to a running fvs2py FVS instance.

        Args:
            fvs_instance: An initialized fvs2py.FVS object
            draw: Parameter draw from get_draw() or sample_draw()
            default_config: The default config dictionary

        Returns:
            Dictionary indicating which attributes were applied
        """
        applied = {}
        multipliers = self.compute_multipliers_for_draw(draw, default_config)
        dims = fvs_instance.dims
        maxsp = dims.get("maxspecies", 0)

        if "growth_mult" in multipliers:
            try:
                arr = self._pad_array(multipliers["growth_mult"], maxsp)
                fvs_instance.set_species_attr("baimult", arr)
                applied["baimult"] = True
            except Exception as e:
                logger.warning(f"Could not apply growth multiplier: {e}")
                applied["baimult"] = False

        if "mort_mult" in multipliers:
            try:
                arr = self._pad_array(multipliers["mort_mult"], maxsp)
                fvs_instance.set_species_attr("mortmult", arr)
                applied["mortmult"] = True
            except Exception as e:
                logger.warning(f"Could not apply mortality multiplier: {e}")
                applied["mortmult"] = False

        if "hg_mult" in multipliers:
            try:
                arr = self._pad_array(multipliers["hg_mult"], maxsp)
                fvs_instance.set_species_attr("htgmult", arr)
                applied["htgmult"] = True
            except Exception as e:
                logger.warning(f"Could not apply height growth multiplier: {e}")
                applied["htgmult"] = False

        if "sdi_values" in multipliers:
            try:
                arr = self._pad_array(multipliers["sdi_values"], maxsp, fill=0)
                fvs_instance.set_species_attr("spsdi", arr)
                applied["spsdi"] = True
            except Exception as e:
                logger.warning(f"Could not apply SDI values: {e}")
                applied["spsdi"] = False

        return applied

    def generate_keywords_for_draw(
        self,
        draw: dict[str, dict[str, float]],
        default_config: dict[str, Any],
        draw_idx: int = 0,
    ) -> str:
        """Generate FVS keywords for a specific posterior draw.

        Produces a keyword block suitable for insertion into any FVS keyfile.
        This is the microfvs/standalone pathway for uncertainty propagation.

        Args:
            draw: Parameter draw from get_draw()
            default_config: Default config dictionary
            draw_idx: Draw index (for comment annotation)

        Returns:
            FVS keyword string
        """
        lines = [
            f"!! Uncertainty draw {draw_idx} for variant {self.variant.upper()}",
            "!!",
        ]

        multipliers = self.compute_multipliers_for_draw(draw, default_config)

        # SDI max
        if "sdi_values" in multipliers:
            lines.append("!! SDI max (posterior draw)")
            for i, val in enumerate(multipliers["sdi_values"]):
                if val > 0:
                    lines.append(f"SDIMAX          {i + 1:10d}{val:10.1f}")

        # Mortality multipliers
        if "mort_mult" in multipliers:
            lines.append("!! Mortality multipliers (posterior draw)")
            for i, mult in enumerate(multipliers["mort_mult"]):
                if abs(mult - 1.0) > 0.01:
                    lines.append(
                        f"MORTMULT        {i + 1:10d}{mult:10.4f}       0.0     999.0"
                    )

        # Growth multipliers
        if "growth_mult" in multipliers:
            lines.append("!! Growth multipliers (posterior draw)")
            for i, mult in enumerate(multipliers["growth_mult"]):
                if abs(mult - 1.0) > 0.01:
                    lines.append(f"BAIMULT         {i + 1:10d}{mult:10.4f}")

        # Height growth multipliers
        if "hg_mult" in multipliers:
            lines.append("!! Height growth multipliers (posterior draw)")
            for i, mult in enumerate(multipliers["hg_mult"]):
                if abs(mult - 1.0) > 0.01:
                    lines.append(f"HTGMULT         {i + 1:10d}{mult:10.4f}")

        return "\n".join(lines)

    # =========================================================================
    # Ensemble Summarization
    # =========================================================================

    @staticmethod
    def summarize_ensemble(
        results: list[pd.DataFrame],
        quantiles: list[float] | None = None,
    ) -> pd.DataFrame:
        """Summarize an ensemble of FVS summary tables into credible intervals.

        Given a list of FVS summary DataFrames (one per posterior draw), computes
        the mean, standard deviation, and quantile bands for every numeric column
        at each projection year.

        Args:
            results: List of FVS summary DataFrames (from fvs.summary)
            quantiles: Quantile levels to compute (default: 2.5%, 50%, 97.5%)

        Returns:
            DataFrame with multi level columns: (variable, statistic)
        """
        if quantiles is None:
            quantiles = [0.025, 0.25, 0.50, 0.75, 0.975]

        if not results:
            raise ValueError("No results to summarize")

        # Stack all results with a draw identifier
        stacked = pd.concat(
            [df.assign(_draw=i) for i, df in enumerate(results)],
            ignore_index=True,
        )

        # Identify the year/cycle column
        year_col = None
        for candidate in ("Year", "year", "YEAR", "IY"):
            if candidate in stacked.columns:
                year_col = candidate
                break

        if year_col is None:
            # Use first column as grouping key
            year_col = stacked.columns[0]

        # Numeric columns to summarize (excluding _draw and year_col)
        numeric_cols = [
            c for c in stacked.select_dtypes(include=[np.number]).columns
            if c not in ("_draw", year_col)
        ]

        summary_pieces = []

        for col in numeric_cols:
            grouped = stacked.groupby(year_col)[col]

            stats = grouped.agg(["mean", "std", "count"]).copy()
            stats.columns = pd.MultiIndex.from_tuples(
                [(col, "mean"), (col, "std"), (col, "n_draws")]
            )

            for q in quantiles:
                q_vals = grouped.quantile(q)
                q_label = f"q{q:.3f}".replace(".", "")
                q_df = q_vals.to_frame()
                q_df.columns = pd.MultiIndex.from_tuples([(col, q_label)])
                stats = stats.join(q_df)

            summary_pieces.append(stats)

        summary = pd.concat(summary_pieces, axis=1)
        summary.index.name = year_col

        return summary

    @staticmethod
    def ensemble_to_long(
        summary: pd.DataFrame,
        variables: list[str] | None = None,
    ) -> pd.DataFrame:
        """Convert ensemble summary to long format for plotting.

        Produces a tidy DataFrame with columns: year, variable, mean,
        lower (q0.025), upper (q0.975), and median (q0.50).

        Args:
            summary: Output from summarize_ensemble()
            variables: Subset of variables to include (None = all)

        Returns:
            Long format DataFrame suitable for ggplot or matplotlib
        """
        year_col = summary.index.name or "Year"
        records = []

        # Get unique top level variable names
        top_vars = summary.columns.get_level_values(0).unique()

        if variables is not None:
            top_vars = [v for v in top_vars if v in variables]

        for var in top_vars:
            var_data = summary[var].copy()
            var_data[year_col] = summary.index

            row = {year_col: var_data[year_col].values}
            row["variable"] = var

            for stat in var_data.columns:
                if stat == year_col:
                    continue
                row[stat] = var_data[stat].values

            # Build records row by row
            n = len(var_data)
            for i in range(n):
                rec = {year_col: var_data[year_col].values[i], "variable": var}
                for stat in var_data.columns:
                    if stat == year_col:
                        continue
                    rec[stat] = var_data[stat].values[i]
                records.append(rec)

        return pd.DataFrame(records)

    # =========================================================================
    # Internal Helpers
    # =========================================================================

    def _extract_species_vector(
        self,
        params: dict[str, float],
        base_name: str,
    ) -> np.ndarray | None:
        """Extract a species indexed parameter vector from a draw.

        Handles both indexed form (b0[1], b0[2], ...) and scalar form.

        Args:
            params: Parameter dictionary from a single draw
            base_name: Base parameter name to search for

        Returns:
            numpy array of species values, or None if not found
        """
        # Try indexed form first: base_name[1], base_name[2], ...
        indexed = {}
        for key, val in params.items():
            if key.startswith(base_name + "["):
                try:
                    idx = int(key.split("[")[1].rstrip("]"))
                    indexed[idx] = float(val) if not isinstance(val, (int, float)) else val
                except (ValueError, IndexError):
                    continue

        if indexed:
            max_idx = max(indexed.keys())
            arr = np.ones(max_idx, dtype=np.float64)
            for idx, val in indexed.items():
                arr[idx - 1] = val  # 1 indexed to 0 indexed
            return arr

        # Try scalar form
        if base_name in params:
            val = params[base_name]
            if isinstance(val, (int, float)):
                return np.array([val], dtype=np.float64)

        return None

    @staticmethod
    def _pad_array(arr: np.ndarray, target_len: int, fill: float = 1.0) -> np.ndarray:
        """Pad or trim array to target length."""
        if len(arr) < target_len:
            return np.pad(arr, (0, target_len - len(arr)), constant_values=fill)
        elif len(arr) > target_len:
            return arr[:target_len]
        return arr
