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
        calibration_dir: Optional[str | Path] = None,
    ):
        """Initialize the uncertainty engine.

        Args:
            variant: FVS variant code (e.g., 'ne', 'ca', 'sn')
            config_dir: Path to config directory containing calibrated/ subdir
            draws_path: Direct path to a draws JSON file (overrides config_dir)
            seed: Random seed for reproducible draw sequences
            calibration_dir: calibration/ root holding output/variants/<v>/ and
                data/equation_availability_full.csv (default: the sibling of
                config_dir), used to put draws on the MAP multiplier basis
        """
        self.variant = variant.lower()
        self._config_dir = (Path(config_dir) if config_dir is not None
                            else Path(__file__).parent)
        self._calibration_dir_override = calibration_dir
        self._cal_config: dict[str, Any] | None = None
        self._ln_sdi_ref: np.ndarray | None = None
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

    def _calibration_dir(self) -> Path:
        """calibration/ root next to the config directory (multipliers.R inputs)."""
        if self._calibration_dir_override is not None:
            return Path(self._calibration_dir_override)
        return self._config_dir.parent / "calibration"

    def _calibrated_config(self) -> dict[str, Any]:
        """The calibrated (MAP) config: SDICON, FIAJSP, BAMAX, maxsp."""
        if self._cal_config is None:
            from config.config_loader import FvsConfigLoader
            self._cal_config = FvsConfigLoader(
                self.variant, version="calibrated", config_dir=self._config_dir,
            ).config
        return self._cal_config

    @staticmethod
    def _ln_sdimax_draw(sd_params: dict[str, Any], fia: list, maxsp: int) -> np.ndarray:
        """Per-slot ln SDImax of one stand_density draw at the centred qmd:
        Intercept + r_SPCD[code,Intercept]; population Intercept where the
        species has no random effect."""
        from config.calibration_multipliers import spcd_effects
        icpt = sd_params.get("Intercept")
        out = np.full(maxsp, np.nan)
        if icpt is None:
            return out
        eff, other = spcd_effects(sd_params)
        for i in range(min(maxsp, len(fia))):
            code = fia[i]
            if code is not None and code in eff:
                out[i] = icpt + eff[code]
            elif np.isfinite(other):
                out[i] = icpt + other
            else:
                out[i] = icpt
        return out

    def _ln_sdimax_reference(self, fia: list, maxsp: int) -> np.ndarray:
        """Posterior median of the per-slot ln SDImax over every stored draw."""
        if self._ln_sdi_ref is None:
            if self._draws_data is None:
                self._load_draws()
            sd = self._draws_data.get("draws", {}).get("stand_density", [])
            if not sd:
                self._ln_sdi_ref = np.full(maxsp, np.nan)
            else:
                mat = np.vstack([
                    self._ln_sdimax_draw(
                        {k: (v[0] if isinstance(v, list) and len(v) == 1 else v)
                         for k, v in d.items()}, fia, maxsp)
                    for d in sd[: self._n_draws]
                ])
                self._ln_sdi_ref = np.nanmedian(mat, axis=0)
        return self._ln_sdi_ref

    def compute_multipliers_for_draw(
        self,
        draw: dict[str, dict[str, float]],
        default_config: Optional[dict[str, Any]] = None,
    ) -> dict[str, np.ndarray]:
        """Per-species FVS adjustments for one posterior draw, on the MAP basis.

        BAIMULT, MORTMULT and HTGMULT come from
        config.calibration_multipliers.multipliers_from_params, the Python port
        of calibration/R/multipliers.R that also reproduces the shipped MAP
        calibration_multipliers block from the posterior medians, so a draw and
        MAP use one function and one basis (D4, 2026-10-30). BAIMULT is the
        diameter-growth scale sqrt(dds), exactly as generate_keywords applies to
        the MAP dds_multiplier.

        SDImax: the MAP SDICON (or R5SDI/R4SDI/FMSDI/SDIDEF) value is a blend
        of quantile-regression and Bayesian estimates, so it cannot be rebuilt
        from one draw. Each draw instead perturbs MAP multiplicatively by the
        posterior deviation of that species' ln SDImax from its posterior
        median, sdi = MAP * exp(ln_sdi_draw - median(ln_sdi)), which centres the
        draws on MAP and keeps the posterior spread. Species whose MAP value is
        missing (the NA entries) are NaN here and are never written, so FVS keeps
        its own default for them, as it does under the MAP keyword block.

        default_config is accepted for backward compatibility and not used.

        Returns:
            Dict with growth_mult (BAIMULT), mort_mult, hg_mult and, when the
            calibrated config carries an SDI vector, sdi_values (NaN = skip).
        """
        from config.calibration_multipliers import (
            baimult_from_dds, config_maxsp, fia_codes, multipliers_from_params,
        )
        from config.config_loader import FvsConfigLoader

        cal = self._calibrated_config()
        fia = fia_codes(cal)
        maxsp = config_maxsp(cal)
        m = multipliers_from_params(draw, cal, self.variant, self._calibration_dir())
        result = {
            "growth_mult": baimult_from_dds(m["dds_multiplier"]),
            "mort_mult": m["mort_multiplier"],
            "hg_mult": m["htg_multiplier"],
        }

        map_sdi = FvsConfigLoader._find_sdi_param(cal.get("categories", {}))
        if map_sdi is not None and "stand_density" in draw:
            base = np.array(
                [np.nan if (v is None or isinstance(v, str)) else float(v) for v in map_sdi],
                dtype=np.float64,
            )
            base = np.where(base > 0, base, np.nan)
            n = min(base.size, maxsp)
            ln_draw = self._ln_sdimax_draw(draw["stand_density"], fia, maxsp)[:n]
            ln_ref = self._ln_sdimax_reference(fia, maxsp)[:n]
            dev = ln_draw - ln_ref
            dev = np.where(np.isfinite(dev), dev, 0.0)
            result["sdi_values"] = base[:n] * np.exp(dev)
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
                # Never write 0 (or NaN) for a species without a value: keep
                # whatever FVS already holds for it (D2, 2026-10-30).
                cur = np.asarray(fvs_instance.get_species_attr("spsdi"),
                                 dtype=np.float64).copy()
                vals = np.asarray(multipliers["sdi_values"], dtype=np.float64)
                n = min(cur.size, vals.size)
                ok = np.isfinite(vals[:n]) & (vals[:n] > 0)
                cur[:n][ok] = vals[:n][ok]
                fvs_instance.set_species_attr("spsdi", cur)
                applied["spsdi"] = True
            except Exception as e:
                logger.warning(f"Could not apply SDI values: {e}")
                applied["spsdi"] = False

        return applied

    def generate_keywords_for_draw(
        self,
        draw: dict[str, dict[str, float]],
        default_config: Optional[dict[str, Any]] = None,
        draw_idx: int = 0,
        include_comments: bool = True,
    ) -> str:
        """Generate FVS keywords for a specific posterior draw.

        Records are written by the same formatters FvsConfigLoader.
        generate_keywords uses for the MAP block (10-column keyword field, then
        10-column numeric fields; SDIMAX species in field 1 and value in field
        2; MORTMULT, BAIMULT and HTGMULT as scheduled activities with a blank
        date in field 1). Before 2026-10-30 this method padded every keyword to
        16 columns, so FVS read field 1 of SDIMAX as blank (all species), each
        species index landed in the SDI field, and the last record (108) won for
        every species in every draw (D1).

        Stand-level keywords that are not posterior quantities (BAMAX) are
        copied from the MAP block so that a draw differs from MAP only in its
        drawn parameters.

        Args:
            draw: Parameter draw from get_draw()
            default_config: accepted for backward compatibility, not used
            draw_idx: Draw index (for comment annotation)
            include_comments: write '!!' comment lines

        Returns:
            FVS keyword string
        """
        import warnings
        from config.config_loader import FvsConfigLoader

        c = include_comments
        lines = []
        if c:
            lines += [f"!! Uncertainty draw {draw_idx} for variant {self.variant.upper()}", "!!"]

        mult = self.compute_multipliers_for_draw(draw)
        cal = self._calibrated_config()

        sdi = mult.get("sdi_values")
        emit_sdimax = bool(cal.get("_emit_sdimax", True))
        if sdi is not None and emit_sdimax:
            vals = [None if not np.isfinite(v) else float(v) for v in sdi]
            lines.append(FvsConfigLoader._format_sdimax_keywords(vals, c))

        bamax = FvsConfigLoader._find_bamax_param(cal.get("categories", {}))
        if bamax is not None:
            try:
                lines.append(FvsConfigLoader._format_bamax_keywords(bamax, c))
            except ValueError as exc:
                warnings.warn(f"BAMAX skipped for variant {self.variant}: {exc}")

        lines.append(FvsConfigLoader._format_mortmult_keywords(mult["mort_mult"], c))
        lines.append(FvsConfigLoader._format_baimult_keywords(mult["growth_mult"], c))
        if np.any(np.abs(mult["hg_mult"] - 1.0) > 0.01):
            lines.append(FvsConfigLoader._format_htgmult_keywords(mult["hg_mult"], c))

        return "\n".join(x for x in lines if x)


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
