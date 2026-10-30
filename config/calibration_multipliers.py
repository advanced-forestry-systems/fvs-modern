"""
Per-species calibration multipliers on the MAP serializer basis.

Python port of calibration/R/multipliers.R (compute_calibration_multipliers)
for the three families the FVS keyword emitter consumes:

  mort_multiplier  mortality-rate ratio from the survival logit including the
                   species random effect:
                   (1 - plogis(Intercept + r_SPCD[s])) / (1 - plogis(Intercept))
  dds_multiplier   diameter-growth intercept shift exp(b0[i] - mean(b0)) on the
                   dense Stan species index, mapped to FVS slots through
                   <variant>/diameter_growth_species_index.csv when present,
                   else placed in slot i (the R serializer's approximate path)
  htg_multiplier   same as dds, for height increment

Each family is gated by calibration/data/equation_availability_full.csv and
clipped to [0.1, 10], exactly as in the R serializer.

The same function is applied to (a) the posterior medians in the summary CSVs,
which reproduces the shipped calibration_multipliers block (the MAP keyword
set), and (b) a single posterior draw, which gives that draw's multipliers.
Using one function for both is what keeps the posterior draws centred on MAP.
Before 2026-10-30 the draws went through UncertaintyEngine.
compute_multipliers_for_draw, which compared the calibrated b0 to the default
FVS B1 coefficient (a different basis), emitted no MORTMULT, and so produced
draw BAIMULT of 1.11 to 2.30 on all 108 NE species against MAP 0.86 to 1.15
on 58 species.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Mapping, Optional

import numpy as np
import pandas as pd

CLIP_LO = 0.1
CLIP_HI = 10.0

_RE_SPCD = re.compile(r'^"?r_SPCD\[([^,\]]+),Intercept\]')
_RE_B0 = re.compile(r"^b0\[(\d+)\]$")


def _plogis(x):
    return 1.0 / (1.0 + np.exp(-np.asarray(x, dtype=np.float64)))


def fia_codes(config: Mapping[str, Any]) -> list[Optional[int]]:
    """FIA SPCD per FVS slot (index 0 = FVS species 1), None where blank."""
    raw = (config.get("categories", {}).get("species_definitions", {})
           .get("FIAJSP", []))
    out: list[Optional[int]] = []
    for x in raw:
        try:
            out.append(int(str(x).strip()))
        except (TypeError, ValueError):
            out.append(None)
    return out


def config_maxsp(config: Mapping[str, Any]) -> int:
    m = config.get("maxsp")
    try:
        m = int(m)
    except (TypeError, ValueError):
        m = 0
    return m if m > 0 else max(len(fia_codes(config)), 1)


def availability(calibration_dir: Optional[Path], variant: str) -> dict[str, bool]:
    """Row of equation_availability_full.csv; fail open, as the R code does."""
    default = dict(HD=True, MORT=True, CR=True, DG=True, SDI=True, HI=True)
    if calibration_dir is None:
        return default
    f = Path(calibration_dir) / "data" / "equation_availability_full.csv"
    if not f.exists():
        return default
    tab = pd.read_csv(f)
    row = tab[tab["variant"].astype(str).str.upper() == variant.upper()]
    if row.empty:
        return default
    r = row.iloc[0]
    return {k: str(r.get(k, "TRUE")).strip().upper() == "TRUE" for k in default}


def read_crosswalk(path: Optional[Path]) -> Optional[pd.DataFrame]:
    if path is None or not Path(path).exists():
        return None
    ix = pd.read_csv(path)
    if not {"species_idx", "SPCD"}.issubset(ix.columns):
        return None
    return ix


def spcd_effects(params: Mapping[str, Any]) -> tuple[dict[int, float], float]:
    """r_SPCD[<SPCD>,Intercept] -> {SPCD: value}, plus the OTHER group (NaN if absent)."""
    vals: dict[int, float] = {}
    other = np.nan
    for k, v in params.items():
        m = _RE_SPCD.match(str(k))
        if not m:
            continue
        key = m.group(1).strip()
        if key.upper() == "OTHER":
            other = float(v)
            continue
        try:
            vals[int(key)] = float(v)
        except ValueError:
            continue
    return vals, other


def dense_b0(params: Mapping[str, Any]) -> np.ndarray:
    """b0[i] dense-index intercepts as a vector (NaN where absent)."""
    idx = {}
    for k, v in params.items():
        m = _RE_B0.match(str(k))
        if m:
            idx[int(m.group(1))] = float(v)
    if not idx:
        return np.array([], dtype=np.float64)
    out = np.full(max(idx), np.nan)
    for i, v in idx.items():
        out[i - 1] = v
    return out


def mort_multiplier(intercept: float, effects: dict[int, float], other: float,
                    fia: list, maxsp: int) -> np.ndarray:
    arr = np.ones(maxsp)
    if intercept is None or not np.isfinite(intercept):
        return arr
    base = 1.0 - _plogis(intercept)
    for i in range(min(maxsp, len(fia))):
        code = fia[i]
        if code is not None and code in effects:
            arr[i] = (1.0 - _plogis(intercept + effects[code])) / base
        elif np.isfinite(other):
            arr[i] = (1.0 - _plogis(intercept + other)) / base
    return np.clip(arr, CLIP_LO, CLIP_HI)


def dense_factor(b0: np.ndarray, fia: list, maxsp: int,
                 xwalk: Optional[pd.DataFrame]) -> np.ndarray:
    arr = np.ones(maxsp)
    if b0.size == 0 or not np.any(np.isfinite(b0)):
        return arr
    m = np.exp(b0 - np.nanmean(b0))
    m[~np.isfinite(m)] = 1.0
    if xwalk is not None:
        for di, sp in zip(xwalk["species_idx"], xwalk["SPCD"]):
            try:
                di = int(di); sp = int(sp)
            except (TypeError, ValueError):
                continue
            if 1 <= di <= m.size and sp in fia:
                slot = fia.index(sp)
                if slot < maxsp:
                    arr[slot] = m[di - 1]
    else:
        n = min(m.size, maxsp)
        arr[:n] = m[:n]
    return np.clip(arr, CLIP_LO, CLIP_HI)


def multipliers_from_params(
    params: Mapping[str, Mapping[str, Any]],
    config: Mapping[str, Any],
    variant: str,
    calibration_dir: Optional[Path] = None,
) -> dict[str, np.ndarray]:
    """MAP-serializer multipliers from one parameter set.

    params maps component name ('mortality', 'diameter_growth',
    'height_increment') to a flat {parameter_name: value} dict, which is the
    shape of both a posterior draw (UncertaintyEngine.get_draw) and the
    medians read by params_from_summaries().
    """
    fia = fia_codes(config)
    maxsp = config_maxsp(config)
    avail = availability(calibration_dir, variant)
    vdir = (Path(calibration_dir) / "output" / "variants" / variant.lower()
            if calibration_dir is not None else None)

    mort = np.ones(maxsp)
    mp = params.get("mortality")
    if avail["MORT"] and mp:
        icpt = mp.get("Intercept")
        eff, oth = spcd_effects(mp)
        if icpt is not None:
            mort = mort_multiplier(float(icpt), eff, oth, fia, maxsp)

    dds = np.ones(maxsp)
    dg = params.get("diameter_growth")
    if avail["DG"] and dg:
        xw = read_crosswalk(vdir / "diameter_growth_species_index.csv") if vdir else None
        dds = dense_factor(dense_b0(dg), fia, maxsp, xw)

    htg = np.ones(maxsp)
    hi = params.get("height_increment")
    if avail["HI"] and hi:
        xw = read_crosswalk(vdir / "height_increment_species_index.csv") if vdir else None
        htg = dense_factor(dense_b0(hi), fia, maxsp, xw)

    return {"mort_multiplier": mort, "dds_multiplier": dds, "htg_multiplier": htg}


def _median_col(df: pd.DataFrame) -> pd.Series:
    for nm in ("p50", "median", "Estimate", "estimate", "mean"):
        if nm in df.columns:
            return df[nm]
    return df.iloc[:, 1]


def params_from_summaries(variant_output_dir: Path) -> dict[str, dict[str, float]]:
    """Posterior medians from the summary CSVs, keyed like a draw (MAP inputs)."""
    d = Path(variant_output_dir)

    def load(*names):
        out: dict[str, float] = {}
        for n in names:
            f = d / n
            if f.exists():
                df = pd.read_csv(f)
                for k, v in zip(df.iloc[:, 0].astype(str), _median_col(df)):
                    out.setdefault(k.strip('"'), float(v))
        return out

    return {
        "mortality": load("mortality_summary.csv", "mortality_posterior.csv"),
        "diameter_growth": load("diameter_growth_summary.csv"),
        "height_increment": load("height_increment_summary.csv"),
    }


def baimult_from_dds(dds: np.ndarray) -> np.ndarray:
    """DDS-scale factor -> diameter-growth BAIMULT, as config_loader.generate_keywords does."""
    return np.sqrt(np.clip(np.asarray(dds, dtype=np.float64), 0.01, 100.0))
