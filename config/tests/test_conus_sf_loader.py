"""
Tests for the CONUS species-free (Leg B) path in config_loader.py.

Builds a synthetic variant JSON whose categories_conus_sf block matches exactly
what calibration/R/63_conus_sf_to_variant_json.R writes, then checks that the
loader decodes it, computes the trait effect for in- and out-of-bundle species,
evaluates the species-free linear predictor, and resolves the Leg A / Leg B
source correctly under version='conus_hybrid'.

Runnable two ways:
    pytest config/tests/test_conus_sf_loader.py
    python  config/tests/test_conus_sf_loader.py     # prints PASS/FAIL
"""
from __future__ import annotations

import importlib.util
import json
import math
from pathlib import Path

_HERE = Path(__file__).resolve()
_CONFIG_DIR = _HERE.parent.parent            # .../config
_LOADER_PATH = _CONFIG_DIR / "config_loader.py"

# Import config_loader by path (no package assumptions).
_spec = importlib.util.spec_from_file_location("config_loader", _LOADER_PATH)
config_loader = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(config_loader)
FvsConfigLoader = config_loader.FvsConfigLoader


# ---------------------------------------------------------------------------
# Synthetic fixture: two traits, two species in-bundle, small L1/L2 RE tables.
# ---------------------------------------------------------------------------

# Trait standardization constants (as 63_...R stores them).
SCALE_MEAN = [0.50, 3.0]          # wood_specific_gravity, shade_tolerance_num
SCALE_SD = [0.10, 1.0]
GAMMA = [0.20, -0.05]             # trait coefficients on standardized traits

# Two in-bundle species with known raw traits.
SP = [12, 316]
RAW_WSG = [0.40, 0.60]            # -> std: (0.40-0.50)/0.10 = -1.0 ; +1.0
RAW_SHT = [2.0, 5.0]             # -> std: (2-3)/1 = -1.0 ; (5-3)/1 = +2.0
STD_WSG = [(RAW_WSG[i] - SCALE_MEAN[0]) / SCALE_SD[0] for i in range(2)]
STD_SHT = [(RAW_SHT[i] - SCALE_MEAN[1]) / SCALE_SD[1] for i in range(2)]
# Stored per-species trait effect = std_traits . gamma
TRAIT_EFFECT = [STD_WSG[i] * GAMMA[0] + STD_SHT[i] * GAMMA[1] for i in range(2)]

B0 = 0.75                         # intercept
B1 = 0.30                         # ln(dbh) coefficient (a covariate)
RE_L1 = {"level": ["8", "9"], "mean": [0.10, -0.10]}
RE_L2 = {"level": ["8.1", "9.4"], "mean": [0.02, -0.03]}


def _build_sf_component() -> dict:
    return {
        "model": "dg_kuehne_cspi_traits1_b1",
        "fixed_effects": {
            "param": ["b0", "b1", "sigma_L1", "sigma"],
            "mean": [B0, B1, 0.4, 0.9],
            "sd": [0.01, 0.01, 0.01, 0.01],
        },
        "trait_gamma": {
            "trait_col": ["wood_specific_gravity", "shade_tolerance_num"],
            "gamma_mean": GAMMA,
            "scale_mean": SCALE_MEAN,
            "scale_sd": SCALE_SD,
        },
        "species": {
            "SPCD": SP,
            "trait_effect_mean": TRAIT_EFFECT,
            "raw_wood_specific_gravity": RAW_WSG,
            "std_wood_specific_gravity": STD_WSG,
            "raw_shade_tolerance_num": RAW_SHT,
            "std_shade_tolerance_num": STD_SHT,
        },
        "re_L1": RE_L1,
        "re_L2": RE_L2,
        "re_L3": {"level": [], "mean": []},   # jsonlite empty -> [] handled by loader
        "re_FT": {"level": [], "mean": []},
        "hybrid_source_map": {"SPCD": [12, 316], "source": ["leg_a", "leg_b"]},
    }


def _write_fixture(tmp_dir: Path) -> None:
    calibrated = tmp_dir / "calibrated"
    calibrated.mkdir(parents=True, exist_ok=True)

    cfg = {
        "variant": "tt",
        "variant_name": "Synthetic Test",
        "maxsp": 3,
        "categories": {"bark_ratio": {"BKRAT": [0.9, 0.9, 0.9]}},
        # Leg A present for SPCD 12 only, so the hybrid source map should
        # resolve 12 -> leg_a and 316 -> leg_b.
        "categories_conus": {
            "diameter_growth": {
                "model": "dg_kuehne_cspi_traits1",
                "species_intercepts": {"SPCD": [12], "mean": [0.05], "sd": [0.01]},
            }
        },
        "categories_conus_sf": {"diameter_growth": _build_sf_component()},
    }
    # Both the default and the calibrated copy exist (loader reads calibrated/).
    (calibrated / "tt.json").write_text(json.dumps(cfg, indent=2))
    (tmp_dir / "tt.json").write_text(json.dumps({"variant": "tt", "maxsp": 3,
                                                 "categories": cfg["categories"]}, indent=2))


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------

def _loader(tmp_dir: Path, version: str) -> "FvsConfigLoader":
    return FvsConfigLoader("tt", version=version, config_dir=tmp_dir)


def test_sf_block_present_and_decoded(tmp_path):
    _write_fixture(tmp_path)
    ld = _loader(tmp_path, "conus_sf")
    assert ld.has_conus_sf_block
    assert "diameter_growth" in ld.conus_sf_components_present()

    rt = ld.get_conus_sf_runtime_block("diameter_growth")
    assert rt["intercept_name"] == "b0"
    assert math.isclose(rt["intercept"], B0)
    # b1 is a covariate; sigma* excluded.
    assert set(rt["covariate"]) == {"b1"}
    assert math.isclose(rt["covariate"]["b1"], B1)
    assert rt["gamma"]["wood_specific_gravity"] == GAMMA[0]
    assert rt["scale"]["shade_tolerance_num"] == (SCALE_MEAN[1], SCALE_SD[1])
    assert rt["re_L1"]["8"] == 0.10 and rt["re_L2"]["9.4"] == -0.03


def test_trait_effect_in_and_out_of_bundle(tmp_path):
    _write_fixture(tmp_path)
    rt = _loader(tmp_path, "conus_sf").get_conus_sf_runtime_block("diameter_growth")

    # In-bundle: returns the stored value.
    assert math.isclose(FvsConfigLoader.sf_trait_effect(rt, 12), TRAIT_EFFECT[0])

    # Out-of-bundle: supply raw traits, expect standardize-then-dot-gamma.
    raw = {999: {"wood_specific_gravity": 0.70, "shade_tolerance_num": 4.0}}
    expected = ((0.70 - SCALE_MEAN[0]) / SCALE_SD[0]) * GAMMA[0] \
        + ((4.0 - SCALE_MEAN[1]) / SCALE_SD[1]) * GAMMA[1]
    got = FvsConfigLoader.sf_trait_effect(rt, 999, raw_traits=raw)
    assert math.isclose(got, expected)

    # Unknown species with no traits -> 0.0 (engine then falls back).
    assert FvsConfigLoader.sf_trait_effect(rt, 8888) == 0.0


def test_sf_linear_predictor(tmp_path):
    _write_fixture(tmp_path)
    ld = _loader(tmp_path, "conus_sf")
    eta = ld.sf_linear_predictor(
        "diameter_growth", spcd=12,
        eco_codes={"L1": "8", "L2": "8.1"},
        covariates={"b1": 2.0},
    )
    expected = B0 + TRAIT_EFFECT[0] + RE_L1["mean"][0] + RE_L2["mean"][0] + B1 * 2.0
    assert math.isclose(eta, expected)


def test_hybrid_source_resolution(tmp_path):
    _write_fixture(tmp_path)
    ld = _loader(tmp_path, "conus_hybrid")
    assert ld.resolve_species_source("diameter_growth", 12) == "leg_a"
    assert ld.resolve_species_source("diameter_growth", 316) == "leg_b"
    # Species not in the map at all -> leg_b.
    assert ld.resolve_species_source("diameter_growth", 4040) == "leg_b"


# ---------------------------------------------------------------------------
# Plain-python runner (no pytest dependency)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    import tempfile
    import traceback

    tests = [
        test_sf_block_present_and_decoded,
        test_trait_effect_in_and_out_of_bundle,
        test_sf_linear_predictor,
        test_hybrid_source_resolution,
    ]
    passed = 0
    for t in tests:
        with tempfile.TemporaryDirectory() as d:
            try:
                t(Path(d))
                print(f"PASS  {t.__name__}")
                passed += 1
            except Exception:
                print(f"FAIL  {t.__name__}")
                traceback.print_exc()
    print(f"\n{passed}/{len(tests)} tests passed")
    raise SystemExit(0 if passed == len(tests) else 1)
