"""Check the fixed-effect modifier equation against its Stan design."""
import math
import unittest
from pathlib import Path

from config.config_loader import FvsConfigLoader


TYPES = ("plant", "fire", "insect", "disease", "wind", "harvest", "cutting", "siteprep")


class SpeciesFreeModifierEtaTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.loader = FvsConfigLoader("ne", version="conus_sf",
                                     config_dir=Path(__file__).resolve().parents[1])
        # Exercise the formula in memory only. Current production configs
        # deliberately lack this provenance gate and cannot be evaluated.
        cls.loader.config["categories_conus_sf_modifier_validated"] = True
        cls.blocks = cls.loader.config["categories_conus_sf_modifier"]

    def flags(self, **active):
        return {name: active.get(name, 0) for name in TYPES}

    def test_production_coefficients_are_blocked_until_refit(self):
        loader = FvsConfigLoader("ne", version="conus_sf",
                                 config_dir=Path(__file__).resolve().parents[1])
        with self.assertRaisesRegex(RuntimeError, "event-code mapping"):
            loader.sf_modifier_fixed_eta("diameter_growth", self.flags(), 0, 0)

    def test_common_form_uses_intercept_and_separate_decay_envelopes(self):
        b = self.blocks["height_growth"]
        actual = self.loader.sf_modifier_fixed_eta(
            "height_growth", self.flags(fire=1, cutting=1), 0.5, 0.25)
        expected = b["alpha_0"] + b["alpha"]["fire"] * 0.5 + b["alpha"]["cutting"] * 0.25
        self.assertAlmostEqual(actual, expected, places=12)
        self.assertAlmostEqual(self.loader.sf_modifier_fixed_eta(
            "diameter_growth", self.flags(), 0, 0),
            self.blocks["diameter_growth"]["alpha_0"], places=12)

    def test_trait_effect_is_standardized_and_applies_to_active_type(self):
        b = self.blocks["height_crown_base"]
        scale = b["trait_standardization"]
        traits = dict(zip(scale["trait"], scale["scale_mean"]))
        at_mean = self.loader.sf_modifier_fixed_eta(
            "height_crown_base", self.flags(plant=1), 0, 0, traits)
        self.assertAlmostEqual(at_mean, b["alpha_0"] + b["alpha"]["plant"], places=12)
        traits[scale["trait"][0]] += scale["scale_sd"][0]
        shifted = self.loader.sf_modifier_fixed_eta(
            "height_crown_base", self.flags(plant=1), 0, 0, traits)
        self.assertAlmostEqual(shifted - at_mean,
                               b["trait_gamma_alpha"]["plant"][0], places=12)

    def test_missing_event_or_trait_inputs_fail_closed(self):
        with self.assertRaises(ValueError):
            self.loader.sf_modifier_fixed_eta("diameter_growth", {"plant": 0}, 0, 0)
        with self.assertRaises(ValueError):
            self.loader.sf_modifier_fixed_eta("diameter_growth", self.flags(fire=1), math.nan, 0)
        with self.assertRaises(ValueError):
            self.loader.sf_modifier_fixed_eta("height_crown_base", self.flags(plant=1), 0, 0)

    def test_no_implicit_output_multiplier(self):
        eta = self.loader.sf_modifier_fixed_eta("mortality", self.flags(insect=1), 1, 0)
        b = self.blocks["mortality"]
        self.assertAlmostEqual(eta, b["alpha_0"] + b["alpha"]["insect"], places=12)
        self.assertNotAlmostEqual(eta, math.exp(eta), places=6)

    def test_all_variant_blocks_have_finite_stress_outputs(self):
        root = Path(__file__).resolve().parents[1]
        configs = [p for p in (root / "calibrated").glob("*.json")
                   if ".sf_preview." not in p.name and "_draws" not in p.name]
        self.assertEqual(len(configs), 25)
        checked = 0
        for path in configs:
            loader = FvsConfigLoader(path.stem, version="conus_sf", config_dir=root)
            loader.config["categories_conus_sf_modifier_validated"] = True
            for component, block in loader.config["categories_conus_sf_modifier"].items():
                traits = None
                if block["form"] == "traitmed":
                    scale = block["trait_standardization"]
                    traits = {name: mean + 2 * sd for name, mean, sd in zip(
                        scale["trait"], scale["scale_mean"], scale["scale_sd"])}
                for active, decay in ((0, 0.0), (1, 1.0)):
                    eta = loader.sf_modifier_fixed_eta(
                        component, self.flags(**{name: active for name in TYPES}),
                        decay, decay, traits)
                    self.assertTrue(math.isfinite(eta))
                    checked += 1
        self.assertEqual(checked, 300)


if __name__ == "__main__":
    unittest.main()
