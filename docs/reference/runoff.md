# Model run-off (October 2026)

The run-off is a pre-registered comparison of three model families on the same FIA stands, scored by the same code. Everything is in [`validation/runoff_20261030/`](https://github.com/advanced-forestry-systems/fvs-modern/tree/main/validation/runoff_20261030) (merged in PR 139): the pre-registration with amendments A1 to A3, the technical report `runoff_report_20261030.md`, the scripts, the Slurm jobs and the result tables. This page summarises the findings that matter to users of fvs-modern; read the report for the methods and the full numbers.

!!! note "The run-off does not change the catalog"
    `fvs2py.catalog.select_config` still reads the configuration catalog described on [Validation](validation.md). The run-off uses a different decision rule (long-term gates first, then the smallest near-term basal area bias, then simplicity among tied arms), so its winner for a variant can differ from the catalog's recommended arm. For NE, for example, the run-off picks `calibrated_native` while the catalog recommends `default_gompit`.

## Design

- **Contestants.** (A) the fvs-modern engine in four arms (default or calibrated parameters with native or GOMPIT mortality); (B) the fvs-conus projector, v2 in four arms and v6 as a diagnostic only, with ingrowth off; (C) unmodified USDA FVS built from the official repository at commit b183c10f, the fork base of fvs-modern. USDA FVS cannot build the Acadian (ACD) variant, so C runs ACD stands on its NE library.
- **Stands.** The 4,748 stands of the catalog sample across 18 variants, each started from its own FIA condition trees.
- **Long term.** 100 years on a 5-year step with no regeneration or management, scored with the Bakuzis stress test and the catalog gates.
- **Near term.** One cycle equal to each plot's own remeasurement interval, scored against observed time-2 survivor basal area, trees per hectare, QMD and volume.

No file in the folder carries a plot identifier or coordinates; the per-stand joined files stay on Cardinal.

## Findings

**fvs-modern default native is USDA FVS.** On every stand and every year of the 100-year runs, fvs-modern with default parameters and native mortality gave the same basal area and stem counts as unmodified USDA FVS in 15 of 18 variants (`results/c_vs_a_default_native_full.csv`). The exceptions are WS (99.98 percent of stand-year rows identical), SN (56.0 percent, maximum difference 2.53 m² ha⁻¹, still to be traced) and ACD (5.6 percent, expected because USDA FVS runs the NE library there). Differences between fvs-modern and USDA FVS therefore come from the calibrated parameters and GOMPIT.

**Every engine arm over-predicts near-term survivor basal area**, by +1.7 percent (TT default GOMPIT) to +51.5 percent (EM default native, and USDA FVS). Part of this is the target: no contestant simulates the disturbance and harvest losses that survivor basal area includes.

**GOMPIT mortality lowers near-term bias** in all 18 variants with default parameters and in 16 of 18 with calibrated parameters (the exceptions are CI and SN).

**fvs-conus v2 overgrows** (+21.2 to +109.7 percent near-term basal area bias) and overshoots the FIA envelope in the eastern variants. **fvs-conus v6 collapses basal area** (its b1 logit arm under-predicts by 20.8 to 67.6 percent) and was never admissible; the report attributes this to a stand basal area cap that appears miscalibrated.

**Decisions** (`results/runoff_decision_by_variant.csv`):

| Outcome | Variants |
|:---|:---|
| fvs-modern arm wins | ACD (`default_native`), BM (`default_gompit`), CS (`calibrated_gompit`), LS (`calibrated_native`), NE (`calibrated_native`), SN (`default_gompit`), WS (`calibrated_native`) |
| fvs-conus v2 arm wins | TT; and CA, EC, EM, SO, where a v2 arm is the only admissible arm |
| No admissible arm | CI, CR, IE, PN, UT, WC |

The report reads the CA, EC, EM and SO results as least bad rather than good, because the only admissible arm there has a near-term bias above 50 percent.

**The October 30 calibration fix changed no admissibility verdict** in the run-off (`results/fix_effect_calibrated_v2.csv`); results from the pre-fix calibrated path are kept in `results/sensitivity/`.

## Limitations

The report lists them: the near-term target is survivor basal area, so every loss to disturbance counts against every model; common stand sets range from 35 stands (PN) to 290; the sample has one site and one density level, so the Bakuzis site and density relations could not be tested; intervals are rounded to whole years; and the USDA FVS result for ACD is an NE substitution.
