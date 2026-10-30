# Variants

fvs-modern builds the 22 US variants of the USDA Forest Vegetation Simulator, the two Canadian variants (BC and ON) and the Acadian variant (ACD). Each is a separate shared library, `lib/FVS<code>.so`, and `fvs2py` infers the variant from that file name. The `fvs2py.enums.FvsVariant` enumeration lists the same 25 codes.

## Coverage by variant

The table combines five sources: the variant names in `deployment/fvs2py/fvs2py/enums.py`; the species slot count `maxsp` in the default configuration `config/<code>.json` (the engine reports its own count in `FVS.dims["maxspecies"]`; for NE both are 108); whether `base/gompmort.f90` appears in the variant's source list `src-converted/bin/FVS<code>_sourceList.txt`, which decides whether GOMPIT survival is available; the list of variants with non-unity per-species diameter growth multipliers in `KNOWN_ISSUES.md`; and the `status` column of `config/catalog/catalog_recommendation_by_variant.csv`.

| Code | Region (fvs2py `FvsVariant`) | Species slots (`maxsp`) | GOMPIT compiled | Non-unity DG multipliers | Catalog status |
|:---|:---|---:|:---|:---|:---|
| ACD | Acadian (Maine, New Hampshire, Vermont) | 108 | yes | yes | RECOMMENDED |
| AK | Southeast Alaska and Coastal British Columbia | 23 | yes | no | NOT EVALUATED |
| BC | British Columbia | 15 | no | no | not in catalog |
| BM | Blue Mountains | 18 | yes | no | RECOMMENDED |
| CA | Inland California and Southern Cascades (ICASCA) | 50 | yes | yes | NOT RECOMMENDED |
| CI | Central Idaho | 19 | yes | no | NOT RECOMMENDED |
| CR | Central Rockies | 38 | yes | no | NOT RECOMMENDED |
| CS | Central States | 96 | yes | yes | RECOMMENDED |
| EC | East Cascades | 32 | yes | no | NOT RECOMMENDED |
| EM | Eastern Montana | 19 | yes | no | NOT RECOMMENDED |
| IE | Inland Empire | 23 | yes | no | NOT RECOMMENDED |
| KT | Kootenai, Kaniksu, and Tally Lake (KooKanTL) | 11 | yes | yes | not in catalog |
| LS | Lake States | 68 | yes | yes | RECOMMENDED |
| NC | Klamath Mountains (and northern California) | 12 | yes | yes | NOT EVALUATED |
| NE | Northeastern US | 108 | yes | no | RECOMMENDED |
| OC | Organon Southwest | 50 | yes | no | not in catalog |
| ON | Ontario | 72 | yes | yes | not in catalog |
| OP | Organon Pacific Northwest | 39 | yes | no | not in catalog |
| PN | Pacific Northwest Coast | 39 | yes | no | NOT RECOMMENDED |
| SN | Southern US | 90 | yes | no | RECOMMENDED |
| SO | South Central Oregon and Northeast California (SORNEC) | 33 | yes | no | NOT RECOMMENDED |
| TT | Tetons | 18 | yes | no | NOT RECOMMENDED |
| UT | Utah | 24 | yes | no | NOT RECOMMENDED |
| WC | Westside Cascades | 39 | yes | no | NOT RECOMMENDED |
| WS | Western Sierra Nevada | 43 | yes | no | RECOMMENDED |

"Not in catalog" means the variant has no row in the catalog file, so `fvs2py.catalog.select_config` raises `KeyError` for it, as it does for the two NOT EVALUATED variants. [Choosing a configuration](choosing-a-configuration.md) gives the full catalog with the arm recommended for each evaluated variant.

## Calibrated configurations

`CALIBRATION.md` and the repository README describe Bayesian calibrated parameter sets for all 25 variants under `config/calibrated/<code>.json`, with posterior draws in `config/calibrated/<code>_draws.json`. The calibrated components per `CALIBRATION.md` are diameter growth, height-diameter, height increment (six variants), mortality, crown ratio change, SDI maximum and the self-thinning slope. `KNOWN_ISSUES.md` narrows this for the shipped files: the per-species diameter growth multipliers are non-unity in only the seven variants marked in the table, and height increment was adopted for six (`bc`, `ci`, `em`, `ie`, `kt`, `ws`).

!!! warning "Unverified"
    The per-variant calibrated JSON and draw files are not in the source snapshot used for this page. Its `config/calibrated/` directory holds only `<code>.sf_preview.json` files for the same 25 codes, and `config/calibrated_sdifix/` holds an alternative set for AK, CR, CS, LS, NE, SN, TT, UT and WS. Which variants actually have a loadable `config/calibrated/<code>.json` could not be confirmed; the catalog's `calibrated_distinct` column is TRUE for every engine arm, meaning the calibrated configuration loaded for all 18 evaluated variants when the catalog was built.

The Acadian variant does not respond to FVS keyword multipliers (`CALIBRATION.md`: a 5x `MORTMULT` changed 100 year biomass by only about 6 percent), so keyword-based calibration has little effect on ACD.

## Other codes in the repository

- `adk`: a source list `src-converted/bin/FVSadk_sourceList.txt` exists, but `adk` is not in the build script's default list or in `FvsVariant`. Build it explicitly if you need it; it is not covered by this guide.
- `cn`: `config/cn.json` is a CONUS unified configuration ("CONUS unified (Johnson, Marshall, Weiskittel; constrained nlsLM, 2026-04-25)", 115 species), not a variant library.

## Canadian variants

!!! warning "Unverified"
    The README reports all 25 shared libraries building and loading (25/25 under RTLD_LAZY), while `KNOWN_ISSUES.md` (section "PN/SN/IE keyrdr.f90 EOF blocker", resolved May 9, 2026) reports 23 of 25 loading, with BC and ON failing on separate compilation gaps. Which statement describes commit fb4cf10 was not checked. BC is also the only variant built without GOMPIT.

## Variant-specific notes

- **NE** switches FIA volume, biomass and carbon (`FIAVBC`) on by default, so stands without an ecodivision log an `FVS43` warning ([Troubleshooting](troubleshooting.md)).
- **NE** carries the frozen v1 mortality limitations listed in `KNOWN_ISSUES.md` (background mortality halved in source, a density gate that most plots never reach, and a basal area clamp at BAMAX), which show up as a basal area plateau in long runs.
- **ACD** ignores keyword multipliers, as noted above.
- **TT** is NOT RECOMMENDED for long horizons, and the engine arm `select_config("tt")` returns is chosen by a tie-break; see [Choosing a configuration](choosing-a-configuration.md).
