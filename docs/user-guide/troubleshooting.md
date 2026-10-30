# Troubleshooting

This page collects the problems a user is likely to meet, drawn from `KNOWN_ISSUES.md`, `STOP_CODES.md` and the code at commit fb4cf10, and keeps only items the code still shows. Each entry gives the symptom, the cause and what to do.

## Exit codes and messages

The engine keeps a single severity code that only ever rises during a run: each message routine in `src-converted/base/errgro.f90` raises it to its own level if it is lower. `fvs2py` exposes the code as `FVS.exit_code`; a standalone executable converts it to its STOP code in `src-converted/base/main.f90`.

| `exit_code` | Standalone STOP | Meaning | Messages that set it (errgro.f90) |
|:---|:---|:---|:---|
| 0 | 0 | No notes, warnings or errors | |
| 1 | 10 | Notes and warnings only; the run completed | FVS03 forest code outside the model's range, FVS07, FVS08 too few projectable trees, FVS09 plot count mismatch, FVS14 habitat or ecoregion, FVS43 FIAVBC ecodivision |
| 2 | 20 | Non-fatal errors; usually a keyword was ignored | FVS01 invalid keyword, FVS02 no STOP record, FVS04 parameter missing or incorrect, FVS05, FVS06, FVS10, FVS12, FVS13, FVS15, FVS16, FVS17 |
| 3 | 30 | Serious error; output may be partial | FVS11 requested extension not in this library, FVS18 |
| 4 | 40 | Fatal input error | FVS19 incorrect record type or end of file |

To find the cause, search the main output listing (`<keyfile>.out`) for `FVS` followed by two digits and `ERROR` or `WARNING`.

!!! note "Differences from older documentation"
    `STOP_CODES.md` lists a STOP 50 for a user-requested abort; no code at this commit sets severity 5, so it does not occur. It also groups FVS03 with severity 2; in `errgro.f90` FVS03 raises severity 1. The `exit_code` property docstring in `fvs2py` labels the codes "input data error", "keyword or expression error" and so on, which does not match the message assignments above.

An exit code of 2 does not mean the projection failed. `STOP_CODES.md` records that running a Northeast keyword file on CS, LS, SN, KT or EM gives STOP 20 because of the NE forest code, NE species codes and NE plot design, while the projection still completes. Use the variant's own forest codes and species codes.

## FVS43 warning on Northeast stands

**Symptom.** `FVS43 WARNING: FIAVBC REQUESTED BUT ECODIVISION CODE NOT RECOGNIZED OR MISSING. USING "0000" BY DEFAULT`, exit code 1.

**Cause.** NE turns FIA volume, biomass and carbon on by default (`src-converted/ne/grinit.f90`, line 166) and needs an ecodivision, which a keyword-only stand or a database stand without `ECOREGION` does not supply.

**Fix.** Supply `ECOREGION` in `FVS_StandInit`. The warning does not stop growth or mortality; the [Quickstart](../getting-started/quickstart.md) run carries it.

## Species silently replaced

**Symptom.** Exit code 1 and `NOTE: INPUT SPECIES CODE (...) WAS SET TO (...) FOR THIS PROJECTION.` in the listing.

**Cause.** A code the variant does not carry is mapped to a catch-all species (`src-converted/base/intree.f90`, line 421). `tools/fvs_input_guards.py` records five different invalid codes all producing identical NE output.

**Fix.** Check species before running; with `inventory_keyfile.make_inventory_keyfile` an unmapped FIA code raises `KeyError`, and `filter_mapped_species` lists what would be dropped.

## A calibrated run gives default results

**Symptom.** `config_version="calibrated"` (or `"custom"`) produces the same numbers as a default run.

**Causes and fixes.**

1. The configuration could not be loaded. `apply_calibrated_config()` logs a warning (`config_loader not found...`, `Config file not found: ...` or `Could not apply calibrated config: ...`) and returns `{}`, and the run continues with defaults. Put the repository root on `PYTHONPATH`, check that `config/calibrated/<variant>.json` exists, and check the returned dictionary.
2. It was a second keyword file on the same `FVS` object. The overlay is applied automatically only on the first run of an object. Stop at point 7 and call `apply_calibrated_config()` for each keyword file ([Running projections](running-projections.md)).
3. The variant's calibration does not touch the quantity you are looking at. The in-process overlay sets only the attributes it can derive (for NE only `spsdi`), and 18 of 25 variants carry unity diameter growth multipliers ([Choosing a configuration](choosing-a-configuration.md)). With the keyword path, a family in which every multiplier is within 0.01 of 1 produces no records and a message on standard error ending `the arm runs on the uncalibrated growth model. This is a calibration content question, not a keyword formatting failure.`

## Keyword path warnings

`BAMAX skipped for variant <code>: ...` is emitted by `FvsConfigLoader.generate_keywords()` when a variant stores a per-species basal area maximum (the code comment names ACD and IE). `BAMAX` is a stand-level keyword with no species field, so the vector is skipped and the rest of the calibration block is still written.

## GOMPIT does not take effect

**Symptom.** Results with `FVS_GOMPIT=1` match native survival, or the listing shows `GOMPMORT: FVS_GOMPIT set but FVS_GOMPIT_COEF empty; gompit mortality NOT enabled.` or `GOMPMORT: cannot open coeff file ...`.

**Causes and fixes.**

- `FVS_GOMPIT_COEF` is unset or points at a missing file. Use the value from `select_config(...).env` or the absolute path of `config/greg_mortality_coefficients.csv`.
- The variables were set after the first run in the process. The engine resolves GOMPIT once per process; start a new process.
- The variant is BC, which is built without GOMPIT.
- A value such as `0`, `no` or `false` (anything starting with `0`, `n`, `N`, `f` or `F`) disables it.

When GOMPIT is active the listing shows `GOMPMORT enabled: <n> fitted species read, <m> matched to this variant.`; species without a match keep native mortality.

## Crashes (segmentation faults)

**Tree attributes.** `FVS.get_tree_attr`, `FVS.set_tree_attr` and `FVS.tree_table` crash because `fvs2py` passes five arguments to an engine routine that takes six. Do not call them at this commit; write tree lists with the `TREELIST` keyword instead ([Running projections](running-projections.md#tree-lists)).

**Two FVS objects on one library.** A second `FVS` object built on the same library file in the same process shares the first one's engine state, and the reload performed by `load_keyfile()` cannot reset a library that another handle still holds. Use one object per library per process.

!!! warning "Unverified"
    In the documentation build, ordinary runs of the quickstart stand also ended with a segmentation fault inside the engine's main routine in roughly a third to a half of fresh processes, with default, calibrated and custom settings alike; runs that completed gave identical results. That build lacked the SQLite, FOFEM and `apisubsc.c` C sources and the NVEL `.for` files and stubbed 39 unresolved symbols, so the crashes may be an artefact of the incomplete build. Whether a complete build shows them was not checked. If you see intermittent crashes, run each stand in its own process and retry failures.

## Library will not load

- `ImportError: ... are needed routines that are not available in library`: the library was not built by `build_fvs_libraries.sh`, or was built with a compiler that changes symbol names. Only `gfortran` is supported; the script reverts `FC` to `gfortran` and warns otherwise.
- `OSError` from ctypes about undefined symbols: open the library lazily (`fvs2py` does this). The NVEL symbols are intentionally left unresolved.
- `OSError` naming `libfvs_stubs.so` or `libfvs_stubs_final.so`: keep the stub libraries in the same directory as `FVS<code>.so`; on Linux the variant library finds them through an `$ORIGIN` run path.
- A library that loads but reads zero trees: rebuild with the current script, which shadows the variant parameter header and fails the build with `[CONCHR] ASSERTION FAILED` when the `/CONCHR/` common block sizes disagree.
- Native Windows: not supported for building; use WSL2 or Docker (README).

## Catalog errors

- `KeyError: 'variant AK: NOT EVALUATED (no FIA thresholds for this variant)'` (also NC): the variant was not gated.
- `KeyError: 'variant KT is not in the catalog (no fitting plots assigned or not run)'` (also BC, OC, ON, OP).
- `ValueError: ... is a CONUS projector arm; it is not runnable through fvs2py.`: a forced `survival` or `parameter_set` led to a projector arm. Drop the constraint or choose an engine arm.
- `FileNotFoundError` for the catalog CSVs: `fvs2py` is not running from a repository checkout. Set `FVS_CATALOG_DIR` or pass `catalog_dir`.

## Long-horizon behaviour that is not a bug in your inputs

`KNOWN_ISSUES.md` lists defects of the Northeast v1 release that were diagnosed and deliberately left in place. They change long projections in recognisable ways:

- **Basal area plateau.** NE inflates mortality until residual basal area is within 1 ft² ac⁻¹ of BAMAX (`src-converted/vls/morts.f90`, the loop that starts at line 805; `KNOWN_ISSUES.md` cites lines 750 to 786 of an earlier revision). Long runs flatten at that ceiling, as the quickstart stand does at 242 ft² ac⁻¹. Do not read a stocking-guide pass as evidence of fidelity.
- **Background mortality** is a diameter-only logistic halved in source (`vls/morts.f90`, lines 561 and 567 at this commit), and the density-dependent mortality gate (`PMSDIL = 55`) is crossed by only about 22 percent of the benchmark cohort. The `FVS_MORT_HALVE` and `FVS_MORT_PMSDIL` variables expose these constants but no setting has been adopted.
- **Reineke slopes** of the native arms fall outside the accepted band; the configuration catalog's Reineke gate is the reason most western variants are NOT RECOMMENDED for long horizons.

`KNOWN_ISSUES.md` also notes, outside the NE release table, that the crown ratio change model has the weakest fit of the calibrated components (R² 0.13 to 0.38 across variants) and that ingrowth uses observed FIA recruitment rates rather than a regeneration model, which it describes as adequate for 5 to 15 year projections but not for long-term regeneration dynamics.

## Results made before October 30, 2026

PR 135 (merged October 30, 2026) fixed four defects that change results produced earlier. Regenerate affected outputs:

- posterior draw keyword blocks from `UncertaintyEngine.generate_keywords_for_draw()` used a 16-column layout that FVS misread ([Uncertainty](uncertainty.md#the-keyword-path-for-draws));
- the in-process calibrated overlay wrote an SDI maximum of 0 for species without a calibrated value ([Running projections](running-projections.md#what-the-in-process-overlay-changes));
- the FIA stand generator gave a condition the trees of the whole plot without `CONDPROP_UNADJ` expansion ([Stand inputs](stand-inputs.md#condition-proportion));
- posterior draws were converted to multipliers on a different basis from the MAP configuration ([Uncertainty](uncertainty.md#running-an-ensemble)).

## Issues listed in KNOWN_ISSUES.md that no longer apply

- The six-column offset in the `BAMAX`, `BAIMULT` and `HTGMULT` writers of `config/config_loader.py` (Table 1 of `KNOWN_ISSUES.md`) has been fixed in `FvsConfigLoader`: all four calibration keywords now left-justify the keyword in 10 columns, and the scheduled multipliers carry the date field. Since PR 135 the draw writer in `config/uncertainty.py` uses the same formatters.
- The PN, SN and IE keyword reader EOF blocker and the ACD build problem on Cardinal are recorded as resolved.
