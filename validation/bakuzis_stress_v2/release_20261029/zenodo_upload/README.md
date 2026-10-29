# A configuration catalog for fvs-modern across the FVS variants of the conterminous United States

Aaron R. Weiskittel, Center for Research on Sustainable Forests, University of Maine. October 29, 2026. Cite the concept DOI https://doi.org/10.5281/zenodo.23777723 (this version: 10.5281/zenodo.23777724).

## What this record holds

The fvs-modern configuration catalog answers one question per Forest Vegetation Simulator (FVS) variant: which survival model (the variant's native mortality or GOMPIT crown closure survival) and which parameter set (fvs-modern default or Bayesian calibrated maximum a posteriori) stays inside the FIA fitting envelope and obeys Reineke self-thinning over 50 and 100 year projections. For each of 18 variants, up to 300 FIA plots were projected 100 years at a 5 year step under all four arms (Cardinal arrays 15976824 and 15984133), gated by a rule pre-registered before any run was scored and repaired after a blind red team, and then stress tested against 144 alternative rules and 500 stand bootstraps (job 15990793). Seven variants have at least one admissible fvs-modern arm (ACD, BM, CS, LS, NE, SN, WS), eleven do not (CA, CI, CR, EC, EM, IE, PN, SO, TT, UT, WC), and no verdict changes under the stress test; eight admissible and eight not admissible arms are labelled fragile.

Lead file: `00_fvs-modern_catalog_graphical_abstract.png` (three panels: variant verdicts with bootstrap and grid stability, the 72 engine arms in the plane of the two deciding gates with Wilson 95 percent intervals, and per arm verdict stability). The Wilson intervals are the only sampling uncertainty shown; the bootstrap probabilities are themselves the uncertainty summary of the verdicts and carry no second interval; the rule grid is a complete enumeration and has none.

## Files

- `00_fvs-modern_catalog_graphical_abstract.png`: the graphical abstract (300 dpi).
- `2026-10-29_fvs-modern-catalog-report.pdf` and `.docx`, with the Markdown source: the technical report (method, rule, red team repairs, results, stress test, selector, limitations).
- `2026-10-29_fvs-modern-catalog-deck.pdf` and `.pptx`: the presentation (14 slides).
- `data/catalog_recommendation_by_variant.csv`, `data/catalog_gates_by_arm.csv`: the catalog as shipped in fvs-modern `config/catalog/` (primary rule verdicts in `recommended`, `status`, `admissible`; alternative rule in `admissible_alt`, `rule_sensitive`). Produced by `code/catalog_gates.R`.
- `data/stress_by_arm.csv`, `data/stress_by_variant.csv`: the stress test. Produced by `code/stress_catalog.R`.
- `data/scorecard_cells.csv`, `data/scorecard_index.md`: every Bakuzis relation cell for every variant and arm, and the index. Produced by `code/bakuzis_scorecards.R`.
- `data/fia_variant_thresholds_plot_acd.csv`: the FIA plot level thresholds by variant. A compiled reference table: plot level percentiles of the fvs-conus fitting data (`conus_remeasurement_pairs_metric_cond_v2_cspiv6.rds`, not deposited, see below) with ACD set equal to NE; the script that computed it from the restricted pairs file is not deposited.
- `data/2026-10-29_bakuzis-assessment-memo_sections-1-16.md`: the running memo of the whole Bakuzis assessment (Sections 1 to 16); Section 15 is the catalog, Section 16 the stress test.
- `data/catalog_README.md`, `data/catalog_implementation-notes.md`: method documentation and the thirteen numbered deviations from the pre-registered plan plus the shipping note.
- `data_dictionary.csv`: every column of every CSV.
- `code/conus_catalog_runs.py`: the engine runner (stand sample, CONDPROP expansion, species relabelling, NOAUTOES, checkpoints). `code/bakuzis_stress_v2.R`: the Bakuzis scorer. `code/catalog_gates.R`: the gates. `code/prep_scorecard_long.R`, `code/bakuzis_scorecards.R`: per variant scorecards. `code/stress_catalog.R`: the stress test. `code/catalog.py`: the selector (`fvs2py/catalog.py` in the repository). `code/release_figures.py`, `code/build_report_docx.py`, `code/build_deck_pptx.py`: this record's figures, report and deck.
- `fig/`: the report and deck figures, including the self-thinning limb by variant (default native arm) and the full Bakuzis matrix for NE default GOMPIT.
- `check_deposit.py`, `files_to_upload.txt`, `zenodo_metadata.json`, `RUNBOOK.md`: the deposit's own gates and manifest.

## Repository

The same files are on the fvs-modern main branch of the advanced-forestry-systems organization at commit 1a520f5e (pull requests 118 to 126), under `config/catalog/` and `validation/bakuzis_stress_v2/`; the report, deck and this record's scripts land in `validation/bakuzis_stress_v2/release_20261029/` in the pull request that accompanies this record. Trajectory level engine output (4,485 stands in total, up to 300 per variant, x 4 arms x 21 steps) stays on Cardinal scratch and firebreather and is reproducible from `conus_catalog_runs.py` with the fvs-modern libraries.

## Not deposited

The fvs-conus remeasurement pairs file and the stand sample frame (`fia_variant_assignment.csv`), the stand initialization tables and the seed gate file: they carry FIA plot identifiers. The catalog is described by variant and stand count only, and no file in this record holds a plot identifier or coordinate (restricted data audit in RUNBOOK.md). Trajectory level output: size, and reproducible from the deposited runner.

## Known limitations

Per variant near term accuracy exists only for the Northeast cohort of Paper 2. The calibrated set is a point estimate (MAP). No arm establishes regeneration. The site based Bakuzis relations are confounded by stand selection on the unstratified sample and are reported, not gated. 263 stands (256 of them eastern) fail the year zero gate, nearly all entering the engine above their FIA basal area with verified input, and are excluded, unexplained. Thresholds are the sample maximum in variants with fewer than 1,000 fitting plots. The projector b1 GOMPIT arm has no NE or SN rows. AK and NC are not evaluated.

## License

CC BY 4.0.
