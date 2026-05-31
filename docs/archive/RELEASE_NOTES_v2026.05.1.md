# Release notes — v2026.05.1 (2026-04-26)

## What's new

This release applies the FIA-stand Bakuzis evaluation to the
manuscript and ships the FIA stratified sampling pipeline as part
of the open-source Bakuzis runner.

### Headline scientific finding

When the Bakuzis matrix is evaluated on real FIA stand combinations
sampled by site index and basal area band rather than synthetic
algorithmic stands, the calibrated posterior credible band achieves
50 percent compliance with Eichhorn's rule for the Acadian (ACD)
variant compared with 20 percent for both the default parameters
and the calibrated MAP point estimate. This 2.5-fold uplift was
observed previously on synthetic stands and is now corroborated
on real stands. It is the strongest single empirical argument that
posterior propagation carries biological information beyond what
the MAP point estimate alone captures.

### Northeast Sukachev finding

NE Sukachev compliance remains 0 percent under all configurations
on real FIA stands. The prior synthetic-stand version of this
analysis suggested this was a synthetic stand artifact; the FIA
result confirms it is a real property of the NE variant's self-
thinning behavior. The manuscript revision explicitly corrects the
prior framing.

## Software changes

- **FIA stratified sampling pipeline**: new
  `calibration/python/fia_stand_generator.py` draws up to five real
  FIA condition pairs per Bakuzis cell from each variant's state
  coverage, with site index and basal area band stratification.
- **Replicate-aware aggregation**:
  `calibration/python/bakuzis_uncertainty_aggregate.py` extended for
  the replicate dimension, with horizon alignment and 5-year
  quantization that handles plots with different FIA inventory years
  cleanly.
- **Per-variant figures**: the matplotlib figures script now
  auto-renders trajectory, divergence, and band-growth panels for
  every variant present in the summary, not just NE.
- **Marshall format adapter**: new `marshall_to_fia_csv.py`
  converts curated Marshall-style FIA CSVs to the standard FIA
  DataMart format the FIA stand generator expects. Ready for PN/SN/IE
  Bakuzis extension once the FVS-PN/SN library issues are addressed.
- **PN/SN library diagnosis**:
  `calibration/python/PN_SN_LIBRARY_DIAGNOSIS.md` documents the
  undefined `morcon_` symbol that prevents PN and SN libraries from
  loading via fvs2py. PN is missing morts.f90 source entirely; SN
  has the source but the library was built without including it.

## Manuscript changes

- Section 4.6 rewritten with FIA-derived numbers (33 of 36 cells
  populated per variant).
- Section 3.3 methods describes the FIA stratified sampling path.
- Section 5.4 limitations updated to reflect 33-of-36 cell coverage
  and the open path to PN/SN/IE extension.
- Internal v2 self-review docx
  (`manuscript/20260426_ems_self-review_fvs-modern_v2.docx`) compiled.
  Indicative recommendation moved from Minor Revision to Accept,
  pending Zenodo DOI.

## Operational notes

- Zenodo GitHub webhook for `holoros/fvs-modern` is still pending
  enablement. Once enabled, this tag (v2026.05.1) will mint a DOI
  that should replace the placeholder in the manuscript software
  availability block before submission.

## Backwards compatibility

The aggregator and figures scripts handle both synthetic-mode
(no replicate column) and FIA-mode (replicate column present)
ensemble outputs without configuration. Existing pipelines that
produce synthetic-stand outputs continue to work unchanged.

## Verification

- `python3 calibration/python/bakuzis_uncertainty_aggregate.py
  --input-dir calibration/output/bakuzis` reproduces the laws
  compliance table shipped in
  `calibration/output/bakuzis/bakuzis_laws_compliance.csv`.
- `python3 calibration/python/bakuzis_uncertainty_figures.py`
  regenerates the seven Bakuzis figures from the shipped summary.
- `pandoc manuscript/fvs_combined_draft.md -o
  manuscript/fvs_combined_draft.docx` regenerates the manuscript
  docx (7,255 words, 11 images).
