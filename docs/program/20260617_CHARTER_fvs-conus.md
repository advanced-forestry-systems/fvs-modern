# fvs-conus — program charter and roadmap
2026-06-17. Track 2 of the two-track FVS modernization program. Companion: CHARTER_fvs-modern and
INTEGRATION_roadmap.

## One-sentence scope

Develop and implement CONUS-wide, FIA-aligned tree- and stand-level equations — in both species-dependent
and species-independent (trait-driven) forms — that fix the structural shortcomings of the legacy FVS
equations, with core modeling properties of annualization, compatibility, tractability, and robustness,
built on non-traditional site-productivity measures (BGI, CSPI, ESI) that avoid the age- and
species-dependence of site index.

## Why it matters / positioning

fvs-conus is the rebuild track. Where fvs-modern re-estimates the existing equation forms, fvs-conus
replaces the forms themselves to remove known defects: inconsistent species/variant equations (e.g.
Douglas-fir has 4+ different crown-width and growth forms across 17 variants), non-annualized periodic
models, incompatibility between components, and reliance on site index (which is age- and
species-dependent and undefined for mixed/uneven-aged stands). It is the longer-horizon scientific
contribution; fvs-modern is its operational delivery vehicle and benchmark.

## Core design principles

- ANNUALIZATION: 1-yr time step so any projection length/cycle composes correctly.
- COMPATIBILITY: components (diameter growth, height, crown, mortality, ingrowth, volume) share state and
  are mutually consistent (e.g. HT-DBH consistent with height increment; CCF/crown-width consistent
  across variants).
- TRACTABILITY: closed-form / efficient evaluation; deployable inside or alongside the FVS engine.
- ROBUSTNESS: trait-driven species-free fallback so the model predicts for unmeasured/rare species
  (471 FIA species vs 245 FVS species; ~232 to map).
- SITE PRODUCTIVITY WITHOUT SITE INDEX: BGI (biophysical growth index), CSPI (climate site productivity
  index), ESI (ecological/environmental site index) as continuous, age- and species-independent drivers.

## What exists (banked on Cardinal, fvs-conus/output)

- Per-variant Bayesian fits: diameter growth (Wykoff/Kuehne), height-diameter, height increment, crown
  ratio, with posteriors and uncertainty (RDS draws, summaries, maps).
- CONUS trait-driven diameter growth (Kuehne, CSPI traits): species random effects replaced by traits;
  validated that traits reproduce species RE across two base models (HG ORGANON and DG Kuehne), all 8
  trait coefficients excluding zero, sigma within ~4-5% of the species-RE model.
- Existing FVS-vs-fvs-conus benchmark (allvar_calibration.csv, comparisons_overstory/): the species-free
  equations REDUCE scatter (RMSE) and improve height/species consistency, but on the current
  (not disturbance-stratified) benchmark do NOT reduce stand-level BA/volume bias (BA +13->+17, volume
  +23->+29 vs default). This is the key result to re-examine on the disturbance-clean basis.
- Companion productivity manuscript (BGI/CSPI 30 m CONUS surface) — pipeline #3.

## What's left

1. Re-benchmark fvs-conus equations on the DISTURBANCE-CLEAN (COND-undisturbed) basis (shared with
   fvs-modern) so its bias is not harvest-inflated; compare to default and to the keyword-calibrated FVS.
2. Combine fvs-conus equations with the stand-level density levers (brms max SDI + ingrowth) — test the
   hypothesis that the equations fix scatter/height/species and the density layer fixes stand bias.
3. Crown-width / CCF unification (Marshall MCW recovery): one CONUS-consistent per-species MCW from the
   variant CCF coefficients (+ Russell-Weiskittel Acadian), feeding the compatible crown/competition term.
4. Finish the remaining species-free component refits (mortality, ingrowth, volume) to full CONUS
   coverage; resolve species mapping (genus / softwood-hardwood) for the ~232 unmapped FIA species.
5. Site-productivity integration: confirm BGI/CSPI/ESI as the productivity driver across all components;
   quantify the gain over site index.
6. Engine integration (Product B): fvs2py in-process path to run the equations INSIDE the FVS engine
   (maintainer-level, David Diaz/MicroFVS) OR a full-dynamics standalone projector with ingrowth.
7. Full-data production refits where pilots used subsets (e.g. B1 HG refit).

## Deliverables / artifacts

- Repo: holoros/fvs-modern (sf_integration_dev, conus/) and the fvs-conus working tree on Cardinal.
- Manuscript v2 (trait-driven species-free framework): Section 3.7 (DG Kuehne), Discussion, figures —
  scope in 20260513_manuscript_v2_scope.md. Plus the productivity (BGI/CSPI) manuscript (pipeline #3).
- OLI/GMUG hierarchical Bayesian FVS manuscript (#8, with Greg Johnson, David Marshall) overlaps both
  tracks — the Bayesian framing is shared.

## Recommended next session (fvs-conus)

Re-benchmark the species-free equations on the disturbance-clean basis and test the combined
equations+density product; implement the Marshall MCW unification; advance manuscript v2 (Section 3.7 +
Discussion + Figure 1).
