# fvs-modern x fvs-conus — integration roadmap
2026-06-17. How the two tracks connect and advance together.

## The relationship in one picture

fvs-modern and fvs-conus are a modernization pincer on the same target (FVS predictions over CONUS,
aligned to FIA, with uncertainty):

- fvs-modern keeps the FVS ENGINE and FORMS, swaps the PARAMETERS (Bayesian, contemporary FIA) and a few
  structural levers (max SDI). Deployable now; agency-compatible.
- fvs-conus rebuilds the EQUATION FORMS (annualized, compatible, trait-driven, site-index-free) to fix
  structural defects. Longer horizon; the eventual content.

They are not alternatives. fvs-conus equations are the future CONTENT; fvs-modern is the DELIVERY VEHICLE
and the VALIDATION/BENCHMARK backbone. The end state is: fvs-conus equations running under the fvs-modern
calibration layer (max SDI + density), inside the FVS engine, validated on one disturbance-aware FIA
benchmark with uncertainty.

## The shared backbone (built this session)

1. DISTURBANCE-AWARE FIA BENCHMARK (COND treatment/disturbance stratification). Both tracks must report
   bias on the disturbance-clean stratum (or simulate removals). Today fvs-modern uses it; fvs-conus does
   not yet — its benchmark is harvest-inflated, which is why the species-free equations appear to not
   reduce stand bias. Harmonizing both onto this one benchmark is the single highest-value integration.
2. brms SITE-SPECIFIC MAX SDI + SITE PRODUCTIVITY (BGI/CSPI/ESI). The plot-level max SDI is the long-term
   density governor for fvs-modern AND the productivity driver for fvs-conus's site-index-free equations.
   One productivity surface serves both.
3. UNCERTAINTY. Both tracks are Bayesian; both carry posterior intervals. The predictive-interval
   machinery (verified ~93% coverage) is shared.
4. CROWN WIDTH / CCF. Marshall's MCW recovery gives a CONUS-consistent crown-width per species that feeds
   the competition term in both the calibrated engine and the fvs-conus crown component.

## The decisive integration experiment (next)

Run, on ONE disturbance-clean (COND-undisturbed) plot set, with held-out validation, all four arms and
report BA/TPH/QMD/volume with CIs:
  A. default FVS
  B. fvs-modern keyword-calibrated (brms SDImax + sign-aware ingrowth + BAIMULT + HT-DBH)
  C. fvs-conus species-free equations (standalone projector)
  D. COMBINED: fvs-conus equations + brms SDImax + ingrowth (the hypothesized best)
Hypothesis: C improves scatter/height/species consistency; B improves stand bias; D gets both. This
figure is the spine of both manuscripts and resolves whether the equations or the density layer (or both)
are doing the work.

## Manuscript architecture

- Paper 1 (fvs-modern, pipeline #1): National Bayesian recalibration + disturbance-aware benchmark +
  four-lever adjustment, with uncertainty. Operational, agency-facing.
- Paper 2 (fvs-conus, manuscript v2): CONUS-wide trait-driven, site-index-free, annualized/compatible
  equations; the species-free framework.
- Paper 3 (BGI/CSPI productivity surface, pipeline #3): the productivity backbone both rely on.
- Paper 4 (OLI/GMUG, pipeline #8, with Johnson/Marshall): the hierarchical Bayesian FVS framework that
  unifies the statistical approach.
Sequence: shared benchmark + productivity (3) underpin 1 and 2; 4 frames the methodology. The four-arm
integration figure appears in both 1 and 2 (from each paper's angle).

## Ordered joint next steps

1. Held-out spatial validation (both tracks) — executed this session for fvs-modern (ne/sn/kt/pn).
2. Disturbance-clean re-benchmark of fvs-conus equations (port the COND stratification into
   19_fia_benchmark_engine.R / the OSM projector).
3. The four-arm integration experiment (A/B/C/D) on one basis, with CIs.
4. Crown-width/MCW unification (Marshall + R&W Acadian) feeding both.
5. Removal-simulation converse test (validates the disturbance-artifact headline for both).
6. Productivity (BGI/CSPI/ESI) confirmed as the common driver; quantify gain vs site index.
7. Engine integration (fvs2py) so fvs-conus equations run inside FVS under the fvs-modern layer.
8. Two manuscripts advanced in parallel off the shared benchmark; Zenodo deposit (new version of the
   fvs_perseus_conus concept DOI) at submission.

## Bottom line

Keep them as two repos/manuscripts but ONE benchmark and ONE productivity backbone. The next concrete
move that advances both is the disturbance-clean, held-out, four-arm comparison — it tells us how much of
the improvement is equations (fvs-conus) vs density calibration (fvs-modern) vs both, and it is the
shared headline figure.
