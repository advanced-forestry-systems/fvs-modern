# Paper 1 (fvs-modern/manuscript/fvs_combined_draft.md) hold inserts, 27 October 2026

Apply to the 19 October 2026 combined draft, which is on the firebreather clone and on GitHub, not in the FVS_Modern folder (that folder's `Paper1_fvs-modern_calibration_EMS_2026-08-04.md` carries no Bakuzis section and needs no change). Back the draft up with `_PREV_20261027` first. Nothing below removes text; each insert is additive and placed as marked.

## Insert A, at the top of Section 3.3 "Bakuzis matrix factorial", before "The Bakuzis matrix (Bakuzis 1969) is a classical framework"

*[ON HOLD, 27 October 2026. The committed posterior-median ensemble behind this section (calibration/output/bakuzis/bakuzis_uncertainty_summary_long.csv) starts every cell at 1.1 to 2.4 m² ha⁻¹ of basal area regardless of its Low, Medium or High density label, and no cell reaches a relative density above 0.47 in 100 years. The density axis described below therefore did not reach the engine, most plausibly because the FIA tree records were written to TREEINIT without their per-acre expansion factors; this is being checked. Figures 8 to 10 and Section 4.6 stay as written until the ensemble is re-seeded with a year-0 gate that refuses any cell whose initial basal area departs from its label by more than 20 percent.]*

## Insert B, replacing nothing, appended to law 2 in the numbered list in Section 3.3 after "...across Low, Medium, and High site classes at matched species and density."

(This is a site-ordering test of basal area. Eichhorn's rule proper, that stand volume is a near-unique function of dominant height independent of site, is tested separately in the versioned scorer as the coefficient of variation of volume across site classes at matched dominant height; both are reported in the re-run.)

## Insert C, appended to the end of Section 4.6 "Bakuzis matrix: outcomes"

A re-scoring of the same ensemble with a scorer that reports, for each relation, the number of trajectories it was evaluated on (fvs-modern `validation/bakuzis_stress_v2`) found that the Reineke self-thinning relation was not applicable on any of the 114 default and calibrated trajectories, because none entered the self-thinning zone, and that the Sukachev and yield-density relations fail in every arm for the same reason: the density labels carry no density. The Northeast Sukachev compliance of zero in `bakuzis_laws_compliance.csv` is the same finding seen through the earlier scorer. The near-identity of default and calibrated trajectories reported above is consistent with parameters that cannot express themselves on near-empty stands, and should not be read as evidence that recalibration leaves long-horizon central tendencies unchanged until the re-seeded ensemble has been run.

## Insert D, limitations paragraph in the Discussion, one sentence

The Bakuzis ensemble reported here is held pending a re-seeded run; its density axis was not realised in the committed data, so its law-compliance fractions and the default-versus-calibrated divergence in Figure 9 are not yet interpretable.

## Content inventory

Added: four inserts (A to D). Changed: nothing. Removed: nothing. Figures 8 to 10 remain in place and are labelled on hold by Insert A.
