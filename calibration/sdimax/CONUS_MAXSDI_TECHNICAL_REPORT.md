# Localizing maximum stand density index for the Forest Vegetation Simulator: evidence that the species-weighted maximum is biased and a model-agnostic replacement

**Prepared for:** the USDA Forest Service Forest Vegetation Simulator (FVS) staff
**Prepared by:** A. Weiskittel, Center for Research on Sustainable Forests, University of Maine, with collaborators
**Date:** 14 June 2026
**Status:** technical brief for discussion. All analyses are on FIA remeasurement data; no FVS source or default was modified.

---

## Executive summary

FVS, like most growth-and-yield systems, sets a stand's maximum stand density index (maximum SDI) by weighting fixed per-species maximum-SDI constants by each species' share of the stand. This brief presents evidence, from FIA remeasurement data across the conterminous United States, that this species-weighted maximum has two quantifiable problems and proposes a model-agnostic replacement.

First, the species-weighted maximum is biased high by about 28 percent relative to an FIA-estimated maximum, and it has effectively no skill at the plot level: even after removing the mean bias it explains about 2 percent of the plot-to-plot variation in the maximum (a negative raw coefficient of determination). The reason is structural. Per-species maxima are conservative pure-stand upper bounds, and basal-area-weighting them across a real, mixed stand systematically overstates the achievable maximum while carrying no information about location.

Second, and decisively, maximum SDI is not directly observable, so the case cannot rest on matching any one estimate. The criterion that is not circular is whether the choice of maximum improves the prediction of observed stand dynamics. Relative density (SDI divided by the maximum) is the variable that drives self-thinning, so the better maximum is the one whose relative density better predicts observed density loss in remeasurement. On 82,130 remeasured plots, a localized FIA-derived maximum predicts observed self-thinning about 85 percent better than the species-weighted maximum (deviance explained 0.107 versus 0.058), and it wins in every region, East and West.

The recommendation is to treat the maximum SDI as a localized, data-derived stand attribute, looked up by location and composition from an FIA-based surface, rather than computed from per-species constants inside the model. In FVS this is the per-stand SDIMAX keyword; in any other engine it is the stand's maximum-SDI input. This decouples the density limit from the model's internal species table and is reusable by any growth-and-yield model.

---

## 1. Background: how FVS sets the density limit

FVS constrains density-dependent mortality and competition through stand density index (Reineke 1933) and its maximum. Each species carries a maximum-SDI constant in the variant's species table, and a stand's maximum is formed by weighting those constants by composition (basal-area share by default), which can be overridden per stand with the SDIMAX keyword. Relative density, SDI divided by this maximum, then governs the onset and rate of density-dependent mortality. The density limit is therefore one of the most consequential single quantities in any long-term projection: it sets the ceiling the stand self-thins toward, and small errors in it compound over a rotation.

This species-weighting method is not unique to FVS. ORGANON, the regional FVS variants, and essentially any model that builds a stand maximum from fixed per-species constants share the same construction, so the findings below are about the method, not about one model.

## 2. The problem in two parts

### 2.1 The maximum is not observable

A stand's maximum SDI is a latent quantity. It is the density the stand would carry at the self-thinning limit, which most stands are not at, so it must be estimated under assumptions rather than measured. This has a methodological consequence that shapes everything below: comparing two estimates of the maximum by how well one reproduces the other is partly circular and cannot settle which is better. The only non-circular test is predictive, against observed data, and we use that test in Section 4.

### 2.2 Species-weighting is biased high and carries little location signal

Using an FIA-estimated maximum SDI as the reference (described in Section 3), on 95,206 remeasured plots we compared the FVS variant-specific species-weighted maximum (computed by basal-area-weighting the per-species constants from the FVS variant species tables over each plot's species; 88 percent of stems matched a species constant) against the FIA reference:

| approach | mean (trees/ha) | bias vs FIA reference | raw R-squared | bias-corrected R-squared |
|---|---:|---:|---:|---:|
| FIA-estimated maximum (reference) | 869 | — | — | — |
| FVS species-weighted (variant-specific) | 1,110 | +28% | -0.61 | 0.02 |

Two findings. The species-weighted maximum overstates the maximum by about 28 percent, and its negative raw R-squared means it predicts the plot maximum worse than simply assigning every plot the overall mean. Removing the mean bias does not rescue it: the bias-corrected R-squared is 0.02, so it carries almost no information about where the real maximum is from plot to plot. The upward bias follows directly from averaging conservative pure-stand upper bounds across mixed stands; the absence of location signal follows from the maximum being a function of species table entries that do not vary in space.

### 2.3 What does localize the maximum

On 113,270 plots with a valid FIA-estimated maximum, the maximum averages 862 trees/ha with a standard deviation of 336 (coefficient of variation 0.39), so there is large, real variation to capture. The variance explained by candidate localizers:

| localizer | R-squared |
|---|---:|
| Forest type (FORTYPCD) | 0.21 |
| EPA Level-3 ecoregion | 0.15 |
| Geographic smooth s(lon, lat) | 0.16 |
| Site class (SICOND) | 0.02 |
| Forest type + geography | 0.25 |
| Forest type + geography + site class | 0.22 |

Maximum SDI is driven by what grows somewhere (composition) and where it grows (geography), and not by how productive the site is: site class is essentially irrelevant (0.02) and adds nothing on top of forest type and geography. This is worth flagging on its own, because it is intuitive to summarize the maximum by site quality, and the data say not to. Forest type plus a spatial smooth explains 0.25, decisively more than species-weighting, but still only a quarter of the plot-level variation; the remaining structure is local and is recovered only by a per-stand value from a wall-to-wall surface.

## 3. Data and methods

**FIA remeasurement.** The analysis uses conterminous-US FIA remeasured plots, paired by control number to their prior measurement, with stand density index at both visits, the remeasurement interval, composition, EPA ecoregion codes, coordinates, forest type, and site class. Density metrics are in metric units (trees per hectare; SDI on the metric convention).

**FIA-estimated maximum (the reference and the localized surface).** The reference maximum is a Bayesian (brms) estimate of plot-level maximum SDI fit to the FIA data (173,700 plots), and its wall-to-wall form is the TreeMap 2022-based 30 m CONUS raster of SDI, maximum SDI, and relative density (Chivhenge, Weiskittel, Woodall, D'Amato, and Daigneault; Zenodo 10.5281/zenodo.19509367). The brms plot values are the FIA-plot form of the same surface and join to the remeasurement data by plot key. These estimates are independent of any FVS species table.

**The non-circular validation.** Because the maximum is not observable, the test of record is predictive. For each plot we computed relative density at the first visit as SDI divided by the maximum, under each candidate maximum, and the observed annual density change from the remeasurement as the negative log ratio of trees per hectare divided by the interval in years. We then fit, separately for each candidate maximum, a smooth model of observed density change on relative density and recorded the deviance explained. The maximum whose relative density better predicts the observed density loss is the better maximum on the only ground that is not circular. We report this nationally and split East and West at longitude 103 W.

## 4. Results: the localized maximum predicts observed self-thinning better

On 82,130 remeasured plots, predicting observed annual density change from relative density:

| region | n | deviance explained, RD from FIA-localized maximum | deviance explained, RD from FVS species-weighted |
|---|---:|---:|---:|
| All | 82,130 | 0.107 | 0.058 |
| East | 62,062 | 0.124 | 0.072 |
| West | 20,068 | 0.057 | 0.034 |

The localized maximum predicts observed self-thinning about 85 percent better than species-weighting overall, and it wins in every region. The correlation of relative density with observed mortality is also higher for the localized maximum everywhere (0.33 versus 0.24 nationally). Because this test uses observed density change as truth, it is not subject to the circularity of comparing two estimates of an unobservable quantity, and it points the same way as the bias and skill findings in Section 2.

A mechanism note ties the two FVS-side density problems together. The native species-weighted maximum is too high (+28 percent), so its relative density is too low and it under-predicts self-thinning. Separately, when the maximum is recalibrated downward inside an engine configuration it can swing the other way and over-predict self-thinning. In a controlled FVS run on Northeast plots, varying only the maximum, a low calibrated SDIMAX over-thinned the stand (trees per hectare biased -26.5 percent) while the native species-weighted and the FIA-localized maxima both returned density to near-unbiased (-0.4 and +3.5 percent). The two FVS-side sources err in opposite directions; the FIA-localized maximum sits between them and predicts observed self-thinning best.

## 5. Recommendation: a localized, data-derived, model-agnostic maximum

Set the density limit from a localized, data-derived maximum SDI, looked up by location and composition, rather than from composition-weighted per-species constants. In order of fidelity:

1. **Per-stand value from a wall-to-wall surface (best).** Assign each stand the maximum SDI at its coordinates from an FIA-derived surface. The TreeMap-based 30 m CONUS maximum-SDI raster (Zenodo 10.5281/zenodo.19509367) is exactly such a product, and the brms FIA plot values are its plot-level form. This captures the full local and regional structure, not just the quarter a summary recovers, and it is reusable by any model because it is a number attached to a location.
2. **Composition-and-geography model (compact, portable).** Where a closed form is wanted, the maximum as forest type plus a spatial smooth (R-squared 0.25) already beats species-weighting and travels as a small table or function. Do not include site class; it carries no signal for the maximum.
3. **Avoid the species-weighted constant approach** in any model. It is biased high and uninformative about the real maximum.

This is deliberately model-agnostic. The maximum becomes a per-stand input rather than an internal computation, which decouples the density limit from any one model's species table and lets FVS, ORGANON, and other engines consume the same value.

## 6. Implementation in FVS

The change is small and does not touch FVS source or species tables. For each stand, look up the localized maximum (raster value at the stand coordinates, or the plot's FIA value where it is an FIA plot, or the forest-type-plus-geography fallback) and set it through the SDIMAX keyword per stand. Setting every species in the stand to the localized stand maximum makes the basal-area-weighted stand maximum equal that localized value regardless of composition, which is the intended behavior. We have implemented this as a small per-stand lookup-and-keyword module that resolves the maximum from the raster, the FIA plot table, or the fallback, and emits the SDIMAX block. Units convert from trees per hectare to the FVS internal trees-per-acre convention. The same per-stand value serves as the maximum input to any non-FVS engine.

## 7. Caveats and the proposed confirming test

The 28 percent bias and the predictive-skill advantage are national results on FIA remeasurement; the magnitude varies regionally. In the Northeast specifically, the native species-weighted maximum is already close to adequate for density (the FIA value there is only marginally better), so the payoff of localization is concentrated where species-weighting is most wrong: the West and structurally complex mixed stands. The natural confirming test, which we recommend and intend to run, is the same density benchmark on a Western variant (for example PN or CR), comparing the species-weighted SDIMAX against the localized value, where the localized surface should reduce both the density bias and any over-thinning. We would welcome the FVS staff's view on the most useful variant and FIA stratification for that test.

A second caveat is definitional consistency. The FIA-estimated maximum is on a specific SDI convention (metric, summation method); adopting it operationally requires matching the convention FVS uses internally so the relative density that drives self-thinning is computed consistently. This is a units-and-definitions check, not a modeling obstacle.

## 8. Side-by-side FVS projections (Pacific Northwest demonstration)

To show the operational consequence inside FVS rather than only in the statistical self-thinning test,
we ran paired FVS Pacific Northwest (PN) projections on 101 remeasured Oregon and Washington plots,
identical in every respect except the maximum SDI: the FVS default (species-weighted, internal) versus
the localized FIA-derived value supplied per stand through the SDIMAX keyword. The PN region is the
intended hard case, high-density Douglas-fir and western hemlock where species-weighting should be
most consequential. We measured each projection two ways: against the observed remeasurement (does the
localized maximum predict the observed density better) and over a 100-year horizon (how the choice
compounds).

**The localized maximum improves the density prediction, and the improvement scales with how dense the
stand is, exactly as the mechanism predicts.** Error in projected density (trees per hectare) against
the observed remeasurement:

| stand set | n | default (species-weighted) % RMSE | localized (FIA-derived) % RMSE |
|---|---:|---:|---:|
| all stands | 101 | 46 | 42 |
| binding (relative density > 0.45) | 68 | 39 | 34 |
| dense (relative density > 0.6) | 30 | 35 | 26 |

In the densest stands, where the density limit actually governs the projection, localizing the maximum
cuts the density error by about a quarter (35 to 26 percent RMSE). Basal area and quadratic mean
diameter are essentially unchanged (within a point), which is the correct and reassuring signature: the
maximum SDI governs density and self-thinning, not tree size, so a correct maximum should move density
and leave size alone. The improvement is concentrated where it should be and absent where it should be.

**A candid limitation that is itself a useful finding for the FVS staff.** FVS's projection sensitivity
to the maximum SDI is bounded by design: the density limit drives the density-dependent component of
mortality once a stand approaches it, so in stands well below the limit the choice of maximum changes
little, and over a 100-year horizon the per-stand density difference between the two maxima is modest
on average (a few trees per hectare) but real, signed, and growing with relative density (panel B of
the figure). The practical reading is that correcting the maximum matters most for dense-stand and
long-horizon applications (carbon, fuels, density management), and matters little for young or open
stands. This bounds the claim honestly: the strong statistical evidence for the localized maximum (the
85 percent better self-thinning prediction) translates inside FVS into a meaningful density improvement
concentrated in the binding regime, not a wholesale change to every projection.

![Default versus localized maximum SDI in FVS PN projections](pn_maxSDI_demo.png)

*Figure. Paired FVS PN projections, default versus localized maximum SDI, 101 remeasured OR/WA plots.
(A) Density error against the observed remeasurement falls with the localized maximum, most in dense,
binding stands. (B) Over a 100-year projection the per-stand density difference is real, signed
(localized thins slightly more), and grows with relative density.*

This demonstration is reproducible (`calibration/python/pn_sdimax_sidebyside.py` in fvs-modern). The
natural extension, which we recommend, is to repeat it on an interior-West variant where the native
species-weighted maximum and the FIA value diverge most, and to pair it with the same demonstration on
the regional variants where species-weighting is known to be furthest off.

## 9. Bottom line

Maximum SDI is not observable, so the case rests on predictive skill, and on that ground the result is clear: a localized, FIA-derived maximum predicts observed self-thinning about 85 percent better than the species-weighted maximum, in every region, while species-weighting is also biased about 28 percent high and explains almost none of the plot-level variation. The recommended change is to set the density limit from a localized, data-derived maximum SDI per stand, supplied to FVS through the SDIMAX keyword and to any other model as its maximum-SDI input, decoupling the density limit from the internal species table. It is a small operational change with a measurable improvement in the density dynamics the model exists to reproduce, and it generalizes beyond FVS.

---

### Selected references

Chivhenge, E., A. Weiskittel, C. Woodall, A. D'Amato, and A. Daigneault. A 30 m wall-to-wall raster of stand density index, maximum stand density index, and relative density for the conterminous United States from TreeMap 2022 and FIA. Zenodo. https://doi.org/10.5281/zenodo.19509367

Houtman, R., et al. TreeMap 2022: a tree-level model of the forests of the conterminous United States.

Reineke, L. H. 1933. Perfecting a stand-density index for even-aged forests. Journal of Agricultural Research 46:627-638.

Weiskittel, A. R., D. W. Hann, J. A. Kershaw, and J. K. Vanclay. 2011. Forest Growth and Yield Modeling. Wiley.
