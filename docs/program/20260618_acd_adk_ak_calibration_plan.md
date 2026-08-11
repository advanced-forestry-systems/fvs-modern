# ACD, ADK, AK calibration via Canadian NFI (MAGPlot), and the four-way comparison
2026-06-18. Plan and status for adding the Acadian (ACD), Adirondack (ADK), and Alaska (AK) variants to the
calibration and stress test, and the full comparison across variants, species, ecoregions, and landowners.

## Data source: MAGPlot (Multi-Agency Ground Plot) database

Natural Resources Canada, Open Government Licence - Canada (attribution, permissive). Pan-Canadian
remeasured forest ground plots: NFI plus 12 jurisdictions, 1949 to present, with five site tables
(location, treatment, disturbance) and three tree tables (DBH, height, crown class, species, status, stem
and health condition). Remeasurement is explicit (subsequent measurements per site), so growth and
mortality pairs are derivable the same way as from FIA. Open release v1.0.0/v1.1.0 excludes NL and SK;
exact coordinates and restricted jurisdictions need the NFI Data Request Form.

Download package (English): MAGPlot_latest_data_package.zip, plus MAGPlot_DB_Data_Dictionary_En.pdf, from
the dataset's Data and Resources section
(https://open.canada.ca/data/en/dataset/1a73441d-cbe9-485b-bc91-84e07f97193f).

## Variant-to-data mapping

| variant | meaning | calibration data |
|---|---|---|
| ADK | Adirondacks (NY spruce-fir / northern hardwood) | NY FIA remeasurement, Adirondack counties (US data, already on Cardinal); scaffolded from Acadian |
| ACD | Acadian (Maine + Maritime Canada) | US: ME/NH/VT FIA (already used); Canada: MAGPlot New Brunswick, Nova Scotia, Prince Edward Island, plus regional NFI plots |
| AK | SE Alaska and coastal British Columbia | MAGPlot coastal British Columbia NFI plots as the calibration source for the coastal-BC half of the variant; SE Alaska proper still needs a US Alaska inventory source |

ADK does not need MAGPlot; it calibrates against existing NY FIA. ACD and AK are the two that need MAGPlot.

## Status

- ADK variant scaffolded: src-converted/adk/ created from acd/ with a calibration species-map template
  (add_variant.sh --name adk --base acd --fvs-src src-converted). Next: edit the species map, set MAXSP,
  build FVSadk.so, then calibrate against NY Adirondack FIA. No external data dependency.
- ACD and AK calibration is gated on the MAGPlot data being staged on Cardinal.

## The one dependency I cannot fulfill directly

I cannot download the MAGPlot data package myself: this environment's content rules prohibit fetching URLs
through bash (curl, wget) or other programmatic means. The dataset page is readable (I used it for this
plan), but staging the binary data package to Cardinal is a step you (or a one-line wget on Cardinal) need
to run. Once it lands at, for example, /fs/scratch/PUOM0008/crsfaaron/MAGPlot/, I proceed automatically:

    # on Cardinal (you run this; I cannot fetch URLs here)
    mkdir -p /fs/scratch/PUOM0008/crsfaaron/MAGPlot && cd $_
    wget -O MAGPlot_latest_data_package.zip "https://ca.nfis.org/fss/fss?command=retrieveByName&fileName=MAGPlot_latest_data_package.zip&fileNameSpace=magplot&format=xml&promptToSave=true"
    unzip MAGPlot_latest_data_package.zip

## Calibration pipeline once MAGPlot is staged (ACD, AK)

1. Read the MAGPlot site and tree tables; build remeasurement pairs per site (t1, t2) with DBH, height,
   species, status, and the treatment/disturbance flags (for the disturbance-clean stratification).
2. Map MAGPlot species codes to FVS species via the data dictionary; map jurisdictions to variants
   (NB/NS/PE to ACD, coastal BC to AK).
3. Compute observed diameter growth, height growth, and mortality from the pairs; fit the same
   forest-type + ecoregion component models (the cspi_traits / forest_eco family) on the Canadian pairs,
   or apply the existing CONUS species-free equations and benchmark, exactly as for the US variants.
4. Add ACD and AK to the species-specific-vs-species-free stress test and the four-way comparison.

## The four-way comparison (variants x species x ecoregion x landowner, disturbance-clean)

Buildable now on the US FIA data, then extended with ACD/ADK/AK:

- variant: the FVS variant.
- species: FIA SPCD (US) and MAGPlot species (Canada), harmonized to the FVS species list.
- ecoregion: EPA L1/L2/L3 (US); the Canadian ecozone/ecoregion fields (Canada).
- landowner: FIA OWNGRPCD (national forest, other federal, state and local, private) for the US; the
  MAGPlot ownership/jurisdiction field for Canada.
- basis: COND-undisturbed (disturbance-clean), default versus calibrated bias and RMSE for BA, TPH, QMD,
  and merch volume, faceted by the four dimensions.

The US FIA condition-level predictions carry variant and ecoregion already; OWNGRPCD joins from FIA COND
on PLT_CN. The species dimension uses the per-species observed and predicted growth (the holdout-by-species
tables) crossed with ecoregion and landowner where sample allows. I will build this on the US data while
the MAGPlot staging is resolved, then add the Canadian variants.

## Recommended order

1. Build and calibrate ADK against NY Adirondack FIA (no external dependency).
2. Build the four-way comparison on the US variants now.
3. Once MAGPlot is staged: build ACD and AK Canadian remeasurement pairs, calibrate, and fold into the
   stress test and the comparison.
