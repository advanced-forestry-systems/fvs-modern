FVS Variant: ADK
Base Variant: ACD
Description: Adirondack variant seeded from Acadian; calibration via MAGPlot NFI
Created: 2026-06-18

Calibration Parameters
======================
This directory holds species-specific calibration parameters for the
ADK variant. Modify these files to adjust growth, mortality,
and volume equations for your region.

Files:
  * species_map.csv     Species code mapping (FIA code -> local index)
  * height_dbh.csv      Height-diameter model coefficients
  * site_index.csv      Site index curve parameters
  * mortality.csv       Background mortality rates
  * crown_ratio.csv     Crown ratio model coefficients
  * volume.csv          Volume equation assignments

To apply calibrations, update the corresponding Fortran source files
in src-converted/adk/ and rebuild with:
  cd src-converted && OSTYPE=linux-gnu make FVSadk.so
