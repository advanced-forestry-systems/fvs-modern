#!/usr/bin/env Rscript
# test_conus_sf_parity.R
# Assert the Fortran engine and the R reference (conus_sf_predict.R) agree per tree
# to tolerance, for the MAP draw and a few posterior draws, for each leg and the blend.
# STATUS: scaffold.
TOL <- 1e-4
# TODO: load a sample of FIA plots, run engine (fvs2py) and R reference, compare per component
stop("TODO: implement parity harness against fvs2py shared library")
