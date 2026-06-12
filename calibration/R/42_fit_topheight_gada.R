#!/usr/bin/env Rscript
# 42_fit_topheight_gada.R  --  Stand-level top-height constraint (Garcia tradition)
# Age-independent, annualized, path-invariant Chapman-Richards GADA difference equation:
#   H2 = A * (1 - (1 - (H1/A)^(1/b2)) * exp(-b1*dt))^b2 ,  A = a0 + a1*bgi_z   (GADA on asymptote)
# First of three stand states (H, then N tied to SDImax, then G). Fit to plot-level top-height
# pairs aggregated from FIA remeasurement survivors. Reference implementation; see
# fvs-conus/stand_level/fit_topheight_gada.R for the run harness on Cardinal.
suppressPackageStartupMessages({ library(data.table) })
gada <- function(H1, dt, A, b1, b2) A*(1 - (1-(H1/A)^(1/b2))*exp(-b1*dt))^b2
# data prep: per plot, top height = mean HT of trees in the top DBH quintile at t1 and t2,
# keep plausible increments, standardize the productivity index, then nls Model 1 (fixed A)
# and Model 2 (A = a0 + a1*bgi_z). Path-invariance unit test: 10x1yr == 1x10yr.
# TODO: wire as Phase A diagnostic against 36_conus_benchmark.R before any constraint.
