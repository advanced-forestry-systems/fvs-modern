#!/usr/bin/env Rscript
# 43_fit_ingrowth_composition.R  --  Trait-driven ingrowth species composition
# The count model (35_fit_ingrowth_negbinom / 35b hurdle) says HOW MANY recruit; this model
# says WHO. Hierarchical multinomial logit (or Dirichlet-multinomial) for species/group shares
# of recruits per plot, conditioned on:
#   - standing overstory composition (seed source, BA share by species/group at t1)
#   - site and climate (BGI or CSPI)
#   - disturbance/treatment state
#   - shade tolerance and seed-dispersal traits (the trait block; this is where traits earn
#     their place for ingrowth even though the count modifier stays common)
# Output: per-plot species shares; product with the count gives recruits-by-species that seed
# the tree list and supply the +N term the stand-level density equation expects.
# STATUS: scaffold. Settle hurdle vs negbinom base by held-out ELPD before production.
suppressPackageStartupMessages({ library(data.table) })
# TODO: build recruit-by-species response from TREESTATUS1 absent & TREESTATUS2==1 at DBH>=2.54
# TODO: fit hierarchical multinomial (brms/cmdstanr) with EPA L1-L3 REs + trait block
stop("scaffold: implement composition response build + hierarchical multinomial fit")
