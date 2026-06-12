#!/usr/bin/env Rscript
# conus_sf_predict.R
# Reference predictor the Fortran engine must match (parity-tested by test_conus_sf_parity.R).
# Implements the per-species blend of Leg A and Leg B and the Layer 2 modifier.

# Blended intercept for species s, draw j.
#   legB = standardize(traits_s) %*% gamma_j         (always defined)
#   legA = b0_j[s] if s is in species_dependent, else legB
#   w_s  = n_s / (n_s + kappa)   (0 if s not in Leg A)
blend_intercept <- function(comp, spcd, draw = NULL) {
  legB <- as.numeric(((comp$traits[[as.character(spcd)]] - comp$species_independent$scale_mean) /
                       comp$species_independent$scale_sd) %*% comp$species_independent$gamma)
  i <- match(spcd, comp$species_dependent$spcd)
  if (is.na(i)) return(legB)
  n   <- comp$species_dependent$n_obs[i]
  w   <- n / (n + comp$blend$kappa)
  legA <- comp$species_dependent$b0[i]
  w * legA + (1 - w) * legB
}

# apply Layer 2 disturbance/treatment modifier to a base prediction
apply_modifier <- function(base, comp, dstrb) {
  m <- comp$modifier
  if (is.null(m) || identical(m$form, "none")) return(base)
  # common: base * exp(alpha_0 + sum alpha[type]); trait_mediated: add gamma %*% traits for subset types
  stop("TODO: implement common + trait_mediated modifier application")
}

# full component prediction for one tree, one draw
predict_component <- function(comp, tree, draw = NULL) {
  b0 <- blend_intercept(comp, tree$SPCD, draw)
  # eta = b0 + shared slopes %*% tree covariates + site term  (TODO per form)
  stop("TODO: assemble linear predictor per component form, then apply_modifier and link")
}
