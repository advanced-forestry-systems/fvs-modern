#!/usr/bin/env Rscript
# 44_stand_state_transitions.R -- Garcia-tradition stand-level transitions (reference forms).
# Age-independent, annualized, path-invariant. Validated on ~110k FIA plot pairs.
#
# DENSITY (self-thinning limb), power form on relative density RD = SDI/SDImax:
#   N2 = N1 * exp( -a * RD^c * dt )        a~0.45, c~2.73  (R2 on lnN ~0.60)
#   net density = ingrowth (negbinom count) MINUS this self-thinning loss
density_selfthin <- function(N1, RD, dt, a, c) N1 * exp(-a * RD^c * dt)
#
# BASAL AREA (monomolecular toward an SDImax-driven carrying capacity):
#   Gmax = g1 * SDImax ;  G2 = G1 + (Gmax - G1)*(1 - exp(-k*dt))   g1~0.036 k~0.053/yr
basal_area_transition <- function(G1, SDImax, dt, g1, k){ Gmax <- g1*SDImax; G1 + (Gmax - G1)*(1 - exp(-k*dt)) }
#
# TOP HEIGHT (Chapman-Richards / Bertalanffy-Richards GADA difference form):
#   H2 = A*(1 - (1-(H1/A)^(1/b2))*exp(-b1*dt))^b2 ,  A = local/site parameter (CSPI-driven)
# NEXT: refit via resde (reducible SDE) with measurement error + local site parameter,
#   to fix the negative increment skill caused by treating noisy top height as exact.
topht_gada <- function(H1, dt, A, b1, b2){ r <- pmin(pmax(H1/A,1e-6),0.999); A*(1-(1-r^(1/b2))*exp(-b1*dt))^b2 }
