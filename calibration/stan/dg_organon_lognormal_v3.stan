// ORGANON-form Diameter Growth Model for CONUS-wide FVS — v3 (2026-05-08)
// Based on Hann et al. (2002, 2006) SWO diameter growth equation
//
// CHANGES FROM v2 (dg_organon_lognormal.stan):
//   1. K1 and K2 are now FIXED CONSTANTS instead of estimated parameters.
//      v2 chains diverged catastrophically (rhat 1.86 to 3.79 across 16
//      parameters, ESS 4.4 to 5.8) because K1*b1, K2*b2, and the b3*ln(CR)
//      term collectively defined trade-off ridges along which many parameter
//      combinations gave nearly identical likelihoods. The four chains found
//      different points on the ridges and never agreed.
//      Fix: K1 = 1.0 cm (Hann SWO metric form), K2 = 1.0 (linear DBH term,
//      eliminating the power-law identifiability with b2). With K1, K2
//      fixed, b1 multiplies ln(DBH+1) and b2 multiplies DBH directly, both
//      on interpretable scales with well-defined regression coefficients.
//   2. Tightened priors that the v2 fit revealed as miscalibrated:
//        mu_b0:    N(-5.0, 2.0)  -> N(-2.0, 1.5)   (v2 fit 95% CI was [-1.36, -1.31])
//        b2:       N(-0.04, 0.03) -> N(0, 0.02)    (rescaled for K2=1; was K2=2-scale)
//        sigma_eco: exp(2)        -> exp(1)        (data wanted larger sigma_eco)
//   3. Submit script (v3) bumps warmup 1000 -> 2000 and adapt_delta 0.9 -> 0.95
//      to give the sampler more adaptation time on what should now be a
//      much better-conditioned posterior.
//
// Model form (all metric units: DBH in cm, BA in m2/ha, BAL in m2/ha, site in m):
//   ln(DG_annual) = b0[species, ecodiv] + b1*ln(DBH + 1.0) + b2*DBH
//                 + b3*ln((CR + 0.2)/1.2) + b4*ln(SITE_PROD)
//                 + b5*BAL/ln(DBH + 2.7) + b6*sqrt(BA)
//                 + b7*CLIM1 + b8*CLIM2 + b9*RD + b10*RD*ln(BAL)
//
// Hierarchy:
//   b0[s,d] = mu_b0 + sigma_sp * z_sp[s] + sigma_eco * z_eco[d]
//   (Non-centered parameterization for efficient HMC sampling)

data {
  int<lower=1> N;                      // Number of tree-period observations
  int<lower=1> N_species;              // Number of species
  int<lower=1> N_ecodiv;               // Number of ecodivisions

  // Observed periodic diameter growth (DBH2 - DBH1)
  vector[N] dg_obs;

  // Tree-level covariates at measurement start (t1)
  // All inputs in metric units: DBH in cm, BA in m2/ha, BAL in m2/ha, site in m
  vector[N] dbh;                      // DBH1 in cm
  vector[N] ln_cr_adj;                // ln((CR1 + 0.2) / 1.2)
  vector[N] ln_site_prod;             // ln(ClimateSI), ln(BGI), or ln(climate PC)
  vector[N] bal_ratio;                // BAL1 / ln(DBH1 + 2.7), precomputed with K4=2.7
  vector[N] ln_bal;                   // ln(BAL1 + 5.0), precomputed
  vector[N] sqrt_ba;                  // sqrt(BA1) in m2/ha
  vector[N] clim1;                    // Climate variable 1 (standardized)
  vector[N] clim2;                    // Climate variable 2 (standardized)
  vector[N] rd;                       // Relative density (0 to ~3.0), precomputed in R

  // Measurement interval
  vector<lower=0>[N] years;           // Remeasurement interval (years)

  // Group indices
  array[N] int<lower=1, upper=N_species> species_id;
  array[N] int<lower=1, upper=N_ecodiv> ecodiv_id;
}

transformed data {
  // Fixed structural constants (was: estimated in v2, caused identifiability collapse)
  real K1 = 1.0;                      // Additive DBH offset (Hann SWO metric form)
  real K2 = 1.0;                      // Power on DBH (linear; simplest identifiable form)

  // Pre-compute the transformed DBH terms once
  vector[N] ln_dbh_k1 = log(dbh + K1);
  vector[N] dbh_k2;
  for (n in 1:N) {
    dbh_k2[n] = pow(dbh[n], K2);     // = dbh[n] when K2 == 1, kept generic for clarity
  }
}

parameters {
  // Hierarchical intercept: species x ecodivision
  real mu_b0;                          // Grand mean intercept
  real<lower=0.01> sigma_sp;           // SD of species intercepts
  real<lower=0.01> sigma_eco;          // SD of ecodivision intercepts
  vector[N_species] z_sp;             // Standardized species effects
  vector[N_ecodiv] z_eco;             // Standardized ecodivision effects

  // Fixed slope coefficients (shared across species)
  real b1;                             // ln(DBH + K1): size effect
  real b2;                             // DBH^K2 (= DBH when K2=1): linear size effect
  real b3;                             // ln((CR+0.2)/1.2): crown ratio
  real b4;                             // ln(SITE_PROD): site productivity
  real b5;                             // BAL/ln(DBH+K4): competition
  real b6;                             // sqrt(BA): stand density
  real b7;                             // CLIM1: climate variable 1
  real b8;                             // CLIM2: climate variable 2
  real b9;                             // RD: relative density (standalone)
  real b10;                            // RD * ln(BAL): relative density interaction

  // Observation error
  real<lower=0.01> sigma;
}

transformed parameters {
  // Non-centered species and ecodivision effects
  vector[N_species] b0_sp = sigma_sp * z_sp;
  vector[N_ecodiv] b0_eco = sigma_eco * z_eco;
}

model {
  // ---- Priors (tightened from v2) ----

  // Hierarchical intercept
  mu_b0 ~ normal(-2.0, 1.5);          // v2 fit pulled hard from N(-5, 2) prior; recenter
  sigma_sp ~ exponential(1);
  sigma_eco ~ exponential(1);          // Loosened from exp(2): v2 wanted ~1.4
  z_sp ~ std_normal();
  z_eco ~ std_normal();

  // Slope priors (informed by ORGANON SWO parameter ranges; b2 rescaled for K2=1)
  b1 ~ normal(0.5, 0.4);              // Positive size effect on ln(DBH+1)
  b2 ~ normal(0, 0.02);               // Linear DBH effect, sign uncertain (small magnitude)
  b3 ~ normal(1.0, 0.5);              // Positive CR effect (~1.0-1.16)
  b4 ~ normal(0.5, 0.4);              // Positive site productivity effect
  b5 ~ normal(-0.008, 0.005);         // Negative competition effect
  b6 ~ normal(-0.04, 0.03);           // Negative stand density effect
  b7 ~ normal(0, 1);                  // Climate 1 (weakly informative)
  b8 ~ normal(0, 1);                  // Climate 2 (weakly informative)
  b9 ~ normal(-0.5, 0.5);             // Relative density (higher RD reduces growth)
  b10 ~ normal(0, 0.3);               // RD x ln(BAL) interaction (weakly informative)

  // Observation error
  sigma ~ exponential(1);

  // ---- Likelihood ----
  // Direct annual model: predict annual DG, multiply by period length
  // (This is the PAI approximation; full Cao/Weiskittel iteration is in R)
  {
    vector[N] ln_dg_annual;
    vector[N] dg_pred;

    for (n in 1:N) {
      // Linear predictor (uses precomputed transformed DBH terms)
      ln_dg_annual[n] = mu_b0 + b0_sp[species_id[n]] + b0_eco[ecodiv_id[n]]
                       + b1 * ln_dbh_k1[n]
                       + b2 * dbh_k2[n]
                       + b3 * ln_cr_adj[n]
                       + b4 * ln_site_prod[n]
                       + b5 * bal_ratio[n]
                       + b6 * sqrt_ba[n]
                       + b7 * clim1[n]
                       + b8 * clim2[n]
                       + b9 * rd[n]
                       + b10 * rd[n] * ln_bal[n];

      // Clamp linear predictor to a safe range so exp() never overflows.
      real ln_dg_safe = fmin(fmax(ln_dg_annual[n], -30.0), 20.0);

      // Scale to predicted periodic growth
      dg_pred[n] = exp(ln_dg_safe) * years[n];
    }

    // Log-normal likelihood
    for (n in 1:N) {
      if (dg_pred[n] > 0.001 && dg_obs[n] > 0.001) {
        target += lognormal_lpdf(dg_obs[n] | log(dg_pred[n]), sigma);
      }
    }
  }
}

generated quantities {
}
