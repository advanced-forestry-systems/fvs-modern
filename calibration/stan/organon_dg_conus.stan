// ORGANON-form Diameter Growth Model for CONUS-wide FVS
// Based on Hann et al. (2002, 2006) SWO diameter growth equation
// Extended with species x ecodivision hierarchy, climate covariates, and relative density
// Annualization handled in R via Cao (2000) / Weiskittel et al. (2007) method;
// this Stan model fits to annualized predicted periodic growth vs observed.
//
// Model form (all metric units: DBH in cm, BA in m2/ha, BAL in m2/ha, site in m):
//   ln(DG_annual) = b0[species, ecodiv] + b1*ln(DBH + K1) + b2*DBH^K2
//                 + b3*ln((CR + 0.2)/1.2) + b4*ln(SITE_PROD)
//                 + b5*BAL/ln(DBH + K4) + b6*sqrt(BA)
//                 + b7*CLIM1 + b8*CLIM2 + b9*RD + b10*RD*ln(BAL)
//
// Estimated parameters: K1 (additive DBH offset), K2 (power on DBH)
// Fixed constants: K4 = 2.7 (in metric units)
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
  vector[N] dbh;                      // DBH1 in cm (raw; K1 applied in model)
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

parameters {
  // Hierarchical intercept: species x ecodivision
  real mu_b0;                          // Grand mean intercept
  real<lower=0.01> sigma_sp;           // SD of species intercepts
  real<lower=0.01> sigma_eco;          // SD of ecodivision intercepts
  vector[N_species] z_sp;             // Standardized species effects
  vector[N_ecodiv] z_eco;             // Standardized ecodivision effects

  // Estimated structural parameters
  real<lower=0.1> K1;                 // Additive DBH offset for ln(DBH + K1)
  real<lower=0.01> K2;                // Power on DBH in b2*DBH^K2 term

  // Fixed slope coefficients (shared across species)
  real b1;                             // ln(DBH + K1): size effect
  real b2;                             // DBH^K2: quadratic-like size effect
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
  // ---- Priors ----

  // Hierarchical intercept
  mu_b0 ~ normal(-5.0, 2.0);          // ORGANON SWO DF ~ -5.35
  sigma_sp ~ exponential(1);
  sigma_eco ~ exponential(2);          // Tighter: less variation across regions
  z_sp ~ std_normal();
  z_eco ~ std_normal();

  // Estimated structural parameters
  K1 ~ normal(6.0, 2.0);              // Additive DBH offset, centered at ~6.0
  K2 ~ normal(1.0, 0.5);              // Power on DBH, centered at 1.0

  // Slope priors (informed by ORGANON SWO parameter ranges)
  b1 ~ normal(0.8, 0.5);              // Positive size effect (~0.84 in SWO)
  b2 ~ normal(-0.04, 0.03);           // Negative coefficient on DBH^K2 term
  b3 ~ normal(1.0, 0.5);              // Positive CR effect (~1.0-1.16)
  b4 ~ normal(0.8, 0.5);              // Positive site productivity effect
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
      // Compute ln(DBH + K1) inline using estimated K1
      real ln_dbh_k1_n = log(dbh[n] + K1);

      // Compute DBH^K2 term
      real dbh_k2_n = pow(dbh[n], K2);

      // Linear predictor
      ln_dg_annual[n] = mu_b0 + b0_sp[species_id[n]] + b0_eco[ecodiv_id[n]]
                       + b1 * ln_dbh_k1_n
                       + b2 * dbh_k2_n
                       + b3 * ln_cr_adj[n]
                       + b4 * ln_site_prod[n]
                       + b5 * bal_ratio[n]
                       + b6 * sqrt_ba[n]
                       + b7 * clim1[n]
                       + b8 * clim2[n]
                       + b9 * rd[n]
                       + b10 * rd[n] * ln_bal[n];

      // Clamp linear predictor to a safe range so exp() never overflows.
      // Mirrors the HG-style hardening that resolved repeated 'is inf' warnings
      // during warmup. The bound is wide enough not to bias inference but tight
      // enough that exp(20)*years stays finite even for years=20.
      real ln_dg_safe = fmin(fmax(ln_dg_annual[n], -30.0), 20.0);

      // Annual DG with low-CR adjustment (precomputed as part of ln_cr_adj)
      // Scale by measurement period for predicted periodic growth
      dg_pred[n] = exp(ln_dg_safe) * years[n];
    }

    // Weighted normal likelihood (weight by sqrt of predicted)
    // Equivalent to: (dg_obs - dg_pred) / sqrt(dg_pred) ~ Normal(0, sigma)
    for (n in 1:N) {
      if (dg_pred[n] > 0.001) {
        target += normal_lpdf(dg_obs[n] | dg_pred[n],
                              sigma * sqrt(fmax(dg_pred[n], 0.01)));
      }
    }
  }
}

generated quantities {
}
