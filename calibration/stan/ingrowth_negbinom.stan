// ============================================================================
// ingrowth_negbinom.stan
//
// FVS-CONUS ingrowth (recruitment) model. Per-plot count of new trees that
// enter the inventory at the second remeasurement period, conditional on
// overstory state at t1 and climate/site covariates.
//
// Negative binomial likelihood with log link. Allows over-dispersion relative
// to Poisson, which is the norm for recruitment count data.
//
// Linear predictor (log of expected new-trees per plot per year):
//   ln(lambda) = b0 + trait_dom_effect + z_L1 + z_L2 + z_L3
//              + b1 ln(BA + 1)         // overstory density
//              + b2 ln(BAL_mean + 1)   // suppression
//              + b3 SDI / SDImax       // relative density
//              + b4 ln_csi             // site productivity
//              + b5 stand_age          // age effect
//              + b6 clim_pca1          // climate component 1
//              + b7 ln(years)          // exposure / annualization
//
// trait_dom_effect uses the dominant trait composition of the existing
// overstory as a fixed effect, not species-specific RE. This is a
// plot-level model so per-species random intercepts do not apply directly.
// ============================================================================
data {
  int<lower=1> N_plots;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=0> P_trait;

  array[N_plots] int<lower=0> n_recruits;   // observed count of new trees
  vector[N_plots] log_years;                 // log of measurement interval (offset)

  vector[N_plots] ln_ba;
  vector[N_plots] ln_bal;
  vector[N_plots] rd;          // relative density
  vector[N_plots] ln_csi;
  vector[N_plots] stand_age;
  vector[N_plots] clim_pca1;

  array[N_plots] int<lower=1, upper=N_L1> L1_idx;
  array[N_plots] int<lower=1, upper=N_L2> L2_idx;
  array[N_plots] int<lower=1, upper=N_L3> L3_idx;

  matrix[N_plots, P_trait > 0 ? P_trait : 1] W_dom;  // dominant-species trait composition
}
parameters {
  real b0;
  real b1;
  real b2;
  real b3;
  real b4;
  real b5;
  real b6;

  vector[P_trait] gamma;

  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;

  real<lower=0> phi;             // negative binomial dispersion (inverse)
}
transformed parameters {
  vector[N_plots] trait_dom_effect;
  if (P_trait > 0) {
    trait_dom_effect = W_dom * gamma;
  } else {
    trait_dom_effect = rep_vector(0.0, N_plots);
  }
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
  vector[N_L2] z_L2 = sigma_L2 * z_L2_raw;
  vector[N_L3] z_L3 = sigma_L3 * z_L3_raw;
}
model {
  b0 ~ normal(0.0, 3.0);
  b1 ~ normal(-0.3, 0.5);
  b2 ~ normal(-0.2, 0.5);
  b3 ~ normal(-0.5, 1.0);
  b4 ~ normal(0.5, 0.5);
  b5 ~ normal(-0.005, 0.02);
  b6 ~ normal(0.0, 0.5);

  gamma ~ normal(0, 0.5);
  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  sigma_L1 ~ normal(0, 1.0);
  sigma_L2 ~ normal(0, 0.5);
  sigma_L3 ~ normal(0, 0.3);

  phi ~ gamma(2.0, 0.5);

  vector[N_plots] eta =
      b0
    + trait_dom_effect
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + b1 * ln_ba
    + b2 * ln_bal
    + b3 * rd
    + b4 * ln_csi
    + b5 * stand_age
    + b6 * clim_pca1
    + log_years;          // offset for annualization

  vector[N_plots] eta_safe;
  for (i in 1:N_plots) eta_safe[i] = fmin(fmax(eta[i], -20.0), 12.0);

  n_recruits ~ neg_binomial_2_log(eta_safe, phi);
}
generated quantities {
  vector[N_plots] log_lik;
  vector[N_plots] eta_gq =
      b0
    + trait_dom_effect
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + b1 * ln_ba
    + b2 * ln_bal
    + b3 * rd
    + b4 * ln_csi
    + b5 * stand_age
    + b6 * clim_pca1
    + log_years;
  for (i in 1:N_plots) {
    log_lik[i] = neg_binomial_2_log_lpmf(n_recruits[i] |
      fmin(fmax(eta_gq[i], -20.0), 12.0), phi);
  }
}
