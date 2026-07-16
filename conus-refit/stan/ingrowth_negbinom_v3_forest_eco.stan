// ============================================================================
// ingrowth_negbinom_v3_forest_eco.stan
// v3: adds a FIA forest-type-group random intercept (z_FT) on top of the v2
//     ecoregion levels (z_L1/z_L2/z_L3), so recruitment depends on both forest
//     type and ecoregion. Non-centered. New data: N_FT, FT_idx; new params:
//     z_FT_raw, sigma_FT. (derived from ingrowth_negbinom_v2.stan)
//
// FVS-CONUS ingrowth (recruitment) model. Plot-level count of new trees,
// negative binomial likelihood, log link, log(years) offset.
//
// v2 (2026-05-15): replaces stand_age with top height (HT40_1, m) and
// switches from Curtis RD to SDImax-relative RD (rd_sdimax).
//
// Linear predictor:
//   ln(lambda) = b0 + trait_dom_effect + z_L1 + z_L2 + z_L3
//              + b1 ln(BA + 1)
//              + b2 ln(BAL_mean + 1)
//              + b3 rd_sdimax       <-- relative density vs SDImax
//              + b4 ln(CSI)         <-- climate site index
//              + b5 ln(HT40 + 1)    <-- top height (replaces stand_age)
//              + b6 clim_pca1
//              + log(years)
// ============================================================================
data {
  int<lower=1> N_plots;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=1> N_FT;             // number of FIA forest-type groups
  int<lower=0> P_trait;

  array[N_plots] int<lower=0> n_recruits;
  vector[N_plots] log_years;

  vector[N_plots] ln_ba;
  vector[N_plots] ln_bal;
  vector[N_plots] rd_sdimax;     // SDI / SDImax_brms (~0 to ~1.4)
  vector[N_plots] ln_csi;
  vector[N_plots] ln_ht40;       // log(top height + 1)
  vector[N_plots] clim_pca1;

  array[N_plots] int<lower=1, upper=N_L1> L1_idx;
  array[N_plots] int<lower=1, upper=N_L2> L2_idx;
  array[N_plots] int<lower=1, upper=N_L3> L3_idx;
  array[N_plots] int<lower=1, upper=N_FT> FT_idx;   // FIA forest-type group index

  matrix[N_plots, P_trait > 0 ? P_trait : 1] W_dom;
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
  vector[N_FT] z_FT_raw;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
  real<lower=0> sigma_FT;

  real<lower=0> phi;
}
transformed parameters {
  vector[N_plots] trait_dom_effect;
  if (P_trait > 0) trait_dom_effect = W_dom * gamma;
  else trait_dom_effect = rep_vector(0.0, N_plots);
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
  vector[N_L2] z_L2 = sigma_L2 * z_L2_raw;
  vector[N_L3] z_L3 = sigma_L3 * z_L3_raw;
  vector[N_FT] z_FT = sigma_FT * z_FT_raw;
}
model {
  // Priors reflect biological expectations on growth-and-yield ingrowth dynamics
  b0 ~ normal(0.0, 3.0);
  b1 ~ normal(-0.3, 0.5);          // ln(BA): high BA -> mature canopy -> less ingrowth
  b2 ~ normal(-0.2, 0.5);          // ln(BAL): suppression -> less ingrowth
  b3 ~ normal(-1.0, 1.0);          // rd_sdimax: dense stands have less ingrowth
  b4 ~ normal(0.5, 0.5);           // ln(CSI): better sites support more ingrowth
  b5 ~ normal(-0.3, 0.5);          // ln(HT40): tall (mature) canopy -> less ingrowth
  b6 ~ normal(0.0, 0.5);           // clim_pca1: agnostic

  gamma ~ normal(0, 0.5);
  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  z_FT_raw ~ std_normal();
  sigma_L1 ~ normal(0, 1.0);
  sigma_L2 ~ normal(0, 0.5);
  sigma_L3 ~ normal(0, 0.3);
  sigma_FT ~ normal(0, 0.3);

  phi ~ gamma(2.0, 0.5);

  vector[N_plots] eta =
      b0
    + trait_dom_effect
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + b1 * ln_ba
    + b2 * ln_bal
    + b3 * rd_sdimax
    + b4 * ln_csi
    + b5 * ln_ht40
    + b6 * clim_pca1
    + log_years;

  vector[N_plots] eta_safe;
  for (i in 1:N_plots) eta_safe[i] = fmin(fmax(eta[i], -20.0), 12.0);

  n_recruits ~ neg_binomial_2_log(eta_safe, phi);
}
generated quantities {
  vector[N_plots] log_lik;
  vector[N_plots] eta_gq =
      b0
    + trait_dom_effect
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + b1 * ln_ba + b2 * ln_bal + b3 * rd_sdimax
    + b4 * ln_csi + b5 * ln_ht40 + b6 * clim_pca1
    + log_years;
  for (i in 1:N_plots) {
    log_lik[i] = neg_binomial_2_log_lpmf(n_recruits[i] |
      fmin(fmax(eta_gq[i], -20.0), 12.0), phi);
  }
}
