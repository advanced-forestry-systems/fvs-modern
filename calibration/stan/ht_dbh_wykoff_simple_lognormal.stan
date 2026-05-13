// ============================================================================
// ht_dbh_wykoff_simple_lognormal.stan
//
// FVS-CONUS static height-diameter model, Wykoff (1986) form with
// simplified hierarchy: species + EPA L1 random intercepts only.
// Log-normal likelihood to eliminate the variance bimodality that caused
// rhat ~1.53 in the normal-likelihood version.
//
// v2: Trait architecture refactored. Species functional traits (W * gamma)
//     now enter as FIXED EFFECTS in the linear predictor, separate from the
//     species random intercept. The species RE (z_sp) captures only residual
//     species variation not explained by traits, with a tighter prior on
//     sigma_sp to encourage shrinkage. This gives interpretable gamma
//     coefficients and a model that generalizes across species via traits.
//
//   log(HT - 1.37) = b0 + b1/(DBH+1) + W[sp]*gamma + z_sp[sp] + z_L1 + covs
//
// Back-transform: HT = 1.37 + exp(eta)
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;

  int<lower=0> P_trait;

  vector<lower=0>[N_obs] ht_obs;     // total height (m), must be > 1.37
  vector<lower=0>[N_obs] dbh;

  vector[N_obs] bal;
  vector[N_obs] sqrt_ba;
  vector[N_obs] ln_cspi_shift;
  vector[N_obs] ba_x_rd;
  vector[N_obs] bal_x_rd;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
transformed data {
  // Pre-compute the log-transformed response
  vector[N_obs] log_ht_above_bh;
  for (i in 1:N_obs)
    log_ht_above_bh[i] = log(fmax(ht_obs[i] - 1.37, 0.01));
}
parameters {
  real b0;
  real<upper=0> b1;             // classical Wykoff has b1 < 0

  real a_bal;
  real a_ba;
  real a_cspi;
  real a_bard;
  real a_blrd;

  vector[P_trait] gamma;
  vector[N_sp] z_sp_raw;
  vector[N_L1] z_L1_raw;
  real<lower=0> sigma_sp;
  real<lower=0> sigma_L1;

  real<lower=0> sigma;           // single residual SD on log scale
}
transformed parameters {
  // Traits as fixed effects (separate from species RE)
  vector[N_sp] trait_effect;
  if (P_trait > 0) {
    trait_effect = W * gamma;
  } else {
    trait_effect = rep_vector(0.0, N_sp);
  }
  // Species RE: purely residual, zero-mean (captures what traits miss)
  vector[N_sp] z_sp = sigma_sp * z_sp_raw;
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
}
model {
  b0 ~ normal(3.2, 0.8);
  b1 ~ normal(-6.0, 3.0);

  a_bal  ~ normal( 0.0, 0.02);
  a_ba   ~ normal( 0.0, 0.05);
  a_cspi ~ normal( 0.2, 0.3);
  a_bard ~ normal( 0.0, 0.02);
  a_blrd ~ normal( 0.0, 0.02);

  gamma ~ normal(0, 0.5);
  z_sp_raw ~ std_normal();
  z_L1_raw ~ std_normal();
  sigma_sp ~ normal(0, 0.2);    // tighter: traits should explain most species variation
  sigma_L1 ~ normal(0, 0.3);

  sigma ~ normal(0, 0.5);

  // Linear predictor: eta = log(HT - 1.37) in expectation
  // trait_effect[sp] = fixed species effect from traits
  // z_sp[sp] = residual species RE (shrunk toward zero)
  vector[N_obs] eta =
      b0
    + trait_effect[sp_idx] + z_sp[sp_idx] + z_L1[L1_idx]
    + a_bal  * bal
    + a_ba   * sqrt_ba
    + a_cspi * ln_cspi_shift
    + a_bard * ba_x_rd
    + a_blrd * bal_x_rd
    + b1 ./ (dbh + 1.0);

  // Log-normal likelihood: no exp(), no heteroscedastic variance
  log_ht_above_bh ~ normal(eta, sigma);
}
generated quantities {
  vector[N_obs] mu_pred;
  vector[N_obs] log_lik;
  {
    vector[N_obs] eta =
        b0
      + trait_effect[sp_idx] + z_sp[sp_idx] + z_L1[L1_idx]
      + a_bal*bal + a_ba*sqrt_ba + a_cspi*ln_cspi_shift
      + a_bard*ba_x_rd + a_blrd*bal_x_rd
      + b1 ./ (dbh + 1.0);
    // Back-transform to natural scale (with log-normal bias correction)
    for (i in 1:N_obs)
      mu_pred[i] = 1.37 + exp(eta[i] + 0.5 * sigma * sigma);
    for (i in 1:N_obs)
      log_lik[i] = normal_lpdf(log_ht_above_bh[i] | eta[i], sigma);
  }
}
