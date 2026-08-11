// ============================================================================
// hcb_unified.stan
//
// Unified HCB model: species-free B1 hierarchy PLUS a species random effect.
// Same as hcb_organon_speciesfree.stan (traits W*gamma + nested L1/L2/L3
// ecoregion + forest-type REs), with a per-species random intercept z_sp added
// back. z_sp shrinks toward 0 (the trait-predicted mean) via a half-normal
// sigma_sp, so EVERY species gets an estimate: data-rich species pull away from
// the trait prediction, data-poor species shrink to it. No rare-species pooling.
//
//   eta = h0 + trait_effect[sp] + z_sp[sp]
//       + z_L1 + z_L2 + z_L3 + z_FT
//       + h1 ln(HT) + h2 ln(DBH) + h3 BAL/(HT+1) + h4 sqrt(BA) + h5 ln(CSPI)
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=1> N_FT;
  int<lower=0> P_trait;

  vector<lower=0, upper=1>[N_obs] ratio;
  vector[N_obs] ln_ht;
  vector[N_obs] ln_dbh;
  vector[N_obs] bal_over_ht;
  vector[N_obs] sqrt_ba;
  vector[N_obs] ln_cspi_shift;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;
  array[N_obs] int<lower=1, upper=N_FT> FT_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
parameters {
  real h0;
  real h1; real h2; real h3; real h4; real h5;

  vector[P_trait] gamma;

  vector[N_sp] z_sp_raw;            // NEW: species random intercept (non-centered)
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_FT] z_FT_raw;

  real<lower=0> sigma_sp;           // NEW: species RE scale (shrinkage toward traits)
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
  real<lower=0> sigma_FT;

  real<lower=1> phi;
}
transformed parameters {
  vector[N_sp] trait_effect;
  if (P_trait > 0) {
    trait_effect = W * gamma;
  } else {
    trait_effect = rep_vector(0.0, N_sp);
  }
  vector[N_sp] z_sp = sigma_sp * z_sp_raw;
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
  vector[N_L2] z_L2 = sigma_L2 * z_L2_raw;
  vector[N_L3] z_L3 = sigma_L3 * z_L3_raw;
  vector[N_FT] z_FT = sigma_FT * z_FT_raw;
}
model {
  h0 ~ normal( 0.0, 2.0);
  h1 ~ normal( 0.5, 0.5);
  h2 ~ normal(-0.3, 0.5);
  h3 ~ normal( 0.0, 0.5);
  h4 ~ normal( 0.0, 0.5);
  h5 ~ normal(-0.2, 0.5);

  gamma ~ normal(0, 0.5);

  z_sp_raw ~ std_normal();
  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  z_FT_raw ~ std_normal();
  sigma_sp ~ normal(0, 0.5);       // residual species variation beyond traits
  sigma_L1 ~ normal(0, 0.5);
  sigma_L2 ~ normal(0, 0.3);
  sigma_L3 ~ normal(0, 0.3);
  sigma_FT ~ normal(0, 0.3);

  phi ~ gamma(2, 0.05);

  vector[N_obs] eta =
      h0 + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + h1 * ln_ht + h2 * ln_dbh
    + h3 * bal_over_ht + h4 * sqrt_ba
    + h5 * ln_cspi_shift;

  vector[N_obs] mu_raw = inv_logit(eta);
  vector[N_obs] mu = 0.001 + 0.998 * mu_raw;

  ratio ~ beta(mu * phi, (1 - mu) * phi);
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] mu_pred;
  {
    vector[N_obs] eta =
        h0 + trait_effect[sp_idx] + z_sp[sp_idx]
      + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
      + h1 * ln_ht + h2 * ln_dbh
      + h3 * bal_over_ht + h4 * sqrt_ba
      + h5 * ln_cspi_shift;
    mu_pred = 0.001 + 0.998 * inv_logit(eta);
  }
  for (i in 1:N_obs) {
    log_lik[i] = beta_lpdf(ratio[i] | mu_pred[i] * phi, (1 - mu_pred[i]) * phi);
  }
}
