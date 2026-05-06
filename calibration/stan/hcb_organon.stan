// ============================================================================
// hcb_organon.stan
//
// FVS-CONUS height to crown base, static (cross-sectional) logistic model
// with nested EPA L1 / L2 / L3 + species random intercepts to match DG / HG.
//
// ORGANON formulation:
//   HCB = HT / (1 + exp(eta))
//   HCB / HT is on (0, 1); eta linear in predictors.
//
// Linear predictor for logit(HCB / HT):
//   eta = h0 + z_sp + z_L1 + z_L2 + z_L3
//       + h1 ln(HT)
//       + h2 ln(DBH)
//       + h3 BAL / (HT + 1)
//       + h4 sqrt(BA)
//       + h5 ln(CSPI_shift)
//
// Fitting: binomial-style on the HCB / HT ratio is awkward because the
// ratio is real on (0, 1) not a count. We use a Beta likelihood on the
// ratio with a shared precision phi; the cross-sectional form has no
// annualization.
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;

  int<lower=0> P_trait;                   // species-level trait columns (0 = disable prior)

  vector<lower=0, upper=1>[N_obs] ratio;  // HCB / HT
  vector[N_obs] ln_ht;
  vector[N_obs] ln_dbh;
  vector[N_obs] bal_over_ht;
  vector[N_obs] sqrt_ba;
  vector[N_obs] ln_cspi_shift;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;  // species x trait design (centered, scaled)
}
parameters {
  real h0;
  real h1;
  real h2;
  real h3;
  real h4;
  real h5;

  vector[P_trait] gamma;                // trait loadings mapping W -> species RE mean

  vector[N_sp] z_sp_raw;
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  real<lower=0> sigma_sp;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;

  real<lower=1> phi;     // Beta precision
}
transformed parameters {
  vector[N_sp] mu_sp;
  if (P_trait > 0) {
    mu_sp = W * gamma;
  } else {
    mu_sp = rep_vector(0.0, N_sp);
  }
  vector[N_sp] z_sp = mu_sp + sigma_sp * z_sp_raw;
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
  vector[N_L2] z_L2 = sigma_L2 * z_L2_raw;
  vector[N_L3] z_L3 = sigma_L3 * z_L3_raw;
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
  sigma_sp ~ normal(0, 0.5);
  sigma_L1 ~ normal(0, 0.5);
  sigma_L2 ~ normal(0, 0.3);
  sigma_L3 ~ normal(0, 0.3);

  phi ~ gamma(2, 0.05);   // mean 40, generous

  vector[N_obs] eta =
      h0
    + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + h1 * ln_ht
    + h2 * ln_dbh
    + h3 * bal_over_ht
    + h4 * sqrt_ba
    + h5 * ln_cspi_shift;

  // Beta mean parameterization: mu = inv_logit(eta); shape1 = mu*phi, shape2 = (1-mu)*phi
  // Clamp mu away from 0/1 to prevent beta shape parameters hitting zero,
  // which causes chains to get stuck with beta_lpdf exceptions at full data scale.
  vector[N_obs] mu_raw = inv_logit(eta);
  vector[N_obs] mu = 0.001 + 0.998 * mu_raw;

  ratio ~ beta(mu * phi, (1 - mu) * phi);
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] mu_pred;
  {
    vector[N_obs] eta =
        h0
      + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
      + h1 * ln_ht + h2 * ln_dbh + h3 * bal_over_ht
      + h4 * sqrt_ba + h5 * ln_cspi_shift;
    mu_pred = 0.001 + 0.998 * inv_logit(eta);
  }
  for (i in 1:N_obs) {
    log_lik[i] = beta_lpdf(ratio[i] | mu_pred[i] * phi, (1 - mu_pred[i]) * phi);
  }
}
