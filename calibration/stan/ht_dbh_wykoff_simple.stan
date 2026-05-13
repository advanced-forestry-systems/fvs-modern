// ============================================================================
// ht_dbh_wykoff_simple.stan  (v2 — eta cap)
//
// FVS-CONUS static height-diameter model, Wykoff (1986) form with
// simplified hierarchy: species + EPA L1 random intercepts only.
//
// Motivation: ht_dbh_wykoff.stan with full L1/L2/L3 hierarchy produced
// extreme multimodality (rhat ~1.52, ESS ~7) across all parameters.
// sigma_L3 was bimodal (median 0.17 vs mean 1.93), creating ridges where
// L2/L3 random effects trade off against fixed effects. Dropping L2/L3
// eliminates this confounding while L1 (~12 ecoregions) captures the
// major PNW vs Southeast vs Lake States height-diameter differences.
//
// v2: added fmin(eta, 20) before exp() to prevent overflow during warmup.
//     exp(20) ~ 485M m — far beyond any real tree, so no posterior impact.
//
//   HT = 1.37 + exp(b0 + b1 / (DBH + 1))
//
// b0 carries the L1 + species RE with optional trait-informed prior,
// plus stand covariates. b1 is a single free scalar (no RE).
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;

  int<lower=0> P_trait;

  vector<lower=0>[N_obs] ht_obs;
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

  real<lower=0> s0;
  real<lower=-0.5, upper=1.0> s1;
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
  sigma_sp ~ normal(0, 0.5);
  sigma_L1 ~ normal(0, 0.3);

  s0 ~ normal(1.0, 0.5);
  s1 ~ normal(0.3, 0.2);

  vector[N_obs] eta =
      b0
    + z_sp[sp_idx] + z_L1[L1_idx]
    + a_bal  * bal
    + a_ba   * sqrt_ba
    + a_cspi * ln_cspi_shift
    + a_bard * ba_x_rd
    + a_blrd * bal_x_rd
    + b1 ./ (dbh + 1.0);

  // Cap eta to prevent exp() overflow during warmup exploration
  vector[N_obs] eta_safe;
  for (i in 1:N_obs) eta_safe[i] = fmin(eta[i], 20.0);
  vector[N_obs] mu_ht = 1.37 + exp(eta_safe);
  vector[N_obs] sigma_i;
  for (i in 1:N_obs) sigma_i[i] = s0 * pow(dbh[i] + 1.0, s1);

  ht_obs ~ normal(mu_ht, sigma_i);
}
generated quantities {
  vector[N_obs] mu_pred;
  vector[N_obs] log_lik;
  {
    vector[N_obs] eta =
        b0
      + z_sp[sp_idx] + z_L1[L1_idx]
      + a_bal*bal + a_ba*sqrt_ba + a_cspi*ln_cspi_shift
      + a_bard*ba_x_rd + a_blrd*bal_x_rd
      + b1 ./ (dbh + 1.0);
    vector[N_obs] eta_safe;
    for (i in 1:N_obs) eta_safe[i] = fmin(eta[i], 20.0);
    mu_pred = 1.37 + exp(eta_safe);
    for (i in 1:N_obs) {
      real si = s0 * pow(dbh[i] + 1.0, s1);
      log_lik[i] = normal_lpdf(ht_obs[i] | mu_pred[i], fmax(si, 1e-6));
    }
  }
}
