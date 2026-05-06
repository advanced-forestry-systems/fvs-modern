// ============================================================================
// ht_dbh_wykoff.stan
//
// FVS-CONUS static height-diameter model, Wykoff (1986) form used as the
// FVS default across most variants:
//
//   HT = 1.37 + exp(b0 + b1 / (DBH + 1))
//
// b0 carries the nested EPA L1 / L2 / L3 + species RE with optional
// trait-informed prior, plus stand covariates (BAL, sqrt_BA, CSPI,
// BA . RD, BAL . RD).  b1 is a single free scalar (no RE) to keep
// identifiability of the curvature term.
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;

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
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;

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
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  real<lower=0> sigma_sp;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;

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
  vector[N_L2] z_L2 = sigma_L2 * z_L2_raw;
  vector[N_L3] z_L3 = sigma_L3 * z_L3_raw;
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
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  sigma_sp ~ normal(0, 0.5);
  sigma_L1 ~ normal(0, 0.3);
  sigma_L2 ~ normal(0, 0.2);
  sigma_L3 ~ normal(0, 0.2);

  s0 ~ normal(1.0, 0.5);
  s1 ~ normal(0.3, 0.2);

  vector[N_obs] eta =
      b0
    + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + a_bal  * bal
    + a_ba   * sqrt_ba
    + a_cspi * ln_cspi_shift
    + a_bard * ba_x_rd
    + a_blrd * bal_x_rd
    + b1 ./ (dbh + 1.0);

  vector[N_obs] mu_ht = 1.37 + exp(eta);
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
      + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
      + a_bal*bal + a_ba*sqrt_ba + a_cspi*ln_cspi_shift
      + a_bard*ba_x_rd + a_blrd*bal_x_rd
      + b1 ./ (dbh + 1.0);
    mu_pred = 1.37 + exp(eta);
    for (i in 1:N_obs) {
      real si = s0 * pow(dbh[i] + 1.0, s1);
      log_lik[i] = normal_lpdf(ht_obs[i] | mu_pred[i], si);
    }
  }
}
