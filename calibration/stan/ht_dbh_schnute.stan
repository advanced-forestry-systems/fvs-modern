// ============================================================================
// ht_dbh_schnute.stan
//
// FVS-CONUS static height-diameter model, Schnute (1981) generalized form
// restricted to the two-parameter (a, b) case used by Kershaw et al. (2008)
// for H-D:
//
//   HT = 1.37 + (HT_ref - 1.37) * [(1 - exp(-a*DBH)) / (1 - exp(-a*DBH_ref))]^b
//
// HT_ref is a free positive scalar (asymptotic reference), DBH_ref is a
// user-supplied anchor (default 25 cm). The flexible a, b pair lets this
// form mimic both Chapman-Richards (b > 1) and saturating (b < 1)
// behavior.
//
// Random effects on log(HT_ref) only: nested EPA + species with optional
// trait-informed prior.  Stand covariates modulate HT_ref just as the
// asymptote in the Chapman form.
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
  real<lower=1> dbh_ref;       // anchor DBH (cm), e.g. 25

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
  real a0;                     // log(HT_ref - 1.37) intercept; HT_ref ~ 20 m -> log(18.6) ~ 2.9
  real<lower=0.001> a_rate;    // Schnute a  (>0)
  real<lower=0.2, upper=4.0> b_shape; // Schnute b

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
  a0      ~ normal(2.9, 0.7);
  a_rate  ~ normal(0.05, 0.03);
  b_shape ~ normal(1.1, 0.4);

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

  vector[N_obs] log_ref =
      a0
    + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + a_bal  * bal
    + a_ba   * sqrt_ba
    + a_cspi * ln_cspi_shift
    + a_bard * ba_x_rd
    + a_blrd * bal_x_rd;

  vector[N_obs] href = exp(log_ref);             // HT_ref - 1.37 (m)
  real denom = 1.0 - exp(-a_rate * dbh_ref);
  if (denom < 1e-8) denom = 1e-8;

  vector[N_obs] mu_ht;
  for (i in 1:N_obs) {
    real num = 1.0 - exp(-a_rate * dbh[i]);
    if (num < 1e-8) num = 1e-8;
    mu_ht[i] = 1.37 + href[i] * pow(num / denom, b_shape);
  }

  vector[N_obs] sigma_i;
  for (i in 1:N_obs) sigma_i[i] = s0 * pow(dbh[i] + 1.0, s1);

  ht_obs ~ normal(mu_ht, sigma_i);
}
generated quantities {
  vector[N_obs] mu_pred;
  vector[N_obs] log_lik;
  {
    vector[N_obs] log_ref =
        a0
      + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
      + a_bal*bal + a_ba*sqrt_ba + a_cspi*ln_cspi_shift
      + a_bard*ba_x_rd + a_blrd*bal_x_rd;
    vector[N_obs] href = exp(log_ref);
    real denom = 1.0 - exp(-a_rate * dbh_ref);
    if (denom < 1e-8) denom = 1e-8;
    for (i in 1:N_obs) {
      real num = 1.0 - exp(-a_rate * dbh[i]);
      if (num < 1e-8) num = 1e-8;
      mu_pred[i] = 1.37 + href[i] * pow(num / denom, b_shape);
      real si = s0 * pow(dbh[i] + 1.0, s1);
      log_lik[i] = normal_lpdf(ht_obs[i] | mu_pred[i], si);
    }
  }
}
