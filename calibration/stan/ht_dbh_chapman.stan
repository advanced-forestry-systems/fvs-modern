// ============================================================================
// ht_dbh_chapman.stan
//
// FVS-CONUS static height-diameter model, Chapman-Richards form.
//
//   HT = 1.37 + A * (1 - exp(-b * DBH))^c
//
// A, b, c are each log-/logit-transformed linear predictors with
//   (i)  nested EPA L1 / L2 / L3 + species random intercepts on A
//   (ii) optional trait-informed species prior on the species RE
//            z_sp ~ Normal( W . traits_sp, sigma_sp )
//   (iii) stand-level covariates on log(A) to let competition and
//         productivity modulate the asymptote:
//              log(A) = a0 + z_sp + z_L1 + z_L2 + z_L3
//                     + a_bal  * bal
//                     + a_ba   * sqrt_ba
//                     + a_cspi * ln_cspi_shift
//                     + a_bard * ba_x_rd
//                     + a_blrd * bal_x_rd
//
// b and c are held as single free scalars (no RE) for identifiability,
// matching Weiskittel, Hann, Kershaw, Vanclay (2011, Ch. 10).
//
// Likelihood: Normal on (HT - 1.37), heteroscedastic SD proportional to
// a power of DBH (sigma = s0 * DBH^s1) after Hann-Hanus (2002).
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;

  int<lower=0> P_trait;                 // number of species-level trait columns (0 = disable)

  vector<lower=0>[N_obs] ht_obs;        // total height, m
  vector<lower=0>[N_obs] dbh;           // cm

  vector[N_obs] bal;
  vector[N_obs] sqrt_ba;
  vector[N_obs] ln_cspi_shift;
  vector[N_obs] ba_x_rd;
  vector[N_obs] bal_x_rd;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W; // species trait design matrix (centered, scaled)
}
parameters {
  real a0;                     // log asymptote intercept
  real<lower=0> b_rate;        // Chapman rate
  real<lower=0.3, upper=4.0> c_shape; // Chapman shape

  real a_bal;
  real a_ba;
  real a_cspi;
  real a_bard;
  real a_blrd;

  vector[P_trait] gamma;       // trait coefficients mapping traits -> species RE mean
  vector[N_sp] z_sp_raw;
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  real<lower=0> sigma_sp;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;

  real<lower=0> s0;            // residual SD scale
  real<lower=-0.5, upper=1.0> s1; // DBH power on residual SD
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
  a0      ~ normal(3.2, 0.8);       // exp(3.2) ~ 24 m plausible asymptote
  b_rate  ~ normal(0.04, 0.03);
  c_shape ~ normal(1.1, 0.4);

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

  vector[N_obs] log_A =
      a0
    + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + a_bal  * bal
    + a_ba   * sqrt_ba
    + a_cspi * ln_cspi_shift
    + a_bard * ba_x_rd
    + a_blrd * bal_x_rd;

  vector[N_obs] A = exp(log_A);
  vector[N_obs] mu_ht;
  for (i in 1:N_obs) {
    real term = 1.0 - exp(-b_rate * dbh[i]);
    if (term < 1e-8) term = 1e-8;
    mu_ht[i] = 1.37 + A[i] * pow(term, c_shape);
  }

  vector[N_obs] sigma_i;
  for (i in 1:N_obs) sigma_i[i] = s0 * pow(dbh[i] + 1.0, s1);

  ht_obs ~ normal(mu_ht, sigma_i);
}
generated quantities {
  vector[N_obs] mu_pred;
  vector[N_obs] log_lik;
  {
    vector[N_obs] log_A =
        a0
      + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
      + a_bal*bal + a_ba*sqrt_ba + a_cspi*ln_cspi_shift
      + a_bard*ba_x_rd + a_blrd*bal_x_rd;
    vector[N_obs] A = exp(log_A);
    for (i in 1:N_obs) {
      real term = 1.0 - exp(-b_rate * dbh[i]);
      if (term < 1e-8) term = 1e-8;
      mu_pred[i] = 1.37 + A[i] * pow(term, c_shape);
      real si = s0 * pow(dbh[i] + 1.0, s1);
      log_lik[i] = normal_lpdf(ht_obs[i] | mu_pred[i], si);
    }
  }
}
