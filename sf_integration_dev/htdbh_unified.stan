// ============================================================================
// htdbh_unified.stan
//
// Unified height-diameter model: species-free B1 Wykoff log-normal hierarchy
// (traits W*gamma + nested L1/L2/L3 ecoregion + forest-type REs) PLUS a
// per-species random intercept z_sp that shrinks toward the trait-predicted
// mean (half-normal sigma_sp). Every species gets an estimate; no pooling.
//
//   log(HT - 1.37) = b0 + trait_effect[sp] + z_sp[sp]
//                  + z_L1 + z_L2 + z_L3 + z_FT
//                  + a_bal*BAL + a_ba*sqrt(BA) + a_cspi*ln(CSPI_shift)
//                  + a_bard*(BA*RD) + a_blrd*(BAL*RD) + b1/(DBH+1)
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=1> N_FT;
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
  array[N_obs] int<lower=1, upper=N_FT> FT_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
transformed data {
  vector[N_obs] log_ht_above_bh;
  for (i in 1:N_obs)
    log_ht_above_bh[i] = log(fmax(ht_obs[i] - 1.37, 0.01));
}
parameters {
  real b0;
  real<upper=0> b1;

  real a_bal;
  real a_ba;
  real a_cspi;
  real a_bard;
  real a_blrd;

  vector[P_trait] gamma;

  vector[N_sp] z_sp_raw;            // NEW species random intercept
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_FT] z_FT_raw;

  real<lower=0> sigma_sp;           // NEW species RE scale
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
  real<lower=0> sigma_FT;

  real<lower=0> sigma;
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
  z_FT_raw ~ std_normal();
  sigma_sp ~ normal(0, 0.5);
  sigma_L1 ~ normal(0, 0.5);
  sigma_L2 ~ normal(0, 0.3);
  sigma_L3 ~ normal(0, 0.3);
  sigma_FT ~ normal(0, 0.3);

  sigma ~ normal(0, 0.5);

  vector[N_obs] eta =
      b0
    + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + a_bal  * bal
    + a_ba   * sqrt_ba
    + a_cspi * ln_cspi_shift
    + a_bard * ba_x_rd
    + a_blrd * bal_x_rd
    + b1 ./ (dbh + 1.0);

  log_ht_above_bh ~ normal(eta, sigma);
}
generated quantities {
  vector[N_obs] mu_pred;
  vector[N_obs] log_lik;
  {
    vector[N_obs] eta =
        b0
      + trait_effect[sp_idx] + z_sp[sp_idx]
      + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
      + a_bal*bal + a_ba*sqrt_ba + a_cspi*ln_cspi_shift
      + a_bard*ba_x_rd + a_blrd*bal_x_rd
      + b1 ./ (dbh + 1.0);
    for (i in 1:N_obs)
      mu_pred[i] = 1.37 + exp(eta[i] + 0.5 * sigma * sigma);
    for (i in 1:N_obs)
      log_lik[i] = normal_lpdf(log_ht_above_bh[i] | eta[i], sigma);
  }
}
