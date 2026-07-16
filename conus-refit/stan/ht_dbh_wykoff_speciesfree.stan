// ============================================================================
// ht_dbh_wykoff_speciesfree.stan
//
// Species-free B1 architecture for the static height-diameter relationship,
// Wykoff (1986) form with a log-normal likelihood. Built from
// ht_dbh_wykoff_simple_lognormal.stan (species-specific) by replacing the
// species random intercept with the species-free trait projection and the
// nested ecoregion + forest-type random effects, matching the other B1
// components (CR, HCB, mortality).
//
//   log(HT - 1.37) = b0 + trait_effect[sp]
//                  + z_L1 + z_L2 + z_L3 + z_FT
//                  + a_bal*BAL + a_ba*sqrt(BA) + a_cspi*ln(CSPI_shift)
//                  + a_bard*(BA*RD) + a_blrd*(BAL*RD)
//                  + b1 / (DBH + 1)
//
//   trait_effect = W * gamma   (carries all species variation; no z_sp)
//   Back-transform: HT = 1.37 + exp(eta + 0.5*sigma^2)
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=1> N_FT;
  int<lower=0> P_trait;

  vector<lower=0>[N_obs] ht_obs;     // total height (m), must be > 1.37
  vector<lower=0>[N_obs] dbh;        // DBH (cm)

  vector[N_obs] bal;                 // BAL metric (BAL_SW + BAL_HW)
  vector[N_obs] sqrt_ba;             // sqrt(BA metric)
  vector[N_obs] ln_cspi_shift;       // ln(climate site index, shifted)
  vector[N_obs] ba_x_rd;             // BA metric * relative density
  vector[N_obs] bal_x_rd;            // BAL metric * relative density

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
  real<upper=0> b1;             // classical Wykoff has b1 < 0

  real a_bal;
  real a_ba;
  real a_cspi;
  real a_bard;
  real a_blrd;

  vector[P_trait] gamma;

  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_FT] z_FT_raw;

  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
  real<lower=0> sigma_FT;

  real<lower=0> sigma;           // residual SD on log scale
}
transformed parameters {
  vector[N_sp] trait_effect;
  if (P_trait > 0) {
    trait_effect = W * gamma;
  } else {
    trait_effect = rep_vector(0.0, N_sp);
  }
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

  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  z_FT_raw ~ std_normal();
  sigma_L1 ~ normal(0, 0.5);
  sigma_L2 ~ normal(0, 0.3);
  sigma_L3 ~ normal(0, 0.3);
  sigma_FT ~ normal(0, 0.3);

  sigma ~ normal(0, 0.5);

  vector[N_obs] eta =
      b0
    + trait_effect[sp_idx]
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
      + trait_effect[sp_idx]
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
