// ============================================================================
// ht_dbh_chapman_speciesfree.stan
//
// FVS-CONUS static height-diameter, Chapman-Richards form, SPECIES-FREE (B1).
// Built from ht_dbh_chapman.stan by replacing the (trait-informed) species
// random intercept z_sp with a pure trait projection trait_effect = W * gamma,
// matching the B1 architecture of the other species-free components (DG, HG,
// CR, HCB, survival): all species variation is carried by traits, with nested
// EPA L1/L2/L3 random effects on the log asymptote. No per-species free
// deviation (no z_sp_raw / sigma_sp).
//
//   HT = 1.37 + A * (1 - exp(-b_rate * DBH))^c_shape
//   log(A) = a0 + trait_effect[sp] + z_L1 + z_L2 + z_L3
//          + a_bal*BAL + a_ba*sqrt(BA) + a_cspi*ln(CSPI_shift)
//          + a_bard*(BA*RD) + a_blrd*(BAL*RD)
//   trait_effect = W * gamma
//
// Monotone non-decreasing in DBH by construction (A>0, 0<term<1, c_shape>0);
// asymptotes at 1.37 + A, which tames large-tree over-tall behaviour relative
// to the Wykoff reciprocal form. Heteroscedastic SD sigma = s0*(DBH+1)^s1.
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;

  int<lower=0> P_trait;

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

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
parameters {
  real a0;                     // log asymptote intercept
  real<lower=0> b_rate;        // Chapman rate
  real<lower=0.6, upper=1.6> c_shape; // Chapman shape

  real a_bal;
  real a_ba;
  real a_cspi;
  real a_bard;
  real a_blrd;

  vector[P_trait] gamma;       // trait coefficients (carry all species variation)

  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;

  real<lower=0> s0;            // residual SD scale
  real<lower=-0.5, upper=1.0> s1; // DBH power on residual SD
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
}
model {
  a0      ~ normal(3.2, 0.8);
  b_rate  ~ normal(0.04, 0.015);
  c_shape ~ normal(1.0, 0.1);

  a_bal  ~ normal( 0.0, 0.02);
  a_ba   ~ normal( 0.0, 0.05);
  a_cspi ~ normal( 0.2, 0.3);
  a_bard ~ normal( 0.0, 0.02);
  a_blrd ~ normal( 0.0, 0.02);

  gamma ~ normal(0, 0.5);
  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  sigma_L1 ~ normal(0, 0.5);
  sigma_L2 ~ normal(0, 0.3);
  sigma_L3 ~ normal(0, 0.2);

  s0 ~ normal(1.0, 0.5);
  s1 ~ normal(0.3, 0.2);

  vector[N_obs] log_A =
      a0
    + trait_effect[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + a_bal  * bal
    + a_ba   * sqrt_ba
    + a_cspi * ln_cspi_shift
    + a_bard * ba_x_rd
    + a_blrd * bal_x_rd;

  vector[N_obs] A;
  for (gi_ in 1:N_obs) A[gi_] = exp(fmin(log_A[gi_], 4.2));  // cap log-asymptote to prevent inf overflow
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
      + trait_effect[sp_idx]
      + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
      + a_bal*bal + a_ba*sqrt_ba + a_cspi*ln_cspi_shift
      + a_bard*ba_x_rd + a_blrd*bal_x_rd;
    vector[N_obs] A;
  for (gi_ in 1:N_obs) A[gi_] = exp(fmin(log_A[gi_], 4.2));  // cap log-asymptote to prevent inf overflow
    for (i in 1:N_obs) {
      real term = 1.0 - exp(-b_rate * dbh[i]);
      if (term < 1e-8) term = 1e-8;
      mu_pred[i] = 1.37 + A[i] * pow(term, c_shape);
      real si = s0 * pow(dbh[i] + 1.0, s1);
      log_lik[i] = normal_lpdf(ht_obs[i] | mu_pred[i], si);
    }
  }
}
