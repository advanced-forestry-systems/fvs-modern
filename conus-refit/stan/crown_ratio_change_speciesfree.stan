// ============================================================================
// crown_ratio_change_speciesfree.stan
//
// Species-free B1 architecture for annual crown ratio change. Mirrors the
// DG_Kuehne v4_full hierarchical pattern that emerged from the May 2026
// LOO comparison: trait_effect = W * gamma (no z_sp random intercept),
// nested EPA L1/L2/L3 ecoregion REs, plus forest type RE (FORTYPCD).
//
// Linear predictor:
//   delta_CR = b0 + trait_effect[sp]
//            + z_L1 + z_L2 + z_L3 + z_FT
//            + b1 DBH + b2 DBH^2
//            + b3 BA + b4 BAL + b5 CR_init
//            + b6 ln_csi
//
// Notes:
//   - Crown ratio change is treated as already annualized (delta_CR / years).
//   - CR_init regression-to-mean is preserved as a main effect (b5).
//   - ln_csi is the climate-derived site index. Will be replaceable with
//     bgi piecewise or mapdd5 once the DG winning architecture is fixed.
//   - No species-specific intercept; trait_effect carries all species
//     variation through W * gamma.
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=1> N_FT;
  int<lower=0> P_trait;

  vector[N_obs] delta_CR_a;        // annualized change in crown ratio
  vector[N_obs] dbh;
  vector[N_obs] dbh_sq;
  vector[N_obs] ba_metric;
  vector[N_obs] bal_metric;
  vector[N_obs] cr_init;
  vector[N_obs] ln_csi;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;
  array[N_obs] int<lower=1, upper=N_FT> FT_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
parameters {
  real b0;
  real b1; real b2; real b3; real b4; real b5; real b6;

  vector[P_trait] gamma;

  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_FT] z_FT_raw;

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
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
  vector[N_L2] z_L2 = sigma_L2 * z_L2_raw;
  vector[N_L3] z_L3 = sigma_L3 * z_L3_raw;
  vector[N_FT] z_FT = sigma_FT * z_FT_raw;
}
model {
  b0 ~ normal( 0.0, 1.0);     // intercept (looser for divergence-prone smoke)
  b1 ~ normal( 0.0, 0.1);
  b2 ~ normal( 0.0, 0.01);
  b3 ~ normal( 0.0, 0.05);
  b4 ~ normal( 0.0, 0.05);
  b5 ~ normal(-0.1, 0.5);
  b6 ~ normal( 0.0, 0.5);

  gamma ~ normal(0, 0.5);

  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  z_FT_raw ~ std_normal();
  sigma_L1 ~ normal(0, 0.3);
  sigma_L2 ~ normal(0, 0.2);
  sigma_L3 ~ normal(0, 0.2);
  sigma_FT ~ normal(0, 0.2);

  sigma ~ normal(0, 0.5);

  vector[N_obs] mu =
      b0 + trait_effect[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + b1 * dbh + b2 * dbh_sq
    + b3 * ba_metric + b4 * bal_metric
    + b5 * cr_init
    + b6 * ln_csi;

  delta_CR_a ~ normal(mu, sigma);
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] mu_pred =
      b0 + trait_effect[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + b1 * dbh + b2 * dbh_sq
    + b3 * ba_metric + b4 * bal_metric
    + b5 * cr_init
    + b6 * ln_csi;
  for (i in 1:N_obs) {
    log_lik[i] = normal_lpdf(delta_CR_a[i] | mu_pred[i], sigma);
  }
}
