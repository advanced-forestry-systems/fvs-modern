// ============================================================================
// gompit_survival_speciesfree.stan
//
// Survival-FRAMED species-free tree survival model. Same species-free
// hierarchical architecture as gompit_mortality_speciesfree.stan but the
// linear predictor eta is on the SURVIVAL scale (higher eta = higher
// survival), matching Greg Johnson's framing in the CONUS mortality document.
//
// Modeling survival directly (rather than mortality hazard) is the user's
// preferred framing and is more interpretable: a positive coefficient means
// the covariate PROMOTES survival.
//
// Link (cloglog on the annual survival probability, with exposure offset for
// variable measurement interval T):
//   annual_hazard      = exp(-eta)               // higher eta -> lower hazard
//   P(survive T years) = exp(-exp(-eta) * T)
//   alive ~ Bernoulli(P_surv_T)
//
// This is mathematically identical to the mortality-hazard formulation with
// all linear-predictor signs flipped, so LOO and predictions match the
// mortality model exactly. The value is interpretive (survival-positive
// coefficients) plus the new predictors below.
//
// Predictors (informed by Greg's residual analysis):
//   eta_surv = b0 + trait_effect[sp]
//            + z_L1 + z_L2 + z_L3 + z_FT
//            + b1 DBH + b2 DBH^2
//            + b3 CR + b3b CR^2          (vigorous trees survive more; nonlinear)
//            + b4 ln_csi                 (better sites survive more)
//            + b5 BAL + b6 sqrt_ba_rd    (competition lowers survival)
//            + b7 cch + b7b cch^2        (crown closure at tree tip, U-shaped)
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=1> N_FT;
  int<lower=0> P_trait;

  array[N_obs] int<lower=0, upper=1> alive;   // survived (1) or died (0)
  vector<lower=0>[N_obs] T_years;

  vector[N_obs] dbh;
  vector[N_obs] dbh_sq;
  vector[N_obs] cr_init;
  vector[N_obs] ln_csi;
  vector[N_obs] bal_metric;
  vector[N_obs] sqrt_ba_rd;
  vector[N_obs] cch;

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;
  array[N_obs] int<lower=1, upper=N_FT> FT_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;
}
transformed data {
  real cch_mean = mean(cch);
  real cch_sd   = sd(cch) > 0 ? sd(cch) : 1.0;
  vector[N_obs] cch_z    = (cch - cch_mean) / cch_sd;
  vector[N_obs] cch_z_sq = cch_z .* cch_z;
  vector[N_obs] cr_sq    = cr_init .* cr_init;
}
parameters {
  real b0;
  real b1; real b2; real b3; real b4; real b5; real b6;
  real b3b;
  real b7; real b7b;

  vector[P_trait] gamma;

  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  vector[N_FT] z_FT_raw;
  vector[N_sp] z_sp_raw;            // NEW species RE

  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;
  real<lower=0> sigma_FT;
  real<lower=0> sigma_sp;           // NEW species RE scale
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
  vector[N_sp] z_sp = sigma_sp * z_sp_raw;
}
model {
  // Priors are on the SURVIVAL scale. Note the sign reversals relative to the
  // mortality model: covariates that raise mortality now LOWER survival.
  b0 ~ normal( 5.0, 2.0);     // high baseline annual survival
  b1 ~ normal( 0.0, 0.05);    // DBH
  b2 ~ normal( 0.0, 0.001);   // DBH^2
  b3 ~ normal( 1.0, 1.0);     // CR: vigorous trees survive more (positive now)
  b3b ~ normal( 0.0, 0.5);    // CR curvature
  b4 ~ normal( 0.3, 0.3);     // ln_csi: better sites survive more (positive now)
  b5 ~ normal(-0.01, 0.05);   // BAL: competition lowers survival (negative now)
  b6 ~ normal(-0.1, 0.3);     // density: lowers survival (negative now)
  b7 ~ normal( 0.0, 0.5);     // cch linear
  b7b ~ normal( 0.0, 0.5);    // cch quadratic (U-shape)

  gamma ~ normal(0, 0.5);

  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  z_FT_raw ~ std_normal();
  z_sp_raw ~ std_normal();
  sigma_L1 ~ normal(0, 0.5);
  sigma_L2 ~ normal(0, 0.3);
  sigma_L3 ~ normal(0, 0.3);
  sigma_FT ~ normal(0, 0.3);
  sigma_sp ~ normal(0, 0.5);

  vector[N_obs] eta =
      b0 + trait_effect[sp_idx] + z_sp[sp_idx]
    + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
    + b1 * dbh + b2 * dbh_sq
    + b3 * cr_init + b3b * cr_sq
    + b4 * ln_csi
    + b5 * bal_metric + b6 * sqrt_ba_rd
    + b7 * cch_z + b7b * cch_z_sq;

  vector[N_obs] eta_safe;
  for (i in 1:N_obs) eta_safe[i] = fmin(fmax(eta[i], -5.0), 20.0);

  // Survival framing: annual hazard = exp(-eta), survival over T = exp(-hazard*T)
  for (i in 1:N_obs) {
    real log_p_surv = -exp(-eta_safe[i]) * T_years[i];
    real log_p_mort = log1m_exp(log_p_surv);
    target += alive[i] == 1 ? log_p_surv : log_p_mort;
  }
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] p_surv_annual;
  {
    vector[N_obs] eta =
        b0 + trait_effect[sp_idx] + z_sp[sp_idx]
      + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
      + b1 * dbh + b2 * dbh_sq
      + b3 * cr_init + b3b * cr_sq
      + b4 * ln_csi
      + b5 * bal_metric + b6 * sqrt_ba_rd
      + b7 * cch_z + b7b * cch_z_sq;
    for (i in 1:N_obs) {
      real eta_s = fmin(fmax(eta[i], -5.0), 20.0);
      real log_p_surv = -exp(-eta_s) * T_years[i];
      real log_p_mort = log1m_exp(log_p_surv);
      log_lik[i] = alive[i] == 1 ? log_p_surv : log_p_mort;
      p_surv_annual[i] = exp(-exp(-eta_s));
    }
  }
}
