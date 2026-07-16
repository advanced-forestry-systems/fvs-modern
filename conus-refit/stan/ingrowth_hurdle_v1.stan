// ============================================================================
// ingrowth_hurdle_v1.stan
//
// Stage 1 + Stage 2 of the four-stage ingrowth model: hurdle negative
// binomial. Models occurrence and conditional count separately, with shared
// covariate structure but independent posteriors.
//
//   Stage 1: P(n_recruits > 0 | plot) = logit^-1(eta_p)
//   Stage 2: n_recruits | n > 0 ~ NegBin2(mu_c, phi), truncated at 0
//            log(mu_c) = eta_c + log(years)
//
// Both eta_p and eta_c use the same plot covariates and dominant-species
// trait vector W; coefficients are independent (subscripts _p and _c).
//
// Stages 3 (DBH distribution) and 4 (species composition) are separate
// models — see ingrowth_species_composition_v1.stan for stage 4.
// ============================================================================
data {
  int<lower=1> N_plots;
  int<lower=1> N_sp;             // dominant species lookup
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;
  int<lower=0> P_trait;

  array[N_plots] int<lower=0> n_recruits;
  vector[N_plots] log_years;     // measurement interval offset (count side only)

  vector[N_plots] ln_ba;
  vector[N_plots] ln_bal;
  vector[N_plots] rd;
  vector[N_plots] ht40;
  vector[N_plots] ln_csi;
  vector[N_plots] clim_pca1;

  array[N_plots] int<lower=1, upper=N_sp> dom_sp_idx;
  array[N_plots] int<lower=1, upper=N_L1> L1_idx;
  array[N_plots] int<lower=1, upper=N_L2> L2_idx;
  array[N_plots] int<lower=1, upper=N_L3> L3_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W_dom;
}
parameters {
  // --- Occurrence (stage 1) ---
  real a0_p;
  vector[P_trait] gamma_p;
  real b1_p; real b2_p; real b3_p; real b4_p; real b5_p; real b6_p;
  vector[N_L1] z_L1_p_raw;
  vector[N_L2] z_L2_p_raw;
  vector[N_L3] z_L3_p_raw;
  real<lower=0> sigma_L1_p;
  real<lower=0> sigma_L2_p;
  real<lower=0> sigma_L3_p;

  // --- Count given occurrence (stage 2) ---
  real a0_c;
  vector[P_trait] gamma_c;
  real b1_c; real b2_c; real b3_c; real b4_c; real b5_c; real b6_c;
  vector[N_L1] z_L1_c_raw;
  vector[N_L2] z_L2_c_raw;
  vector[N_L3] z_L3_c_raw;
  real<lower=0> sigma_L1_c;
  real<lower=0> sigma_L2_c;
  real<lower=0> sigma_L3_c;

  real<lower=0> phi;             // NB dispersion (count side)
}
transformed parameters {
  vector[N_sp] trait_p;
  vector[N_sp] trait_c;
  if (P_trait > 0) {
    trait_p = W_dom * gamma_p;
    trait_c = W_dom * gamma_c;
  } else {
    trait_p = rep_vector(0.0, N_sp);
    trait_c = rep_vector(0.0, N_sp);
  }
  vector[N_L1] z_L1_p = sigma_L1_p * z_L1_p_raw;
  vector[N_L2] z_L2_p = sigma_L2_p * z_L2_p_raw;
  vector[N_L3] z_L3_p = sigma_L3_p * z_L3_p_raw;
  vector[N_L1] z_L1_c = sigma_L1_c * z_L1_c_raw;
  vector[N_L2] z_L2_c = sigma_L2_c * z_L2_c_raw;
  vector[N_L3] z_L3_c = sigma_L3_c * z_L3_c_raw;
}
model {
  // Priors
  a0_p ~ normal(0.0, 2.0);
  a0_c ~ normal(0.0, 2.0);
  gamma_p ~ normal(0, 0.5);
  gamma_c ~ normal(0, 0.5);

  b1_p ~ normal(-0.2, 0.5); b1_c ~ normal(-0.3, 0.5);   // ln_ba
  b2_p ~ normal(-0.2, 0.5); b2_c ~ normal(-0.2, 0.5);   // ln_bal
  b3_p ~ normal(-0.3, 1.0); b3_c ~ normal(-0.5, 1.0);   // rd
  b4_p ~ normal( 0.0, 0.5); b4_c ~ normal( 0.0, 0.5);   // ht40
  b5_p ~ normal( 0.3, 0.5); b5_c ~ normal( 0.3, 0.5);   // ln_csi (site, expect +)
  b6_p ~ normal( 0.0, 0.5); b6_c ~ normal( 0.0, 0.5);   // clim_pca1

  z_L1_p_raw ~ std_normal(); z_L2_p_raw ~ std_normal(); z_L3_p_raw ~ std_normal();
  z_L1_c_raw ~ std_normal(); z_L2_c_raw ~ std_normal(); z_L3_c_raw ~ std_normal();
  sigma_L1_p ~ normal(0, 1.0); sigma_L1_c ~ normal(0, 1.0);
  sigma_L2_p ~ normal(0, 0.5); sigma_L2_c ~ normal(0, 0.5);
  sigma_L3_p ~ normal(0, 0.3); sigma_L3_c ~ normal(0, 0.3);
  phi ~ gamma(2.0, 0.5);

  // Linear predictors
  vector[N_plots] eta_p =
      a0_p + trait_p[dom_sp_idx]
    + z_L1_p[L1_idx] + z_L2_p[L2_idx] + z_L3_p[L3_idx]
    + b1_p * ln_ba + b2_p * ln_bal + b3_p * rd
    + b4_p * ht40 + b5_p * ln_csi + b6_p * clim_pca1;

  vector[N_plots] eta_c =
      a0_c + trait_c[dom_sp_idx]
    + z_L1_c[L1_idx] + z_L2_c[L2_idx] + z_L3_c[L3_idx]
    + b1_c * ln_ba + b2_c * ln_bal + b3_c * rd
    + b4_c * ht40 + b5_c * ln_csi + b6_c * clim_pca1
    + log_years;

  // Bound eta_c to avoid numerical pathologies
  vector[N_plots] eta_c_safe;
  for (i in 1:N_plots) eta_c_safe[i] = fmin(fmax(eta_c[i], -20.0), 12.0);

  // Hurdle likelihood
  for (i in 1:N_plots) {
    if (n_recruits[i] == 0) {
      target += bernoulli_logit_lpmf(0 | eta_p[i]);
    } else {
      target += bernoulli_logit_lpmf(1 | eta_p[i])
              + neg_binomial_2_log_lpmf(n_recruits[i] | eta_c_safe[i], phi)
              - log1m_exp(neg_binomial_2_log_lpmf(0 | eta_c_safe[i], phi));
    }
  }
}
generated quantities {
  vector[N_plots] log_lik;
  vector[N_plots] eta_p_gq =
      a0_p + trait_p[dom_sp_idx]
    + z_L1_p[L1_idx] + z_L2_p[L2_idx] + z_L3_p[L3_idx]
    + b1_p * ln_ba + b2_p * ln_bal + b3_p * rd
    + b4_p * ht40 + b5_p * ln_csi + b6_p * clim_pca1;
  vector[N_plots] eta_c_gq =
      a0_c + trait_c[dom_sp_idx]
    + z_L1_c[L1_idx] + z_L2_c[L2_idx] + z_L3_c[L3_idx]
    + b1_c * ln_ba + b2_c * ln_bal + b3_c * rd
    + b4_c * ht40 + b5_c * ln_csi + b6_c * clim_pca1
    + log_years;
  for (i in 1:N_plots) {
    real ec = fmin(fmax(eta_c_gq[i], -20.0), 12.0);
    if (n_recruits[i] == 0) {
      log_lik[i] = bernoulli_logit_lpmf(0 | eta_p_gq[i]);
    } else {
      log_lik[i] = bernoulli_logit_lpmf(1 | eta_p_gq[i])
                 + neg_binomial_2_log_lpmf(n_recruits[i] | ec, phi)
                 - log1m_exp(neg_binomial_2_log_lpmf(0 | ec, phi));
    }
  }
}
