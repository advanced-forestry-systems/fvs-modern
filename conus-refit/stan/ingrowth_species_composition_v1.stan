// ============================================================================
// ingrowth_species_composition_v1.stan
//
// Stage 4 of the four-stage ingrowth model: species composition. For each
// plot with at least one recruit, model the species mix of those recruits
// using multinomial likelihood with trait-driven covariate effects.
//
// Linear predictor (log-alpha) for plot p, species s:
//
//   log_alpha[p, s] = alpha_0[s]                            // species intercept
//                  + sum_c b_c * x_c[p]                     // shared covariate
//                  + sum_c (W_sp[s,:] %*% gamma_c) * x_c[p] // trait-modulated
//
// Probability of species s at plot p:  p[p,s] = softmax(log_alpha[p,:])[s]
// Likelihood:                          y[p,:] ~ multinomial(N_p, p[p,:])
//
// alpha_0 is hierarchical with mean from a trait projection (so trait
// effects partially share strength across species).
//
// Stage 1+2 (occurrence + conditional count) is in ingrowth_hurdle_v1.stan.
// Stage 3 (DBH dist of recruits) is a separate empirical/Weibull model.
// ============================================================================
data {
  int<lower=1> N_plots;          // plots with at least one recruit
  int<lower=1> N_sp;             // recruit species (top-N + OTHER bucket)
  int<lower=1> P_cov;            // plot covariates
  int<lower=0> P_trait;          // species traits

  array[N_plots] int<lower=1> total_recruits;        // sum across species per plot
  array[N_plots, N_sp] int<lower=0> y_count;         // recruit counts per (plot, species)

  matrix[N_plots, P_cov] X;      // standardized plot covariates
  matrix[N_sp, P_trait > 0 ? P_trait : 1] W_sp;     // standardized species traits
}
parameters {
  vector[N_sp] alpha_0_raw;                          // species intercepts (centered)
  real<lower=0> sigma_alpha;
  vector[P_trait] gamma_int;                         // trait-driven mean of alpha_0

  vector[P_cov] b;                                   // shared covariate slopes
  matrix[P_trait, P_cov] gamma_cov;                  // trait-modulated covariate slopes
}
transformed parameters {
  // Hierarchical species intercept: mean = W_sp * gamma_int, scale = sigma_alpha
  vector[N_sp] alpha_0;
  if (P_trait > 0) {
    alpha_0 = W_sp * gamma_int + sigma_alpha * alpha_0_raw;
  } else {
    alpha_0 = sigma_alpha * alpha_0_raw;
  }

  // Species-specific covariate slopes: beta[s, c] = b[c] + W_sp[s,:] %*% gamma_cov[:, c]
  // For efficiency, compute trait-modulated covariate contribution X * gamma_cov^T * W_sp^T
  // log_alpha[p, s] = alpha_0[s] + X[p,:] * b + X[p,:] * gamma_cov^T * W_sp[s,:]^T
  matrix[N_plots, N_sp] log_alpha;
  {
    vector[N_plots] X_b = X * b;                     // shared part [N_plots]
    matrix[N_plots, N_sp] trait_modul = (X * gamma_cov') * W_sp';  // [N_plots, N_sp]
    for (i in 1:N_plots) {
      for (s in 1:N_sp) {
        log_alpha[i, s] = alpha_0[s] + X_b[i] + trait_modul[i, s];
      }
    }
  }
}
model {
  // Priors
  alpha_0_raw ~ std_normal();
  sigma_alpha ~ normal(0, 1);
  gamma_int ~ normal(0, 0.5);
  b ~ normal(0, 0.5);
  to_vector(gamma_cov) ~ normal(0, 0.3);

  // Multinomial likelihood per plot
  for (i in 1:N_plots) {
    target += multinomial_logit_lpmf(y_count[i] | to_vector(log_alpha[i]));
  }
}
generated quantities {
  vector[N_plots] log_lik;
  for (i in 1:N_plots) {
    log_lik[i] = multinomial_logit_lpmf(y_count[i] | to_vector(log_alpha[i]));
  }
}
