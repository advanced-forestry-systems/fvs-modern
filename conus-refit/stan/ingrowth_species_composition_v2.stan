// ingrowth_species_composition_v2.stan
// Memory-fixed variant of v1: identical model, but the huge N_plots x N_sp
// log_alpha matrix is now a LOCAL in model{} and generated quantities{} rather
// than a transformed parameter, so it is NOT stored in every draw. v1 OOM'd in
// save_object because log_alpha (N_plots*N_sp per draw) ballooned the draws
// array. Only the actual params (alpha_0, b, gamma_*, sigma_alpha) are saved now.
data {
  int<lower=1> N_plots; int<lower=1> N_sp; int<lower=1> P_cov; int<lower=0> P_trait;
  array[N_plots] int<lower=1> total_recruits;
  array[N_plots, N_sp] int<lower=0> y_count;
  matrix[N_plots, P_cov] X;
  matrix[N_sp, P_trait > 0 ? P_trait : 1] W_sp;
}
parameters {
  vector[N_sp] alpha_0_raw;
  real<lower=0> sigma_alpha;
  vector[P_trait] gamma_int;
  vector[P_cov] b;
  matrix[P_trait, P_cov] gamma_cov;
}
transformed parameters {
  vector[N_sp] alpha_0;
  if (P_trait > 0) alpha_0 = W_sp * gamma_int + sigma_alpha * alpha_0_raw;
  else             alpha_0 = sigma_alpha * alpha_0_raw;
}
model {
  alpha_0_raw ~ std_normal();
  sigma_alpha ~ normal(0, 1);
  gamma_int ~ normal(0, 0.5);
  b ~ normal(0, 0.5);
  to_vector(gamma_cov) ~ normal(0, 0.3);
  {
    vector[N_plots] X_b = X * b;
    matrix[N_plots, N_sp] trait_modul = (X * gamma_cov') * W_sp';
    for (i in 1:N_plots) {
      row_vector[N_sp] la;
      for (s in 1:N_sp) la[s] = alpha_0[s] + X_b[i] + trait_modul[i, s];
      target += multinomial_logit_lpmf(y_count[i] | to_vector(la));
    }
  }
}
generated quantities {
  vector[N_plots] log_lik;
  {
    vector[N_plots] X_b = X * b;
    matrix[N_plots, N_sp] trait_modul = (X * gamma_cov') * W_sp';
    for (i in 1:N_plots) {
      row_vector[N_sp] la;
      for (s in 1:N_sp) la[s] = alpha_0[s] + X_b[i] + trait_modul[i, s];
      log_lik[i] = multinomial_logit_lpmf(y_count[i] | to_vector(la));
    }
  }
}
