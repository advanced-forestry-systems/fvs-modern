// ingrowth_hurdle_v2_z.stan
// Standardized variant of ingrowth_hurdle_v1.stan. Identical hurdle NB model,
// but z-scores the continuous predictors (ln_ba, ln_bal, rd, ht40, ln_csi,
// clim_pca1) in transformed data to fix the ill-conditioned geometry that
// caused 100% max-treedepth saturation and 48h timeouts. log_years stays a
// raw offset. REs were already non-centered. Coeffs are on the standardized
// scale (back-transform for engine use).
data {
  int<lower=1> N_plots; int<lower=1> N_sp; int<lower=1> N_L1; int<lower=1> N_L2; int<lower=1> N_L3;
  int<lower=0> P_trait;
  array[N_plots] int<lower=0> n_recruits;
  vector[N_plots] log_years;
  vector[N_plots] ln_ba; vector[N_plots] ln_bal; vector[N_plots] rd;
  vector[N_plots] ht40; vector[N_plots] ln_csi; vector[N_plots] clim_pca1;
  array[N_plots] int<lower=1, upper=N_sp> dom_sp_idx;
  array[N_plots] int<lower=1, upper=N_L1> L1_idx;
  array[N_plots] int<lower=1, upper=N_L2> L2_idx;
  array[N_plots] int<lower=1, upper=N_L3> L3_idx;
  matrix[N_sp, P_trait > 0 ? P_trait : 1] W_dom;
}
transformed data {
  vector[N_plots] z_ba; vector[N_plots] z_bal; vector[N_plots] z_rd;
  vector[N_plots] z_ht40; vector[N_plots] z_csi; vector[N_plots] z_pca1;
  { real s;
    s=sd(ln_ba);     z_ba  =(ln_ba     -mean(ln_ba))    /(s>0?s:1);
    s=sd(ln_bal);    z_bal =(ln_bal    -mean(ln_bal))   /(s>0?s:1);
    s=sd(rd);        z_rd  =(rd         -mean(rd))       /(s>0?s:1);
    s=sd(ht40);      z_ht40=(ht40       -mean(ht40))     /(s>0?s:1);
    s=sd(ln_csi);    z_csi =(ln_csi     -mean(ln_csi))   /(s>0?s:1);
    s=sd(clim_pca1); z_pca1=(clim_pca1  -mean(clim_pca1))/(s>0?s:1);
  }
}
parameters {
  real a0_p; vector[P_trait] gamma_p;
  real b1_p; real b2_p; real b3_p; real b4_p; real b5_p; real b6_p;
  vector[N_L1] z_L1_p_raw; vector[N_L2] z_L2_p_raw; vector[N_L3] z_L3_p_raw;
  real<lower=0> sigma_L1_p; real<lower=0> sigma_L2_p; real<lower=0> sigma_L3_p;
  real a0_c; vector[P_trait] gamma_c;
  real b1_c; real b2_c; real b3_c; real b4_c; real b5_c; real b6_c;
  vector[N_L1] z_L1_c_raw; vector[N_L2] z_L2_c_raw; vector[N_L3] z_L3_c_raw;
  real<lower=0> sigma_L1_c; real<lower=0> sigma_L2_c; real<lower=0> sigma_L3_c;
  real<lower=0> phi;
}
transformed parameters {
  vector[N_sp] trait_p; vector[N_sp] trait_c;
  if (P_trait > 0) { trait_p = W_dom * gamma_p; trait_c = W_dom * gamma_c; }
  else { trait_p = rep_vector(0.0, N_sp); trait_c = rep_vector(0.0, N_sp); }
  vector[N_L1] z_L1_p = sigma_L1_p * z_L1_p_raw;
  vector[N_L2] z_L2_p = sigma_L2_p * z_L2_p_raw;
  vector[N_L3] z_L3_p = sigma_L3_p * z_L3_p_raw;
  vector[N_L1] z_L1_c = sigma_L1_c * z_L1_c_raw;
  vector[N_L2] z_L2_c = sigma_L2_c * z_L2_c_raw;
  vector[N_L3] z_L3_c = sigma_L3_c * z_L3_c_raw;
}
model {
  a0_p ~ normal(0.0,2.0); a0_c ~ normal(0.0,2.0);
  gamma_p ~ normal(0,0.5); gamma_c ~ normal(0,0.5);
  b1_p ~ normal(-0.2,0.5); b1_c ~ normal(-0.3,0.5);
  b2_p ~ normal(-0.2,0.5); b2_c ~ normal(-0.2,0.5);
  b3_p ~ normal(-0.3,1.0); b3_c ~ normal(-0.5,1.0);
  b4_p ~ normal( 0.0,0.5); b4_c ~ normal( 0.0,0.5);
  b5_p ~ normal( 0.3,0.5); b5_c ~ normal( 0.3,0.5);
  b6_p ~ normal( 0.0,0.5); b6_c ~ normal( 0.0,0.5);
  z_L1_p_raw ~ std_normal(); z_L2_p_raw ~ std_normal(); z_L3_p_raw ~ std_normal();
  z_L1_c_raw ~ std_normal(); z_L2_c_raw ~ std_normal(); z_L3_c_raw ~ std_normal();
  sigma_L1_p ~ normal(0,1.0); sigma_L1_c ~ normal(0,1.0);
  sigma_L2_p ~ normal(0,0.5); sigma_L2_c ~ normal(0,0.5);
  sigma_L3_p ~ normal(0,0.3); sigma_L3_c ~ normal(0,0.3);
  phi ~ gamma(2.0,0.5);
  vector[N_plots] eta_p = a0_p + trait_p[dom_sp_idx]
    + z_L1_p[L1_idx] + z_L2_p[L2_idx] + z_L3_p[L3_idx]
    + b1_p*z_ba + b2_p*z_bal + b3_p*z_rd + b4_p*z_ht40 + b5_p*z_csi + b6_p*z_pca1;
  vector[N_plots] eta_c = a0_c + trait_c[dom_sp_idx]
    + z_L1_c[L1_idx] + z_L2_c[L2_idx] + z_L3_c[L3_idx]
    + b1_c*z_ba + b2_c*z_bal + b3_c*z_rd + b4_c*z_ht40 + b5_c*z_csi + b6_c*z_pca1 + log_years;
  vector[N_plots] eta_c_safe;
  for (i in 1:N_plots) eta_c_safe[i] = fmin(fmax(eta_c[i], -20.0), 12.0);
  for (i in 1:N_plots) {
    if (n_recruits[i] == 0) target += bernoulli_logit_lpmf(0 | eta_p[i]);
    else target += bernoulli_logit_lpmf(1 | eta_p[i])
              + neg_binomial_2_log_lpmf(n_recruits[i] | eta_c_safe[i], phi)
              - log1m_exp(neg_binomial_2_log_lpmf(0 | eta_c_safe[i], phi));
  }
}
generated quantities {
  vector[N_plots] log_lik;
  vector[N_plots] eta_p_gq = a0_p + trait_p[dom_sp_idx]
    + z_L1_p[L1_idx] + z_L2_p[L2_idx] + z_L3_p[L3_idx]
    + b1_p*z_ba + b2_p*z_bal + b3_p*z_rd + b4_p*z_ht40 + b5_p*z_csi + b6_p*z_pca1;
  vector[N_plots] eta_c_gq = a0_c + trait_c[dom_sp_idx]
    + z_L1_c[L1_idx] + z_L2_c[L2_idx] + z_L3_c[L3_idx]
    + b1_c*z_ba + b2_c*z_bal + b3_c*z_rd + b4_c*z_ht40 + b5_c*z_csi + b6_c*z_pca1 + log_years;
  for (i in 1:N_plots) {
    real ec = fmin(fmax(eta_c_gq[i], -20.0), 12.0);
    if (n_recruits[i]==0) log_lik[i]=bernoulli_logit_lpmf(0 | eta_p_gq[i]);
    else log_lik[i]=bernoulli_logit_lpmf(1 | eta_p_gq[i])
       + neg_binomial_2_log_lpmf(n_recruits[i] | ec, phi)
       - log1m_exp(neg_binomial_2_log_lpmf(0 | ec, phi));
  }
}
