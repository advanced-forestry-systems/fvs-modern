// crown_ratio_t2_annualized.stan
// Annualized, path-invariant crown-ratio transition (Garcia rate form):
//   logit(CR2) = E + (logit(CR1) - E) * exp(-k * T_years)
//   E = equilibrium logit(CR) = b0 + trait + REs + b1 dbh + b2 dbh^2 + b3 ba + b4 bal + b6 ln_csi
// At T->0, CR2->CR1 (no change over zero interval); k is the annual approach rate.
// Replaces the interval-free logit(CR2)=...+b_cr1 logit(CR1) form.
data {
  int<lower=1> N_obs; int<lower=1> N_sp; int<lower=1> N_L1; int<lower=1> N_L2;
  int<lower=1> N_L3; int<lower=1> N_FT; int<lower=0> P_trait;
  vector[N_obs] cr2_logit; vector[N_obs] cr1_logit;
  vector[N_obs] dbh; vector[N_obs] dbh_sq; vector[N_obs] ba_metric; vector[N_obs] bal_metric; vector[N_obs] ln_csi;
  vector<lower=0>[N_obs] T_years;                         // NEW: measurement interval
  array[N_obs] int<lower=1,upper=N_sp> sp_idx; array[N_obs] int<lower=1,upper=N_L1> L1_idx;
  array[N_obs] int<lower=1,upper=N_L2> L2_idx; array[N_obs] int<lower=1,upper=N_L3> L3_idx;
  array[N_obs] int<lower=1,upper=N_FT> FT_idx; matrix[N_sp, P_trait>0?P_trait:1] W;
}
parameters {
  real b0; real<lower=0> k;                               // k = annual approach rate (NEW, replaces b_cr1)
  real b1; real b2; real b3; real b4; real b6; vector[P_trait] gamma;
  vector[N_sp] z_sp_raw; vector[N_L1] z_L1_raw; vector[N_L2] z_L2_raw; vector[N_L3] z_L3_raw; vector[N_FT] z_FT_raw;
  real<lower=0> sigma_sp; real<lower=0> sigma_L1; real<lower=0> sigma_L2; real<lower=0> sigma_L3; real<lower=0> sigma_FT; real<lower=0> sigma;
}
transformed parameters {
  vector[N_sp] trait_effect; if (P_trait>0) trait_effect = W*gamma; else trait_effect = rep_vector(0.0,N_sp);
  vector[N_sp] z_sp = sigma_sp*z_sp_raw; vector[N_L1] z_L1 = sigma_L1*z_L1_raw; vector[N_L2] z_L2 = sigma_L2*z_L2_raw;
  vector[N_L3] z_L3 = sigma_L3*z_L3_raw; vector[N_FT] z_FT = sigma_FT*z_FT_raw;
}
model {
  b0 ~ normal(0,1.0); k ~ normal(0.05,0.05);
  b1 ~ normal(0,0.1); b2 ~ normal(0,0.01); b3 ~ normal(0,0.05); b4 ~ normal(0,0.05); b6 ~ normal(0,0.5); gamma ~ normal(0,0.5);
  z_sp_raw ~ std_normal(); z_L1_raw ~ std_normal(); z_L2_raw ~ std_normal(); z_L3_raw ~ std_normal(); z_FT_raw ~ std_normal();
  sigma_sp ~ normal(0,0.3); sigma_L1 ~ normal(0,0.3); sigma_L2 ~ normal(0,0.2); sigma_L3 ~ normal(0,0.2); sigma_FT ~ normal(0,0.2); sigma ~ normal(0,1.0);
  vector[N_obs] E = b0 + trait_effect[sp_idx] + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx] + z_FT[FT_idx]
                  + b1*dbh + b2*dbh_sq + b3*ba_metric + b4*bal_metric + b6*ln_csi;
  vector[N_obs] mu = E + (cr1_logit - E) .* exp(-k * T_years);
  cr2_logit ~ normal(mu, sigma);
}
