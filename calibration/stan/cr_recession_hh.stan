// ============================================================================
// cr_recession_hh.stan
//
// FVS-CONUS crown recession (ΔHCB) model following Hann & Hanus (2004,
// Can. J. For. Res. 34: 1993-2003), adapted for CONUS-wide application
// with nested EPA L1/L2/L3 + species random intercepts and optional
// trait-informed species prior.
//
// This is the incremental method: ΔHCB is modeled directly rather than
// differencing two static HCB predictions.  The ORGANON-family dynamic
// crown model updates HCB each growth cycle as:
//     HCB2 = HCB1 + ΔHCB
//     CR2  = 1 - HCB2 / HT2
//
// Hann & Hanus best form (Eq. 9, with modifications from their Eqs. 9-10):
//
//   ΔHCB = (CL_S + ΔH) / (1 + exp(eta))
//
// where CL_S + ΔH = (HT1 - HCB1) + (HT2 - HT1) = HT2 - HCB1 is the
// maximum possible crown recession (the crown cannot recede beyond its
// full length plus height growth).
//
// We define the ratio:
//   r = ΔHCB / (CL_S + ΔH)  on (0, 1)
//
// so that  r = 1 / (1 + exp(eta)) = inv_logit(-eta).
// We use a Beta(mu*phi, (1-mu)*phi) likelihood with mu = inv_logit(-eta).
//
// The ratio r is dimensionless and approximately interval-independent
// because both numerator and denominator scale with the growth period
// length (through ΔH in the denominator).
//
// Linear predictor (adapted from Hann & Hanus 2004 Eqs. 9/10 for CONUS):
//   eta = r0 + z_sp + z_L1 + z_L2 + z_L3
//       + r1 * ln(CR_S)         crown ratio effect (nonlinear via log)
//       + r2 * CR_S             crown ratio linear term
//       + r3 * sqrt(BA)         stand density (BHA proxy)
//       + r4 * ln(BAL/BA + 1)   competition from above (CCF proxy)
//       + r5 * CR_S/RD          relative crown position (CR/CCF proxy)
//       + r6 * ln(CSPI_shift)   site productivity
//
// Expected behavior:
//   ΔHCB increases (higher r → lower eta) with:
//     - higher CCF/BAL (more side shading at crown base)
//     - intermediate CR (peak near 0.6, then decline; captured by
//       ln(CR) + CR interaction)
//     - lower stand density for same tree position
//   ΔHCB decreases with:
//     - very low CR (little crown to lose)
//     - very high CR / open stands (minimal side shading)
//     - higher site quality (Hann & Hanus report negative BHA effect)
//
// Reference:
//   Hann, D.W. and Hanus, M.L. 2004. Evaluation of nonspatial approaches
//   and equation forms used to predict tree crown recession.
//   Can. J. For. Res. 34: 1993-2003. doi:10.1139/X04-076
// ============================================================================
data {
  int<lower=1> N_obs;
  int<lower=1> N_sp;
  int<lower=1> N_L1;
  int<lower=1> N_L2;
  int<lower=1> N_L3;

  int<lower=0> P_trait;                   // species-level trait columns (0 = disable prior)

  vector<lower=0, upper=1>[N_obs] ratio;  // ΔHCB / (CL_S + ΔH), strictly (0,1)

  vector[N_obs] ln_cr;         // log(CR1)
  vector[N_obs] cr;            // CR1 (crown ratio at start)
  vector[N_obs] sqrt_ba;       // sqrt(BA1)
  vector[N_obs] ln_bal_ba;     // log(BAL1 / BA1 + 1)
  vector[N_obs] cr_over_rd;    // CR1 / rd_add (relative crown position)
  vector[N_obs] ln_cspi_shift; // log(CSPI + shift)

  array[N_obs] int<lower=1, upper=N_sp> sp_idx;
  array[N_obs] int<lower=1, upper=N_L1> L1_idx;
  array[N_obs] int<lower=1, upper=N_L2> L2_idx;
  array[N_obs] int<lower=1, upper=N_L3> L3_idx;

  matrix[N_sp, P_trait > 0 ? P_trait : 1] W;  // species x trait design (centered, scaled)
}
parameters {
  real r0;
  real r1;  // ln(CR)
  real r2;  // CR
  real r3;  // sqrt(BA)
  real r4;  // ln(BAL/BA + 1)
  real r5;  // CR / RD
  real r6;  // ln(CSPI_shift)

  vector[P_trait] gamma;                // trait loadings mapping W -> species RE mean

  vector[N_sp] z_sp_raw;
  vector[N_L1] z_L1_raw;
  vector[N_L2] z_L2_raw;
  vector[N_L3] z_L3_raw;
  real<lower=0> sigma_sp;
  real<lower=0> sigma_L1;
  real<lower=0> sigma_L2;
  real<lower=0> sigma_L3;

  real<lower=1> phi;     // Beta precision
}
transformed parameters {
  vector[N_sp] mu_sp;
  if (P_trait > 0) {
    mu_sp = W * gamma;
  } else {
    mu_sp = rep_vector(0.0, N_sp);
  }
  vector[N_sp] z_sp = mu_sp + sigma_sp * z_sp_raw;
  vector[N_L1] z_L1 = sigma_L1 * z_L1_raw;
  vector[N_L2] z_L2 = sigma_L2 * z_L2_raw;
  vector[N_L3] z_L3 = sigma_L3 * z_L3_raw;
}
model {
  // Priors guided by Hann & Hanus (2004) Table 6, Eq. 9
  // Their b5,0 ~ -5, b5,1 ~ -3.8, b5,2 ~ 5.7-6.4, b5,3 ~ -0.3, b5,4 ~ -0.27,
  // b5,5 ~ 121-128. Note: their parameterization is in the denominator
  // (1 + exp(eta)) so positive eta means LESS recession.
  r0 ~ normal(-5.0, 3.0);       // intercept
  r1 ~ normal(-4.0, 2.0);       // ln(CR): negative means higher CR → less recession (larger eta)
  r2 ~ normal( 6.0, 3.0);       // CR: positive offsets ln(CR), creates peaked response
  r3 ~ normal( 0.0, 0.5);       // sqrt(BA): stand density effect
  r4 ~ normal(-0.3, 0.5);       // ln(BAL/BA + 1): higher competition → more recession (lower eta)
  r5 ~ normal( 0.5, 1.0);       // CR/RD: higher relative crown → more recession suppressed
  r6 ~ normal(-0.2, 0.5);       // ln_cspi_shift: site productivity

  gamma ~ normal(0, 0.5);
  z_sp_raw ~ std_normal();
  z_L1_raw ~ std_normal();
  z_L2_raw ~ std_normal();
  z_L3_raw ~ std_normal();
  sigma_sp ~ normal(0, 0.5);
  sigma_L1 ~ normal(0, 0.5);
  sigma_L2 ~ normal(0, 0.3);
  sigma_L3 ~ normal(0, 0.3);

  phi ~ gamma(2, 0.05);   // mean 40, generous for ratio data

  vector[N_obs] eta =
      r0
    + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
    + r1 * ln_cr
    + r2 * cr
    + r3 * sqrt_ba
    + r4 * ln_bal_ba
    + r5 * cr_over_rd
    + r6 * ln_cspi_shift;

  // Hann & Hanus convention: r = 1 / (1 + exp(eta)) = inv_logit(-eta)
  // Clamp mu away from 0/1 to prevent beta shape parameters hitting zero,
  // which causes chains to get stuck with beta_lpdf exceptions at full data scale.
  vector[N_obs] mu_raw = inv_logit(-eta);
  vector[N_obs] mu = 0.001 + 0.998 * mu_raw;

  ratio ~ beta(mu * phi, (1 - mu) * phi);
}
generated quantities {
  vector[N_obs] log_lik;
  vector[N_obs] mu_pred;
  {
    vector[N_obs] eta =
        r0
      + z_sp[sp_idx] + z_L1[L1_idx] + z_L2[L2_idx] + z_L3[L3_idx]
      + r1 * ln_cr + r2 * cr + r3 * sqrt_ba
      + r4 * ln_bal_ba + r5 * cr_over_rd + r6 * ln_cspi_shift;
    mu_pred = 0.001 + 0.998 * inv_logit(-eta);
  }
  for (i in 1:N_obs) {
    log_lik[i] = beta_lpdf(ratio[i] | mu_pred[i] * phi, (1 - mu_pred[i]) * phi);
  }
}
