data {
  int<lower=0> N;                // number of observations
  int<lower=0> Nsp;              // number of species
  int<lower=0> Nsite;            // number of sites
  int<lower=0> sp[N];            // species labels
  int<lower=0> site[N];          // site labels
  matrix[N,29] X;                // environmental values
  vector[N] y;                   // annual growth
}


parameters {
  vector[29] beta_pop;           // population-level parameters
  
  vector[Nsite] g_site_hat;      // standardized site effects
  real<lower=0> tau_site;        // scale of site random effects
  
  matrix[5,Nsp] sp_eff_hat;      // species random effects
  vector<lower=0>[5] tau_sp;     // scale of species random effects
  cholesky_factor_corr[5] L;     // correlation of species random effects
  
  real<lower=0> beta;             // rate parameter of gamma distribution
  real cutoff;                    // cutoff for zero growth
}


transformed parameters{
  
  vector[Nsite] g_site;              // site random effects
  matrix[Nsp,5] sp_eff;              // species random effects

  vector[N] mu_pop;
  vector[N] growth;
  
  vector<lower=0,upper=1>[N] p_zero;

  
  g_site = g_site_hat*tau_site;
  sp_eff = (diag_pre_multiply(tau_sp, L) * sp_eff_hat)';
  
  mu_pop = X*beta_pop; 
  for(i in 1:N) {
    growth[i] = exp(mu_pop[i] + g_site[site[i]] + X[i,1:5]*sp_eff[sp[i],]');
    p_zero[i] = logistic_cdf(cutoff, log(growth[i]), 1);
  }
}


model {
  
  // priors
  beta_pop ~ std_normal();
  g_site_hat ~ std_normal();
  to_vector(sp_eff_hat) ~ std_normal();
  beta ~ lognormal(0.75, 0.75);
  cutoff ~ std_normal(); 

  // hyperpriors
  tau_site ~ normal(0, 0.25);
  tau_sp ~ normal(0, 0.25);
  L ~ lkj_corr_cholesky(2);
  
  // likelihood
  for(i in 1:N) {
    if(y[i] == 0) {
      target += bernoulli_lpmf(1 | p_zero[i]);
      } else {
      target += bernoulli_lpmf(0 | p_zero[i]) + gamma_lpdf(y[i] | growth[i]*beta, beta);
    }
  }
}
