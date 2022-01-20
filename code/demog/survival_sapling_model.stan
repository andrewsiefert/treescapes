data {
  int<lower=0> N;                // number of observations
  int<lower=0> Nsp;              // number of species
  int<lower=0> Nsite;            // number of sites
  int<lower=0> sp[N];            // species labels
  int<lower=0> site[N];          // site labels
  matrix[N,25] X;                // environmental values
  int<lower=0,upper=1> y[N];     // survival
  vector[N] yrs;                 // census interval (yrs/10)
}


parameters {
  vector[25] beta_pop;           // population-level parameters
  
  vector[Nsite] g_site_hat;      // standardized site effects
  real<lower=0> tau_site;        // scale of site random effects
  
  matrix[4,Nsp] sp_eff_hat;      // species random effects
  vector<lower=0>[4] tau_sp;     // scale of species random effects
  cholesky_factor_corr[4] L;     // correlation of species random effects
}


transformed parameters{
  
  vector[Nsite] g_site;              // site random effects
  matrix[Nsp,4] sp_eff;              // species random effects

  vector[N] mu_pop;
  vector[N] theta;
  
  g_site = g_site_hat*tau_site;
  sp_eff = (diag_pre_multiply(tau_sp, L) * sp_eff_hat)';
  
  mu_pop = X*beta_pop; 
  for(i in 1:N) {
    theta[i] = inv_logit(mu_pop[i] + g_site[site[i]] + X[i,1:4]*sp_eff[sp[i],]')^yrs[i];
  }
}


model {
  
  // priors
  beta_pop ~ std_normal();
  g_site_hat ~ std_normal();
  to_vector(sp_eff_hat) ~ std_normal();

  // hyperpriors
  tau_site ~ normal(0, 0.25);
  tau_sp ~ normal(0, 0.25);
  L ~ lkj_corr_cholesky(2);
  
  // likelihood
  y ~ bernoulli(theta);
}


generated quantities {
  vector[N] log_lik;
  
  for(i in 1:N) {
    log_lik[i] = bernoulli_lpmf(y[i] | theta[i]);
  }
}
