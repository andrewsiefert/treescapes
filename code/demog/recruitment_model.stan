data {
  int<lower=0> N;
  int<lower=0> Nsp;
  int<lower=0> Nsite;
  int<lower=1,upper=Nsp> sp[N];
  int<lower=1, upper=Nsite> site[N];
  vector[N] yrs;
  matrix[N,21] X;                    // population-level predictors
  matrix[Nsp,3] trait;               // species trait values
  int y[N];

  // diameter vector lengths
  int c_obs; 
  int s_obs;

  // diameter vectors
  vector[c_obs] c_dia; 
  vector[s_obs] s_dia;

  // diameter segment positions
  int c_pos[N]; 
  int s_pos[N]; 

  // diameter segment lengths
  int c_n[N]; 
  int s_n[N]; 
  
  // repro data
  int<lower=0> N_r;                  // number of observations
  int<lower=0> Nsite_r;              // number of sites
  int<lower=0> Ntree_r;              // number of trees
  int<lower=0> Nsp_r;                // number of species

  int<lower=0> site_r[N_r];          // site index
  int<lower=0> tree_r[N_r];          // tree index
  int<lower=0> sp_r[N_r];            // species index

  vector[N_r] dia_r;                 // diameter
  matrix[Nsp_r,3] trait_r;           // species trait values

  int<lower=0, upper=1> repro[N_r];  // reproductive status
}


parameters {
  
  real d0_pop;                   // size threshold
  real slope_pop;
  real str_pop;                  // potential standardized number of recruits
  real alpha_pop;                // size effect
  
  vector[3] beta_d0;             // trait effects on size-repro threshold
  vector[3] beta_slope;          // trait effects on size-repro slope
  vector[3] beta_alpha;          // trait effects on size-fecundity effect
  
  vector[21] beta_pop;           // population-level parameters
  
  vector[Nsite] g_site_hat;      // standardized site effects
  real<lower=0> tau_site;        // scale of site random effects

  matrix[7,Nsp] sp_eff_hat;      // species random effects
  vector<lower=0>[7] tau_sp;     // scale of species random effects
  cholesky_factor_corr[7] L;     // correlation of species random effects
  
  real<lower=0> disp;
  
  vector[Nsite_r] g_site_hat_r;      // standardized site effects for repro model
  vector[Ntree_r] g_tree_hat_r;      // standardized site effects for repro model
  real<lower=0> tau_site_r;          // scale of site random effects for repro model
  real<lower=0> tau_tree_r;          // scale of site random effects for repro model
  
  matrix[2,Nsp_r] sp_eff_hat_r;        // species random effects
  vector<lower=0>[2] tau_sp_r;         // scale of species random effects
  cholesky_factor_corr[2] L_r;         // correlation of species random effects

}


transformed parameters{
  
  vector[Nsite] g_site;              // site random effects 
  matrix[Nsp,7] sp_eff;              // species random effects
  
  vector[Nsp] d0;
  vector[Nsp] slope;
  vector[Nsp] str;
  vector[Nsp] alpha;

  vector[N] S;                       // potential recruitment
  vector[N] R;                       // expected recruitment
  vector[N] pop_eff;                 // population-level climate, crowding, trait effects 
 
  vector[Nsp_r] d0_r;
  vector[Nsp_r] slope_r; 
  vector[Nsite_r] g_site_r;          // site effects for repro model
  vector[Ntree_r] g_tree_r;          // site effects for repro model
  matrix[Nsp_r,2] sp_eff_r;            // species random effects

  vector[N_r] theta;                 // probability of being reproductive  

  g_site = g_site_hat*tau_site;
  sp_eff = (diag_pre_multiply(tau_sp, L) * sp_eff_hat)';
  
  d0 = d0_pop + trait*beta_d0 + sp_eff[,1];
  slope = exp(slope_pop + trait*beta_slope + sp_eff[,2]);
  str = exp(str_pop + sp_eff[,3]);
  alpha = exp(alpha_pop + sp_eff[,4] + trait*beta_alpha);
  
  pop_eff = X*beta_pop;
  
 // get densities above threshold, effective densities, recruitment rate, expected recruits

  for(i in 1:N) {
    vector[c_n[i]] dia = segment(c_dia, c_pos[i], c_n[i]); 
    vector[s_n[i]] sap_dia = segment(s_dia, s_pos[i], s_n[i]); 
    vector[c_n[i]] f_c;
    vector[s_n[i]] f_s;

      // recruitment for individual canopy trees
      for(j in 1:c_n[i]) {
        if (dia[j]>0)
          f_c[j] = inv_logit(slope[sp[i]] * (log(dia[j]) - d0[sp[i]])) * dia[j]^alpha[sp[i]];
        else
          f_c[j] = 0;
      }
        
      // recruitment for individual saplings
      for(k in 1:s_n[i]) {
        if (sap_dia[k]>0)
          f_s[k] = inv_logit(slope[sp[i]] * (log(sap_dia[k]) - d0[sp[i]])) * sap_dia[k]^alpha[sp[i]];
        else
          f_s[k] = 0;
      }

    S[i] = str[sp[i]]*(sum(f_c)/12.5 + sum(f_s));         // potential annual recruitment
    R[i] = S[i] * exp(pop_eff[i] + g_site[site[i]] + X[i,1:3]*sp_eff[sp[i],5:7]') * yrs[i];
  }
  
  // mastif reproduction
  g_site_r = g_site_hat_r*tau_site_r; 
  g_tree_r = g_tree_hat_r*tau_tree_r; 
  sp_eff_r = (diag_pre_multiply(tau_sp_r, L_r) * sp_eff_hat_r)';
  
  d0_r = d0_pop + trait_r*beta_d0 + sp_eff_r[,1];
  slope_r = exp(slope_pop + trait_r*beta_slope + sp_eff_r[,2]);

  for(i in 1:N_r) {
    theta[i] = inv_logit(slope_r[sp_r[i]] * (log(dia_r[i]) - d0_r[sp_r[i]]) + g_site_r[site_r[i]] + g_tree_r[tree_r[i]]);
  }
}


model {
  
  // priors
  
  d0_pop ~ std_normal();
  slope_pop ~ std_normal();
  str_pop ~ normal(-1, 1);
  alpha_pop ~ normal(-0.5, 0.5);
  
  beta_d0 ~ std_normal();
  beta_slope ~ std_normal();
  beta_alpha ~ normal(0, 0.2);
  
  beta_pop ~ std_normal();
  
  g_site_hat ~ std_normal();
  tau_site ~ normal(0, 0.2);

  to_vector(sp_eff_hat) ~ std_normal();
  tau_sp ~ normal(0, 0.2);
  L ~ lkj_corr_cholesky(2);

  disp ~ std_normal();
  
  g_site_hat_r ~ std_normal();
  tau_site_r ~ normal(0, 0.2);
  
  g_tree_hat_r ~ std_normal();
  tau_tree_r ~ normal(0, 0.2);
  
  to_vector(sp_eff_hat_r) ~ std_normal();
  tau_sp_r ~ normal(0, 0.2);
  L_r ~ lkj_corr_cholesky(2);
  
  
  // likelihood
  repro ~ bernoulli(theta);
  y ~ neg_binomial_2(R, 1/sqrt(disp));
}
