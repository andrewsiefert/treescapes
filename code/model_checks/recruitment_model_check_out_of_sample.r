library(tidyverse)
library(rstan)
library(doParallel)


source("code/transformer.r")
tf <- readRDS("data/transformers.rds")

# Load and prepare test data ----
test <- readRDS("data/recruitment_data_test.rds")

d <- test$df %>% 
  mutate(wd = transform(wood_density, tf$trait$wood_density_log, log = T),
         sla = transform(sla, tf$trait$sla), 
         hmax = transform(hmax, tf$trait$hmax_log, log = T),
         mat = transform(mat, tf$env$mat),
         crowd_c = transform(canopy_nbr_ba, tf$recr$canopy_nbr_ba_log, log = T), 
         crowd_s = transform(sapling_nbr_ba, tf$recr$sapling_nbr_ba),
         wd_sq = wd^2,
         sla_sq = sla^2,
         hmax_sq = hmax^2)

f <- lm(ingrowth ~ mat + crowd_c + crowd_s + 
          wd + sla + hmax + 
          wd_sq + sla_sq + hmax_sq + 
          wd:sla + wd:hmax + sla:hmax + 
          mat:wd + mat:sla + mat:hmax +
          mat:wd_sq + mat:sla_sq + mat:hmax_sq +
          mat:wd:sla + mat:wd:hmax + mat:sla:hmax,
        data = d)

X <- model.matrix(f)[,-1]




# Load model and get posterior samples ----
fit <- readRDS("results/models/recruitment_model.rds") 
samples <- rstan::extract(fit)
rm(fit)


# Get posterior expectations ----
N <- nrow(d)
S <- length(samples$str_pop)
sp <- d$sp

n_ss <- 1e3
ss <- sample(1:S, n_ss)

pop_effects <- X %*% t(samples$beta_pop)


cl <- makeCluster(19)
registerDoParallel(cl)

recr_out <- foreach(s = ss) %dopar% {
  
  d0 <- samples$d0[s,]
  slope <- samples$slope[s,]
  alpha <- samples$alpha[s,]
  str <- samples$str[s,]
  
  pop_eff <- pop_effects[,s]
  sp_eff <- samples$sp_eff[s,,]
  
  exp <- rep(0, N)

  for(i in 1:N) {
    
    dia <- d$canopy_dia[[i]]/30
    sap_dia <- d$sapling_dia[[i]]/30
    
    c_n <- length(dia)
    s_n <- length(sap_dia)
    
    f_c <- rep(0, c_n)
    f_s <- rep(0, s_n)
    
    # recruitment for individual canopy trees
    for(j in 1:c_n) {
      if (dia[j]>0) {
        f_c[j] = plogis(slope[sp[i]] * (log(dia[j]) - d0[sp[i]])) * dia[j]^alpha[sp[i]];
      } else {
        f_c[j] = 0
      }
    }
    
    # recruitment for individual saplings
    for(k in 1:s_n) {
      if (sap_dia[k]>0) {
        f_s[k] = plogis(slope[sp[i]] * (log(sap_dia[k]) - d0[sp[i]])) * sap_dia[k]^alpha[sp[i]];
      } else {
        f_s[k] = 0;
      }
    }
    
    recr <- str[sp[i]]*(sum(f_c)/12.5 + sum(f_s))         # potential annual recruitment
    exp[i] <- recr * exp(pop_eff[i] + X[i,1:3] %*% sp_eff[sp[i],5:7]) * test$yrs[i]
  }
  pred <- rnbinom(N, mu = exp, size = 1/sqrt(samples$disp[s]))
  return(list(exp = exp, pred = pred))
}

stopCluster(cl)


recr_exp <- sapply(recr_out, function(i) i$exp)
recr_pred <- sapply(recr_out, function(i) i$pred)


# Calculate Bayesian R2 ----
get_r2_bayes <- function(obs, pred) {
  var_pred <- mean((pred-mean(pred))^2)
  var_resid <- mean((obs-pred)^2)
  return(var_pred/(var_pred+var_resid))
}

r2_bayes <- apply(recr_exp, 2, function(i) get_r2_bayes(d$ingrowth, i))
r_pred_obs <- apply(recr_exp, 2, function(i) cor(i, d$ingrowth))

# Calculate prediction intervals and coverage ----
df <- tibble(obs = d$ingrowth, 
             exp_mean = rowMeans(recr_exp),
             exp_med = apply(recr_exp, 1, median),
             pred_mean = rowMeans(recr_pred),
             pred_med = apply(recr_pred, 1, median),
             q2.5 = apply(recr_pred, 1, quantile, 0.025),
             q5 = apply(recr_pred, 1, quantile, 0.05),
             q95 = apply(recr_pred, 1, quantile, 0.95),
             q97.5 = apply(recr_pred, 1, quantile, 0.975))
df$cover95 <- (df$obs >= df$q2.5) * (df$obs <= df$q97.5)
df$cover90 <- (df$obs >= df$q5) * (df$obs <= df$q95)

recr_pred_50 <- recr_pred[,1:50]


# Export results
out <- list(df = df, r2 = r2_bayes, r = r_pred_obs, recr_pred = recr_pred_50)
saveRDS(out, "results/model_checks/recruitment_check_out_of_sample.rds")
