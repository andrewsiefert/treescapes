library(tidyverse)
library(rstan)
library(doParallel)


source("code/transformer.r")
tf <- readRDS("data/transformers.rds")

# Load data
data <- read_csv("data/growth_data_train.csv")

d <- data %>% 
  filter(ann_dia_growth < 40) %>%
  mutate(wd = transform(wood_density, tf$trait$wood_density_log, log = T),
         sla = transform(sla, tf$trait$sla), 
         hmax = transform(hmax, tf$trait$hmax_log, log = T),
         mat = transform(mat, tf$env$mat),
         size = transform(prevdia, tf$tree$prevdia),
         log_size = transform(prevdia, tf$tree$log_prevdia, log = T),
         crowd = transform(canopy_nbr_ba, tf$tree$log_canopy_nbr_ba, log = T), 
         wd_sq = wd^2,
         sla_sq = sla^2,
         hmax_sq = hmax^2,
         sp =  as.numeric(factor(sp)),
         site = as.numeric(factor(site)))

f <- lm(ann_dia_growth ~ 1 + mat + log_size + size + crowd + 
          wd + sla + hmax + 
          wd_sq + sla_sq + hmax_sq + 
          wd:sla + wd:hmax + sla:hmax + 
          wd:log_size + sla:log_size + hmax:log_size +
          wd:size + sla:size + hmax:size +
          mat:wd + mat:sla + mat:hmax + 
          mat:wd_sq + mat:sla_sq + mat:hmax_sq +
          mat:wd:sla + mat:wd:hmax + mat:sla:hmax,
        data = d)

X <- model.matrix(f)

d_list <- list(N = nrow(d), 
               Nsp = n_distinct(d$sp),
               Nsite = n_distinct(d$site),
               sp =  as.numeric(factor(d$sp)),
               site = as.numeric(factor(d$site)),
               X = X,
               y = d$ann_dia_growth)

# Load model ----
fit <- readRDS("results/models/growth_model.rds")
samples <- rstan::extract(fit)
rm(fit)

# Get posterior predictions ----
n <- nrow(d)
n_samples <- length(samples$cutoff)
S <- 1000
ps <- sample(1:n_samples, S)


cl <- makeCluster(19)
registerDoParallel(cl)

out <- foreach(s=ps) %dopar% {
  mu_pop <- X %*% samples$beta_pop[s,]
  site_eff <- samples$g_site[s,d$site]
  sp_eff <- rep(0, n)
  for(i in 1:n) {
    sp_eff[i] <- X[i,1:5] %*% samples$sp_eff[s,d$sp[i],]
  }
  y_exp <- exp(mu_pop + site_eff + sp_eff)
  p_zero <-  plogis(samples$cutoff[s], log(y_exp), 1)
  var_pred <- mean((y_exp - mean(y_exp))^2)
  var_resid <- mean((d$ann_dia_growth - y_exp)^2)
  y_pred <- (1 - rbinom(n, 1, p_zero)) * rgamma(n, y_exp*samples$beta[s], samples$beta[s])
  list(exp = y_exp, pred = y_pred, r2 = var_pred/(var_pred + var_resid), r = cor(y_exp, d$ann_dia_growth))
}

stopCluster(cl)

growth_exp <- do.call(cbind, lapply(out, function(i) i$exp))
growth_pred <- do.call(cbind, lapply(out, function(i) i$pred))
r2_bayes <- do.call(c, lapply(out, function(i) i$r2))
r_pred_obs <- do.call(c, lapply(out, function(i) i$r))

hist(r2_bayes)


# Calculate prediction intervals and coverage ----
df <- tibble(obs = d$ann_dia_growth, 
             exp_mean = rowMeans(growth_exp),
             exp_med = apply(growth_exp, 1, median),
             pred_mean = rowMeans(growth_pred),
             pred_med = apply(growth_pred, 1, median),
             q2.5 = apply(growth_pred, 1, quantile, 0.025),
             q5 = apply(growth_pred, 1, quantile, 0.05),
             q95 = apply(growth_pred, 1, quantile, 0.95),
             q97.5 = apply(growth_pred, 1, quantile, 0.975))
df$cover95 <- (df$obs >= df$q2.5) * (df$obs <= df$q97.5)
df$cover90 <- (df$obs >= df$q5) * (df$obs <= df$q95)

# take subset of growth predictions
growth_pred_s <- growth_pred %>% as_tibble() %>% sample_n(2000) 
growth_pred_s <- growth_pred_s[,1:100]

# Export results
out <- list(df = df, r2 = r2_bayes, r = r_pred_obs, growth_pred = growth_pred_s)
saveRDS(out, "results/model_checks/growth_check_in_sample.rds")
