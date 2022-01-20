library(tidyverse)
library(rstan)
library(doParallel)


source("code/transformer.r")
tf <- readRDS("data/transformers.rds")

# Load and prepare test data ----
test <- read_csv("data/survival_canopy_data_test.csv")

d <- test %>% 
  mutate(wd = transform(wood_density, tf$trait$wood_density_log, log = T),
         sla = transform(sla, tf$trait$sla), 
         hmax = transform(hmax, tf$trait$hmax_log, log = T),
         mat = transform(mat, tf$env$mat),
         size = transform(prevdia, tf$canopy$prevdia),
         log_size = transform(prevdia, tf$canopy$log_prevdia, log = T),
         crowd = transform(canopy_nbr_ba, tf$tree$log_canopy_nbr_ba, log = T), 
         wd_sq = wd^2,
         sla_sq = sla^2,
         hmax_sq = hmax^2)

f <- lm(surv ~ 1 + mat + size + log_size + crowd + 
          wd + sla + hmax + 
          wd_sq + sla_sq + hmax_sq + 
          wd:sla + wd:hmax + sla:hmax + 
          wd:size + sla:size + hmax:size +
          wd:log_size + sla:log_size + hmax:log_size +
          mat:wd + mat:sla + mat:hmax + 
          mat:wd_sq + mat:sla_sq + mat:hmax_sq +
          mat:wd:sla + mat:wd:hmax + mat:sla:hmax,
        data = d)

X <- model.matrix(f)


# Load model and get posterior samples ----
fit <- readRDS("results/models/survival_canopy_model.rds")
samples <- rstan::extract(fit)
rm(fit)

n <- nrow(d)
n_samples <- nrow(samples$beta_pop)


# Get posterior expectations ----
S <- 1000
ps <- sample(1:n_samples, S)

cl <- makeCluster(10)
registerDoParallel(cl)

out <- foreach(s=ps) %dopar% {
  mu_pop <- X %*% samples$beta_pop[s,]
  ranef <- rep(0, n)
  for(i in 1:n) {
    ranef[i] <- X[i,1:5] %*% samples$sp_eff[s,d$sp[i],]
  }
  exp <- plogis(mu_pop + ranef)^d$remper
  pred <- rbinom(n, 1, exp)
  list(exp = exp, pred = pred)
}

stopCluster(cl)

surv_exp <- sapply(out, function(i) i$exp)
surv_pred <- sapply(out, function(i) i$pred)

# Get posterior predictions
surv_pred <- apply(surv_exp, 2, function(i) rbinom(n, 1, i))

# Calculate Bayesian AUC
cl <- makeCluster(20)
registerDoParallel(cl)

auc_bayes <- foreach(i=1:ncol(surv_pred), .combine = "c") %dopar% {
  pROC::roc(d$surv, surv_exp[,i])$auc
}

stopCluster(cl)

# Calculate prediction intervals and coverage ----
df <- tibble(obs = d$surv, 
             theta_mean = rowMeans(surv_exp),
             theta_mid = apply(surv_exp, 1, median),
             q2.5 = apply(surv_pred, 1, quantile, 0.025),
             q5 = apply(surv_pred, 1, quantile, 0.05),
             q95 = apply(surv_pred, 1, quantile, 0.95),
             q97.5 = apply(surv_pred, 1, quantile, 0.975))
df$cover95 <- (df$obs >= df$q2.5) * (df$obs <= df$q97.5)
df$cover90 <- (df$obs >= df$q5) * (df$obs <= df$q95)

# Export results
out <- list(df = df, auc = auc_bayes, post_mean = colMeans(surv_pred))
saveRDS(out, "results/model_checks/canopy_survival_check_out_of_sample.rds")
