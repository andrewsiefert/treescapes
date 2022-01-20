library(tidyverse)
library(rstan)
library(doParallel)


source("code/transformer.r")
tf <- readRDS("data/transformers.rds")

# Load and prepare test data ----
d <- readRDS("data/recruitment_data_train.rds")$df

# Load model and get posterior samples ----
fit <- readRDS("results/models/recruitment_model.rds") 
samples <- rstan::extract(fit)
rm(fit)

recr_exp <- t(samples$R)

# Get posterior predictions ----
n <- nrow(d)
n_samples <- 1000
ss <- sample(1:length(samples$disp), n_samples)

recr_pred <- sapply(ss, function(s) {
  rnbinom(n, mu = recr_exp[,s], size = 1/sqrt(samples$disp[s]))
})


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

recr_pred_s <- recr_pred[,1:400]


# Export results
out <- list(df = df, r2 = r2_bayes, r = r_pred_obs, recr_pred = recr_pred_s)
saveRDS(out, "results/model_checks/recruitment_check_in_sample.rds")
