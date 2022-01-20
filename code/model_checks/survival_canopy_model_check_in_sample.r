library(tidyverse)
library(rstan)
library(doParallel)


# Load and prepare test data ----
d <- read_csv("data/survival_canopy_data_train.csv")

# Load model and get posterior samples ----
fit <- readRDS("results/models/survival_canopy_model.rds")

S <- 4000
n_samples <- 1000
ss <- sample(1:S, n_samples)

surv_exp <- rstan::extract(fit)$theta[ss,] %>% t()

rm(fit)


# Get posterior predictions ----
n <- nrow(d)

surv_pred <- apply(surv_exp, 2, function(i) rbinom(n, 1, i))


# Calculate Bayesian AUC
cl <- makeCluster(20)
registerDoParallel(cl)

auc_bayes <- foreach(i=1:ncol(surv_pred), .combine = "c") %dopar% {
  pROC::roc(d$surv, surv_exp[,i])$auc
}

stopCluster(cl)

hist(auc_bayes)

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
saveRDS(out, "results/model_checks/canopy_survival_check_in_sample.rds")
