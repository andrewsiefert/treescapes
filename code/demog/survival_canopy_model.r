library(tidyverse)
library(rstan)

options(mc.cores = parallel::detectCores())

source("code/transformer.r")
tf <- readRDS("data/transformers.rds")

# load and prepare data
data <- read_csv("data/survival_canopy_data_train.csv")

d <- data %>% 
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

d_list <- list(N = nrow(d), 
               Nsp = n_distinct(d$sp),
               Nsite = n_distinct(d$site),
               sp =  as.numeric(factor(d$sp)),
               site = as.numeric(factor(d$site)),
               X = X,
               y = d$surv, 
               yrs = d$remper)

# fit the model
model <- stan_model("code/demog/survival_canopy_model.stan")

fit <- sampling(model, data = d_list, chains = 4, cores = 4, 
                pars = c("tau_site", "tau_sp", "g_site", "sp_eff", "beta_pop", "theta"),
                save_warmup = F, init_r = 1)

# save the model
saveRDS(fit, "results/models/survival_canopy_model.rds")
