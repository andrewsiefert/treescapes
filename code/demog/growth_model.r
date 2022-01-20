library(tidyverse)
library(rstan)

options(mc.cores = parallel::detectCores())

source("code/transformer.r")
tf <- readRDS("data/transformers.rds")

# load and prepare data
data <- read_csv("data/growth_data_train.csv")

d <- data %>% 
  filter(ann_dia_growth < 40) %>%
  mutate(t1 = transform(wood_density, tf$trait$wood_density_log, log = T),
         t2 = transform(sla, tf$trait$sla), 
         t3 = transform(hmax, tf$trait$hmax_log, log = T),
         mat = transform(mat, tf$env$mat),
         size = transform(prevdia, tf$tree$prevdia),
         log_size = transform(prevdia, tf$tree$log_prevdia, log = T),
         crowd = transform(canopy_nbr_ba, tf$tree$log_canopy_nbr_ba, log = T), 
         t1_sq = t1^2,
         t2_sq = t2^2,
         t3_sq = t3^2)

f <- lm(ann_dia_growth ~ 1 + mat + log_size + size + crowd + 
          t1 + t2 + t3 + 
          t1_sq + t2_sq + t3_sq + 
          t1:t2 + t1:t3 + t2:t3 + 
          t1:log_size + t2:log_size + t3:log_size +
          t1:size + t2:size + t3:size +
          mat:t1 + mat:t2 + mat:t3 + 
          mat:t1_sq + mat:t2_sq + mat:t3_sq +
          mat:t1:t2 + mat:t1:t3 + mat:t2:t3,
        data = d)

X <- model.matrix(f)

d_list <- list(N = nrow(d), 
               Nsp = n_distinct(d$sp),
               Nsite = n_distinct(d$site),
               sp =  as.numeric(factor(d$sp)),
               site = as.numeric(factor(d$site)),
               X = X,
               y = d$ann_dia_growth)

# fit the model
model <- stan_model("code/demog/growth_model.stan")

fit <- sampling(model, data = d_list, chains = 4, cores = 4, 
                pars = c("tau_site", "tau_sp", "g_site", "sp_eff", "beta_pop", "beta", "cutoff"),
                save_warmup = F, init_r = 1)

# save the model
saveRDS(fit, "results/models/growth_model.rds")
