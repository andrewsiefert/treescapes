library(tidyverse)
library(rstan)

source("code/transformer.r")
tf <- readRDS("data/transformers.rds")


# Survival, canopy----------------------------------------------------------------

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

X <- model.matrix(f) %>% as_tibble() %>% janitor::clean_names()

fit <- readRDS("results/models/survival_canopy_model.rds")

samples <- rstan::extract(fit)$beta_pop %>% as_tibble()
names(samples) <- names(X)

saveRDS(samples, "results/demog_pars/surv_canopy_pars.rds")



# Survival, sapling----------------------------------------------

data <- read_csv("data/survival_sapling_data_train.csv")

d <- data %>% 
  mutate(wd = transform(wood_density, tf$trait$wood_density_log, log = T),
         sla = transform(sla, tf$trait$sla), 
         hmax = transform(hmax, tf$trait$hmax_log, log = T),
         mat = transform(mat, tf$env$mat),
         logsize = transform(prevdia, tf$sapling$log_prevdia, log = T),
         crowd = transform(canopy_nbr_ba, tf$tree$log_canopy_nbr_ba, log = T), 
         wd_sq = wd^2,
         sla_sq = sla^2,
         hmax_sq = hmax^2)

f <- lm(surv ~ 1 + mat + logsize + crowd + 
          wd + sla + hmax + 
          wd_sq + sla_sq + hmax_sq + 
          wd:sla + wd:hmax + sla:hmax + 
          wd:logsize + sla:logsize + hmax:logsize +
          mat:wd + mat:sla + mat:hmax + 
          mat:wd_sq + mat:sla_sq + mat:hmax_sq +
          mat:wd:sla + mat:wd:hmax + mat:sla:hmax,
        data = d)

X <- model.matrix(f) %>% as_tibble() %>% janitor::clean_names()

fit <- readRDS("results/models/survival_sapling_model.rds")

samples <- rstan::extract(fit)$beta_pop %>% as_tibble()
names(samples) <- names(X)

saveRDS(samples, "results/demog_pars/surv_sapling_pars.rds")


# Growth -------------------------------------------------------

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

X <- model.matrix(f) %>% as_tibble() %>% janitor::clean_names()

fit <- readRDS("results/models/growth_model.rds")

samples <- rstan::extract(fit, pars = "g_site", include = F)
samples$beta_pop <- as_tibble(samples$beta_pop)
names(samples$beta_pop) <- names(X)

saveRDS(samples, "results/demog_pars/growth_pars.rds")


# Recruitment ------------------------------------------------------------

data <- readRDS("data/recruitment_data_train.rds")

d <- data$df %>% 
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

X <- model.matrix(f)[,-1] %>% as_tibble() %>% janitor::clean_names()

fit <- readRDS("results/models/recruitment_model.rds")

samples <- rstan::extract(fit)[c('d0_pop', 'slope_pop', 'str_pop', 'alpha_pop',
                                 'beta_d0', 'beta_slope', 'beta_alpha', 'beta_pop')]
rm(fit)

samples$beta_pop <- as_tibble(samples$beta_pop)
names(samples$beta_pop) <- names(X)

saveRDS(samples, "results/demog_pars/recruitment_pars.rds")
