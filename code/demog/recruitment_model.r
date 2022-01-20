library(tidyverse)
library(rstan)

options(mc.cores = parallel::detectCores())

source("code/transformer.r")
tf <- readRDS("data/transformers.rds")


# load and prepare FIA data -----------------------------------------------

# FIA data
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

X <- model.matrix(f)[,-1]

data$X <- X

data$trait <- d %>%   
  select(sp, wd, sla, hmax) %>% 
  distinct() %>%
  arrange(sp) %>%
  select(-sp) %>%
  as.matrix()

# MASTIF repro data
d2 <- read_csv("data/mastif_data.csv")

data$N_r <- nrow(d2)
data$Nsite_r <- n_distinct(d2$site)
data$Ntree_r <- n_distinct(d2$tree)
data$Nsp_r <- n_distinct(d2$AccSpeciesName)

data$site_r <- as.numeric(factor(d2$site)) 
data$tree_r <- as.numeric(factor(d2$tree)) 
data$sp_r <- as.numeric(factor(d2$AccSpeciesName)) 

data$dia_r <- d2$diam/30
data$repro <- d2$repro
data$trait_r <- d2 %>% 
  select(species, wd, sla, hmax) %>% 
  distinct() %>%
  arrange(species) %>%
  select(-species) %>%
  as.matrix()


# fit the model
model <- stan_model("code/demog/recruitment_model.stan")

fit <- sampling(model, data = data, chains = 4, cores = 4, 
                init_r = 1,
                pars = c("beta_pop", "d0_pop", "slope_pop", "alpha_pop", "str_pop", 
                         "beta_d0", "beta_slope", "beta_alpha",
                         "d0", "slope", "str", "alpha", "d0_r", "slope_r",
                         "tau_site", "tau_sp", "tau_site_r", "tau_tree_r", "tau_sp_r",
                         "R", "theta", "disp", "sp_eff", "sp_eff_r"), 
                control = list(adapt_delta = 0.9),
                save_warmup = F)

# save the model
saveRDS(fit, "results/models/recruitment_model.rds")
