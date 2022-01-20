library(tidyverse)
library(plyr)


source("code/transformer.r")
tf <- readRDS("data/transformers.rds")

# Load demographic model parameters -------------------------------------------

surv_canopy_pars <- readRDS("results/demog_pars/surv_canopy_pars.rds") 
surv_sapling_pars <- readRDS("results/demog_pars/surv_sapling_pars.rds") 
growth_pars <- readRDS("results/demog_pars/growth_pars.rds")
recr_pars <- readRDS("results/demog_pars/recruitment_pars.rds")

# Demography functions ----

## Survival, canopy --------------------------------------------------------

s_z_canopy <- function(z, traits = rep(0, 3), env = 0, pars) {
  
  df <- tibble(intercept = 1,
               mat = env,
               size = transform(z, tf$canopy$prevdia),
               log_size = transform(z, tf$canopy$log_prevdia, log = T),
               crowd = 0,
               wd = traits[1],
               sla = traits[2], 
               hmax = traits[3],
               wd_sq = wd^2,
               sla_sq = sla^2,
               hmax_sq = hmax^2, 
               wd_sla = wd*sla,
               wd_hmax = wd*hmax,
               sla_hmax = sla*hmax,
               size_wd = wd*size,
               size_sla = sla*size,
               size_hmax = size*hmax,
               log_size_wd = wd*log_size,
               log_size_sla = sla*log_size,
               log_size_hmax = hmax*log_size,
               mat_wd = mat*wd,
               mat_sla = mat*sla,
               mat_hmax = mat*hmax,
               mat_wd_sq = mat*wd_sq,
               mat_sla_sq = mat*sla_sq,
               mat_hmax_sq = mat*hmax_sq,
               mat_wd_sla = mat*wd*sla,
               mat_wd_hmax = mat*wd*hmax,
               mat_sla_hmax = mat*sla*hmax)
  df <- df[,names(surv_canopy_pars)]
  
  surv <- plogis((as.matrix(df) %*% t(pars)))
  return(surv)
}

s_z_canopy(z = 20, pars = surv_canopy_pars) %>% hist()

## Survival, sapling -------------------------------------------------------

s_z_sapling <- function(z, traits = rep(0, 3), env = 0, pars) {
  
  df <- tibble(intercept = 1,
               mat = env,
               logsize = transform(z, tf$sapling$log_prevdia, log = T),
               crowd = 0,
               wd = traits[1],
               sla = traits[2], 
               hmax = traits[3],
               wd_sq = wd^2,
               sla_sq = sla^2,
               hmax_sq = hmax^2, 
               wd_sla = wd*sla,
               wd_hmax = wd*hmax,
               sla_hmax = sla*hmax,
               logsize_wd = wd*logsize,
               logsize_sla = sla*logsize,
               logsize_hmax = hmax*logsize,
               mat_wd = mat*wd,
               mat_sla = mat*sla,
               mat_hmax = mat*hmax,
               mat_wd_sq = mat*wd_sq,
               mat_sla_sq = mat*sla_sq,
               mat_hmax_sq = mat*hmax_sq,
               mat_wd_sla = mat*wd*sla,
               mat_wd_hmax = mat*wd*hmax,
               mat_sla_hmax = mat*sla*hmax)
  df <- df[,names(surv_sapling_pars)]
  
  surv <- plogis((as.matrix(df) %*% t(pars)))
  return(surv)
}

s_z_sapling(z = 5, pars = surv_sapling_pars) %>% hist()


## Growth ------------------------------------------------------------------

g_z <- function(z, traits = rep(0, 3), env = 0, pars) {
  
  beta <- pars$beta_pop
  
  df <- tibble(intercept = 1,
               mat = env,
               log_size = transform(z, tf$tree$log_prevdia, log = T),
               size = transform(z, tf$tree$prevdia),
               crowd = 0,
               t1 = traits[1],
               t2 = traits[2], 
               t3 = traits[3],
               t1_sq = t1^2,
               t2_sq = t2^2,
               t3_sq = t3^2, 
               t1_t2 = t1*t2,
               t1_t3 = t1*t3,
               t2_t3 = t2*t3,
               log_size_t1 = t1*log_size,
               log_size_t2 = t2*log_size,
               log_size_t3 = t3*log_size,
               size_t1 = t1*size,
               size_t2 = t2*size,
               size_t3 = size*t3,
               mat_t1 = mat*t1,
               mat_t2 = mat*t2,
               mat_t3 = mat*t3,
               mat_t1_sq = mat*t1_sq,
               mat_t2_sq = mat*t2_sq,
               mat_t3_sq = mat*t3_sq,
               mat_t1_t2 = mat*t1*t2,
               mat_t1_t3 = mat*t1*t3,
               mat_t2_t3 = mat*t2*t3)
  
  df <- df[,names(beta)]
  
  growth <- exp(as.matrix(df) %*% t(beta))

  return(growth)
}

g_z(z = 20, pars = growth_pars) %>% hist()


## Recruitment -------------------------------------------------------------

f_z <- function(z, traits = rep(0, 3), env = 0, pars) {
  
  beta_pop <- pars$beta_pop
  
  df <- tibble(mat = env,
               crowd_c = 0, 
               crowd_s = 0, 
               wd = traits[1],
               sla = traits[2], 
               hmax = traits[3],
               wd_sq = wd^2,
               sla_sq = sla^2,
               hmax_sq = hmax^2, 
               wd_sla = wd*sla,
               wd_hmax = wd*hmax,
               sla_hmax = sla*hmax,
               mat_wd = mat*wd,
               mat_sla = mat*sla,
               mat_hmax = mat*hmax,
               mat_wd_sq = mat*wd_sq,
               mat_sla_sq = mat*sla_sq,
               mat_hmax_sq = mat*hmax_sq,
               mat_wd_sla = mat*wd*sla,
               mat_wd_hmax = mat*wd*hmax,
               mat_sla_hmax = mat*sla*hmax)
  
  df <- df[,names(beta_pop)]
  
  size <- z/30
  
  traits <- matrix(traits, nrow = 1)
  
  d0 <- pars$d0_pop + as.vector(traits %*% t(pars$beta_d0))
  slope <- exp(pars$slope_pop + as.vector(traits %*% t(pars$beta_slope)))
  str <- exp(pars$str_pop)
  alpha <- exp(pars$alpha_pop + as.vector(traits %*% t(pars$beta_alpha)))
  
  pop_eff <- exp(as.matrix(df) %*% t(beta_pop))
  
  recr_n <- c(str)*plogis(c(slope)*(log(size) - c(d0)))*size^c(alpha) * c(pop_eff)
  
  return(recr_n)
  
}

hist(f_z(z = 20, pars = recr_pars))


# Demographic tradeoffs ----

## prepare data to construct grids
d <- read_csv("data/growth_data_train.csv")

seq_trim <- function(x, n = 100, q = 0.99) {seq(quantile(x, (1-q)/2), quantile(x, 1-(1-q)/2), length.out = n)}

wd_seq <- seq_trim(d$wood_density %>% transform(tf$trait$wood_density_log, log = T), 7)
sla_seq <- seq_trim(d$sla %>% transform(tf$trait$sla), 7)
hmax_seq <- seq_trim(d$hmax %>% transform(tf$trait$hmax_log, log = T), 7)

mat_levs <- c(5, 10, 15) %>% transform(tf$env$mat)

## Growth-survival tradeoff ----

sm <- 5

### Wood density ----
gs_wd <- expand.grid(wd = wd_seq,
                     sla = 0, 
                     hmax = 0, 
                     mat = mat_levs) %>%
  as_tibble()

gs_wd_s <- apply(gs_wd, 1, function(i) s_z_sapling(z = sm, traits = i[1:3], env = i[4], pars = surv_sapling_pars))
gs_wd_g <- apply(gs_wd, 1, function(i) g_z(z = sm, traits = i[1:3], env = i[4], pars = growth_pars))

gs_wd2 <- gs_wd %>%
  mutate(surv_est = colMeans(gs_wd_s),
         surv_lo = apply(gs_wd_s, 2, quantile, 0.05),
         surv_hi = apply(gs_wd_s, 2, quantile, 0.95),
         growth_est = colMeans(gs_wd_g),
         growth_lo = apply(gs_wd_g, 2, quantile, 0.05),
         growth_hi = apply(gs_wd_g, 2, quantile, 0.95), 
         mat = backtransform(mat, tf$env$mat) %>% paste0("°C") %>% fct_reorder(mat), 
         wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         tradeoff = "gs", 
         trait = "wd",  
         select(tradeoff, trait, everything()))

### SLA ----
gs_sla <- expand.grid(wd = 0,
                     sla = sla_seq, 
                     hmax = 0, 
                     mat = mat_levs) %>%
  as_tibble()

gs_sla_s <- apply(gs_sla, 1, function(i) s_z_sapling(z = sm, traits = i[1:3], env = i[4], pars = surv_sapling_pars))
gs_sla_g <- apply(gs_sla, 1, function(i) g_z(z = sm, traits = i[1:3], env = i[4], pars = growth_pars))

gs_sla2 <- gs_sla %>%
  mutate(surv_est = colMeans(gs_sla_s),
         surv_lo = apply(gs_sla_s, 2, quantile, 0.05),
         surv_hi = apply(gs_sla_s, 2, quantile, 0.95),
         growth_est = colMeans(gs_sla_g),
         growth_lo = apply(gs_sla_g, 2, quantile, 0.05),
         growth_hi = apply(gs_sla_g, 2, quantile, 0.95), 
         mat = backtransform(mat, tf$env$mat) %>% paste0("°C") %>% fct_reorder(mat), 
         sla = backtransform(sla, tf$trait$sla),
         tradeoff = "gs", 
         trait = "sla",  
         select(tradeoff, trait, everything()))

### Maximum height ----
gs_hmax <- expand.grid(wd = 0,
                     sla = 0, 
                     hmax = hmax_seq, 
                     mat = mat_levs) %>%
  as_tibble()

gs_hmax_s <- apply(gs_hmax, 1, function(i) s_z_sapling(z = sm, traits = i[1:3], env = i[4], pars = surv_sapling_pars))
gs_hmax_g <- apply(gs_hmax, 1, function(i) g_z(z = sm, traits = i[1:3], env = i[4], pars = growth_pars))

gs_hmax2 <- gs_hmax %>%
  mutate(surv_est = colMeans(gs_hmax_s),
         surv_lo = apply(gs_hmax_s, 2, quantile, 0.05),
         surv_hi = apply(gs_hmax_s, 2, quantile, 0.95),
         growth_est = colMeans(gs_hmax_g),
         growth_lo = apply(gs_hmax_g, 2, quantile, 0.05),
         growth_hi = apply(gs_hmax_g, 2, quantile, 0.95), 
         mat = backtransform(mat, tf$env$mat) %>% paste0("°C") %>% fct_reorder(mat), 
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T),
         tradeoff = "gs", 
         trait = "hmax",  
         select(tradeoff, trait, everything()))

## Stature-recruitment tradeoff ----

sm <- 8
lg <- 60

### Wood density ----
sr_wd <- expand.grid(wd = wd_seq,
                     sla = 0, 
                     hmax = 0, 
                     mat = mat_levs) %>%
  as_tibble()

sr_wd_s <- apply(sr_wd, 1, function(i) s_z_canopy(z = lg, traits = i[1:3], env = i[4], pars = surv_canopy_pars))
sr_wd_g <- apply(sr_wd, 1, function(i) g_z(z = lg, traits = i[1:3], env = i[4], pars = growth_pars))
sr_wd_r <- apply(sr_wd, 1, function(i) f_z(z = sm, traits = i[1:3], env = i[4], pars = recr_pars))

sr_wd2 <- sr_wd %>%
  mutate(surv_est = colMeans(sr_wd_s),
         surv_lo = apply(sr_wd_s, 2, quantile, 0.05),
         surv_hi = apply(sr_wd_s, 2, quantile, 0.95),
         growth_est = colMeans(sr_wd_g),
         growth_lo = apply(sr_wd_g, 2, quantile, 0.05),
         growth_hi = apply(sr_wd_g, 2, quantile, 0.95), 
         recr_est = colMeans(sr_wd_r),
         recr_lo = apply(sr_wd_r, 2, quantile, 0.05),
         recr_hi = apply(sr_wd_r, 2, quantile, 0.95), 
         mat = backtransform(mat, tf$env$mat) %>% paste0("°C") %>% fct_reorder(mat), 
         wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         tradeoff = "sr", 
         trait = "wd",  
         select(tradeoff, trait, everything()))

### SLA ----
sr_sla <- expand.grid(wd = 0,
                     sla = sla_seq, 
                     hmax = 0, 
                     mat = mat_levs) %>%
  as_tibble()

sr_sla_s <- apply(sr_sla, 1, function(i) s_z_canopy(z = lg, traits = i[1:3], env = i[4], pars = surv_canopy_pars))
sr_sla_g <- apply(sr_sla, 1, function(i) g_z(z = lg, traits = i[1:3], env = i[4], pars = growth_pars))
sr_sla_r <- apply(sr_sla, 1, function(i) f_z(z = sm, traits = i[1:3], env = i[4], pars = recr_pars))

sr_sla2 <- sr_sla %>%
  mutate(surv_est = colMeans(sr_sla_s),
         surv_lo = apply(sr_sla_s, 2, quantile, 0.05),
         surv_hi = apply(sr_sla_s, 2, quantile, 0.95),
         growth_est = colMeans(sr_sla_g),
         growth_lo = apply(sr_sla_g, 2, quantile, 0.05),
         growth_hi = apply(sr_sla_g, 2, quantile, 0.95), 
         recr_est = colMeans(sr_sla_r),
         recr_lo = apply(sr_sla_r, 2, quantile, 0.05),
         recr_hi = apply(sr_sla_r, 2, quantile, 0.95), 
         mat = backtransform(mat, tf$env$mat) %>% paste0("°C") %>% fct_reorder(mat), 
         sla = backtransform(sla, tf$trait$sla),
         tradeoff = "sr", 
         trait = "sla",  
         select(tradeoff, trait, everything()))

### Maximum height ----
sr_hmax <- expand.grid(wd = 0,
                     sla = 0, 
                     hmax = hmax_seq, 
                     mat = mat_levs) %>%
  as_tibble()

sr_hmax_s <- apply(sr_hmax, 1, function(i) s_z_canopy(z = lg, traits = i[1:3], env = i[4], pars = surv_canopy_pars))
sr_hmax_g <- apply(sr_hmax, 1, function(i) g_z(z = lg, traits = i[1:3], env = i[4], pars = growth_pars))
sr_hmax_r <- apply(sr_hmax, 1, function(i) f_z(z = sm, traits = i[1:3], env = i[4], pars = recr_pars))

sr_hmax2 <- sr_hmax %>%
  mutate(surv_est = colMeans(sr_hmax_s),
         surv_lo = apply(sr_hmax_s, 2, quantile, 0.05),
         surv_hi = apply(sr_hmax_s, 2, quantile, 0.95),
         growth_est = colMeans(sr_hmax_g),
         growth_lo = apply(sr_hmax_g, 2, quantile, 0.05),
         growth_hi = apply(sr_hmax_g, 2, quantile, 0.95), 
         recr_est = colMeans(sr_hmax_r),
         recr_lo = apply(sr_hmax_r, 2, quantile, 0.05),
         recr_hi = apply(sr_hmax_r, 2, quantile, 0.95), 
         mat = backtransform(mat, tf$env$mat) %>% paste0("°C") %>% fct_reorder(mat), 
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T),
         tradeoff = "sr", 
         trait = "hmax",  
         select(tradeoff, trait, everything()))

#### Plot ----
ggplot(sr_hmax2, aes(x = growth_est, xmin = growth_lo, xmax = growth_hi, 
                   y = recr_est, ymin = recr_lo, ymax = recr_hi, 
                   color = hmax)) + 
  geom_pointrange() +
  geom_errorbarh() +
  geom_point(size = 1.8) +
  facet_wrap(~mat)  + 
  theme(panel.grid.major = element_blank()) +
  scale_color_viridis_c(option = "magma") +
  scale_y_log10() +
  labs(x = "Growth, 60 cm dbh",
       y = "Recruitment, 5 dbh", 
       color = "Maximum height")

# Combine data

df <- bind_rows(gs_wd2, gs_sla2, gs_hmax2, sr_wd2, sr_sla2, sr_hmax2) %>%
  select(tradeoff, trait, everything())

write_csv(df, "results/demographic_tradeoffs/demographic_tradeoffs.csv")
