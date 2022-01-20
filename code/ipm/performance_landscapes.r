library(tidyverse)
library(plyr)


source("code/transformer.r")
tf <- readRDS("data/transformers.rds")

# load models

surv_canopy_pars <- readRDS("results/demog_pars/surv_canopy_pars.rds") %>% apply(2, mean)
surv_sapling_pars <- readRDS("results/demog_pars/surv_sapling_pars.rds") %>% apply(2, mean)

growth_pars <- readRDS("results/demog_pars/growth_pars.rds")
growth_pars$beta_pop <- apply(growth_pars$beta_pop, 2, mean)
growth_pars$beta <- mean(growth_pars$beta)
growth_pars$cutoff <- mean(growth_pars$cutoff)

recr_pars <- readRDS("results/demog_pars/recruitment_pars.rds")
recr_pars$d0_pop <- mean(recr_pars$d0_pop)
recr_pars$slope_pop <- mean(recr_pars$slope_pop)
recr_pars$str_pop <- mean(recr_pars$str_pop)
recr_pars$alpha_pop <- mean(recr_pars$alpha_pop)
recr_pars$beta_d0 <- apply(recr_pars$beta_d0, 2, mean)
recr_pars$beta_slope <- apply(recr_pars$beta_slope, 2, mean)
recr_pars$beta_alpha <- apply(recr_pars$beta_alpha, 2, mean)
recr_pars$beta_pop <- apply(recr_pars$beta_pop, 2, mean)


# survival, canopy --------------------------------------------------------

s_z_canopy <- function(z, traits = rep(0, 3), env = 0) {
  
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
  
  surv <- plogis((as.matrix(df) %*% surv_canopy_pars))
  return(surv)
}


# survival, sapling -------------------------------------------------------

s_z_sapling <- function(z, traits = rep(0, 3), env = 0) {
  
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
  
  surv <- plogis((as.matrix(df) %*% surv_sapling_pars))
  return(surv)
}


# growth ------------------------------------------------------------------

g_z <- function(z, traits = rep(0, 3), env = 0) {
  
  beta <- growth_pars$beta_pop
  
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
  
  growth <- exp(as.matrix(df) %*% beta)

  return(growth)
}



# recruitment -------------------------------------------------------------

f_z <- function(z, traits = rep(0, 3), env = 0) {
  
  beta_pop <- recr_pars$beta_pop
  
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
  
  d0 <- recr_pars$d0_pop + traits %*% recr_pars$beta_d0
  slope <- exp(recr_pars$slope_pop + traits %*% recr_pars$beta_slope)
  str <- exp(recr_pars$str_pop)
  alpha <- exp(recr_pars$alpha_pop + traits %*% recr_pars$beta_alpha)
  
  pop_eff <- exp(as.matrix(df) %*% beta_pop)
  
  recr_n <- c(str)*plogis(c(slope)*(log(size) - c(d0)))*size^c(alpha) * c(pop_eff)
  
  return(recr_n)
  
}


# get demographic rates for trait/environment combinations -------------------------------

## prepare data to construct grids
d <- read_csv("data/growth_data_train.csv")

wd_r <- range(d$wood_density) %>% transform(tf$trait$wood_density_log, log = T)
sla_r <- range(d$sla) %>% transform(tf$trait$sla)
hmax_r <- range(d$hmax) %>% transform(tf$trait$hmax_log, log = T)

mat_levs <- c(5, 10, 15) %>% transform(tf$env$mat)
wd_seq <-seq(wd_r[1], wd_r[2], length.out = 40)
sla_seq <-seq(sla_r[1], sla_r[2], length.out = 40)
hmax_seq <-seq(hmax_r[1], hmax_r[2], length.out = 40)


# wood density & SLA
wd_sla <- expand.grid(size = c(5, 10, 20, 30, 40, 60), 
                      wd = wd_seq, 
                      sla = sla_seq,
                      hmax = 0, 
                      mat = mat_levs) %>%
  as_tibble()

# wood density and max height
wd_hmax <- expand.grid(size = c(5, 10, 20, 30, 40, 60), 
                       wd = wd_seq, 
                       sla = 0,
                       hmax = hmax_seq, 
                       mat = mat_levs) %>%
  as_tibble()

# SLA and max height
sla_hmax <- expand.grid(size = c(5, 10, 20, 30, 40, 60), 
                        wd = 0, 
                        sla = sla_seq,
                        hmax = hmax_seq, 
                        mat = mat_levs) %>%
  as_tibble()

grid <- bind_rows(wd_sla, wd_hmax, sla_hmax) 


## calculate demographic rates
out <- grid %>%
  mutate(growth = apply(grid, 1, function(i) g_z(i[1], traits = i[2:4], env = i[5])),
         surv_canopy = apply(grid, 1, function(i) s_z_canopy(i[1], traits = i[2:4], env = i[5])),
         surv_sapling = apply(grid, 1, function(i) s_z_sapling(i[1], traits = i[2:4], env = i[5])),
         recruitment = apply(grid, 1, function(i) f_z(i[1], traits = i[2:4], env = i[5])),
         wd_rs = backtransform(wd, tf$trait$wood_density_log, log = T),
         sla_rs = backtransform(sla, tf$trait$sla),
         hmax_rs = backtransform(hmax, tf$trait$hmax_log, log = T),
         mat = backtransform(mat, tf$env$mat))


# save results ------------------------------------------------------------

write_csv(out, "results/ipm/performance_landscapes.csv")
