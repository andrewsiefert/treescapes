library(tidyverse)
library(plyr)


source("code/transformer.r")
tf <- readRDS("data/transformers.rds")


# Load demographic parameters ------------------------------------------------

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

recr_size_pars <- read_csv("results/demog_pars/recr_size_pars.csv")
recr_pars$hu <- mean(plogis(recr_size_pars$b_hu_Intercept))
recr_pars$shape <- mean(recr_size_pars$shape)
recr_pars$mu <- mean(exp(recr_size_pars$b_Intercept))
recr_pars$scale <- recr_pars$mu/recr_pars$shape


# Demographic functions ---------------------------------------------------

## Survival, canopy --------------------------------------------------------

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


## Survival, sapling -------------------------------------------------------

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


## Combined survival ------------------------------------------------------
s_z <- function(z, traits = rep(0, 3), env = 0) {
  out <- NA
  if (z < 12.7) {
    out <- s_z_sapling(z, traits = traits, env = env)
  } else {
    out <- s_z_canopy(z, traits = traits, env = env)
  }
  return(out)
}


## Growth ------------------------------------------------------------------

g_z1z <- function(z1, z, upper = 100, traits = rep(0, 3), env = 0) {
  
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
  
  linpred <- as.matrix(df) %*% beta
  mu <- exp(linpred)
  p_zero <- plogis(growth_pars$cutoff, linpred)
  
  g_inc <- (z1-z)*10
  
  g <- ifelse(g_inc == 0, p_zero, 
              ifelse(z1 == upper, (1-p_zero) * (1-pgamma(g_inc, mu*growth_pars$beta, growth_pars$beta)),
                     ifelse(z1 < upper, (1-p_zero) * dgamma(g_inc, mu*growth_pars$beta, growth_pars$beta), 
                            0)))
  
  return(g)
}



## Recruitment -------------------------------------------------------------

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


# IPMs -------------------------------------------------------------------------

# function to construct IPM and extract lambda for trait/environment combination

get_lambda <- function(x) {
  
  traits <- as.numeric(x[1:3])
  env <- as.numeric(x[4])
  
  # create the mesh points
  L <- 2.5
  U <- 100
  m <- 1000
  h <- (U-L)/m
  y <- L + (1:m)*h - h/2
  
  # growth kernel
  G <- outer(y, y, g_z1z, traits = traits, env = env, upper = max(y)) # growth kernel
  for(i in 1:(m-2)) {
    G[(i+1):(m-1),i] <- G[(i+1):(m-1),i] / sum(G[(i+1):(m-1),i]) * (1-G[i,i]-G[m,i])
  }
  G[m,m-1] <- 1-G[m-1,m-1]
  G[m,m] <- 1
  
  # survival kernel
  S <- sapply(y, function(i) s_z(i, traits = traits, env = env))
  
  # growth/survival kernel
  matP <- sapply(1:nrow(G), function(i) G[,i]*S[i])
  
  # recruitment kernel
  matF <- matrix(0, nrow = m, ncol = m)
  recr_n <- f_z(y, traits = traits, env = env)
  size_dist <- dgamma(y-2.54, shape = recr_pars$shape, scale = recr_pars$scale)
  size_dist[1] <- 0
  size_dist <- (1-recr_pars$hu) * size_dist/sum(size_dist) 
  size_dist[1] <- recr_pars$hu
  for(i in 1:m) {
    matF[,i] <- size_dist*recr_n[i]
  }
  
  # full kernel
  K <- matP + matF
  
  # extract lambda
  lambda <- Re(eigen(K)$values[1])
  return(lambda)
}


## Get lambdas for trait/environment combinations -------------------------------

### Create grid of trait/environment values ------------
d <- read_csv("data/growth_data_train.csv")

wd_r <- range(d$wood_density) %>% transform(tf$trait$wood_density_log, log = T)
sla_r <- range(d$sla) %>% transform(tf$trait$sla)
hmax_r <- range(d$hmax) %>% transform(tf$trait$hmax_log, log = T)

mat_levs <- c(5, 10, 15) %>% transform(tf$env$mat)
wd_seq <-seq(wd_r[1], wd_r[2], length.out = 40)
sla_seq <-seq(sla_r[1], sla_r[2], length.out = 40)
hmax_seq <-seq(hmax_r[1], hmax_r[2], length.out = 40)


# wood density & sla
wd_sla <- expand.grid(wd = wd_seq, 
                      sla = sla_seq,
                      hmax = 0, 
                      mat = mat_levs) %>%
  as_tibble()

# wood density and max height
wd_hmax <- expand.grid(wd = wd_seq, 
                       sla = 0,
                       hmax = hmax_seq, 
                       mat = mat_levs) %>%
  as_tibble()

# sla and max height
sla_hmax <- expand.grid(wd = 0, 
                        sla = sla_seq,
                        hmax = hmax_seq, 
                        mat = mat_levs) %>%
  as_tibble()

grid <- bind_rows(wd_sla, wd_hmax, sla_hmax)


### Calculate lambdas --------

library(doFuture)
registerDoFuture()
plan("multicore")

lambda <- alply(grid, 1, .fun = get_lambda, 
                .parallel = T, 
                .paropts = list(.export = c('grid', 'get_lambda', 'g_z1z', 's_z_sapling', 's_z_canopy', 's_z', 'f_z'),
                                .packages = 'dplyr'))

grid <- grid %>%
  mutate(lambda = unlist(lambda), 
         wd_rs = backtransform(wd, tf$trait$wood_density_log, log = T),
         sla_rs = backtransform(sla, tf$trait$sla),
         hmax_rs = backtransform(hmax, tf$trait$hmax_log, log = T),
         mat = backtransform(mat, tf$env$mat))


# save results ------------------------------------------------------------

write_csv(grid, "results/ipm/fitness_landscapes.csv")