library(plyr)
library(tidyverse)

source("code/transformer.r")
tf <- readRDS("data/transformers.rds")


# Load demographic parameters ------------------------------------------------

surv_canopy_pars <- readRDS("results/demog_pars/surv_canopy_pars.rds")
surv_sapling_pars <- readRDS("results/demog_pars/surv_sapling_pars.rds")

growth_pars <- readRDS("results/demog_pars/growth_pars.rds")

recr_pars <- readRDS("results/demog_pars/recruitment_pars.rds")

recr_size_pars <- read_csv("results/demog_pars/recr_size_pars.csv")
recr_pars$hu <- plogis(recr_size_pars$b_hu_Intercept)
recr_pars$shape <- recr_size_pars$shape
recr_pars$mu <- exp(recr_size_pars$b_Intercept)
recr_pars$scale <- recr_pars$mu/recr_pars$shape


# Demographic functions ---------------------------------------------------

## Survival, canopy --------------------------------------------------------

s_z_canopy <- function(z, traits = rep(0, 3), env = 0, pars, s) {
  
  beta_pop <- pars[s,] 
  
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
  
  df <- df[,names(beta_pop)]
  
  surv <- plogis((as.matrix(df) %*% t(beta_pop)))
  return(surv)
}


## Survival, sapling -------------------------------------------------------

s_z_sapling <- function(z, traits = rep(0, 3), env = 0, pars, s) {
  
  beta_pop <- pars[s,] 
  
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
  
  df <- df[,names(beta_pop)]
  
  surv <- plogis((as.matrix(df) %*% t(beta_pop)))
  return(surv)
}


## Combined survival ---------------
s_z <- function(z, traits = rep(0, 3), env = 0, sapling_pars, canopy_pars, s) {
  out <- NA
  if (z < 12.7) {
    out <- s_z_sapling(z, traits = traits, env = env, sapling_pars, s)
  } else {
    out <- s_z_canopy(z, traits = traits, env = env, canopy_pars, s)
  }
  return(out)
}


## Growth ------------------------------------------------------------------

g_z1z <- function(z1, z, upper = 100, traits = rep(0, 3), env = 0, pars, s) {
  
  beta_pop <- pars$beta_pop[s,]
  beta <- pars$beta[s]
  cutoff <- pars$cutoff[s]
  
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
  
  df <- df[,names(beta_pop)]
  
  linpred <- as.matrix(df) %*% t(beta_pop)
  mu <- exp(linpred)
  p_zero <- plogis(cutoff, linpred)
  
  g_inc <- (z1-z)*10
  
  g <- ifelse(g_inc == 0, p_zero, 
              ifelse(z1 == upper, (1-p_zero) * (1-pgamma(g_inc, mu*beta, beta)),
                     ifelse(z1 < upper, (1-p_zero) * dgamma(g_inc, mu*beta, beta), 
                            0)))
  
  return(g)
}



# Recruitment -------------------------------------------------------------

f_z <- function(z, traits = rep(0, 3), env = 0, pars, s) {
  
  d0_pop <- pars$d0_pop[s]
  slope_pop<- pars$slope_pop[s]
  str_pop <- pars$str_pop[s]
  alpha_pop <- pars$alpha_pop[s]
  beta_d0 <- pars$beta_d0[s,]
  beta_slope <- pars$beta_slope[s,]
  beta_alpha <- pars$beta_alpha[s,]
  beta_pop <- pars$beta_pop[s,]
  
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
  
  d0 <- d0_pop + traits %*% beta_d0
  slope <- exp(slope_pop + traits %*% beta_slope)
  str <- exp(str_pop)
  alpha <- exp(alpha_pop + traits %*% beta_alpha)
  
  pop_eff <- exp(as.matrix(df) %*% t(beta_pop))
  
  recr_n <- c(str)*plogis(c(slope)*(log(size) - c(d0)))*size^c(alpha) * c(pop_eff)
  
  return(recr_n)
  
}


# IPMs -------------------------------------------------------------------------

# function to construct IPM and extract lambda for trait/environment combination

get_lambda <- function(x, sc_pars, ss_pars, g_pars, r_pars, s) {
  
  traits <- as.numeric(x[1:3])
  env <- as.numeric(x[4])
  
  # create the mesh points
  L <- 2.5
  U <- 100
  m <- 1000
  h <- (U-L)/m
  y <- L + (1:m)*h - h/2
  
  # growth kernel
  G <- outer(y, y, g_z1z, traits = traits, env = env, upper = max(y), pars = g_pars, s = s)
  for(i in 1:(m-2)) {
    G[(i+1):(m-1),i] <- G[(i+1):(m-1),i] / sum(G[(i+1):(m-1),i]) * (1-G[i,i]-G[m,i])
  }
  G[m,m-1] <- 1-G[m-1,m-1]
  G[m,m] <- 1
  
  # survival kernel
  S <- sapply(y, function(i) s_z(i, traits = traits, env = env, 
                                 sapling_pars = ss_pars, 
                                 canopy_pars = sc_pars, s = s))
  
  # growth/survival kernel
  matP <- sapply(1:nrow(G), function(i) G[,i]*S[i])
  
  # recruitment kernel
  shape <- r_pars$shape[s]
  scale <- r_pars$scale[s]
  hu <- r_pars$hu[s]
  
  matF <- matrix(0, nrow = m, ncol = m)
  recr_n <- f_z(y, traits = traits, env = env, pars = r_pars, s = s)
  size_dist <- dgamma(y-2.54, shape = shape, scale = scale)
  size_dist[1] <- 0
  size_dist <- (1-hu) * size_dist/sum(size_dist) 
  size_dist[1] <- hu
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
vals <- c(-2, -sqrt(2), 0, sqrt(2), 2)
grid <- expand.grid(wd = vals, 
                    sla = vals, 
                    hmax = vals, 
                    mat = transform(c(5, 10, 15), tf$env$mat)) %>%
  as_tibble() %>%
  filter(sqrt(wd^2 + sla^2 + hmax^2) %in% c(0, 2))


## Calculate lambdas ---------------------------------------------------------

library(doParallel)
cl <- makeCluster(32)
registerDoParallel(cl)

n <- 800        # number of samples
set.seed(42)
samples <- sample(1:4000, n)

start <- Sys.time()

lambda <- foreach(i=1:nrow(grid)) %:%
  foreach(s=1:n, 
          .export = c('grid', 'get_lambda', 'g_z1z', 's_z_sapling', 's_z_canopy', 's_z', 'f_z'), 
          .packages=c('plyr', 'dplyr'),
          .combine = "c") %dopar% {
    get_lambda(grid[i,],
               sc_pars = surv_canopy_pars,
               ss_pars = surv_sapling_pars,
               g_pars = growth_pars,
               r_pars = recr_pars, 
               s = samples[s])
  }

Sys.time() - start

stopCluster(cl)

lambda <- do.call(rbind, lambda)


# save results ------------------------------------------------------------

saveRDS(lambda, "results/ipm/fitness_uncertainty.rds")
