library(tidyverse)
library(ggpubr)

source("code/transformer.r")
tf <- readRDS("data/transformers.rds")


# Setup -------------------------------------------------------------------

## Date prep ----

### Load posteriors ----
surv_canopy_pars <- readRDS("results/demog_pars/surv_canopy_pars.rds")
surv_sapling_pars <- readRDS("results/demog_pars/surv_sapling_pars.rds") 
growth_pars <- readRDS("results/demog_pars/growth_pars.rds")
recr_pars <- readRDS("results/demog_pars/recruitment_pars.rds")


### Load data ----
data <- read_csv("data/growth_data_train.csv")

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

### Predictors ----
mat_levs <- c(5, 15) %>% transform(tf$env$mat)

wd_r <- quantile(d$wd, c(0.01, 0.99))
wd_seq <- seq(wd_r[1], wd_r[2], length.out = 40)
wd <- backtransform(wd_seq, tf$trait$wood_density_log, log =  T)

sla_r <- quantile(d$sla, c(0.01, 0.99))
sla_seq <- seq(sla_r[1], sla_r[2], length.out = 40)
sla <- backtransform(sla_seq, tf$trait$sla)

hmax_r <- quantile(d$hmax, c(0.01, 0.99))
hmax_seq <- seq(hmax_r[1], hmax_r[2], length.out = 40)
hmax <- backtransform(hmax_seq, tf$trait$hmax_log, log =  T)

size_raw <- c(5, 15, 45)
size_levs <- size_raw %>% transform(tf$tree$prevdia)
logsize_levs <- size_raw %>% transform(tf$tree$log_prevdia, log = T)

size_levs_c <- size_raw[2:3] %>% transform(tf$canopy$prevdia)
logsize_levs_c <- size_raw[2:3] %>% transform(tf$canopy$log_prevdia, log = T)

logsize_levs_s <- size_raw[1] %>% transform(tf$sapling$log_prevdia, log = T)

## Plot setup----
### Labels----
wd_lab <- expression(paste("Wood density (g cm"^{-3}*")"))
sla_lab <- expression(paste("Specific leaf area (mm"^{2}*" mg"^{-1}*")"))
hmax_lab <- "Maximum height (m)"

surv_lab <- expression(paste("Survival\nprobability (yr"^{-1}*")"))
growth_lab <- expression(paste("Growth (mm yr"^{-1}*")"))
recruit_lab <- expression(paste("\nRecruitment (recruits ind"^{-1}*" yr"^{-1}*")"))

title_size <- 10
label_size <- 8

### Theme----

set_theme <- theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = title_size),
        axis.text.x = element_text(size = label_size),
        axis.title.y = element_text(size = title_size),
        axis.text.y = element_text(size = label_size),
        strip.text.x = element_text(size = title_size, margin = margin(t=2, b=2, unit="pt")), 
        legend.justification = "left",
        legend.title = element_text(size = title_size),
        legend.text = element_text(size = label_size), 
        legend.key.height = unit(0.4, "cm"))

# Survival ----

## Posterior predictions ----

yrs <- 5

### Canopy ----

# canopy, wood density
df_sc_wd <- expand.grid(intercept = 1, 
                  mat = mat_levs, 
                  size = size_levs_c,
                  wd = wd_seq) %>%
  as_tibble() %>%
  mutate(log_size = size %>% backtransform(tf$canopy$prevdia) %>% transform(tf$canopy$log_prevdia, log = T)) %>%
  mutate(wd_sq = wd^2, 
         mat_wd = mat*wd, 
         mat_wd_sq = mat*wd_sq,
         size_wd = wd*size,
         log_size_wd = wd*log_size, 
         trait = "wd")

df_sc_sla <- expand.grid(intercept = 1, 
                        mat = mat_levs, 
                        size = size_levs_c,
                        sla = sla_seq) %>%
  as_tibble() %>%
  mutate(log_size = size %>% backtransform(tf$canopy$prevdia) %>% transform(tf$canopy$log_prevdia, log = T)) %>%
  mutate(sla_sq = sla^2, 
         mat_sla = mat*sla, 
         mat_sla_sq = mat*sla_sq,
         size_sla = sla*size,
         log_size_sla = sla*log_size, 
         trait = "sla")

df_sc_hmax <- expand.grid(intercept = 1, 
                         mat = mat_levs, 
                         size = size_levs_c,
                         hmax = hmax_seq) %>%
  as_tibble() %>%
  mutate(log_size = size %>% backtransform(tf$canopy$prevdia) %>% transform(tf$canopy$log_prevdia, log = T)) %>%
  mutate(hmax_sq = hmax^2, 
         mat_hmax = mat*hmax, 
         mat_hmax_sq = mat*hmax_sq,
         size_hmax = hmax*size,
         log_size_hmax = hmax*log_size, 
         trait = "hmax")

df_sc <- bind_rows(df_sc_wd, df_sc_sla, df_sc_hmax) %>% 
  mutate_all(~replace_na(., 0)) %>%
  dplyr::select(trait, everything())

beta <- surv_canopy_pars[,names(df_sc)[-1]] %>% as.matrix()

sc <- sapply(1:nrow(beta), function(i) plogis(as.matrix(df_sc[,-1]) %*% beta[i,]))
sc <- sc^yrs

df_sc2 <- df_sc %>%
mutate(survival = apply(sc, 1, median),
       lo = apply(sc, 1, quantile, 0.05),
       hi = apply(sc, 1, quantile, 0.95),
       wd = backtransform(wd, tf$trait$wood_density_log, log = T), 
       sla = backtransform(sla, tf$trait$sla),
       hmax = backtransform(hmax, tf$trait$hmax_log, log = T), 
       mat = round(backtransform(mat, tf$env$mat), 1)) %>%
  arrange(wd) %>%
  mutate(size = backtransform(size, tf$canopy$prevdia)) %>%
  dplyr::select(trait, wd, sla, hmax, mat, size, survival, lo, hi)


### Saplings ----
df_ss_wd <- expand.grid(intercept = 1, 
                     mat = mat_levs, 
                     logsize = logsize_levs_s,
                     wd = wd_seq) %>%
  as_tibble() %>%
  mutate(wd_sq = wd^2, 
         mat_wd = mat*wd, 
         mat_wd_sq = mat*wd_sq,
         logsize_wd = wd*logsize, 
         trait = "wd")

df_ss_sla <- expand.grid(intercept = 1, 
                        mat = mat_levs, 
                        logsize = logsize_levs_s,
                        sla = sla_seq) %>%
  as_tibble() %>%
  mutate(sla_sq = sla^2, 
         mat_sla = mat*sla, 
         mat_sla_sq = mat*sla_sq,
         logsize_sla = sla*logsize, 
         trait = "sla")

df_ss_hmax <- expand.grid(intercept = 1, 
                        mat = mat_levs, 
                        logsize = logsize_levs_s,
                        hmax = hmax_seq) %>%
  as_tibble() %>%
  mutate(hmax_sq = hmax^2, 
         mat_hmax = mat*hmax, 
         mat_hmax_sq = mat*hmax_sq,
         logsize_hmax = hmax*logsize, 
         trait = "hmax")

df_ss <- bind_rows(df_ss_wd, df_ss_sla, df_ss_hmax) %>% 
  mutate_all(~replace_na(., 0)) %>%
  dplyr::select(trait, everything())

beta <- surv_sapling_pars[,names(df_ss)[-1]] %>% as.matrix()

ss <- sapply(1:nrow(beta), function(i) plogis(as.matrix(df_ss[,-1]) %*% beta[i,]))
ss <- ss^yrs

df_ss2 <- df_ss %>%
  mutate(survival = apply(ss, 1, median),
         lo = apply(ss, 1, quantile, 0.05),
         hi = apply(ss, 1, quantile, 0.95),
         wd = backtransform(wd, tf$trait$wood_density_log, log = T), 
         sla = backtransform(sla, tf$trait$sla),
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T), 
         mat = round(backtransform(mat, tf$env$mat), 1)) %>%
  arrange(wd) %>%
  mutate(size = backtransform(logsize, tf$sapling$log_prevdia, log = T)) %>%
  dplyr::select(trait, wd, sla, hmax, mat, size, survival, lo, hi)

### Combine ---
df_surv <- rbind(df_sc2, df_ss2) %>% 
  mutate(size = paste(size, "cm") %>% fct_reorder(size), 
         mat = paste0(mat, "°C") %>% fct_reorder(mat))

## Plots ----
s1 <- df_surv %>%
  filter(trait == "wd") %>%
  ggplot(aes(x = wd, y = survival, ymin = lo, ymax = hi, group = size)) +
  geom_path(aes(color = factor(size))) +
  geom_ribbon(aes(fill = factor(size)), alpha = 0.2) +
  facet_wrap(~mat) + 
  scale_color_viridis_d() + scale_fill_viridis_d() +
  set_theme +
  labs(y = "\nSurvival probability (5 yr)", x = wd_lab, color = "Diameter", fill ="Diameter")

s2 <- df_surv %>%
  filter(trait == "sla") %>%
  ggplot(aes(x = sla, y = survival, ymin = lo, ymax = hi, group = size)) +
  geom_path(aes(color = factor(size))) +
  geom_ribbon(aes(fill = factor(size)), alpha = 0.2) +
  facet_wrap(~mat) + 
  scale_color_viridis_d() + scale_fill_viridis_d() +
  set_theme +
  labs(y = "\nSurvival probability (5 yr)", x = sla_lab, color = "Diameter (cm)", fill ="Diameter (cm)")

s3 <- df_surv %>%
  filter(trait == "hmax") %>%
  ggplot(aes(x = hmax, y = survival, ymin = lo, ymax = hi, group = size)) +
  geom_path(aes(color = factor(size))) +
  geom_ribbon(aes(fill = factor(size)), alpha = 0.2) +
  facet_wrap(~mat) + 
  scale_color_viridis_d() + scale_fill_viridis_d() +
  set_theme +
  labs(y = "\nSurvival probability (5 yr)", x = hmax_lab, color = "Diameter (cm)", fill ="Diameter (cm)")



# Growth ----

## Posterior predictions ----

df_g1 <- expand.grid(intercept = 1, 
                     mat = mat_levs, 
                     size = size_levs,
                     t1 = wd_seq) %>%
  as_tibble() %>%
  mutate(log_size = size %>% backtransform(tf$tree$prevdia) %>% transform(tf$tree$log_prevdia, log = T)) %>%
  mutate(t1_sq = t1^2, 
         mat_t1 = mat*t1, 
         mat_t1_sq = mat*t1_sq,
         size_t1 = t1*size,
         log_size_t1 = t1*log_size, 
         trait = "wd")

df_g2 <- expand.grid(intercept = 1, 
                     mat = mat_levs, 
                     size = size_levs,
                     t2 = sla_seq) %>%
  as_tibble() %>%
  mutate(log_size = size %>% backtransform(tf$tree$prevdia) %>% transform(tf$tree$log_prevdia, log = T)) %>%
  mutate(t2_sq = t2^2, 
         mat_t2 = mat*t2, 
         mat_t2_sq = mat*t2_sq,
         size_t2 = t2*size,
         log_size_t2 = t2*log_size, 
         trait = "sla")

df_g3 <- expand.grid(intercept = 1, 
                     mat = mat_levs, 
                     size = size_levs,
                     t3 = hmax_seq) %>%
  as_tibble() %>%
  mutate(log_size = size %>% backtransform(tf$tree$prevdia) %>% transform(tf$tree$log_prevdia, log = T)) %>%
  mutate(t3_sq = t3^2, 
         mat_t3 = mat*t3, 
         mat_t3_sq = mat*t3_sq,
         size_t3 = t3*size,
         log_size_t3 = t3*log_size, 
         trait = "hmax")

df_g <- bind_rows(df_g1, df_g2, df_g3) %>% 
  mutate_all(~replace_na(., 0)) %>%
  dplyr::select(trait, everything())

beta <- growth_pars$beta_pop[,names(df_g)[-1]] %>% as.matrix()

gr <- sapply(1:nrow(beta), function(i) exp(as.matrix(df_g[,-1]) %*% beta[i,]))

df_gr <- df_g %>%
  mutate(growth = apply(gr, 1, mean),
         lo = apply(gr, 1, quantile, 0.05),
         hi = apply(gr, 1, quantile, 0.95),
         wd = backtransform(t1, tf$trait$wood_density_log, log = T), 
         sla = backtransform(t2, tf$trait$sla),
         hmax = backtransform(t3, tf$trait$hmax_log, log = T), 
         mat = round(backtransform(mat, tf$env$mat), 1) %>% paste0("°C") %>% fct_reorder(mat),
         size = backtransform(size, tf$tree$prevdia)) %>%
  dplyr::select(trait, wd, sla, hmax, mat, size, growth, lo, hi)

## Plots ----

g1 <- df_gr %>%
  filter(trait == "wd") %>%
  ggplot(aes(x = wd, y = growth, ymin = lo, ymax = hi, group = size)) +
  geom_path(aes(color = factor(size))) +
  geom_ribbon(aes(fill = factor(size)), alpha = 0.2) +
  facet_wrap(~mat) + 
  scale_color_viridis_d() + scale_fill_viridis_d() +
  set_theme +
  labs(y = growth_lab, x = wd_lab, color = "Diameter (cm)", fill ="Diameter (cm)")

g2 <- df_gr %>%
  filter(trait == "sla") %>%
  ggplot(aes(x = sla, y = growth, ymin = lo, ymax = hi, group = size)) +
  geom_path(aes(color = factor(size))) +
  geom_ribbon(aes(fill = factor(size)), alpha = 0.2) +
  facet_wrap(~mat) + 
  scale_color_viridis_d() + scale_fill_viridis_d() +
  set_theme +
  labs(y = growth_lab, x = sla_lab)

g3 <- df_gr %>%
  filter(trait == "hmax") %>%
  ggplot(aes(x = hmax, y = growth, ymin = lo, ymax = hi, group = size)) +
  geom_path(aes(color = factor(size))) +
  geom_ribbon(aes(fill = factor(size)), alpha = 0.2) +
  facet_wrap(~mat) + 
  scale_color_viridis_d() + scale_fill_viridis_d() +
  set_theme +
  labs(y = growth_lab, x = hmax_lab)

# Growth and survival figure ----
ggarrange(s1, g1, s2, g2, s3, g3, 
          common.legend = T, legend = "bottom", nrow = 3, ncol = 2,
          labels = c("a", "d", "b", "e", "c", "f"))

ggsave("results/figures/ed_fig2.pdf", width = 8, height = 7, units = "in")


# Recruitment ----

## Posterior predictions ----

size_r <- c(5, 15, 45)

df_r1 <- expand.grid(size = size_r, 
                     mat = mat_levs, 
                     wd = wd_seq) %>%
  as_tibble() %>%
  mutate(wd_sq = wd^2, 
         mat_wd = mat*wd, 
         mat_wd_sq = mat*wd_sq,
         trait = "wd")

df_r2 <- expand.grid(size = size_r, 
                     mat = mat_levs, 
                     sla = sla_seq) %>%
  as_tibble() %>%
  mutate(sla_sq = sla^2, 
         mat_sla = mat*sla, 
         mat_sla_sq = mat*sla_sq,
         trait = "sla")

df_r3 <- expand.grid(size = size_r, 
                     mat = mat_levs, 
                     hmax = hmax_seq) %>%
  as_tibble() %>%
  mutate(hmax_sq = hmax^2, 
         mat_hmax = mat*hmax, 
         mat_hmax_sq = mat*hmax_sq,
         trait = "hmax")

df_r <- bind_rows(df_r1, df_r2, df_r3) %>% 
  mutate_all(~replace_na(., 0)) %>%
  dplyr::select(trait, everything())

traits <- df_r %>% dplyr::select(wd, sla, hmax) %>% as.matrix()

beta_pop <- recr_pars$beta_pop[,names(df_r)[-(1:2)]] %>% as.matrix()

samples <- length(recr_pars$d0_pop)

recr <- matrix(NA, nrow = 720, ncol = samples)
for(s in 1:samples) {
  d0 <- recr_pars$d0_pop[s] + traits %*% recr_pars$beta_d0[s,]
  slope <- exp(recr_pars$slope_pop[s] + traits %*% recr_pars$beta_slope[s,])
  str <- exp(recr_pars$str_pop[s])
  alpha <- exp(recr_pars$alpha_pop[s] + traits %*% recr_pars$beta_alpha[s,])
  
  pop_eff <- exp(as.matrix(df_r[,-(1:2)]) %*% beta_pop[s,])
  
  size <- df_r$size/30
  recr[,s] <- str * plogis(slope*(log(size) - d0)) * size^alpha * pop_eff
}

df_recr <- df_r %>%
  mutate(recruitment = apply(recr, 1, mean),
         lo = apply(recr, 1, quantile, 0.05),
         hi = apply(recr, 1, quantile, 0.95),
         wd = backtransform(wd, tf$trait$wood_density_log, log = T), 
         sla = backtransform(sla, tf$trait$sla),
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T), 
         mat = round(backtransform(mat, tf$env$mat), 1) %>% paste0("°C") %>% fct_reorder(mat), 
         size = paste(size, "cm") %>% fct_reorder(size)) %>%
  dplyr::select(trait, wd, sla, hmax, mat, size, recruitment, lo, hi)

## Plots ----

r1 <- df_recr %>%
  filter(trait == "wd") %>%
  ggplot(aes(x = wd, y = recruitment, ymin = lo, ymax = hi, group = size)) +
  geom_path(aes(color = factor(size))) +
  geom_ribbon(aes(fill = factor(size)), alpha = 0.2) +
  facet_wrap(~mat) + 
  scale_color_viridis_d() + scale_fill_viridis_d() +
  set_theme +
  labs(y = recruit_lab, x = wd_lab, color = "Diameter", fill ="Diameter") +
  scale_y_log10()

r2 <- df_recr %>%
  filter(trait == "sla") %>%
  ggplot(aes(x = sla, y = recruitment, ymin = lo, ymax = hi, group = size)) +
  geom_path(aes(color = factor(size))) +
  geom_ribbon(aes(fill = factor(size)), alpha = 0.2) +
  facet_wrap(~mat) + 
  scale_color_viridis_d() + scale_fill_viridis_d() +
  set_theme +
  labs(y = recruit_lab, x = sla_lab) +
  scale_y_log10()

r3 <- df_recr %>%
  filter(trait == "hmax") %>%
  ggplot(aes(x = hmax, y = recruitment, ymin = lo, ymax = hi, group = size)) +
  geom_path(aes(color = factor(size))) +
  geom_ribbon(aes(fill = factor(size)), alpha = 0.2) +
  facet_wrap(~mat) + 
  scale_color_viridis_d() + scale_fill_viridis_d() +
  set_theme +
  labs(y = recruit_lab, x = hmax_lab) +
  scale_y_log10()

ggarrange(r1, r2, r3, common.legend = T, legend = "bottom", ncol = 1,
          labels = "auto")

ggsave("results/figures/ed_fig3.pdf", width = 5, height = 8, 
       units = "in")

