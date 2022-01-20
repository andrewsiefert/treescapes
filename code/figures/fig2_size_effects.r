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
mat_levs <- 10 %>% transform(tf$env$mat)

wd_levs <- c(-2, 0, 2)
sla_levs <- c(-2, 0, 2)
hmax_levs <- c(-2, 0, 2)

size_levs <- transform(3:60, tf$tree$prevdia)
logsize_levs <- transform(3:60, tf$tree$log_prevdia, log = T)

size_levs_c <- transform(13:60, tf$canopy$prevdia)
logsize_levs_c <- transform(13:60, tf$canopy$log_prevdia, log = T)

logsize_levs_s <- transform(3:12, tf$sapling$log_prevdia, log = T)

## Plot setup----
### Labels----
wd_lab <- expression(paste("Wood density (g cm"^{-3}*")"))
sla_lab <- expression(paste("Specific leaf area (mm"^{2}*" mg"^{-1}*")"))
hmax_lab <- "Maximum height (m)"

s_surv_lab <- expression(atop(Sapling~survival, probability~(yr^-1)))
c_surv_lab <- expression(atop(Canopy~tree~survival, probability~(yr^-1)))
growth_lab <- expression(atop(Diameter~growth, (mm~yr^-1)))
recruit_lab <- expression(atop(Recruitment, (recruits~ind^-1~yr^-1)))

title_size <- 10
label_size <- 8

### Theme----

theme_top <- theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = label_size),
        axis.title.y = element_text(size = title_size),
        axis.text.y = element_text(size = label_size),
        strip.text.x = element_text(size = title_size, margin = margin(t=2, b=2, unit="pt")), 
        legend.justification = "left",
        legend.title = element_text(size = title_size),
        legend.text = element_text(size = label_size), 
        legend.key.height = unit(0.4, "cm"))

theme_middle <- theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = label_size),
        axis.title.y = element_text(size = title_size),
        axis.text.y = element_text(size = label_size),
        strip.background = element_blank(),
        strip.text.x = element_blank(),         
        legend.justification = "left",
        legend.title = element_text(size = title_size),
        legend.text = element_text(size = label_size), 
        legend.key.height = unit(0.4, "cm"))

theme_bottom <- theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = title_size),
        axis.text.x = element_text(size = label_size),
        axis.title.y = element_text(size = title_size),
        axis.text.y = element_text(size = label_size),
        strip.background = element_blank(),
        strip.text.x = element_blank(),    
        legend.justification = "left",
        legend.title = element_text(size = title_size),
        legend.text = element_text(size = label_size), 
        legend.key.height = unit(0.4, "cm"))

# Canopy survival ----

## Posterior predictions ----

yrs <- 5

# canopy, wood density
df_sc_wd <- expand.grid(intercept = 1, 
                  mat = mat_levs, 
                  size = size_levs_c,
                  wd = wd_levs) %>%
  as_tibble() %>%
  mutate(log_size = size %>% backtransform(tf$canopy$prevdia) %>% transform(tf$canopy$log_prevdia, log = T)) %>%
  mutate(wd_sq = wd^2, 
         mat_wd = mat*wd, 
         mat_wd_sq = mat*wd_sq,
         size_wd = wd*size,
         log_size_wd = wd*log_size, 
         trait = "Wood density")

df_sc_sla <- expand.grid(intercept = 1, 
                        mat = mat_levs, 
                        size = size_levs_c,
                        sla = sla_levs) %>%
  as_tibble() %>%
  mutate(log_size = size %>% backtransform(tf$canopy$prevdia) %>% transform(tf$canopy$log_prevdia, log = T)) %>%
  mutate(sla_sq = sla^2, 
         mat_sla = mat*sla, 
         mat_sla_sq = mat*sla_sq,
         size_sla = sla*size,
         log_size_sla = sla*log_size, 
         trait = "Specific leaf area")

df_sc_hmax <- expand.grid(intercept = 1, 
                         mat = mat_levs, 
                         size = size_levs_c,
                         hmax = hmax_levs) %>%
  as_tibble() %>%
  mutate(log_size = size %>% backtransform(tf$canopy$prevdia) %>% transform(tf$canopy$log_prevdia, log = T)) %>%
  mutate(hmax_sq = hmax^2, 
         mat_hmax = mat*hmax, 
         mat_hmax_sq = mat*hmax_sq,
         size_hmax = hmax*size,
         log_size_hmax = hmax*log_size, 
         trait = "Maximum height")

df_sc <- bind_rows(df_sc_wd, df_sc_sla, df_sc_hmax) %>% 
  mutate_all(~replace_na(., 0)) %>%
  select(trait, everything()) %>%
  mutate(trait = fct_relevel(trait, "Wood density", "Specific leaf area"))

beta <- surv_canopy_pars[,names(df_sc)[-1]] %>% as.matrix()

sc <- sapply(1:nrow(beta), function(i) plogis(as.matrix(df_sc[,-1]) %*% beta[i,]))
sc <- sc^yrs

df_sc2 <- df_sc %>%
mutate(survival = apply(sc, 1, median),
       lo = apply(sc, 1, quantile, 0.05),
       hi = apply(sc, 1, quantile, 0.95),
       #wd = backtransform(wd, tf$trait$wood_density_log, log = T), 
       #sla = backtransform(sla, tf$trait$sla),
       #hmax = backtransform(hmax, tf$trait$hmax_log, log = T), 
       mat = round(backtransform(mat, tf$env$mat), 1)) %>%
  arrange(wd) %>%
  mutate(size = backtransform(size, tf$canopy$prevdia)) %>%
  select(trait, wd, sla, hmax, mat, size, survival, lo, hi) %>%
  mutate(trait_val = ifelse(trait == "Wood density", wd, 
                            ifelse(trait == "Specific leaf area", sla, hmax)), 
         trait_val = ifelse(trait_val == -2, "Low", 
                            ifelse(trait_val == 0, "Average", "High")) %>% 
           fct_reorder(trait_val))%>%
  arrange(size)

### Plot ----
surv_canopy <- df_sc2 %>%
  ggplot(aes(x = size, y = survival, ymin = lo, ymax = hi, group = trait_val)) +
  geom_path(aes(color = factor(trait_val))) +
  geom_ribbon(aes(fill = factor(trait_val)), alpha = 0.2) +
  facet_wrap(~trait) + 
  scale_color_viridis_d() + scale_fill_viridis_d() +
  theme_middle +
  labs(y = c_surv_lab, x = "Diameter (cm)", color = "Trait value", fill = "Trait value")


# Sapling survival ----

### Posterior predictions ----
df_ss_wd <- expand.grid(intercept = 1, 
                     mat = mat_levs, 
                     logsize = logsize_levs_s,
                     wd = wd_levs) %>%
  as_tibble() %>%
  mutate(wd_sq = wd^2, 
         mat_wd = mat*wd, 
         mat_wd_sq = mat*wd_sq,
         logsize_wd = wd*logsize, 
         trait = "Wood density")

df_ss_sla <- expand.grid(intercept = 1, 
                        mat = mat_levs, 
                        logsize = logsize_levs_s,
                        sla = sla_levs) %>%
  as_tibble() %>%
  mutate(sla_sq = sla^2, 
         mat_sla = mat*sla, 
         mat_sla_sq = mat*sla_sq,
         logsize_sla = sla*logsize, 
         trait = "Specific leaf area")

df_ss_hmax <- expand.grid(intercept = 1, 
                        mat = mat_levs, 
                        logsize = logsize_levs_s,
                        hmax = hmax_levs) %>%
  as_tibble() %>%
  mutate(hmax_sq = hmax^2, 
         mat_hmax = mat*hmax, 
         mat_hmax_sq = mat*hmax_sq,
         logsize_hmax = hmax*logsize, 
         trait = "Maximum height")

df_ss <- bind_rows(df_ss_wd, df_ss_sla, df_ss_hmax) %>% 
  mutate_all(~replace_na(., 0)) %>%
  select(trait, everything()) %>%
  mutate(trait = fct_relevel(trait, "Wood density", "Specific leaf area"))

beta <- surv_sapling_pars[,names(df_ss)[-1]] %>% as.matrix()

ss <- sapply(1:nrow(beta), function(i) plogis(as.matrix(df_ss[,-1]) %*% beta[i,]))
ss <- ss^yrs

df_ss2 <- df_ss %>%
  mutate(survival = apply(ss, 1, median),
         lo = apply(ss, 1, quantile, 0.05),
         hi = apply(ss, 1, quantile, 0.95),
         #wd = backtransform(wd, tf$trait$wood_density_log, log = T), 
         #sla = backtransform(sla, tf$trait$sla),
         #hmax = backtransform(hmax, tf$trait$hmax_log, log = T), 
         mat = round(backtransform(mat, tf$env$mat), 1)) %>%
  arrange(wd) %>%
  mutate(size = backtransform(logsize, tf$sapling$log_prevdia, log = T)) %>%
  select(trait, wd, sla, hmax, mat, size, survival, lo, hi) %>%
  mutate(trait_val = ifelse(trait == "Wood density", wd, 
                            ifelse(trait == "Specific leaf area", sla, hmax)), 
         trait_val = ifelse(trait_val == -2, "Low", 
                            ifelse(trait_val == 0, "Average", "High")) %>% 
           fct_reorder(trait_val))%>%
  arrange(size)

## Plot ----
surv_sapling <- df_ss2 %>%
  ggplot(aes(x = size, y = survival, ymin = lo, ymax = hi, group = trait_val)) +
  geom_path(aes(color = factor(trait_val))) +
  geom_ribbon(aes(fill = factor(trait_val)), alpha = 0.2) +
  facet_wrap(~trait) + 
  scale_color_viridis_d() + scale_fill_viridis_d() +
  theme_top +
  labs(y = s_surv_lab, x = "Diameter (cm)", color = "Trait value", fill = "Trait value")


# Growth ----

## Posterior predictions ----

df_g1 <- expand.grid(intercept = 1, 
                     mat = mat_levs, 
                     size = size_levs,
                     t1 = wd_levs) %>%
  as_tibble() %>%
  mutate(log_size = size %>% backtransform(tf$tree$prevdia) %>% transform(tf$tree$log_prevdia, log = T)) %>%
  mutate(t1_sq = t1^2, 
         mat_t1 = mat*t1, 
         mat_t1_sq = mat*t1_sq,
         size_t1 = t1*size,
         log_size_t1 = t1*log_size, 
         trait = "Wood density")

df_g2 <- expand.grid(intercept = 1, 
                     mat = mat_levs, 
                     size = size_levs,
                     t2 = sla_levs) %>%
  as_tibble() %>%
  mutate(log_size = size %>% backtransform(tf$tree$prevdia) %>% transform(tf$tree$log_prevdia, log = T)) %>%
  mutate(t2_sq = t2^2, 
         mat_t2 = mat*t2, 
         mat_t2_sq = mat*t2_sq,
         size_t2 = t2*size,
         log_size_t2 = t2*log_size, 
         trait = "Specific leaf area")

df_g3 <- expand.grid(intercept = 1, 
                     mat = mat_levs, 
                     size = size_levs,
                     t3 = hmax_levs) %>%
  as_tibble() %>%
  mutate(log_size = size %>% backtransform(tf$tree$prevdia) %>% transform(tf$tree$log_prevdia, log = T)) %>%
  mutate(t3_sq = t3^2, 
         mat_t3 = mat*t3, 
         mat_t3_sq = mat*t3_sq,
         size_t3 = t3*size,
         log_size_t3 = t3*log_size, 
         trait = "Maximum height")

df_g <- bind_rows(df_g1, df_g2, df_g3) %>% 
  mutate_all(~replace_na(., 0)) %>%
  select(trait, everything()) %>%
  mutate(trait = fct_relevel(trait, "Wood density", "Specific leaf area"))

beta <- growth_pars$beta_pop[,names(df_g)[-1]] %>% as.matrix()

gr <- sapply(1:nrow(beta), function(i) exp(as.matrix(df_g[,-1]) %*% beta[i,]))

df_gr <- df_g %>%
  mutate(growth = apply(gr, 1, mean),
         lo = apply(gr, 1, quantile, 0.05),
         hi = apply(gr, 1, quantile, 0.95),
         #wd = backtransform(t1, tf$trait$wood_density_log, log = T), 
         #sla = backtransform(t2, tf$trait$sla),
         #hmax = backtransform(t3, tf$trait$hmax_log, log = T), 
         mat = round(backtransform(mat, tf$env$mat), 1) %>% paste0("°C") %>% fct_reorder(mat),
         size = backtransform(size, tf$tree$prevdia)) %>%
  rename(wd = t1, sla = t2, hmax = t3) %>%
  select(trait, wd, sla, hmax, mat, size, growth, lo, hi) %>%
  mutate(trait_val = ifelse(trait == "Wood density", wd, 
                            ifelse(trait == "Specific leaf area", sla, hmax)), 
         trait_val = ifelse(trait_val == -2, "Low", 
                            ifelse(trait_val == 0, "Average", "High")) %>% 
           fct_reorder(trait_val))%>%
  arrange(size)

## Plot ----

growth <- df_gr %>%
  ggplot(aes(x = size, y = growth, ymin = lo, ymax = hi, group = trait_val)) +
  geom_path(aes(color = factor(trait_val))) +
  geom_ribbon(aes(fill = factor(trait_val)), alpha = 0.2) +
  facet_wrap(~trait) + 
  scale_color_viridis_d() + scale_fill_viridis_d() +
  theme_middle +
  labs(y = growth_lab, x = "Diameter (cm)", color = "Trait value", fill = "Trait value")


# Recruitment ----

## Posterior predictions ----

size_r <- 5:60

df_r1 <- expand.grid(size = size_r, 
                     mat = mat_levs, 
                     wd = wd_levs) %>%
  as_tibble() %>%
  mutate(wd_sq = wd^2, 
         mat_wd = mat*wd, 
         mat_wd_sq = mat*wd_sq,
         trait = "Wood density")

df_r2 <- expand.grid(size = size_r, 
                     mat = mat_levs, 
                     sla = sla_levs) %>%
  as_tibble() %>%
  mutate(sla_sq = sla^2, 
         mat_sla = mat*sla, 
         mat_sla_sq = mat*sla_sq,
         trait = "Specific leaf area")

df_r3 <- expand.grid(size = size_r, 
                     mat = mat_levs, 
                     hmax = hmax_levs) %>%
  as_tibble() %>%
  mutate(hmax_sq = hmax^2, 
         mat_hmax = mat*hmax, 
         mat_hmax_sq = mat*hmax_sq,
         trait = "Maximum height")

df_r <- bind_rows(df_r1, df_r2, df_r3) %>% 
  mutate_all(~replace_na(., 0)) %>%
  select(trait, everything()) %>%
  mutate(trait = fct_relevel(trait, "Wood density", "Specific leaf area"))

traits <- df_r %>% select(wd, sla, hmax) %>% as.matrix()

beta_pop <- recr_pars$beta_pop[,names(df_r)[-(1:2)]] %>% as.matrix()

n <- nrow(df_r)
samples <- length(recr_pars$d0_pop)

repro <- matrix(NA, nrow = n, ncol = samples)
recr <- matrix(NA, nrow = n, ncol = samples)
recr_n <- matrix(NA, nrow = n, ncol = samples)
for(s in 1:samples) {
  d0 <- recr_pars$d0_pop[s] + traits %*% recr_pars$beta_d0[s,]
  slope <- exp(recr_pars$slope_pop[s] + traits %*% recr_pars$beta_slope[s,])
  str <- exp(recr_pars$str_pop[s])
  alpha <- exp(recr_pars$alpha_pop[s] + traits %*% recr_pars$beta_alpha[s,])
  
  pop_eff <- exp(as.matrix(df_r[,-(1:2)]) %*% beta_pop[s,])
  
  size <- df_r$size/30
  repro_p <- plogis(slope*(log(size) - d0))
  recr_c <- str * size^alpha * pop_eff
  
  repro[,s] <- repro_p
  recr[,s] <- recr_c
  recr_n[,s] <- repro_p*recr_c
}

df_recr <- df_r %>%
  mutate(recruitment = apply(recr_n, 1, median),
         recr_lo = apply(recr_n, 1, quantile, 0.05),
         recr_hi = apply(recr_n, 1, quantile, 0.95),
         repro_p = apply(repro, 1, median),
         repro_lo = apply(repro, 1, quantile, 0.05),
         repro_hi = apply(repro, 1, quantile, 0.95)) %>%
  select(trait, wd, sla, hmax, mat, size, recruitment, recr_lo, recr_hi, repro_p, repro_lo, repro_hi) %>%
  mutate(trait_val = ifelse(trait == "Wood density", wd, 
                            ifelse(trait == "Specific leaf area", sla, hmax)), 
         trait_val = ifelse(trait_val == -2, "Low", 
                            ifelse(trait_val == 0, "Average", "High")) %>% 
           fct_reorder(trait_val))%>%
  arrange(size)

## Plot ----

### Reproductive status ----
repro_status <- df_recr %>%
  ggplot(aes(x = size, y = repro_p, ymin = repro_lo, ymax = repro_hi, group = trait_val)) +
  geom_path(aes(color = factor(trait_val))) +
  geom_ribbon(aes(fill = factor(trait_val)), alpha = 0.2) +
  facet_wrap(~trait) + 
  scale_color_viridis_d() + scale_fill_viridis_d() +
  theme_middle +
  labs(y = "Probability of\nbeing reproductive", x = "Diameter (cm)", color = "Trait value", fill = "Trait value") 

### Per-capita recruitment ----
recruitment <- df_recr %>%
  ggplot(aes(x = size, y = recruitment, ymin = recr_lo, ymax = recr_hi, group = trait_val)) +
  geom_path(aes(color = factor(trait_val))) +
  geom_ribbon(aes(fill = factor(trait_val)), alpha = 0.2) +
  facet_wrap(~trait) + 
  scale_color_viridis_d() + scale_fill_viridis_d() +
  theme_bottom +
  labs(y = recruit_lab, x = "Diameter (cm)", color = "Trait value", fill = "Trait value") +
  scale_y_log10()


svg("results/figures/fig2.svg", 7, 8)
gridExtra::grid.arrange(rbind(ggplotGrob(surv_sapling), ggplotGrob(surv_canopy), 
                              ggplotGrob(growth), ggplotGrob(repro_status), ggplotGrob(recruitment)))
dev.off()

