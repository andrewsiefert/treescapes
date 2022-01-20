library(tidyverse)
library(ggpubr)

surv_canopy_in <- readRDS("results/model_checks/canopy_survival_check_in_sample.rds")
surv_sapling_in <- readRDS("results/model_checks/sapling_survival_check_in_sample.rds")

surv_canopy_out <- readRDS("results/model_checks/canopy_survival_check_out_of_sample.rds")
surv_sapling_out <- readRDS("results/model_checks/sapling_survival_check_out_of_sample.rds")

growth_in <- readRDS("results/model_checks/growth_check_in_sample.rds")
growth_out <- readRDS("results/model_checks/growth_check_out_of_sample.rds")

recr_in <- readRDS("results/model_checks/recruitment_check_in_sample.rds")
recr_out <- readRDS("results/model_checks/recruitment_check_out_of_sample.rds")

growth_lab <- expression(paste("Growth (mm yr"^{-1}*")"))


# functions to add stats

r2_label <- function(stat){
  eq <- substitute(italic(R)^2~"="~r2~"("*lo*", "*hi*")", 
                   list(r2 = format(mean(stat), digits = 3),
                        lo = format(unname(quantile(stat, 0.05)), digits = 3),
                        hi = format(unname(quantile(stat, 0.95)), digits = 3)))
  as.character(as.expression(eq))
}

auc_label <- function(stat){
  eq <- substitute("AUC ="~auc~"("*lo*", "*hi*")", 
                   list(auc = format(mean(stat), digits = 3),
                        lo = format(unname(quantile(stat, 0.05)), digits = 3),
                        hi = format(unname(quantile(stat, 0.95)), digits = 3)))
  as.character(as.expression(eq))
}

cover_label <- function(coverage) {
  paste("Coverage =", round(mean(coverage), 2))
}

# Growth ----

## Growth PPC ----
g_pred <- growth_in$growth_pred %>%
  gather("sample", "pred", everything()) %>%
  group_by(sample) %>%
  mutate(cd = cume_dist(pred)) %>%
  ungroup()  %>%
  add_row(sample = paste0("V", 1:100), pred = 0, cd = 0) %>%
  arrange(pred)

g_obs <- growth_in$df %>%
  select(obs) %>%
  mutate(cd = cume_dist(obs)) %>%
  add_row(obs = 0, cd = 0) %>%
  arrange(obs, cd)

g_ppc <- ggplot(g_pred) + 
  geom_path(aes(x = pred, y = cd, group = sample), color = "lightblue", alpha = 0.5) +
  geom_path(data = g_obs, aes(x = obs, y = cd)) +
  scale_y_continuous(limits = c(0,1)) +
  theme_classic() +
  labs(title = "Growth", x = growth_lab, y = "Cumulative probability") 

## Growth in-sample ---- 
g_in <- growth_in$df %>% 
  sample_frac(0.1) %>%
  ggplot(aes(x = pred_mean, y = obs)) + 
  geom_point(size = 0.1, alpha = 0.2, color = 3) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  #coord_fixed(xlim = c(0, 30), ylim = c(0, 30)) +
  theme_classic() +
  labs(title = "", x = "Predicted growth", y = "Observed growth") +
  annotate("text", x = 0, y = Inf, label = r2_label(growth_in$r2), parse = T, hjust = 0, vjust = 1.5, size = 3) +
  annotate("text", x = 0, y = Inf, label = cover_label(growth_in$df$cover90), hjust = 0, vjust = 4, size = 3)


## Growth out-of-sample ----
g_out <- growth_out$df %>% 
  sample_frac(0.1) %>%
  ggplot(aes(x = pred_mean, y = obs)) + 
  geom_point(size = 0.1, alpha = 0.2, color = 2) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  #coord_fixed(xlim = c(0, 30), ylim = c(0, 30)) +
  theme_classic() +
  labs(title = "", x = "Predicted growth", y = "Observed growth") +
  annotate("text", x = 0, y = Inf, label = r2_label(growth_out$r2), parse = T, hjust = 0, vjust = 1.5, size = 3) +
  annotate("text", x = 0, y = Inf, label = cover_label(growth_out$df$cover90), hjust = 0, vjust = 4, size = 3)

# Recruitment ----

## Recruitment PPC ----
r_pred <- recr_in$recr_pred %>%
  as_tibble() %>%
  gather("sample", "pred", V1:V400) %>%
  group_by(sample) %>%
  mutate(cd = cume_dist(pred)) %>%
  ungroup()  %>%
  add_row(sample = paste0("V", 1:400), pred = 0, cd = 0) %>%
  distinct() %>%
  arrange(sample, pred, cd)

r_obs <- recr_in$df %>%
  select(obs) %>%
  mutate(cd = cume_dist(obs)) %>%
  add_row(obs = 0, cd = 0) %>%
  arrange(obs, cd) %>%
  distinct()

r_ppc <- ggplot(r_pred) + 
  geom_step(aes(x = pred, y = cd, group = sample), color = "lightblue", alpha = 0.5) +
  geom_step(data = r_obs, aes(x = obs, y = cd)) +
  scale_x_continuous(limits = c(0,20)) +
  theme_classic() +
  labs(title = "Recruitment", x = "Recruitment (number of recruits)", y = "Cumulative probability")

## Recruitment in-sample ---- 
r_in <- recr_in$df %>% 
  ggplot(aes(x = exp_mean, y = obs)) + 
  geom_jitter(width = 0, size = 0.4, alpha = 0.3, color = 3) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  #coord_fixed(xlim = c(0, 30), ylim = c(0, 30)) +
  theme_classic() +
  labs(title = "", x = "Predicted recruitment", y = "Observed recruitment") +
  annotate("text", x = 0, y = Inf, label = r2_label(recr_in$r2), parse = T, hjust = 0, vjust = 1.5, size = 3) +
  annotate("text", x = 0, y = Inf, label = cover_label(recr_in$df$cover90), hjust = 0, vjust = 4, size = 3)

## Recruitment out-of-sample ----
r_out <- recr_out$df %>% 
  ggplot(aes(x = exp_mean, y = obs)) + 
  geom_jitter(width = 0, size = 0.4, alpha = 0.3, color = 2) +
  geom_abline(aes(intercept = 0, slope = 1)) +
  #coord_fixed(xlim = c(0, 30), ylim = c(0, 30)) +
  theme_classic() +
  labs(title = "", x = "Predicted recruitment", y = "Observed recruitment") +
  annotate("text", x = 0, y = Inf, label = r2_label(recr_out$r2), parse = T, hjust = 0, vjust = 1.5, size = 3) +
  annotate("text", x = 0, y = Inf, label = cover_label(recr_out$df$cover90), hjust = 0, vjust = 4, size = 3)

# Sapling survival ----

## PPC ----
ss_ppc <- ggplot(tibble(surv = surv_sapling_in$post_mean), aes(x = surv)) + 
  geom_histogram(color = "lightblue3", fill = "lightblue2") +
  geom_vline(xintercept = mean(surv_sapling_in$df$obs)) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  theme_classic() +
  labs(title = "Sapling survival", x = "Survival probability", y = "Frequency")

## In-sample ---- 

ssi <- surv_sapling_in$df %>%
  mutate(bin = ntile(theta_mean, 20)) %>%
  group_by(bin) %>%
  summarize(obs = mean(obs),
            pred = mean(theta_mean)) %>%
  ungroup()

ss_in <- surv_sapling_in$df %>% 
  sample_frac(0.2) %>%
  ggplot(aes(x = theta_mean, y = obs)) + 
  geom_jitter(width = 0, size = 0.1, alpha = 0.2, color = 3) +
  geom_point(data = ssi, aes(x = pred, y = obs), size = 3, shape = "x") +
  theme_classic() +
  labs(title = "", x = "Predicted survival probability", y = "Observed survival") +
  annotate("text", x = 0, y = Inf, label = auc_label(surv_sapling_in$auc), parse = T, hjust = 0, vjust = 1.5, size = 3) +
  annotate("text", x = 0, y = Inf, label = cover_label(surv_sapling_in$df$cover90), hjust = 0, vjust = 3.5, size = 3)


## Out-of-sample ---- 

sso <- surv_sapling_out$df %>%
  mutate(bin = ntile(theta_mean, 20)) %>%
  group_by(bin) %>%
  summarize(obs = mean(obs),
            pred = mean(theta_mean)) %>%
  ungroup()

ss_out <- surv_sapling_out$df %>% 
  sample_frac(0.2) %>%
  ggplot(aes(x = theta_mean, y = obs)) + 
  geom_jitter(width = 0, size = 0.1, alpha = 0.5, color = 2) +
  geom_point(data = sso, aes(x = pred, y = obs), size = 3, shape = "x") +
  theme_classic() +
  labs(title = "", x = "Predicted survival probability", y = "Observed survival") +
  annotate("text", x = 0, y = Inf, label = auc_label(surv_sapling_out$auc), parse = T, hjust = 0, vjust = 1.5, size = 3) +
  annotate("text", x = 0, y = Inf, label = cover_label(surv_sapling_out$df$cover90), hjust = 0, vjust = 3.5, size = 3) 

# Canopy survival ----

## PPC ----
sc_ppc <- ggplot(tibble(surv = surv_canopy_in$post_mean), aes(x = surv)) + 
  geom_histogram(color = "lightblue3", fill = "lightblue2") +
  geom_vline(xintercept = mean(surv_canopy_in$df$obs)) +
  scale_y_continuous(expand = c(0,0,0.05,0)) +
  theme_classic() +
  labs(title = "Canopy tree survival", x = "Survival probability", y = "Frequency")

## In-sample ---- 

sci <- surv_canopy_in$df %>%
  mutate(bin = ntile(theta_mean, 20)) %>%
  group_by(bin) %>%
  summarize(obs = mean(obs),
            pred = mean(theta_mean)) %>%
  ungroup()

sc_in <- surv_canopy_in$df %>% 
  sample_frac(0.1) %>%
  ggplot(aes(x = theta_mean, y = obs)) + 
  geom_jitter(width = 0, size = 0.1, alpha = 0.2, color = 3) +
  geom_point(data = sci, aes(x = pred, y = obs), size = 3, shape = "x") +
  theme_classic() +
  labs(title = "", x = "Predicted survival probability", y = "Observed survival") +
  annotate("text", x = 0, y = Inf, label = auc_label(surv_canopy_in$auc), parse = T, hjust = 0, vjust = 1.5, size = 3) +
  annotate("text", x = 0, y = Inf, label = cover_label(surv_canopy_in$df$cover90), hjust = 0, vjust = 3.5, size = 3) 


## Out-of-sample ---- 

sco <- surv_canopy_out$df %>%
  mutate(bin = ntile(theta_mean, 20)) %>%
  group_by(bin) %>%
  summarize(obs = mean(obs),
            pred = mean(theta_mean)) %>%
  ungroup()

sc_out <- surv_canopy_out$df %>% 
  sample_frac(0.1) %>%
  ggplot(aes(x = theta_mean, y = obs)) + 
  geom_jitter(width = 0, size = 0.1, alpha = 0.5, color = 2) +
  geom_point(data = sco, aes(x = pred, y = obs), size = 3, shape = "x") +
  theme_classic() +
  labs(title = "", x = "Predicted survival probability", y = "Observed survival") +
  annotate("text", x = 0, y = Inf, label = auc_label(surv_canopy_out$auc), parse = T, hjust = 0, vjust = 1, size = 3) +
  annotate("text", x = 0, y = Inf, label = cover_label(surv_canopy_out$df$cover90), hjust = 0, vjust = 3, size = 3)
  
# Plot ----
ggarrange(ss_ppc, ss_in, ss_out,
          sc_ppc, sc_in, sc_out,
          g_ppc, g_in, g_out,
          r_ppc, r_in, r_out, 
          nrow = 4, ncol = 3, 
          labels = "auto", 
          label.y = 0.95)

ggsave("results/figures/ed_fig10.jpg", width = 10, height = 12, units = "in", dpi = 300)
