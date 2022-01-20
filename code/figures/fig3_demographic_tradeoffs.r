library(tidyverse)


d <- read_csv("results/demographic_tradeoffs/demographic_tradeoffs.csv")


set_theme <- theme_bw() +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

nice_colors <- scale_color_viridis_c(option = "magma", begin = 0.1, end = 0.9)

wd_lab <- expression(paste("Wood density (g cm"^{-3}*")"))
hmax_lab <- "Maximum height (m)"


## Wood density ----

### Growth-survival ----
wd1 <- d %>% filter(mat == "10°C", trait == "wd", tradeoff == "gs") %>%
  ggplot(aes(x = growth_est, xmin = growth_lo, xmax = growth_hi,
             y = surv_est, ymin = surv_lo, ymax = surv_hi,
             color = wd)) +
  geom_pointrange() +
  geom_errorbarh(height = 0) +
  set_theme +
  nice_colors +
  labs(x = "Growth (5 cm dbh)", y = "Survival (5 cm dbh)", color = "trait")


### Growth-recruiment ----
wd2 <- d %>% filter(mat == "10°C", trait == "wd", tradeoff == "sr") %>%
  ggplot(aes(x = growth_est, xmin = growth_lo, xmax = growth_hi,
             y = recr_est*10, ymin = recr_lo*10, ymax = recr_hi*10,
             color = wd)) +
  geom_pointrange() +
  geom_errorbarh(height = 0) +
  scale_y_log10(breaks = c(0.001, 0.001, 0.01, 0.1)) +
  set_theme +
  nice_colors +
  labs(x = "Growth (60 cm dbh)", y = "Recruitment (8 cm dbh)")

### Survival-recruitment ----
wd3 <- d %>% filter(mat == "10°C", trait == "wd", tradeoff == "sr") %>%
  ggplot(aes(x = surv_est, xmin = surv_lo, xmax = surv_hi,
             y = recr_est*10, ymin = recr_lo*10, ymax = recr_hi*10,
             color = wd)) +
  geom_pointrange() +
  geom_errorbarh(height = 0) +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1)) +
  set_theme +
  nice_colors +
  labs(x = "Survival (60 cm dbh)", y = "Recruitment (8 cm dbh)")

## SLA ----

### Growth-survival ----
sla1 <- d %>% filter(mat == "10°C", trait == "sla", tradeoff == "gs") %>%
  ggplot(aes(x = growth_est, xmin = growth_lo, xmax = growth_hi,
             y = surv_est, ymin = surv_lo, ymax = surv_hi,
             color = sla)) +
  geom_pointrange() +
  geom_errorbarh(height = 0) +
  set_theme +
  nice_colors +
  labs(x = "Growth (5 cm dbh)", y = "Survival (5 cm dbh)", color = "trait")


### Growth-recruiment ----
sla2 <- d %>% filter(mat == "10°C", trait == "sla", tradeoff == "sr") %>%
  ggplot(aes(x = growth_est, xmin = growth_lo, xmax = growth_hi,
             y = recr_est*10, ymin = recr_lo*10, ymax = recr_hi*10,
             color = sla)) +
  geom_pointrange() +
  geom_errorbarh(height = 0) +
  scale_y_log10(breaks = c(0.001, 0.001, 0.01, 0.1), limits = c(0.01, 0.25)) +
  set_theme +
  nice_colors +
  labs(x = "Growth (60 cm dbh)", y = "Recruitment (8 cm dbh)")

### Survival-recruitment ----
sla3 <- d %>% filter(mat == "10°C", trait == "sla", tradeoff == "sr") %>%
  ggplot(aes(x = surv_est, xmin = surv_lo, xmax = surv_hi,
             y = recr_est*10, ymin = recr_lo*10, ymax = recr_hi*10,
             color = sla)) +
  geom_pointrange() +
  geom_errorbarh(height = 0) +
  scale_y_log10(breaks = c(0.001, 0.01, 0.1), limits = c(0.01, 0.25)) +
  set_theme +
  nice_colors +
  labs(x = "Survival (60 cm dbh)", y = "Recruitment (8 cm dbh)")


## Height ----
### Growth-survival ----
h1 <- d %>% filter(mat == "10°C", trait == "hmax", tradeoff == "gs") %>%
  ggplot(aes(x = growth_est, xmin = growth_lo, xmax = growth_hi,
             y = surv_est, ymin = surv_lo, ymax = surv_hi,
             color = hmax)) +
  geom_pointrange() +
  geom_errorbarh(height = 0) +
  set_theme +
  nice_colors +
  labs(x = "Growth (5 cm dbh)", y = "Survival (5 cm dbh)", color = "trait")


### Growth-recruiment ----
h2 <- d %>% filter(mat == "10°C", trait == "hmax", tradeoff == "sr") %>%
  ggplot(aes(x = growth_est, xmin = growth_lo, xmax = growth_hi,
             y = recr_est*10, ymin = recr_lo*10, ymax = recr_hi*10,
             color = hmax)) +
  geom_pointrange() +
  geom_errorbarh(height = 0) +
  scale_y_log10(breaks = c(0.001, 0.001, 0.01, 0.1)) +
  set_theme +
  nice_colors +
  labs(x = "Growth (60 cm dbh)", y = "Recruitment (8 cm dbh)")

### Survival-recruitment ----
h3 <- d %>% filter(mat == "10°C", trait == "hmax", tradeoff == "sr") %>%
  ggplot(aes(x = surv_est, xmin = surv_lo, xmax = surv_hi,
             y = recr_est*10, ymin = recr_lo*10, ymax = recr_hi*10,
             color = hmax)) +
  geom_pointrange() +
  geom_errorbarh(height = 0) +
  scale_y_log10(breaks = c(0.001, 0.001, 0.01, 0.1)) +
  set_theme +
  nice_colors +
  labs(x = "Survival (60 cm dbh)", y = "Recruitment (8 cm dbh)")

# Plots ----
r1 <- ggpubr::ggarrange(wd1, wd2, wd3, ncol = 3, common.legend = T, legend = "right", 
                        labels = c("a", "b", "c"))
r2 <- ggpubr::ggarrange(sla1, sla2, sla3, ncol = 3, common.legend = T, legend = "right", 
                        labels = c("d", "e", "f"))
r3 <- ggpubr::ggarrange(h1, h2, h3, ncol = 3, common.legend = T, legend = "right", 
                        labels = c("g", "h", "i"))

pdf("results/figures/fig3.pdf", 9.5, 8)
gridExtra::grid.arrange(rbind(ggplotGrob(r1), ggplotGrob(r2), ggplotGrob(r3)))
dev.off()


   