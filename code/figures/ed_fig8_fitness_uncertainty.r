library(tidyverse)
library(gridExtra)
library(raster)
library(sp)


source("code/transformer.r")
tf <- readRDS("data/transformers.rds")


# Trait "maps" ----

## Create data frame ----

tmap <- expand.grid(t1 = c(-2, -sqrt(2), 0, sqrt(2), 2),
                    t2 = c(-2, -sqrt(2), 0, sqrt(2), 2),
                    traits = c("wd-sla", "wd-hmax", "sla-hmax")) %>%
  as_tibble() %>%
  filter(sqrt(t1^2 + t2^2) %in% c(0, 2)) %>%
  mutate(wd = ifelse(str_detect(traits, "wd"), t1, 0),
         sla = ifelse(traits == "wd-sla", t2, ifelse(traits == "sla-hmax", t1, 0)),
         hmax = ifelse(str_detect(traits, "hmax"), t2, 0),
         label = case_when(
           wd == -2 ~ 2,
           wd == 2 ~ 3,
           sla == -2 ~ 4,
           sla == 2 ~ 5,
           hmax == -2 ~ 6, 
           hmax == 2 ~ 7,
           wd == sqrt(2) & sla == sqrt(2) ~ 8,
           wd == sqrt(2) & sla == -sqrt(2) ~ 9,
           wd == -sqrt(2) & sla == -sqrt(2) ~ 10,
           wd == -sqrt(2) & sla == sqrt(2) ~ 11,
           wd == sqrt(2) & hmax == sqrt(2) ~ 12,
           wd == sqrt(2) & hmax == -sqrt(2) ~ 13,
           wd == -sqrt(2) & hmax == -sqrt(2) ~ 14,
           wd == -sqrt(2) & hmax == sqrt(2) ~ 15,
           sla == sqrt(2) & hmax == sqrt(2) ~ 16,
           sla == sqrt(2) & hmax == -sqrt(2) ~ 17,
           sla == -sqrt(2) & hmax == -sqrt(2) ~ 18,
           sla == -sqrt(2) & hmax == sqrt(2) ~ 19,
           TRUE ~ 1
         ),
         t1 = ifelse(traits == "sla-hmax",
                     backtransform(t1, tf$trait$sla),
                     backtransform(t1, tf$trait$wood_density_log, log = T)),
         t2 = ifelse(traits == "wd-sla",
                     backtransform(t2, tf$trait$sla),
                     backtransform(t2, tf$trait$hmax_log, log = T)))

wd_lab <- expression(paste("Wood density (g cm"^{-3}*")"))
sla_lab <- expression(paste("Specific leaf area (mm"^{2}*" mg"^{-1}*")"))
hmax_lab <- "Maximum height (m)"



## Create masks -----

d <- read_csv("results/ipm/fitness_landscapes.csv")

traits <- read_csv("data/growth_data_train.csv") %>%
  distinct(wood_density, sla, hmax) %>%
  rename(wd = wood_density) %>%
  mutate(wd = transform(wd, tf$trait$wood_density_log, log = T),
         sla = transform(sla, tf$trait$sla), 
         hmax = transform(hmax, tf$trait$hmax_log, log = T))

# WD and SLA
tr <- traits %>% dplyr::select(wd, sla) %>% distinct()
hull_pts <- chull(tr)
wd_sla_hull <- tr[c(hull_pts, hull_pts[1]),] %>% 
  as_tibble() %>%
  mutate(x = backtransform(wd, tf$trait$wood_density_log, log = T),
         y = backtransform(sla, tf$trait$sla),
         theta = acos(sum(wd*sla)/10)) %>%
  arrange(theta)

# WD and max height
tr <- traits %>% dplyr::select(wd, hmax) %>% distinct()
hull_pts <- chull(tr)
wd_hmax_hull <- tr[c(hull_pts, hull_pts[1]),] %>% 
  as_tibble() %>%
  mutate(x = backtransform(wd, tf$trait$wood_density_log, log = T),
         y = backtransform(hmax, tf$trait$hmax_log, log = T),
         theta = acos(sum(wd*hmax)/10)) %>%
  arrange(theta)

# SLA and max height
tr <- traits %>% dplyr::select(sla, hmax) %>% distinct()
hull_pts <- chull(tr)
sla_hmax_hull <- tr[c(hull_pts, hull_pts[1]),] %>% 
  as_tibble() %>%
  mutate(x = backtransform(sla, tf$trait$sla),
         y = backtransform(hmax, tf$trait$hmax_log, log = T),
         theta = acos(sum(sla*hmax)/100)) %>%
  arrange(theta)

## Make plots
m1 <- tmap %>% 
  filter(traits == "wd-sla") %>%
  ggplot() + 
  geom_polygon(data = wd_sla_hull, aes(x = x, y = y), alpha = 0.1) +
  geom_label(aes(x = t1, y = t2, label = label)) + 
  scale_x_log10() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = wd_lab, y = sla_lab)

m2 <- tmap %>% 
  filter(traits == "wd-hmax") %>%
  ggplot() + 
  geom_polygon(data = wd_hmax_hull, aes(x = x, y = y), alpha = 0.1) +
  geom_label(aes(x = t1, y = t2, label = label)) + 
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = wd_lab, y = hmax_lab)

m3 <- tmap %>% 
  filter(traits == "sla-hmax") %>%
  ggplot() + 
  geom_polygon(data = sla_hmax_hull, aes(x = x, y = y), alpha = 0.1) +
  geom_label(aes(x = t1, y = t2, label = label)) + 
  scale_y_log10() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = sla_lab, y = hmax_lab)

maps <- ggpubr::ggarrange(m1, m2, m3, nrow = 1)

# Plot fitness ----

## Prepare data ----
vals <- c(-2, -sqrt(2), 0, sqrt(2), 2)
grid <- expand.grid(wd = vals, 
                    sla = vals, 
                    hmax = vals, 
                    mat = transform(c(5, 10, 15), tf$env$mat)) %>%
  as_tibble() %>%
  filter(sqrt(wd^2 + sla^2 + hmax^2) %in% c(0, 2))

lambda <- readRDS("results/ipm/fitness_uncertainty.rds")

df <- grid %>% 
  mutate_at(vars(wd, sla, hmax), round, 1) %>%
  mutate(traits = paste(wd, sla, hmax, sep = "\n"),
         mean = rowMeans(lambda), 
         median = apply(lambda, 1, median),
         q2.5 = apply(lambda, 1, quantile, 0.025),
         q5 = apply(lambda, 1, quantile, 0.05),
         q95 = apply(lambda, 1, quantile, 0.95),
         q97.5 = apply(lambda, 1, quantile, 0.975)) %>%
  arrange(wd, sla, hmax) %>%
  mutate(traits = fct_inorder(traits),
         mat = backtransform(mat, tf$env$mat) %>% paste0("°C") %>% fct_reorder(mat)) %>%
  left_join(tmap %>% dplyr::select(wd:label) %>% mutate_at(vars(wd:hmax), round, 1)) %>%
  arrange(label, mat)

lambda_lab <- expression("Fitness ("*lambda*")")



## Horizontal ----

p1 <- ggplot(df, aes(x = factor(label), y = median, ymin = q5, ymax = q95, color = mat)) +
  geom_hline(aes(yintercept = 1), color = "gray", linetype = "dashed") +
  geom_pointrange(position = position_dodge(width = 0.6), size = 0.4) + 
  scale_color_viridis_d() +
  labs(y = lambda_lab, color = "Mean annual\ntemperature") +
  theme_bw() +
  theme(axis.text.x = element_blank(), 
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

p2 <- df %>% 
  distinct(label, wd, sla, hmax) %>% 
  rename("WD" = wd, "SLA" = sla, "Height" = hmax) %>%
  gather("trait", "value", WD:Height) %>%
  arrange(label) %>%
  ggplot(aes(x = factor(label), y = trait, fill = value)) +
  geom_tile() +
  scale_fill_distiller(palette = "RdBu", direction = 1) +
  labs(fill = "Trait value") +
  theme(axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_blank(),
        axis.title.y = element_blank()) +
  scale_y_discrete(expand = c(0,0))

ggpubr::ggarrange(maps, p1, p2, nrow = 3)
ggsave("results/figures/ed_fig8.pdf", width = 10, height = 9, units = "in")

