library(tidyverse)
library(raster)
library(sp)
library(scales)
library(ggnewscale)
library(ggpubr)

source("code/transformer.r")
tf <- readRDS("data/transformers.rds")


# Setup ------------------------------------------------------------
## Read in data ----
### Performance landscapes----
d <- read_csv("results/ipm/performance_landscapes.csv") %>%
  mutate(survival = ifelse(size <= 10, surv_sapling, surv_canopy), 
         size_lab = paste(size, "cm dbh") %>% fct_reorder(size))

### Traits----
traits <- read_csv("data/growth_data_train.csv") %>%
  distinct(wood_density, sla, hmax) %>%
  rename(wd = wood_density) %>%
  mutate(wd = transform(wd, tf$trait$wood_density_log, log = T),
         sla = transform(sla, tf$trait$sla), 
         hmax = transform(hmax, tf$trait$hmax_log, log = T))

## Plot setup----
### Labels----
wd_lab <- expression(paste("Wood density (g cm"^{-3}*")"))
sla_lab <- expression(paste("Specific leaf area (mm"^{2}*" mg"^{-1}*")"))
hmax_lab <- "Maximum height (m)"

surv_lab <- expression(paste("Survival\nprobability (yr"^{-1}*")"))
growth_lab <- expression(paste("Growth (mm yr"^{-1}*")"))
recruit_lab <- expression(paste("Recruitment\n(recruits ind"^{-1}*" yr"^{-1}*")"))

title_size <- 9
label_size <- 7

### Theme----

set_theme <- theme_bw() +
  theme(axis.title.x = element_text(size = title_size),
        axis.text.x = element_text(size = label_size),
        axis.title.y = element_text(size = title_size),
        axis.text.y = element_text(size = label_size),
        strip.text.x = element_text(size = title_size, margin = margin(t=2, b=2, unit="pt")), 
        strip.text.y = element_text(size = title_size, margin = margin(l=2, r=2, unit="pt")), 
        legend.justification = "left",
        legend.title = element_text(size = title_size),
        legend.text = element_text(size = label_size), 
        legend.key.height = unit(0.4, "cm"))


### Colors----
t_col <- function(color, percent = 50) {
  rgb.val <- col2rgb(color)
  t.col <- rgb(rgb.val[1], rgb.val[2], rgb.val[3], max = 255, 
               alpha = (100 - percent) * 255 / 100)
  invisible(t.col)
}

clear <- t_col("white", 100)
mid <- t_col("white", 25)

## Create masks -----

# WD and SLA
tr <- traits %>% dplyr::select(wd, sla) %>% distinct()
hull_pts <- chull(tr)
hull <- tr[c(hull_pts, hull_pts[1]),] %>% 
  as_tibble()

p <- Polygon(hull)
ps <- Polygons(list(p),1)
sps <- SpatialPolygons(list(ps))


mask_wd_sla <- expand.grid(x = modelr::seq_range(d$wd, 200), 
                           y = modelr::seq_range(d$sla, 200)) %>%
  mutate(value = 1) %>%
  rasterFromXYZ() %>%
  crop(extent(sps)) %>%
  mask(sps, inverse = T) %>%
  as.data.frame(xy = T) %>%
  as_tibble() %>%
  mutate(x = backtransform(x, tf$trait$wood_density_log, log = T),
         y = backtransform(y, tf$trait$sla))

# WD and max height
tr <- traits %>% dplyr::select(wd, hmax) %>% distinct()
hull_pts <- chull(tr)
hull <- tr[c(hull_pts, hull_pts[1]),] %>% 
  as_tibble()

p <- Polygon(hull)
ps <- Polygons(list(p),1)
sps <- SpatialPolygons(list(ps))


mask_wd_hmax <- expand.grid(x = modelr::seq_range(d$wd, 200), 
                            y = modelr::seq_range(d$hmax, 200)) %>%
  mutate(value = 1) %>%
  rasterFromXYZ() %>%
  crop(extent(sps)) %>%
  mask(sps, inverse = T) %>%
  as.data.frame(xy = T) %>%
  as_tibble() %>%
  mutate(x = backtransform(x, tf$trait$wood_density_log, log = T),
         y = backtransform(y, tf$trait$hmax_log, log = T))

tr <- traits %>% dplyr::select(sla, hmax) %>% distinct()
hull_pts <- chull(tr)
hull <- tr[c(hull_pts, hull_pts[1]),] %>% 
  as_tibble()

p <- Polygon(hull)
ps <- Polygons(list(p),1)
sps <- SpatialPolygons(list(ps))


mask_sla_hmax <- expand.grid(x = modelr::seq_range(d$sla, 200), 
                             y = modelr::seq_range(d$hmax, 200)) %>%
  mutate(value = 1) %>%
  rasterFromXYZ() %>%
  crop(extent(sps)) %>%
  mask(sps, inverse = T) %>%
  as.data.frame(xy = T) %>%
  as_tibble() %>%
  mutate(x = backtransform(x, tf$trait$sla),
         y = backtransform(y, tf$trait$hmax_log, log = T))



# Make plots ------------------------------------------------------------

## Survival ----

### Fill limits

surv_lim <- d %>% filter(size %in% c(10, 30)) %>% pull(survival) %>% range()
surv_breaks <- seq(0, 1, 0.005)

### WD - SLA ----
s1 <- d %>% 
  filter(hmax == 0, size %in% c(5, 20, 60))  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         sla = backtransform(sla, tf$trait$sla),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = wd, y = sla)) + 
  geom_raster(aes(fill = survival)) + 
  geom_contour(aes(z = survival), color = "white", size = 0.1) +
  labs(x = wd_lab, y = sla_lab, fill = surv_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_sla$x), oob = squish_infinite) +
  scale_y_continuous(expand = c(0,0), limits = range(mask_wd_sla$y), oob = squish_infinite) +
  scale_fill_viridis_c(na.value = t_col("white", 100), limits = surv_lim) +
  new_scale_fill() +
  geom_raster(data = mask_wd_sla, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(size_lab~mat) +
  set_theme + 
  theme(legend.position = "none")

### WD - height ----
s2 <- d %>% 
  filter(sla == 0, size %in% c(5, 20, 60))  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = wd, y = hmax)) + 
  geom_raster(aes(fill = survival)) + 
  geom_contour(aes(z = survival), color = "white", size = 0.1) +
  labs(x = wd_lab, y = hmax_lab, fill = surv_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_hmax$x), oob = squish_infinite) +
  scale_y_log10(expand = c(0,0), limits = range(mask_wd_hmax$y), oob = squish_infinite) +
  scale_fill_viridis_c(na.value = t_col("white", 100), limits = surv_lim) +
  new_scale_fill() +
  geom_raster(data = mask_wd_hmax, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(size_lab~mat) +
  set_theme + 
  theme(legend.position = "none")


### SLA - height ----

s3 <- d %>% 
  filter(wd == 0, size %in% c(5, 20, 60))  %>%
  mutate(sla = backtransform(sla, tf$trait$sla),
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = sla, y = hmax)) + 
  geom_raster(aes(fill = survival)) + 
  geom_contour(aes(z = survival), color = "white", size = 0.1) +
  labs(x = sla_lab, y = hmax_lab, fill = surv_lab)  +
  scale_x_continuous(expand = c(0,0), limits = range(mask_sla_hmax$x), oob = squish_infinite) +
  scale_y_log10(expand = c(0,0), limits = range(mask_sla_hmax$y), oob = squish_infinite) +
  scale_fill_viridis_c(na.value = t_col("white", 100), limits = surv_lim) +
  new_scale_fill() +
  geom_raster(data = mask_sla_hmax, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(size_lab~mat) +
  set_theme

### Plot ----

ggarrange(s1, s2, s3, common.legend = T, legend = "bottom",
          labels = "auto")

ggsave("results/figures/ed_fig4.pdf", width = 8, height = 8, units = "in")

## Growth ----

### WD-SLA ----
g1 <- d %>%
  filter(hmax == 0, size %in% c(5, 20, 60))  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         sla = backtransform(sla, tf$trait$sla),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = wd, y = sla)) + 
  geom_raster(aes(fill = growth)) + 
  geom_contour(aes(z = growth), color = "white", size = 0.1) +
  labs(x = wd_lab, y = sla_lab, fill = growth_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_sla$x), oob = squish_infinite) +
  scale_y_continuous(expand = c(0,0), limits = range(mask_wd_sla$y), oob = squish_infinite) +
  scale_fill_viridis_c(trans = "log", na.value = t_col("white", 100), breaks = c(1,2,3)) +
  new_scale_fill() +
  geom_raster(data = mask_wd_sla, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(size_lab~mat) +
  set_theme

### WD - height ----
g2 <- d %>%
  filter(sla == 0, size %in% c(5, 20, 60))  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = wd, y = hmax)) + 
  geom_raster(aes(fill = growth)) + 
  geom_contour(aes(z = growth), color = "white", size = 0.1) +
  labs(x = wd_lab, y = hmax_lab, fill = growth_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_hmax$x), oob = squish_infinite) +
  scale_y_log10(expand = c(0,0), limits = range(mask_wd_hmax$y), oob = squish_infinite) +
  scale_fill_viridis_c(trans = "log", na.value = t_col("white", 100)) +
  new_scale_fill() +
  geom_raster(data = mask_wd_hmax, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(size_lab~mat) +
  set_theme

### SLA - height ----

g3 <- d %>%
  filter(wd == 0, size %in% c(5, 20, 60))  %>%
  mutate(sla = backtransform(sla, tf$trait$sla),
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = sla, y = hmax)) + 
  geom_raster(aes(fill = growth)) + 
  geom_contour(aes(z = growth), color = "white", size = 0.1) +
  labs(x = sla_lab, y = hmax_lab, fill = growth_lab)  +
  scale_x_continuous(expand = c(0,0), limits = range(mask_sla_hmax$x), oob = squish_infinite) +
  scale_y_log10(expand = c(0,0), limits = range(mask_sla_hmax$y), oob = squish_infinite) +
  scale_fill_viridis_c(trans = "log", na.value = t_col("white", 100)) +
  new_scale_fill() +
  geom_raster(data = mask_sla_hmax, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(size_lab~mat) +
  set_theme

### Plot ----

ggarrange(g1, g2, g3, common.legend = T, nrow = 2, ncol = 2, legend = "bottom",
          labels = "auto")

ggsave("results/figures/ed_fig5.pdf", width = 8, height = 8, units = "in")


## Recruitment ----

### WD - SLA ----

r1 <- d %>% 
  filter(hmax == 0, size %in% c(5, 20, 60))  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         sla = backtransform(sla, tf$trait$sla),
         mat = paste0(mat, "°C") %>% fct_reorder(mat),
         recruitment = pmax(recruitment, 0.0005))%>%
  ggplot(aes(x = wd, y = sla)) + 
  geom_raster(aes(fill = recruitment)) + 
  geom_contour(aes(z = recruitment), color = "white", breaks = seq(0, 1, 0.005), size = 0.1) +
  labs(x = wd_lab, y = sla_lab, fill = recruit_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_sla$x), oob = squish) +
  scale_y_continuous(expand = c(0,0), limits = range(mask_wd_sla$y), oob = squish_infinite) +
  scale_fill_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1)) + 
  new_scale_fill() +
  geom_raster(data = mask_wd_sla, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(size_lab~mat) +
  set_theme


### WD - height----------------- 

r2 <- d %>% 
  filter(sla == 0, size %in% c(5, 20, 60))  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T),
         mat = paste0(mat, "°C") %>% fct_reorder(mat),
         recruitment = pmax(recruitment, 0.0005))%>%
  ggplot(aes(x = wd, y = hmax)) + 
  geom_raster(aes(fill = recruitment)) + 
  geom_contour(aes(z = recruitment), color = "white", breaks = seq(0, 1, 0.005), size = 0.1) +
  labs(x = wd_lab, y = hmax_lab, fill = recruit_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_hmax$x), oob = squish) +
  scale_y_log10(expand = c(0,0), limits = range(mask_wd_hmax$y), oob = squish_infinite) +
  scale_fill_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1)) + 
  new_scale_fill() +
  geom_raster(data = mask_wd_hmax, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(size_lab~mat) +
  set_theme


### SLA - height ----
r3 <- d %>% 
  filter(wd == 0, size %in% c(5, 20, 60))  %>%
  mutate(sla = backtransform(sla, tf$trait$sla),
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T),
         mat = paste0(mat, "°C") %>% fct_reorder(mat),
         recruitment = pmax(recruitment, 0.0005))%>%
  ggplot(aes(x = sla, y = hmax)) + 
  geom_raster(aes(fill = recruitment)) + 
  geom_contour(aes(z = recruitment), color = "white", breaks = seq(0, 1, 0.005), size = 0.1) +
  labs(x = sla_lab, y = hmax_lab, fill = recruit_lab)  +
  scale_x_continuous(expand = c(0,0), limits = range(mask_sla_hmax$x)) +
  scale_y_log10(expand = c(0,0), limits = range(mask_sla_hmax$y), oob = squish_infinite) +
  scale_fill_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1)) + 
  new_scale_fill() +
  geom_raster(data = mask_sla_hmax, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(size_lab~mat) +
  set_theme


### Plot ----

ggarrange(r1, r2, r3, common.legend = T, nrow = 2, ncol = 2, legend = "bottom",
          labels = "auto")

ggsave("results/figures/ed_fig6.pdf", width = 8, height = 8, units = "in")
