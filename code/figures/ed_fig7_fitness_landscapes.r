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
### Lambda ----
d <- read_csv("results/ipm/fitness_landscapes.csv")

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

lambda_lab <- expression("Fitness ("*lambda*")")

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

### WD - SLA ----
f1 <- d %>% 
  filter(hmax == 0)  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         sla = backtransform(sla, tf$trait$sla),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = wd, y = sla)) + 
  geom_raster(aes(fill = lambda)) + 
  geom_contour(aes(z = lambda), color = "white", size = 0.1) +
  labs(x = wd_lab, y = sla_lab, fill = lambda_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_sla$x), oob = squish_infinite) +
  scale_y_continuous(expand = c(0,0), limits = range(mask_wd_sla$y), oob = squish_infinite) +
  scale_fill_viridis_c(trans = "log", option = "magma", limits = c(0.9126, 1.055)) + 
  new_scale_fill() +
  geom_raster(data = mask_wd_sla, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(~mat) +
  set_theme + 
  theme(legend.position = "none")

### WD - height ----
f2 <- d %>% 
  filter(sla == 0)  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = wd, y = hmax)) + 
  geom_raster(aes(fill = lambda)) + 
  geom_contour(aes(z = lambda), color = "white", size = 0.1) +
  labs(x = wd_lab, y = hmax_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_hmax$x), oob = squish_infinite) +
  scale_y_log10(expand = c(0,0), limits = range(mask_wd_hmax$y), oob = squish_infinite) +
  scale_fill_viridis_c(trans = "log", option = "magma", limits = c(0.9126, 1.055)) + 
  new_scale_fill() +
  geom_raster(data = mask_wd_hmax, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(~mat) +
  set_theme + 
  theme(legend.position = "none")


### SLA - height ----

f3 <- d %>% 
  filter(wd == 0)  %>%
  mutate(sla = backtransform(sla, tf$trait$sla),
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = sla, y = hmax)) + 
  geom_raster(aes(fill = lambda)) + 
  geom_contour(aes(z = lambda), color = "white", size = 0.1) +
  labs(x = sla_lab, y = hmax_lab, fill = lambda_lab)  +
  scale_x_continuous(expand = c(0,0), limits = range(mask_sla_hmax$x), oob = squish_infinite) +
  scale_y_log10(expand = c(0,0), limits = range(mask_sla_hmax$y), oob = squish_infinite) +
  scale_fill_viridis_c(trans = "log", option = "magma", limits = c(0.9126, 1.055)) + 
  new_scale_fill() +
  geom_raster(data = mask_sla_hmax, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(~mat) +
  set_theme + 
  theme(legend.position = "none")

### Plot ----

ggarrange(f1, f2, f3, common.legend = T, legend = "right",
          labels = "auto", ncol = 1)

ggsave("results/figures/ed_fig7.pdf", width = 8, height = 8, units = "in")

