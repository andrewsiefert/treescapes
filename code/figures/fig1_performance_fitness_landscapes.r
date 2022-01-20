library(tidyverse)
library(raster)
library(sp)
library(scales)
library(ggnewscale)
library(grid)
library(egg)

source("code/transformer.r")
tf <- readRDS("data/transformers.rds")


# Setup ------------------------------------------------------------
## Read in data ----
### Performance landscapes----
d <- read_csv("results/ipm/performance_landscapes.csv") %>%
  mutate(survival = ifelse(size == 10, surv_sapling, surv_canopy), 
         size_lab = paste(size, "cm dbh") %>% fct_reorder(size))

### Fitness landscapes----
f <- read_csv("results/ipm/fitness_landscapes.csv") %>%
  rename(fitness = lambda)

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
lambda_lab <- expression("Fitness ("*lambda*")")

title_size <- 9
label_size <- 7

### Themes----

theme_top <- theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   axis.text.y = element_text(size = label_size),
                   strip.text.x = element_text(size = title_size), 
                   legend.justification = "left",
                   legend.title = element_text(size = title_size),
                   legend.text = element_text(size = label_size), 
                   legend.key.height = unit(0.4, "cm"))

theme_middle <- theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.text.y = element_text(size = label_size),
                      strip.background = element_blank(),
                      strip.text.x = element_blank(), 
                      legend.justification = "left",
                      legend.title = element_text(size = title_size),
                      legend.text = element_text(size = label_size), 
                      legend.key.height = unit(0.4, "cm"))

theme_ylab <- theme(axis.title.x = element_blank(),
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    axis.title.y = element_text(size = title_size),
                    axis.text.y = element_text(size = label_size),
                    strip.background = element_blank(),
                    strip.text.x = element_blank(), 
                    legend.justification = "left",
                    legend.title = element_text(size = title_size),
                    legend.text = element_text(size = label_size), 
                    legend.key.height = unit(0.4, "cm"))

theme_bottom <- theme(axis.title.x = element_text(size = title_size),
                      axis.text.x = element_text(size = label_size),
                      axis.text.y = element_text(size = label_size),
                      strip.background = element_blank(),
                      strip.text.x = element_blank(), 
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
mid <- t_col("white", 30)


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

## WD - SLA ----

### Survival ----
p1 <- d %>% 
  filter(size == 20, hmax == 0, mat != 10)  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         sla = backtransform(sla, tf$trait$sla),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = wd, y = sla)) + 
  geom_raster(aes(fill = survival)) + 
  geom_contour(aes(z = survival), color = "white", breaks = seq(0, 5, 0.005), size = 0.1) +
  labs(x = wd_lab, y = "", fill = surv_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_sla$x), oob = squish_infinite) +
  scale_y_continuous(expand = c(0,0), limits = range(mask_wd_sla$y), oob = squish_infinite) +
  scale_fill_viridis_c(na.value = t_col("white", 100), limits = c(0.89, 1)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "none") +
  new_scale_fill() +
  geom_raster(data = mask_wd_sla, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(~mat) +
  theme_top

### Growth ----
p2 <- d %>% 
  filter(size == 20, hmax == 0, mat != 10)  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         sla = backtransform(sla, tf$trait$sla),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = wd, y = sla)) + 
  geom_raster(aes(fill = growth)) + 
  geom_contour(aes(z = growth), color = "white", breaks = seq(0, 5, 0.2), size = 0.1) +
  labs(x = wd_lab, y = "", fill = growth_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_sla$x), oob = squish_infinite) +
  scale_y_continuous(expand = c(0,0), limits = range(mask_wd_sla$y), oob = squish_infinite) +
  scale_fill_viridis_c(na.value = t_col("white", 100), limits = c(1.4, 4.4)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "none") +
  new_scale_fill() +
  geom_raster(data = mask_wd_sla, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(~mat) +
  theme_middle

### Recruitment ----
p3 <- d %>% 
  filter(size == 20, hmax == 0, mat != 10)  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         sla = backtransform(sla, tf$trait$sla),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = wd, y = sla)) + 
  geom_raster(aes(fill = recruitment)) + 
  geom_contour(aes(z = recruitment), color = "white", breaks = seq(0, 1, 0.005), size = 0.1) +
  labs(x = wd_lab, y = "", fill = recruit_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_sla$x), oob = squish) +
  scale_y_continuous(expand = c(0,0), limits = range(mask_wd_sla$y), oob = squish_infinite) +
  scale_fill_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1), limits = c(0.000159, 0.401)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "none") +
  new_scale_fill() +
  geom_raster(data = mask_wd_sla, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(~mat) +
  theme_middle

### Fitness ----
p4 <- f %>% 
  filter(hmax == 0, mat != 10)  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         sla = backtransform(sla, tf$trait$sla),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = wd, y = sla)) + 
  geom_raster(aes(fill = fitness)) + 
  geom_contour(aes(z = fitness), color = "white", breaks = seq(0, 2, 0.005), size = 0.1) +
  labs(x = wd_lab, y = "", fill = lambda_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_sla$x), oob = squish) +
  scale_y_continuous(expand = c(0,0), limits = range(mask_wd_sla$y), oob = squish_infinite) +
  scale_fill_viridis_c(trans = "log", option = "magma", limits = c(0.9126, 1.055)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.position = "none") +
  new_scale_fill() +
  geom_raster(data = mask_wd_sla, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(~mat) +
  theme_bottom + 
  labs(x = "")


## WD - height----------------- 

### Survival ----
p5 <- d %>% 
  filter(size == 20, sla == 0, mat != 10)  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = wd, y = hmax)) + 
  geom_raster(aes(fill = survival)) + 
  geom_contour(aes(z = survival), color = "white", breaks = seq(0, 5, 0.005), size = 0.1) +
  labs(x = wd_lab, y = "", fill = surv_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_hmax$x), oob = squish_infinite) +
  scale_y_log10(expand = c(0,0), limits = range(mask_wd_hmax$y), oob = squish_infinite) +
  scale_fill_viridis_c(na.value = t_col("white", 100), limits = c(0.89, 1)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  new_scale_fill() +
  geom_raster(data = mask_wd_hmax, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(~mat) +
  theme_top

### Growth ----
p6 <- d %>% 
  filter(size == 20, sla == 0, mat != 10)  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = wd, y = hmax)) + 
  geom_raster(aes(fill = growth)) + 
  geom_contour(aes(z = growth), color = "white", breaks = seq(0, 5, 0.2), size = 0.1) +
  labs(x = wd_lab, y = "", fill = growth_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_hmax$x), oob = squish_infinite) +
  scale_y_log10(expand = c(0,0), limits = range(mask_wd_hmax$y), oob = squish_infinite) +
  scale_fill_viridis_c(na.value = t_col("white", 100), limits = c(1.4, 4.3)) +
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  new_scale_fill() +
  geom_raster(data = mask_wd_hmax, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(~mat) +
  theme_middle

### Recruitment ----
p7 <- d %>% 
  filter(size == 20, sla == 0, mat != 10)  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = wd, y = hmax)) + 
  geom_raster(aes(fill = recruitment)) + 
  geom_contour(aes(z = recruitment), color = "white", breaks = seq(0, 1, 0.005), size = 0.1) +
  labs(x = wd_lab, y = "", fill = recruit_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_hmax$x), oob = squish) +
  scale_y_log10(expand = c(0,0), limits = range(mask_wd_hmax$y), oob = squish_infinite) +
  scale_fill_viridis_c(trans = "log", breaks = c(0.001, 0.01, 0.1), limits = c(0.000159, 0.401)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  new_scale_fill() +
  geom_raster(data = mask_wd_hmax, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(~mat) +
  theme_middle

### Fitness ----
p8 <- f %>% 
  filter(sla == 0, mat != 10)  %>%
  mutate(wd = backtransform(wd, tf$trait$wood_density_log, log = T),
         hmax = backtransform(hmax, tf$trait$hmax_log, log = T),
         mat = paste0(mat, "°C") %>% fct_reorder(mat)) %>%
  ggplot(aes(x = wd, y = hmax)) + 
  geom_raster(aes(fill = fitness)) + 
  geom_contour(aes(z = fitness), color = "white", breaks = seq(0, 2, 0.005), size = 0.1) +
  labs(x = wd_lab, y = "", fill = lambda_lab)  +
  scale_x_log10(expand = c(0,0), limits = range(mask_wd_hmax$x), oob = squish_infinite) +
  scale_y_log10(expand = c(0,0), limits = range(mask_wd_hmax$y), oob = squish_infinite) +
  scale_fill_viridis_c(trans = "log", option = "magma", limits = c(0.9126, 1.055)) + 
  theme_bw() + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  new_scale_fill() +
  geom_raster(data = mask_wd_hmax, aes(x = x, y = y, fill = value), show.legend = F) +
  scale_fill_gradient2(low = 1, mid = mid, high = 2, na.value = clear) + 
  facet_grid(~mat) +
  theme_bottom + 
  labs(x = "")

svg("results/figures/fig1.svg", 8.5, 6.5)
grid.arrange(cbind(rbind(ggplotGrob(p1), ggplotGrob(p2), ggplotGrob(p3), ggplotGrob(p4)),
                   rbind(ggplotGrob(p5), ggplotGrob(p6), ggplotGrob(p7), ggplotGrob(p8))))
grid.text(wd_lab, x = 0.45, y = 0.03, just = "top", gp = gpar(fontsize=title_size))
grid.text(sla_lab, x = 0.005, y = unit(0.5, "npc"), rot = 90, just = "top", 
          gp = gpar(fontsize=title_size))
grid.text(hmax_lab, x = 0.43, y = unit(0.5, "npc"), rot = 90, just = "top", 
          gp = gpar(fontsize=title_size))
dev.off()

