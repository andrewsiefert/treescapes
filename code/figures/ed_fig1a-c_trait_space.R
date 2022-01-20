library(tidyverse)
library(raster)
library(ggnewscale)


data <- read_csv("data/species_traits.csv")

traits <- data %>%
  distinct(species, wd, sla, hmax) %>%
  separate(species, c("genus", "species"), remove = F) %>%
  mutate(spp = paste(str_sub(genus, 1,3), str_sub(species, 1,3), sep = "."))

wd_lab <- expression(paste("Wood density (g cm"^{-3}*")"))
sla_lab <- expression(paste("Specific leaf area (mm"^{2}*" mg"^{-1}*")"))
hmax_lab <- "Maximum height (m)"

# WD and SLA
tr <- traits %>% dplyr::select(wd, sla) %>% distinct()
hull_pts <- chull(tr)
wd_sla_hull <- tr[c(hull_pts, hull_pts[1]),] %>% 
  as_tibble() %>%
  mutate(theta = acos(sum(wd*sla/1000)/10)) %>%
  arrange(theta)

# WD and max height
tr <- traits %>% dplyr::select(wd, hmax) %>% distinct()
hull_pts <- chull(tr)
wd_hmax_hull <- tr[c(hull_pts, hull_pts[1]),] %>% 
  as_tibble() %>%
  mutate(theta = acos(sum(wd*hmax/1000)/10)) %>%
  arrange(theta)

# SLA and max height
tr <- traits %>% dplyr::select(sla, hmax) %>% distinct()
hull_pts <- chull(tr)
sla_hmax_hull <- tr[c(hull_pts, hull_pts[1]),] %>% 
  as_tibble() %>%
  mutate(theta = acos(sum(sla*hmax/1000)/100)) %>%
  arrange(theta)

p1 <- ggplot(traits, aes(x = wd, y = sla)) + 
  #geom_point(size = 0.7) +
  ggrepel::geom_text_repel(aes(label = spp), size = 2.5, fontface = "italic", box.padding = 0.1, segment.color = NA) +
  geom_polygon(data = wd_sla_hull, alpha = 0.05) + 
  scale_x_log10() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = wd_lab, y = sla_lab)
  
p2 <- ggplot(traits, aes(x = wd, y = hmax)) + 
  #geom_point(size = 0.7) +
  ggrepel::geom_text_repel(aes(label = spp), size = 2.5, fontface = "italic", box.padding = 0.1, segment.color = NA) +
  geom_polygon(data = wd_hmax_hull, alpha = 0.05) + 
  scale_x_log10() +
  scale_y_log10() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = wd_lab, y = hmax_lab)

p3 <- ggplot(traits, aes(x = sla, y = hmax)) + 
  #geom_point(size = 0.7) +
  ggrepel::geom_text_repel(aes(label = spp), size = 2.5, fontface = "italic", box.padding = 0.1, segment.color = NA) +
  geom_polygon(data = sla_hmax_hull, alpha = 0.05) + 
  scale_y_log10() +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  labs(x = sla_lab, y = hmax_lab)

ggpubr::ggarrange(p1, p2, p3, ncol = 1)
ggsave("results/figures/ed_fig1a-c.jpg", height = 16, width = 6, dpi = 300)
