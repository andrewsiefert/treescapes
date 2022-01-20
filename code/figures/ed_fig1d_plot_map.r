library(sf)
library(tidyverse)


data <- read_csv("data/plot_data.csv") %>%
  distinct(lon, lat, mat)


world <- map_data("world")
states <- map_data("state")

lakes <- rnaturalearth::ne_download(scale = 50, 
                                    type = 'lakes', 
                                    category = 'physical') %>% 
  sf::st_as_sf(lakes110, crs = 4269) %>%
  filter(name_tr %in% c("Erie", "Ontario", "Superior", "Michigan", "Huron", "Champlain"))  

ggplot(world) + 
  geom_polygon(aes(x = long, y = lat, group = group), fill = "white", color = "darkgray") +
  geom_polygon(data = states, aes(x = long, y = lat,group = group), fill = "white", color = "darkgray") +
  geom_sf(data = lakes,
          mapping = aes(geometry = geometry),
          color = "darkgray",
          fill = "lightgray") +
  geom_point(data = data, aes(x = lon, y = lat, color = mat), size = 0.4) +
  coord_sf(xlim = c(-91, -67), ylim = c(30, 47)) +
  theme(panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "lightgray"),
        panel.border = element_rect(colour = "black", fill=NA)) +
  scale_color_viridis_c() +
  labs(x = "", y = "", color = expression(paste("Mean annual\ntemperature", "(", degree*C, ")", sep = ""))) +
  theme(legend.position = "bottom")

ggsave("results/figures/ed_fig1d.pdf", height = 8, width = 8)
       