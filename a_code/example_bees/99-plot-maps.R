
# Load packages -----------------------------------------------------------
source("a_code/functions/ipak.R")

ipak(c("ggplot2","sf","ggthemes","rnaturalearth","ggspatial"))

# Source code -------------------------------------------------------------

source("a_code/example_bees/06-correlations.R")

# Latin America map -------------------------------------------------------

mundo <- ne_countries(scale = "medium", returnclass = "sf")
america <- subset(mundo, region_wb == "Latin America & Caribbean")

america_map <- ggplot(data = america) +
  geom_sf( fill = "gray100", color = "gray60") +  
  geom_rect(xmin = -58, xmax = -34, ymin = -34, ymax = -2.5, 
            fill = NA, colour = "#3B528B") +
  annotate(geom = "text", x = -48, y = 17, label = "Atlantic\nOcean", 
           fontface = "italic", color = "grey22", angle = 0, family = "serif", size = 6) +
  annotate(geom = "text", x = -95, y = -20, label = "Pacific Ocean", 
           fontface = "italic", color = "grey22", family = "serif", size = 6) + 
  coord_sf(datum = NA) + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.title.y = element_blank(), panel.background = element_rect(fill = "white"),
        panel.border = element_rect(fill = NA))

# Brazil's cutout ---------------------------------------------------------

cutout_brasil <- ggplot(data = america) +
  geom_sf(fill =  "gray100", col = "black") +
  coord_sf(xlim = c(-60, -34), ylim = c(-35, 5), expand = T, 
           crs = st_crs(4326), ) 

cutout_brasil2 <- cutout_brasil + annotation_scale(location = "br", width_hint = 0.5, 
                                          text_family = "serif", text_cex = 1.2) +
  annotation_north_arrow(location = "tr", which_north = "true", 
                         pad_x = unit(0.2, "in"), pad_y = unit(0.3, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotate(geom = "text", x = -45, y = -33, label = "Atlantic Ocean", 
           fontface = "italic", color = "grey22", angle = 0, family = "serif", size = 8) +
  labs(x = "Longitude", y = "Latitude") +
  theme(text = element_text(family = "serif", size = 14, colour = "black"))

# Species richness map ----------------------------------------------------

richness <- combined_df |> 
  group_by(cell) |> 
  summarise(richness = n()) |>
  left_join(combined_df, by = "cell") |>
  select(longitude, latitude, richness) |>
  raster::rasterFromXYZ(res = res(environ_layers_cropped), 
                        crs = raster::crs(environ_layers_cropped)) 

ma_shapefile <- st_as_sf(ma_shapefile)
richness <- as.data.frame(richness, xy =TRUE)

brazil_ma <- cutout_brasil2 +  
  geom_sf(data = ma_shapefile, fill = "#C4E6C3",color = "transparent") +
  geom_tile(data = richness, aes(x = x, y = y, fill = richness)) +
  scale_fill_viridis_c(option = "C", direction = -1, na.value = NA) +
  geom_sf(data = america, fill =  NA, col = "black") +
  coord_sf(xlim = c(-60, -34), ylim = c(-35, 5), expand = T, 
           crs = st_crs(4326), ) +
  labs(fill = "Richness") +
  theme_classic() +
  #Wlabs (fill = "EE", title = "(A)") +
  theme(panel.background = element_rect(fill = "aliceblue"),
        legend.background = element_rect(fill = "transparent", color = NA),
        text = element_text(family = "serif", size = 18),
        axis.ticks = element_blank(),
        plot.subtitle = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA, size = 1),
        axis.title = element_blank(),
        axis.line = element_blank() ,
        legend.position = c(0.9, 0.2)) +
  guides(fill = guide_colorbar(direction = "vertical", 
                               title.position = "top",  
                               title.hjust = 0.05)) 

# Combining maps ----------------------------------------------------------

final_map = brazil_ma + annotation_custom(
  grob = ggplotGrob(america_map),
  xmin = -45,
  xmax = -61,
  ymin = -10,
  ymax = 6)

ggsave("c_outputs/figures/bees/richness.png", 
       plot = final_map,
       dpi = 300,
       width = 8,
       height = 10)

# Environmental map -------------------------------------------------------

environ_dataframe <- as.data.frame(environ_layers_cropped, xy = TRUE)

environ_map <- ggplot() +
  geom_tile(data = environ_dataframe, aes(x = x, y = y, fill = WC_tmean1)) +
  geom_sf(data = ma_shapefile, color = "black", fill = NA) +
  scale_fill_viridis_c(na.value = "white") +
  coord_sf() +
  xlim(c(-55.6,-35)) +
  ylim(c(-30,-6)) +
  theme_classic() +
  labs(fill = "Temperature") +
  theme(legend.position = c(.95,.2),
        axis.title = element_blank(), 
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank())

ggsave("c_outputs/figures/bees/environ_map.png", 
       plot = environ_map,
       dpi = 300,
       width = 3.4,
       height = 3.48)

