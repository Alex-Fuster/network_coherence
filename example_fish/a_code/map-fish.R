# Script to map the study area of the fish example

# NAFO shapefiles from: https://open.canada.ca/data/en/dataset/59af1c96-fc8f-4fa0-b398-d65e953eadaa
# downloaded into /b_data/

library(sf)
library(tidyverse)
library(ggOceanMaps)

canada <- rnaturalearth::ne_countries(scale = 10, returnclass = "sf", country = "Canada") 

polys = sf::read_sf("b_data/NAFO_Divisions_SHP/NAFO_Divisions_2021_poly_clipped.shp")

polys_example = polys |>
  dplyr::filter(Division %in% c("2J", "3K", "3L"))

rast <- terra::rast("b_data/HYP_LR_SR_OB_DR/HYP_LR_SR_OB_DR.tif")
plot(rast)

basemap_sf <- rnaturalearth::ne_countries(
  country = "Canada",
  returnclass = "sf") %>% 
  dplyr::select(1) %>% 
  st_transform(crs = terra::crs(rast))  # make sure crs is same as raster
plot(basemap_sf)


fish_vect = ggplot() + geom_sf(data = polys_example,
          fill = "gold", linewidth = 0, alpha = .8) |>
  terra::vect()

basemap(limits = c(-80, -43, 40, 63), rotate = TRUE,#bathy.style = "rcb",
        land.col = "grey30", land.border.col = "grey30", grid.col = NA) + #, bathy.style = "rcb", rotate = TRUE) +
  geom_sf(data = polys_example,
          fill = "gold", linewidth = 0, alpha = .8) +
  theme_void() +
  theme(legend.position = "none",
        panel.background = element_rect(fill = "#d4ebfc"))
ggsave("c_outputs/fish-example/figures/map.png", width = 5.72, height = 5.83)



