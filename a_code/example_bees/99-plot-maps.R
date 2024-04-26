
# Load packages -----------------------------------------------------------

# Source code -------------------------------------------------------------
source("a_code/example_bees/06-correlations.R")

# Species richness map ----------------------------------------------------
plot(environ_layers_cropped, col = "grey",
     legend = FALSE, add = TRUE)

combined_df |> 
  group_by(cell) |> 
  summarise(richness = n()) |>
  left_join(combined_df, by = "cell") |>
  select(longitude, latitude, richness) |>
  raster::rasterFromXYZ(res = res(environ_layers_cropped), 
                        crs = raster::crs(environ_layers_cropped))  |> 
  plot(col = viridis::viridis(10), legend = FALSE)
