# Get environmental values by cell ----------------------------------------

environ_values_pollinators <- terra::extract(environ_layers_cropped, 
                                             dplyr::select(pollinators_occ_selected, 
                                                           decimalLongitude, 
                                                           decimalLatitude),
                                             cells = TRUE,
                                             xy = TRUE)

environ_values_plants <- terra::extract(environ_layers_cropped, 
                                        dplyr::select(plants_occ_selected, 
                                                      decimalLongitude, 
                                                      decimalLatitude),
                                        cells = TRUE,
                                        xy = TRUE)


combined_df <- bind_rows(
  bind_cols(species = pollinators_occ_selected$species, 
            environ_values = environ_values_pollinators$WC_tmean1,
            group = "pollinators",
            cell = environ_values_pollinators$cell,
            longitude = environ_values_pollinators$x,
            latitude = environ_values_pollinators$y),
  bind_cols(species = plants_occ_selected$species, 
            environ_values = environ_values_plants$WC_tmean1,
            group = "plants",
            cell = environ_values_plants$cell,
            longitude = environ_values_plants$x,
            latitude = environ_values_plants$y)
)


# Get environmental curve -------------------------------------------------

# Just plottin'
ggplot(combined_df, aes(x=environ_values, color=species)) +
  theme_void() +
  geom_density(show.legend = FALSE)

