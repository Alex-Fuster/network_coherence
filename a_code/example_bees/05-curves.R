# Get environmental values by cell ----------------------------------------

environ_values_pollinators <- terra::extract(environ_layers_cropped, dplyr::select(pollinators_occ_selected, decimalLongitude, decimalLatitude))

environ_values_plants <- terra::extract(environ_layers_cropped, dplyr::select(plants_occ_selected, decimalLongitude, decimalLatitude))


combined_df <- bind_rows(
  bind_cols(species = pollinators_occ_selected$species, 
            environ_values = environ_values_pollinators,
            group = "pollinators"),
  bind_cols(species = plants_occ_selected$species, 
            environ_values = environ_values_plants,
            group = "plants")
)


# Get environmental curve -------------------------------------------------

# Just plottin'
ggplot(combined_df, aes(x=environ_values, color=species)) +
  theme_void() +
  geom_density(show.legend = FALSE)

