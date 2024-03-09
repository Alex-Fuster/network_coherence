# Load packages -----------------------------------------------------------
source("a_code/functions/ipak.R")

ipak(c("sdmpredictors", "terra", "tidyverse"))


# Load occurrence data ----------------------------------------------------

load("b_data/gbif_occ_selected.Rdata")

# Getting environmental data ----------------------------------------------

# Listing data sets available on sdmpredictors package

list_env_layers <- list_layers(datasets=c("WorldClim", "ENVIREM"), 
                               terrestrial = TRUE, 
                               marine = FALSE, 
                               freshwater = FALSE)$layer_code

subset_env_layers <- str_subset(list_env_layers, 
                                str_c(".tmean."))


# Random sample of environmental variables

environ_layers <- load_layers(sample(subset_env_layers, 1, replace = FALSE), 
                              rasterstack = FALSE, 
                              datadir = "b_data/environmental")

# Define a boundary

boundary <- terra::ext(ma_shapefile)

# Cropping environmental variables based on the extent of the occurrence data.
# Notice the name of the variable changes with sampling, you might need to change
# the value after the dollar sign for environ_layers.

environ_layers_cropped <- terra::crop(environ_layers$WC_tmean1, boundary@cpp$vector)



# Get environmental values by cell ----------------------------------------

environ_values_pollinators <- terra::extract(environ_layers_cropped, select(pollinators_occ_selected, decimalLongitude, decimalLatitude))

environ_values_plants <- terra::extract(environ_layers_cropped, select(plants_occ_selected, decimalLongitude, decimalLatitude))


combined_df <- bind_rows(
  bind_cols(species = pollinators_occ_selected$scientificName, 
            environ_values = environ_values_pollinators,
            group = "pollinators"),
  bind_cols(species = plants_occ_selected$scientificName, 
            environ_values = environ_values_plants,
            group = "plants")
)

  
# Get environmental curve -------------------------------------------------

# Just plottin'
ggplot(combined_df, aes(x=environ_values, color=group)) +
  geom_density()

# Explorin'
combined_df |> 
  group_by(group) |> 
  drop_na() |> 
  summarise(
    min = min(environ_values),
    max = max(environ_values),
    mean = mean(environ_values),
    sd = sd(environ_values)
    )
