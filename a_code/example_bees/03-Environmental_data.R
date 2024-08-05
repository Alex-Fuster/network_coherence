# Load packages -----------------------------------------------------------
source("a_code/functions/ipak.R")

ipak(c("sdmpredictors", "terra", "tidyverse"))


# Load occurrence data ----------------------------------------------------

load("b_data/gbif_occ_selected.Rdata") # Load plants and pollinators occurrence data

# Getting environmental data ----------------------------------------------

# Listing data sets available on sdmpredictors package

list_env_layers <- list_layers(
  datasets=c("WorldClim", "ENVIREM"), # Selecting datasets WorldClim and ENVIREM
  terrestrial = TRUE,                 # Selecting only terrestrial layers
  marine = FALSE,                     # Exclude marine layers
  freshwater = FALSE)$layer_code      # Exclude freshwater layers and get the layer codes

# Subsetting those layers that are temperature related
subset_env_layers <- str_subset(list_env_layers, 
                                str_c(".tmean.")) # Detecting layers with "tmean" in the name


# Random sample of environmental variables

set.seed(42)
environ_layers <- load_layers(sample(
  subset_env_layers, 1, replace = FALSE), # Selecting one layer without replacement
  rasterstack = FALSE,                    # Return a single raster layer
  datadir = "b_data/environmental")       # Save the data in the folder "b_data/environmental"


# Cropping environmental variables based on the extent of the occurrence data.

environ_layers_cropped <- terra::crop(
  terra::mask(rast(environ_layers[[1]]), ma_shapefile), # Masking the raster with the shapefile
    ma_shapefile)                                       # Cropping the raster with the shapefile

# Get environmental values by cell ----------------------------------------

environ_values_pollinators <- terra::extract(environ_layers_cropped, 
                              # Selects the coordinates from pollinators occurrence data
                                             dplyr::select(pollinators_occ_selected, 
                                                           decimalLongitude, 
                                                           decimalLatitude), 
                                             cells = TRUE, # Keep the cell number
                                             xy = TRUE)    # Keep the coordinates

environ_values_plants <- terra::extract(environ_layers_cropped, 
                         # Selects the coordinates from plants occurrence data
                                        dplyr::select(plants_occ_selected, 
                                                      decimalLongitude, 
                                                      decimalLatitude),
                                        cells = TRUE, # Keep the cell number
                                        xy = TRUE)    # Keep the coordinates
