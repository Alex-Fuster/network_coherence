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

