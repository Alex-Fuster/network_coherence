# Load packages -----------------------------------------------------------
source("a_code/functions/ipak.R")

ipak(c("sdmpredictors", "terra"))


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

boundary = terra::ext(c(xmin = min(c(bees_occ_selected$decimalLongitude, plants_occ_selected$decimalLongitude)), 
                        xmax = max(c(bees_occ_selected$decimalLongitude, plants_occ_selected$decimalLongitude)), 
                        ymin = min(c(bees_occ_selected$decimalLatitude, plants_occ_selected$decimalLatitude)), 
                        ymax = max(c(bees_occ_selected$decimalLatitude, plants_occ_selected$decimalLatitude))))

# Cropping environmental variables based on the extent of the marine biogeographical realm of the study area

environ_layers_cropped <- terra::crop(environ_layers$WC_tmean11, boundary@ptr$vector)
