# Pakcages ----------------------------------------------------------------
source("./a_code/functions/ipak.R")
ipak(c("rgbif", "tidyverse"))


# Source objects ----------------------------------------------------------

# If you're starting from here, run the following lines to get only the objects you need. If you're running the code for the first time, source scripts 01 and 01.1.

load("b_data/pollinization_df.RDS")
source("a_code/example_bees/01.1-get_biome_shapefile.R")

# Get GBIF taxon keys -----------------------------------------------------
# Make sure to have your GBIF credentials in place for the code below. 
# Read more about how to do it here: 
# https://docs.ropensci.org/rgbif/articles/gbif_credentials.html

plants_checklist <- unique(interaction_matrix$spp_plants) # get unique plant species names
pollinators_checklist <- unique(interaction_matrix$spp_pollinators) # get unique pollinator species names

gbif_taxon_keys_pollinators <- pollinators_checklist |> 
  name_backbone_checklist()  |>   # match to backbone 
  filter(!matchType == "NONE") |> # get matched names
  pull(usageKey)                  # get taxon keys

gbif_taxon_keys_plants <- plants_checklist |> 
  name_backbone_checklist()  |>   # match to backbone 
  filter(!matchType == "NONE") |> # get matched names
  pull(usageKey)                  # get taxon keys


# Start query -------------------------------------------------------------
pollinators_occ <- 
occ_download(
  type="and",                                           # all conditions must be met
  pred_in("taxonKey", gbif_taxon_keys_pollinators),     # only pollinator species, search by key
  pred_within(ma_shapefile_wkt),                        # only occurrences within the biome
  pred("hasGeospatialIssue", FALSE),                    # only occurrences without geospatial issues
  pred("hasCoordinate", TRUE),                          # only occurrences with coordinates
  pred("occurrenceStatus","PRESENT"),                   # only present occurrences
  pred_gte("year", 1974),                               # only occurrences from 1974 onwards
  pred_not(pred_in("basisOfRecord","FOSSIL_SPECIMEN")), # only occurrences that are not fossils
  pred("country","BR"),                                 # only occurrences from Brazil
  pred_or(                                              # either the establishment means is not managed or introduced
    pred_not(pred_in("establishmentMeans",c("MANAGED","INTRODUCED"))),
    pred_isnull("establishmentMeans")
  ),
  pred_or(                                            # either the coordinate uncertainty is less than 10,000 meters or null
    pred_lt("coordinateUncertaintyInMeters",10000),
    pred_isnull("coordinateUncertaintyInMeters")
  ),
  format = "SIMPLE_CSV"                               # format of the output
)

plants_occ <- 
  occ_download(
    type="and",                                           # all conditions must be met
    pred_in("taxonKey", gbif_taxon_keys_plants),          # only plant species, search by key
    pred_within(ma_shapefile_wkt),                        # only occurrences within the biome
    pred("hasGeospatialIssue", FALSE),                    # only occurrences without geospatial issues
    pred("hasCoordinate", TRUE),                          # only occurrences with coordinates
    pred("occurrenceStatus","PRESENT"),                   # only present occurrences
    pred_gte("year", 1974),                               # only occurrences from 1974 onwards
    pred_not(pred_in("basisOfRecord","FOSSIL_SPECIMEN")), # only occurrences that are not fossils
    pred("country","BR"),                                 # only occurrences from Brazil
    pred_or(                                              # either the establishment means is not managed or introduced
      pred_not(pred_in("establishmentMeans",c("MANAGED","INTRODUCED"))),
      pred_isnull("establishmentMeans")
    ),
    pred_or(                                              # either the coordinate uncertainty is less than 10,000 meters or null
      pred_lt("coordinateUncertaintyInMeters",10000),
      pred_isnull("coordinateUncertaintyInMeters")
    ),
    format = "SIMPLE_CSV"                                 # format of the output
  )


# Retrieve data -----------------------------------------------------------
dir.create("b_data/gbif")                                      # create a folder to save the data

pollinators_occ <- occ_download_get('0204755-240321170329656', # Download results of the query above
                                    path = "b_data/gbif")  |>  # Save the data in the folder b_data/gbif
  occ_download_import()                                        # Import the data

plants_occ <- occ_download_get('0204757-240321170329656',      # Download results of the query above
                               path = "b_data/gbif") |>        # Save the data in the folder b_data/gbif
  occ_download_import()                                        # Import the data


# Select important columns ------------------------------------------------

pollinators_occ_selected <- pollinators_occ |> 
  dplyr::select(gbifID, taxonKey, species, decimalLatitude, decimalLongitude)

plants_occ_selected <- plants_occ |> 
  dplyr::select(gbifID, taxonKey, species, decimalLatitude, decimalLongitude)

# Save RData ---------------------------------------------------------------

save(pollinators_occ_selected, plants_occ_selected, file = "b_data/gbif_occ_selected.Rdata")
