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

plants_checklist <- unique(complete_df$spp_plants)
pollinators_checklist <- unique(complete_df$spp_pollinators)

gbif_taxon_keys_pollinators <- pollinators_checklist |> 
  name_backbone_checklist()  |>  # match to backbone 
  filter(!matchType == "NONE") |> # get matched names
  pull(usageKey) 

gbif_taxon_keys_plants <- plants_checklist |> 
  name_backbone_checklist()  |>  # match to backbone 
  filter(!matchType == "NONE") |> # get matched names
  pull(usageKey) 


# Start query -------------------------------------------------------------
pollinators_occ <- 
occ_download(
  type="and",
  pred_in("taxonKey", gbif_taxon_keys_pollinators),
  pred_within(ma_shapefile_wkt),
  pred("hasGeospatialIssue", FALSE),
  pred("hasCoordinate", TRUE),
  pred("occurrenceStatus","PRESENT"), 
  pred_gte("year", 1974),
  pred_not(pred_in("basisOfRecord","FOSSIL_SPECIMEN")),
  pred("country","BR"),
  pred_or(
    pred_not(pred_in("establishmentMeans",c("MANAGED","INTRODUCED"))),
    pred_isnull("establishmentMeans")
  ),
  pred_or(  
    pred_lt("coordinateUncertaintyInMeters",10000),
    pred_isnull("coordinateUncertaintyInMeters")
  ),
  format = "SIMPLE_CSV"
)

plants_occ <- 
  occ_download(
    type="and",
    pred_in("taxonKey", gbif_taxon_keys_plants),
    pred_within(ma_shapefile_wkt),
    pred("hasGeospatialIssue", FALSE),
    pred("hasCoordinate", TRUE),
    pred("occurrenceStatus","PRESENT"), 
    pred_gte("year", 1974),
    pred_not(pred_in("basisOfRecord","FOSSIL_SPECIMEN")),
    pred("country","BR"),
    pred_or(
      pred_not(pred_in("establishmentMeans",c("MANAGED","INTRODUCED"))),
      pred_isnull("establishmentMeans")
    ),
    pred_or(  
      pred_lt("coordinateUncertaintyInMeters",10000),
      pred_isnull("coordinateUncertaintyInMeters")
    ),
    format = "SIMPLE_CSV"
  )


# Retrieve data -----------------------------------------------------------
dir.create("b_data/gbif")

pollinators_occ <- occ_download_get('0040176-240229165702484', path = "b_data/gbif") |> 
  occ_download_import()

plants_occ <- occ_download_get('0040177-240229165702484', path = "b_data/gbif") |> 
  occ_download_import()


# Select important columns ------------------------------------------------

pollinators_occ_selected <- pollinators_occ |> 
  dplyr::select(gbifID, taxonKey, species, decimalLatitude, decimalLongitude)

plants_occ_selected <- plants_occ |> 
  dplyr::select(gbifID, taxonKey, species, decimalLatitude, decimalLongitude)

# Save RData ---------------------------------------------------------------

save(pollinators_occ_selected, plants_occ_selected, file = "b_data/gbif_occ_selected.Rdata")
