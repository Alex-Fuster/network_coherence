# Pakcages ----------------------------------------------------------------
source("./a_code/functions/ipak.R")
ipak(c("rgbif", "tidyverse"))


# Source objects ----------------------------------------------------------

# If you're starting from here, run the following lines to get only the objects you need. If you're running the code for the first time, source scripts 01 and 01.1.

load("b_data/pollinization_df.RDS")
source("a_code/example_bees/01.1-get_biome_shapefile.R")

# Get GBIF taxon keys -----------------------------------------------------

plants_checklist <- unique(complete_df$spp_plants)
polinators_checklist <- unique(complete_df$spp_polinators)

gbif_taxon_keys_polinators <- polinators_checklist |> 
  name_backbone_checklist()  |>  # match to backbone 
  filter(!matchType == "NONE") |> # get matched names
  pull(usageKey) 

gbif_taxon_keys_plants <- plants_checklist |> 
  name_backbone_checklist()  |>  # match to backbone 
  filter(!matchType == "NONE") |> # get matched names
  pull(usageKey) 


# Start query -------------------------------------------------------------
polinators_occ <- 
occ_download(
  type="and",
  pred_in("taxonKey", gbif_taxon_keys_polinators),
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

polinators_occ <- occ_download_get('0040176-240229165702484', path = "b_data/gbif") %>%
  occ_download_import()

plants_occ <- occ_download_get('0040177-240229165702484', path = "b_data/gbif") %>%
  occ_download_import()


# Select important columns ------------------------------------------------

polinators_occ_selected <- polinators_occ |> 
  select(gbifID, scientificName, decimalLatitude, decimalLongitude) |> 
  unique(species)

plants_occ_selected <- plants_occ |> 
  select(gbifID, scientificName, decimalLatitude, decimalLongitude)

# Save RData ---------------------------------------------------------------

save(plants_occ_selected, plants_occ_selected, file = "b_data/gbif_occ_selected.Rdata")
