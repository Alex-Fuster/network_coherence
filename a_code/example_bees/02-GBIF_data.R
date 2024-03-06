# Pakcages ----------------------------------------------------------------
source("./a_code/functions/ipak.R")
ipak("rgbif")

# Get GBIF taxon keys -----------------------------------------------------

bees_checklist <- unique(interaction_data_clean$bee_species)
plants_checklist <- unique(interaction_data_clean$plant_species)

gbif_taxon_keys_bees <- bees_checklist |> 
  name_backbone_checklist()  |>  # match to backbone 
  filter(!matchType == "NONE") |> # get matched names
  pull(usageKey) 

gbif_taxon_keys_plants <- plants_checklist |> 
  name_backbone_checklist()  |>  # match to backbone 
  filter(!matchType == "NONE") |> # get matched names
  pull(usageKey) 



# Start query -------------------------------------------------------------
bees_occ <- 
occ_download(
  type="and",
  pred_in("taxonKey", gbif_taxon_keys_bees),
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

bees_occ <- occ_download_get('0024325-240229165702484', path = "b_data/gbif") %>%
  occ_download_import()

plants_occ <- occ_download_get('0024326-240229165702484', path = "b_data/gbif") %>%
  occ_download_import()


# Select important columns ------------------------------------------------

bees_occ_selected <- bees_occ |> 
  select(gbifID, scientificName, decimalLatitude, decimalLongitude) |> 
  unique(species)

plants_occ_selected <- plants_occ |> 
  select(gbifID, scientificName, decimalLatitude, decimalLongitude)

# Save RData ---------------------------------------------------------------

save(bees_occ_selected, plants_occ_selected, file = "b_data/gbif_occ_selected.Rdata")
