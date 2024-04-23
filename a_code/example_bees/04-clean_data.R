
# Load packages -----------------------------------------------------------

library(tidyverse)

# Load data ----------------------------------------------------------------

load("b_data/gbif_occ_selected.Rdata")

# Sampling bias ------------------------------------------------------------

# Sample unique points for each species
pollinators_occ_selected <- pollinators_occ_selected |>
  group_by(species) |> 
  distinct(decimalLatitude, decimalLongitude) |> 
  ungroup()


plants_occ_selected <- plants_occ_selected |>
  group_by(species) |> 
  distinct(decimalLatitude, decimalLongitude) |> 
  ungroup()

# Exclude species with only one point

plants_occ_selected <- plants_occ_selected |> 
  group_by(species) |> 
  filter(n() > 1,
         species != "") |> 
  ungroup()

pollinators_occ_selected <- pollinators_occ_selected |>
  group_by(species) |>
  filter(n() > 1,
         species != "")|> 
  ungroup()
