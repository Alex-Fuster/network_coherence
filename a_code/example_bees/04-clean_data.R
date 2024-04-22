
# Load packages -----------------------------------------------------------

library(tidyverse)
library(bdc)
library(taxadb)
library(rgbif)
library(dismo)

# Load data ----------------------------------------------------------------

load("b_data/gbif_occ_selected.Rdata")

# Sampling bias ------------------------------------------------------------

# sample:
pollinators_occ_unique <- gridSample(data.frame(pollinators_occ_selected[,c("decimalLongitude","decimalLatitude")]),
                   environ_layers_cropped, n=1)

plants_occ_unique <- gridSample(data.frame(plants_occ_selected[,c("decimalLongitude","decimalLatitude")]),
                                     environ_layers_cropped, n=1)

# Merge datasets

plants_occ_selected <- plants_occ_selected |> 
  semi_join(plants_occ_unique, by = c("decimalLatitude", "decimalLongitude"))


pollinators_occ_selected <- pollinators_occ_selected |>
  semi_join(pollinators_occ_unique, by = c("decimalLatitude", "decimalLongitude"))

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
