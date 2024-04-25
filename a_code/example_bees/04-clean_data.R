
# Load packages -----------------------------------------------------------

library(tidyverse)

# Load data ----------------------------------------------------------------

load("b_data/gbif_occ_selected.Rdata")

# Combine datasets --------------------------------------------------------

combined_df <- bind_rows(
  bind_cols(species = pollinators_occ_selected$species, 
            environ_values = environ_values_pollinators$WC_tmean1,
            group = "pollinators",
            cell = environ_values_pollinators$cell,
            longitude = environ_values_pollinators$x,
            latitude = environ_values_pollinators$y),
  bind_cols(species = plants_occ_selected$species, 
            environ_values = environ_values_plants$WC_tmean1,
            group = "plants",
            cell = environ_values_plants$cell,
            longitude = environ_values_plants$x,
            latitude = environ_values_plants$y)
)

# Sampling bias ------------------------------------------------------------

# Sample unique points for each species

combined_df <- combined_df |>
  group_by(species) |> 
  distinct(cell, species, .keep_all = TRUE) |> 
  ungroup()


# Exclude species with only one point

combined_df <- combined_df |> 
  group_by(species) |> 
  filter(n() > 1,
         species != "") |> 
  ungroup()

# Clean interaction matrix ------------------------------------------------

interaction_matrix <- interaction_matrix |>
  filter(spp_plants %in% plants_occ_selected$species,
         spp_pollinators %in% pollinators_occ_selected$species)









