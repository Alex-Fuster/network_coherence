# Packages loading --------------------------------------------------------

source("a_code/functions/ipak.R")

ipak(c("tidyverse", "googlesheets4"))


# Load species data and splitting

interaction_data_raw <- read_sheet("1FdADRdStEK-mYwXDi1bshKAHUib1S6Ftvfx_NXDeIbs", trim_ws = TRUE)

interaction_data_clean <- interaction_data_raw |> 
  select(c(id, genero_abelha_catalogo, especie_abelha_catalogo, subespecie_abelha_catalogo, genero_planta_catalogo, especie_planta_catalogo)) |> 
  na.exclude() |> 
  mutate(especie_abelha_catalogo = as.character(especie_abelha_catalogo))


# Recoding values ---------------------------------------------------------
interaction_data_clean$especie_abelha_catalogo <- interaction_data_clean$especie_abelha_catalogo |> 
  recode("NULL" = "")

interaction_data_clean$genero_planta_catalogo <- interaction_data_clean$genero_planta_catalogo |> 
  recode("Eugenia (3 spp.)" = "Eugenia")

interaction_data_clean$especie_planta_catalogo <- interaction_data_clean$especie_planta_catalogo |> 
  recode("sp. (3 spp.)" = "spp.",
         "sp. (3spp)" = "spp.",
         "spp" = "spp.",
         "2" = "spp.")

# Binding species names columns -------------------------------------------
interaction_data_clean <- 
interaction_data_clean |> 
  unite(col = "bee_species",
        genero_abelha_catalogo:subespecie_abelha_catalogo,
        sep = " ",
        remove = TRUE) |> 
  unite(col = "plant_species",
        genero_planta_catalogo:especie_planta_catalogo,
        sep = " ",
        remove = TRUE)

rm(interaction_data_raw)
