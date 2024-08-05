# Packages loading --------------------------------------------------------

source("a_code/functions/ipak.R")

ipak(c("tidyverse", "openxlsx"))



# Load species data and splitting -----------------------------------------
# Varassin I.G. & Sazima M. (2012) Spatial heterogeneity and the distribution of
# bromeliad pollinators in the Atlantic Forest. Acta Oecologica, 43: 104-112.

interaction_matrix <- read.xlsx("http://www.ecologia.ib.usp.br/iwdb/data/plant_pollinator/excel/varassin_sazima_2012.xlsx",
                                  fillMergedCells = FALSE, # to avoid merged cells
                                  rows = c(2:24),          # select rows with data
                                  colNames = TRUE,         # keep column names
                                  sep.names = " ")         # column names separator

interaction_matrix <- interaction_matrix |> 
  rename("spp_plants" = "X1")                              # renaming column


# Creating long df --------------------------------------------------------

interaction_matrix <- interaction_matrix |> 
  pivot_longer(-1,                                    # Exclude the first column in the pivoting
               names_to = "spp_pollinators",          # New column with pollinators
               values_to = "weighted_interaction") |> # New column with interactions
  # Filter out species with incomplete identification
  filter(!str_detect(spp_plants, " sp"), 
         !str_detect(spp_pollinators, " sp.")) 
