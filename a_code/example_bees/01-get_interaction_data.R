# Packages loading --------------------------------------------------------

source("a_code/functions/ipak.R")

ipak(c("tidyverse", "openxlsx"))



# Load species data and splitting -----------------------------------------
# Varassin I.G. & Sazima M. (2012) Spatial heterogeneity and the distribution of
# bromeliad pollinators in the Atlantic Forest. Acta Oecologica, 43: 104-112.

interaction_matrix <- read.xlsx("http://www.ecologia.ib.usp.br/iwdb/data/plant_pollinator/excel/varassin_sazima_2012.xlsx",
                                  fillMergedCells = FALSE,
                                  rows = c(2:24),
                                  colNames = TRUE,
                                  sep.names = " ")

interaction_matrix <- interaction_matrix |> 
  rename("spp_plants" = "X1")


# Creating long df --------------------------------------------------------

interaction_matrix <- interaction_matrix |> 
  pivot_longer(-1,
               names_to = "spp_pollinators",
               values_to = "weighted_interaction") |> 
  filter(!str_detect(spp_plants, " sp"),
         !str_detect(spp_pollinators, " sp.")) 
