
# Load packages -----------------------------------------------------------
source("a_code/functions/ipak.R")

ipak(c("tidyverse", "corrplot"))

# Load scripts ------------------------------------------------------------
source("a_code/example_bees/06-correlations.R")

# Plot pairwise correlations ----------------------------------------------
# Make correlation data frame as a matrix

correlation_matrix <- 
  interaction_matrix |> 
  select(-weighted_interaction) |> # exclude the correlation column
  pivot_wider(names_from = spp_pollinators, # make the matrix wider
              values_from = correlation) |>  # fill the matrix with the weighted interactions
  column_to_rownames(var = "spp_plants") |> # make the plants the row names
  as.matrix() # make it a matrix

# Plot as a heatmap of pairwise correlations

png(height=1800, width=1800, 
    file="c_outputs/figures/bees/heatmap_species_associations.png")

plot.new()
corrplot::corrplot(correlation_matrix, 
                   type = "lower",
                   method = "color", 
                   tl.cex = 2.5, cl.cex = 3, tl.col = "black"
                   )

dev.off()
