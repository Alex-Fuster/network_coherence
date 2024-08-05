# If you're running this script for the first time, make sure to run the previous 
# scripts first. 

source("a_code/example_bees/01-get_interaction_data.R")
source("a_code/example_bees/01.1-get_biome_shapefile.R")
source("a_code/example_bees/03-Environmental_data.R")
source("a_code/example_bees/04-clean_data.R")
source("a_code/example_bees/05-derivatives.R")


interaction_matrix <- interaction_matrix |>
  mutate(correlation = NA)


# Create loop to calculate the correlation between interacting species --------

for (i in 1:nrow(interaction_matrix)) {
  # if(interaction_matrix$weighted_interaction[i] > 0){
    plant <- filter(combined_df, species == interaction_matrix$spp_plants[i])
    pollinator <- filter(combined_df, species == interaction_matrix$spp_pollinators[i])
    cor_df <- merge(plant, pollinator, by = "cell", suffixes = c("_plant","_pollinator"))
    interaction_matrix$correlation[i] <- cor(cor_df$derivative_plant, 
                                             cor_df$derivative_pollinator, 
                                             method = "spearman")
  # } else { 
  #   interaction_matrix$correlation[i] <- NA
  # }
}


