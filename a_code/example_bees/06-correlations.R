source("a_code/example_bees/01-get_interaction_data.R")


interaction_matrix <- interaction_matrix |>
  mutate(correlation = NA)


for (i in 1:nrow(interaction_matrix)) {
  ifelse(interaction_matrix$weighted_interaction > 0,
    
    plant <- filter(combined_df, species == interaction_matrix$spp_plants[i])
    pollinator <- filter(combined_df, species == interaction_matrix$spp_pollinators[i])
    cor_df <- merge(plant, pollinator, by = "cell", suffixes = c("_plant","_pollinator"))
    interaction_matrix$correlation[i] <- cor(cor_df$derivative_plant, 
                                             cor_df$derivative_pollinator, 
                                             method = "spearman"),

    interaction_matrix$correlation[i] <- NA
}

