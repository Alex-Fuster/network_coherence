# Source previous scripts -------------------------------------------------
source("a_code/example_bees/01.1-get_biome_shapefile.R")
source("a_code/example_bees/03-Environmental_data.R") 
source("a_code/example_bees/04-clean_data.R")
source("a_code/example_bees/05-curves.R")
source("a_code/example_bees/06-derivatives.R")

# Get unique species names and environmental values
# species_names <- unique(combined_df$species)
  
cell_values <- combined_df |> 
  dplyr::select(cell) |>
  unique()

# Remove NAs from environ_values
combined_df <- combined_df |>
  filter(!is.na(environ_values))


# Create new column for derivatives
combined_df <- combined_df |> 
  mutate(derivative = NA)

# Loop through each species and calculate derivatives for environmental values
for (j in 1:nrow(combined_df)){
    
  # Calculate derivatives for environmental values
  combined_df$derivative[j] <- derivative_by_species(combined_df, 
                                       combined_df$species[j], 
                                       combined_df$environ_values[j])
  
}


# Clean NAs and empty data

combined_df <- combined_df |>
  filter(!is.na(derivative)) |> 
  mutate(group = as.factor(group))

# Load interacting data ---------------------------------------------------

load("b_data/pollinization_df.RDS")

# Correlation matrix for interacting species ------------------------------



# Correlation matrix for non-interacting species --------------------------


