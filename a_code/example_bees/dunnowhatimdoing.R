# Source previous scripts -------------------------------------------------
source("a_code/example_bees/01.1-get_biome_shapefile.R")
source("a_code/example_bees/03-Environmental_data.R") # careful with the last line!
source("a_code/example_bees/04-clean_data.R")
source("a_code/example_bees/05-curves.R")
source("a_code/example_bees/06-derivatives.R")

# Get unique species names and environmental values
species_names <- unique(combined_df$species)
  
environ_values <- combined_df |> 
  dplyr::select(environ_values) |>
  drop_na() |>
  unique()

# Initialize an empty list to store derivative matrices for each species
derivative_df <- data.frame(species_name = character(),
                            environ_values = numeric(),
                            derivative = numeric())

# Loop through each species and calculate derivatives for environmental values
for (i in species_names) {
  for (j in environ_values){
    
  # Calculate derivatives for environmental values
  derivatives <- derivative_by_species(combined_df, i, j)
  # Append values to the derivative_matrices data frame
  derivative_df <- rbind(derivative_df, 
                                data.frame(species_name = i, 
                                           environ_values = j, 
                                           derivative = derivatives))
}}


# Merge datasets ----------------------------------------------------------

# Merge the derivative matrices with the combined_df
combined_df <- combined_df |>
  left_join(derivative_df, by = c("species" = "species_name",
                                  "environ_values" = "environ_values"))

# Clean NAs and empty data

combined_df <- combined_df |>
  filter(!is.na(derivative)) |> 
  mutate(group = as.factor(group))

# Load interacting data ---------------------------------------------------

load("b_data/pollinization_df.RDS")

# Correlation matrix for interacting species ------------------------------



# Correlation matrix for non-interacting species --------------------------


