# Calculate ENC -----------------------------------------------------------

# Get unique species names
species_names <- unique(combined_df$species)
environ_values <- combined_df |> 
  dplyr::select(environ_values) |>
  drop_na() |>
  unique()

# Initialize an empty list to store derivative matrices for each species
derivative_matrices <- data.frame(species_name = character(), 
                                  environ_values = numeric(),
                                  derivative = numeric())

# Loop through each species and calculate derivatives for environmental values
for (i in species_names) {
  for (j in environ_values){
    
  # Calculate derivatives for environmental values
  derivatives <- derivative_by_species(combined_df, i, j)
  # Append values to the data frame
  derivative_matrices <- rbind(derivative_matrices, c(i, j, derivatives))
}}

# Combine derivative matrices into a single matrix
combined_matrix <- do.call(rbind, derivative_matrices)
}
# # Set row names to species names
# rownames(combined_matrix) <- i
# 
# # Set column names to environmental values
# colnames(combined_matrix) <- j
# return(combined_matrix)
# }

