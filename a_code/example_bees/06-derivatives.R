
# Load packages -----------------------------------------------------------
library(truncnorm)
library(numDeriv)

# Design function ---------------------------------------------------------

pdf_truncnorm <- function(x, min, max, mean, sd) {
  dtruncnorm(x, a = min, b = max, mean = mean, sd = sd)
}

# Function to calculate the derivative of PDF at a specific point using numerical differentiation
calculate_derivative <- function(x, min, max, mean, sd) {
  grad(function(x) pdf_truncnorm(x, min, max, mean, sd), x)
}


combined_df |> 
  filter(group == "pollinators") |> 
  drop_na() |>
  summarise(
    min = min(environ_values),
    max = max(environ_values),
    mean = mean(environ_values),
    sd = sd(environ_values)
  )


# Derivatives by group ----------------------------------------------------

# df = dataframe with all environmental values and species occurrences
# y = group name
# x = environmental value to which we want to know the derivative

derivative_by_group <- function(df, y, x){
  params <- df |> 
    filter(group == y) |> 
    drop_na() |>
    summarise(
      min = min(environ_values),
      max = max(environ_values),
      mean = mean(environ_values),
      sd = sd(environ_values)
    )

  # Calculate the derivative at a specific point
  derivative <- calculate_derivative(x, params$min, params$max, 
                                     params$mean, params$sd)

  return(derivative)
}


# Derivatives by species --------------------------------------------------
# df = dataframe with all environmental values and species occurrences
# y = species name
# x = environmental value to which we want to know the derivative

derivative_by_species <- function(df, y, x){
  params <- df |> 
    filter(species == y) |> 
    drop_na() |>
    summarise(
      min = min(environ_values),
      max = max(environ_values),
      mean = mean(environ_values),
      sd = sd(environ_values)
    )
  
  # Calculate the derivative at a specific point
  derivative <- calculate_derivative(x, params$min, params$max, 
                                     params$mean, params$sd)
  
  return(derivative)
}


# Calculate derivatives per cell and species ------------------------------


# Get unique cells 

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
