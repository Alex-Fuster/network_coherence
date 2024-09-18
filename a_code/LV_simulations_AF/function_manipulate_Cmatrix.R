# Load necessary library
library(ggplot2)
library(MASS)
library(dplyr)
library(tidyr)
library(radiant.data)

# Function to create a correlation matrix with specified mean and sd
createCorMat <- function(S, type, target_mean = 0, target_sd = 0.1){
  # Define the number of off-diagonal elements
  n_off_diag <- S * (S - 1) / 2
  
  # Generate random values for the upper triangular part
  if (type == "positive") {
    # Generate correlations between 0 and 1
    cor_values <- rbeta(n_off_diag, shape1 = 2, shape2 = 5)  # Beta distribution for positive correlations
  } else if (type == "mixed") {
    # Generate correlations between -1 and 1
    cor_values <- runif(n_off_diag, -1, 1)
  } else {
    stop("type not defined")
  }
  
  # Adjust mean and sd to match the target
  cor_values <- scale(cor_values, center = mean(cor_values), scale = sd(cor_values))  # Standardize values
  cor_values <- cor_values * target_sd + target_mean  # Rescale to desired mean and sd
  
  # Ensure values stay between -1 and 1
  cor_values[cor_values < -1] <- -1
  cor_values[cor_values > 1] <- 1
  
  # Create a symmetric matrix with the correlation values
  corMat <- matrix(1, nrow = S, ncol = S)
  corMat[upper.tri(corMat)] <- cor_values
  corMat <- corMat * t(corMat)  # Make the matrix symmetric
  
  
  return(corMat)
}


# Create a correlation matrix with target mean 0 and standard deviation 0.2
S <- 20  # Number of species
cor_matrix <- createCorMat(S, type = "mixed", target_mean = 0, target_sd = 0.3)



cor_matrix %>%
   as_tibble() %>%
   rownames_to_column("Var1") %>%
   pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
   ggplot(aes(Var1, Var2)) +
   geom_tile(aes(fill = value)) +
   scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-1,1)*max(abs(cor_matrix))) +
   ggtitle("Covariance Matrix")





