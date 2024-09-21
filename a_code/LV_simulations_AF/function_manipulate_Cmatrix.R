
library(ggplot2)
library(MASS)
library(dplyr)
library(tidyr)
library(radiant.data)

# Function to create a correlation matrix with specified mean and sd
createCorMat <- function(S, target_mean = 0, target_sd = 0.1){
  # Define the number of off-diagonal elements
  n_off_diag <- S * (S - 1) / 2
  
  cor_values <- runif(n_off_diag, -1, 1)
  
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


# Create a correlation matrix with target mean and standard deviation
S <- 20  # Number of species
cor_matrix <- createCorMat(S,  target_mean = 0, target_sd = 0.3)


#cor_matrix[which(lower.tri(cor_matrix))] <- NA

cor_matrix <- as.data.frame(cor_matrix)
cor_matrix$Var1 = paste0("V", rownames(cor_matrix))
cor_matrix = cor_matrix %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") #|> na.omit()
cor_matrix$Var1 = factor(cor_matrix$Var1, levels = unique(cor_matrix$Var1))
cor_matrix$Var2 = factor(cor_matrix$Var2, levels = rev(unique(cor_matrix$Var2)))


# plot the matrix

ggplot(data = cor_matrix,
       aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-1,1)*max(abs(cor_matrix$value))) +
  ggtitle("Covariance Matrix")


# recover the matrix format and well adjusted

## Pivot the long-format dataframe back to wide format
wide_cor_matrix <- cor_matrix %>%
  pivot_wider(names_from = Var2, values_from = value)

## Remove the row identifier column (Var1) to convert it back to a matrix
wide_cor_matrix <- wide_cor_matrix %>%
  select(-Var1)

## Convert the dataframe back to a matrix
recovered_matrix <- as.matrix(wide_cor_matrix)

# Print the recovered matrix
print(recovered_matrix)


################################# COVARIANCE MATRIX

sd_X = rep(0.1, 20)

covMat <- as.matrix(Matrix::nearPD(sd_X %*% t(sd_X) * as.matrix(recovered_matrix))$mat)



# Convert the covariance matrix to a dataframe and add row/column labels
covMat <- as.data.frame(covMat)
covMat$Var1 <- paste0("V", rownames(covMat))  # Create row labels

# Reshape the dataframe into long format
covMat_long <- covMat %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")

# Set factor levels for proper ordering of rows and columns
covMat_long$Var1 <- factor(covMat_long$Var1, levels = unique(covMat_long$Var1))  # Row order (Var1)
covMat_long$Var2 <- factor(covMat_long$Var2, levels = rev(unique(covMat_long$Var2)))  # Reverse order for Var2

# Plot the covariance matrix using ggplot2
ggplot(data = covMat_long, aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-1, 1) * max(abs(covMat_long$value), na.rm = TRUE)) +
  ggtitle("Covariance Matrix") +
  theme_minimal()

###########################################


# # Generate random values for the upper triangular part
# if (type == "positive") {
#   # Generate correlations between 0 and 1
#   cor_values <- rbeta(n_off_diag, shape1 = 2, shape2 = 5)  # Beta distribution for positive correlations
# } else if (type == "mixed") {
#   # Generate correlations between -1 and 1
#   cor_values <- runif(n_off_diag, -1, 1)
# } else {
#   stop("type not defined")
#}

# Generate U-shaped distribution using Beta(0.5, 0.5)
n <- 10000
u_shaped_values <- rbeta(n, shape1 = 0.5, shape2 = 0.5)
u_shaped_values <- 2 * (u_shaped_values - 0.5)  # Map to [-1, 1]

# Plot the U-shaped distribution
ggplot(data.frame(x = u_shaped_values), aes(x = x)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "blue", alpha = 0.7) +
  ggtitle("U-shaped Distribution") +
  labs(x = "Value", y = "Density") +
  theme_minimal()


# Function to create a correlation matrix with different distributions
createCorMat <- function(S, distribution_type = "normal", 
                         mean = 0, sd = 1, 
                         beta_shape1 = NULL, beta_shape2 = NULL, 
                         ushape1 = NULL, ushape2 = NULL) {
  # Define the number of off-diagonal elements
  n_off_diag <- S * (S - 1) / 2
  
  # Generate correlation values based on the chosen distribution
  if (distribution_type == "normal") {
    # Normal distribution with specified mean and sd
    cor_values <- rnorm(n_off_diag, mean = mean, sd = sd)
    
  } else if (distribution_type == "beta") {
    # Beta distribution with specified shape parameters
    if (is.null(beta_shape1) || is.null(beta_shape2)) {
      stop("For Beta distribution, 'beta_shape1' and 'beta_shape2' must be specified.")
    }
    cor_values <- rbeta(n_off_diag, shape1 = beta_shape1, shape2 = beta_shape2)
    cor_values <- 2 * (cor_values - 0.5)  # Map Beta values from [0, 1] to [-1, 1]
    
  } else if (distribution_type == "u_shaped") {
    # U-shaped distribution with user-specified shape parameters
    if (is.null(ushape1) || is.null(ushape2)) {
      stop("For U-shaped distribution, 'ushape1' and 'ushape2' must be specified.")
    }
    cor_values <- rbeta(n_off_diag, shape1 = ushape1, shape2 = ushape2)
    cor_values <- 2 * (cor_values - 0.5)  # Map Beta values to [-1, 1]
    
  } else {
    stop("Invalid distribution type specified.")
  }
  
  # Adjust mean and sd to match the target (only for normal distribution)
  if (distribution_type == "normal") {
    cor_values <- scale(cor_values, center = mean(cor_values), scale = sd(cor_values))  # Standardize values
    cor_values <- cor_values * sd + mean  # Rescale to desired mean and sd
  }
  
  # Ensure values stay between -1 and 1 for all distributions
  cor_values[cor_values < -1] <- -1
  cor_values[cor_values > 1] <- 1
  
  # Create a symmetric matrix with the correlation values
  corMat <- matrix(1, nrow = S, ncol = S)
  corMat[upper.tri(corMat)] <- cor_values
  corMat <- corMat * t(corMat)  # Make the matrix symmetric
  
  return(corMat)
}



S <- 20
cor_matrix_normal <- createCorMat(S, distribution_type = "normal", mean = 0, sd = 0.3)


cor_matrix_beta <- createCorMat(S, distribution_type = "beta", beta_shape1 = 2, beta_shape2 = 5)


cor_matrix_u <- createCorMat(S, distribution_type = "u_shaped", ushape1 = 0.5, ushape2 = 0.5)


cor_matrix <- cor_matrix_beta

#cor_matrix[which(lower.tri(cor_matrix))] <- NA

cor_matrix <- as.data.frame(cor_matrix)
cor_matrix$Var1 = paste0("V", rownames(cor_matrix))
cor_matrix = cor_matrix %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") #|> na.omit()
cor_matrix$Var1 = factor(cor_matrix$Var1, levels = unique(cor_matrix$Var1))
cor_matrix$Var2 = factor(cor_matrix$Var2, levels = rev(unique(cor_matrix$Var2)))


# plot the matrix

ggplot(data = cor_matrix,
       aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-1,1)*max(abs(cor_matrix$value))) +
  ggtitle("Covariance Matrix")










############ U distrib

library(ggplot2)

# Generate U-shaped distribution using specified shape parameters
n <- 10000
u_shaped_values <- rbeta(n, shape1 = 0.5, shape2 = 0.5)  # Example values for the U-shape
u_shaped_values <- 2 * (u_shaped_values - 0.5)  # Map to [-1, 1]

# Plot the U-shaped distribution
ggplot(data.frame(x = u_shaped_values), aes(x = x)) +
  geom_histogram(aes(y = ..density..), bins = 50, fill = "blue", alpha = 0.7) +
  ggtitle("U-shaped Distribution") +
  labs(x = "Value", y = "Density") +
  theme_minimal()

################# B distrib

# Function to generate and plot Beta distribution, mapped to [-1, 1]
plot_mapped_beta_distribution <- function(shape1 = 2, shape2 = 5, n = 10000) {
  
  # Generate Beta distribution using specified shape parameters
  beta_values <- rbeta(n, shape1 = shape1, shape2 = shape2)
  
  # Map Beta values from [0, 1] to [-1, 1]
  beta_values_mapped <- 2 * (beta_values - 0.5)
  
  # Plot the Beta distribution mapped to [-1, 1]
  ggplot(data.frame(x = beta_values_mapped), aes(x = x)) +
    geom_histogram(aes(y = ..density..), bins = 50, fill = "blue", alpha = 0.7) +
    ggtitle(paste("Beta Distribution Mapped to [-1, 1] with shape1 =", shape1, "and shape2 =", shape2)) +
    labs(x = "Mapped Value", y = "Density") +
    theme_minimal()
}

# Example usage: Mapped Beta distribution with shape1 = 2, shape2 = 5
plot_mapped_beta_distribution(shape1 = 2, shape2 = 5)

# Example usage: Left-skewed Beta distribution with shape1 = 5, shape2 = 2
plot_mapped_beta_distribution(shape1 = 5, shape2 = 2)

# Example usage: U-shaped Beta distribution with shape1 = 0.5, shape2 = 0.5
plot_mapped_beta_distribution(shape1 = 0.5, shape2 = 0.5)
