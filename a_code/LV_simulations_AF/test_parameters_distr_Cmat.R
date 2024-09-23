library(ggplot2)
library(tidyr)
library(dplyr)
library(Matrix)
library(deSolve)

# Function to create a correlation matrix with different distributions
createCorMat <- function(S, distribution_type = "normal", 
                         mean = 0, sd = 1, 
                         beta_shape1 = NULL, beta_shape2 = NULL, 
                         ushape1 = NULL, ushape2 = NULL) {
  n_off_diag <- S * (S - 1) / 2
  
  # Generate correlation values based on the chosen distribution
  if (distribution_type == "normal") {
    cor_values <- rnorm(n_off_diag, mean = mean, sd = sd)
  } else if (distribution_type == "beta") {
    cor_values <- rbeta(n_off_diag, shape1 = beta_shape1, shape2 = beta_shape2)
    cor_values <- 2 * (cor_values - 0.5)  # Map Beta values to [-1, 1]
  } else if (distribution_type == "u_shaped") {
    cor_values <- rbeta(n_off_diag, shape1 = ushape1, shape2 = ushape2)
    cor_values <- 2 * (cor_values - 0.5)  # Map Beta values to [-1, 1]
  } else {
    stop("Invalid distribution type specified.")
  }
  
  if (distribution_type == "normal") {
    cor_values <- scale(cor_values, center = mean(cor_values), scale = sd(cor_values))
    cor_values <- cor_values * sd + mean
  }
  
  cor_values[cor_values < -1] <- -1
  cor_values[cor_values > 1] <- 1
  
  corMat <- matrix(1, nrow = S, ncol = S)
  corMat[upper.tri(corMat)] <- cor_values
  corMat <- corMat * t(corMat)
  
  return(corMat)
}



#############################################


S <- 20

cor_matrix <- createCorMat(S, distribution_type = "normal", mean = 0, sd = 0.3)

#cor_matrix <- createCorMat(S, distribution_type = "beta", beta_shape1 = 2, beta_shape2 = 5)

#cor_matrix <- createCorMat(S, distribution_type = "u_shaped", ushape1 = 0.5, ushape2 = 0.5)



cor_matrix <- as.data.frame(cor_matrix)
cor_matrix$Var1 = paste0("V", rownames(cor_matrix))
cor_matrix = cor_matrix %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") #|> na.omit()
cor_matrix$Var1 = factor(cor_matrix$Var1, levels = unique(cor_matrix$Var1))
cor_matrix$Var2 = factor(cor_matrix$Var2, levels = rev(unique(cor_matrix$Var2)))


# PLOT CORRELATION MATRIX

ggplot(data = cor_matrix,
       aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-1,1)*max(abs(cor_matrix$value))) +
  ggtitle("Covariance Matrix")



# Pivot the long-format dataframe back to wide format
wide_cor_matrix <- cor_matrix %>%
  pivot_wider(names_from = Var2, values_from = value) %>%
  select(-Var1)

# Convert back to a matrix
recovered_matrix <- as.matrix(wide_cor_matrix)

# COMPUTE COVARIANCE MATRIX

sd_X = rep(0.9, 20)

covMat <- as.matrix(Matrix::nearPD(sd_X %*% t(sd_X) * recovered_matrix)$mat)


# PLOT COVARIANCE MATRIX

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