#################################
##### normal ######################


S = 20

cor_matrix <- createCorMat(S = S, 
                           distribution_type = "normal", 
                           mean = 0, sd = 0.1)

# Adjust C matrix
cor_matrix <- as.data.frame(cor_matrix)
cor_matrix$Var1 <- paste0("V", rownames(cor_matrix))
cor_matrix <- cor_matrix %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")
cor_matrix$Var1 <- factor(cor_matrix$Var1, levels = unique(cor_matrix$Var1))
cor_matrix$Var2 <- factor(cor_matrix$Var2, levels = rev(unique(cor_matrix$Var2)))

# plot the matrix

ggplot(data = cor_matrix,
       aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-1,1)*max(abs(cor_matrix$value))) +
  ggtitle("Correlation Matrix")

# Pivot the long-format dataframe back to wide format
wide_cor_matrix <- cor_matrix %>%
  pivot_wider(names_from = Var2, values_from = value) %>%
  select(-Var1)

# Convert back to a matrix
recovered_matrix <- as.matrix(wide_cor_matrix)

sd_X = rep(0.9, 20)

covMat <- as.matrix(Matrix::nearPD(sd_X %*% t(sd_X) * recovered_matrix)$mat)


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



#################################
##### BETA ######################

S = 20

cor_matrix <- createCorMat(S = S, 
                           distribution_type = "beta",
                           beta_shape1 = 0.05,
                           beta_shape2 = 0.05)

# Adjust C matrix
cor_matrix <- as.data.frame(cor_matrix)
cor_matrix$Var1 <- paste0("V", rownames(cor_matrix))
cor_matrix <- cor_matrix %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value")
cor_matrix$Var1 <- factor(cor_matrix$Var1, levels = unique(cor_matrix$Var1))
cor_matrix$Var2 <- factor(cor_matrix$Var2, levels = rev(unique(cor_matrix$Var2)))

# plot the matrix

ggplot(data = cor_matrix,
       aes(x = Var1, y = Var2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-1,1)*max(abs(cor_matrix$value))) +
  ggtitle("Correlation Matrix")

# Pivot the long-format dataframe back to wide format
wide_cor_matrix <- cor_matrix %>%
  pivot_wider(names_from = Var2, values_from = value) %>%
  select(-Var1)

# Convert back to a matrix
recovered_matrix <- as.matrix(wide_cor_matrix)

sd_X = rep(0.9, 20)

covMat <- as.matrix(Matrix::nearPD(sd_X %*% t(sd_X) * recovered_matrix)$mat)


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