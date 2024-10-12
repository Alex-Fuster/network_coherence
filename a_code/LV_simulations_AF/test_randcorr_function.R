pak::pak("randcorr")
library(randcorr)

cor_matrix <- randcorr(10)

# Extract the off-diagonal values
off_diag_values <- cor_matrix[upper.tri(cor_matrix)]

# Convert to a data frame for plotting
df <- data.frame(values = off_diag_values)

# Plot the distribution of values with a density overlay
ggplot(df, aes(x = values)) +
  geom_histogram(aes(y = ..density..), bins = 20, fill = "skyblue", alpha = 0.7) +
  geom_density(color = "red", size = 1) +
  labs(
    title = "Distribution of Correlation Matrix Values (Off-diagonal)",
    x = "Correlation Values",
    y = "Density"
  )+
  xlim(-1, 1) +
  theme_minimal()


MASS::mvrnorm(n = 1, mu = rep(0, 10), Sigma = cor_matrix)
