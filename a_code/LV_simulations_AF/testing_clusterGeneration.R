pak::pak("clusterGeneration")

library(clusterGeneration)
library(ggplot2)
library(reshape2)

# Define the dimension of the matrix
d <- 40  # Number of species/variables

# Define a range of alphad values to explore
alphad_values <- c(0.1, 0.5, 1, 2, 5, 10)

# Initialize a list to store correlation matrices
corr_matrices <- list()

# Generate correlation matrices for each alphad value
set.seed(123)  # For reproducibility
for (alpha in alphad_values) {
  corr_matrix <- rcorrmatrix(d = d, alphad = alpha)
  corr_matrices[[paste0("alphad_", alpha)]] <- corr_matrix
}


# Function to plot a correlation matrix
plot_corr_matrix <- function(corr_matrix, title) {
  # Convert the matrix to a long format for ggplot
  corr_melt <- melt(corr_matrix)
  colnames(corr_melt) <- c("Var1", "Var2", "value")
  
  # Create the heatmap
  ggplot(corr_melt, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_distiller(palette = "RdBu", limits = c(-1, 1), name = "Correlation") +
    ggtitle(title) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title = element_blank())
}

# Plot the correlation matrices
plot_list <- list()
for (alpha in alphad_values) {
  corr_matrix <- corr_matrices[[paste0("alphad_", alpha)]]
  title <- paste("Correlation Matrix (alphad =", alpha, ")")
  p <- plot_corr_matrix(corr_matrix, title)
  plot_list[[paste0("alphad_", alpha)]] <- p
}

# Display the plots
library(gridExtra)
do.call(grid.arrange, c(plot_list, ncol = 2))


# Function to plot histogram of correlation coefficients
plot_corr_histogram <- function(corr_matrix, title) {
  # Extract the upper triangle of the correlation matrix (excluding the diagonal)
  corr_values <- corr_matrix[upper.tri(corr_matrix)]
  
  # Create the histogram
  ggplot(data.frame(correlation = corr_values), aes(x = correlation)) +
    geom_histogram(bins = 20, color = "black", fill = "skyblue") +
    xlim(-1, 1) +
    ggtitle(title) +
    theme_minimal()
}

# Plot histograms
hist_list <- list()
for (alpha in alphad_values) {
  corr_matrix <- corr_matrices[[paste0("alphad_", alpha)]]
  title <- paste("Histogram of Correlations (alphad =", alpha, ")")
  p <- plot_corr_histogram(corr_matrix, title)
  hist_list[[paste0("alphad_", alpha)]] <- p
}

# Display the histograms
do.call(grid.arrange, c(hist_list, ncol = 2))
