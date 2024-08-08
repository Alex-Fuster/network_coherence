###########################
# Covariance matrices
###########################


S <- 25
set.seed(120)

# Function to generate mixed covariance matrices
generate_mixed_cov_matrix <- function(S, range_neg, range_pos) {
  mat <- matrix(runif(S^2, range_neg[1], range_pos[2]), nrow = S)
  diag(mat) <- runif(S, range_neg[1], range_pos[2])
  return(mat)
}

cov_matrices <- list(
  "Strong Positive" = matrix(runif(S^2, 0.6, 1), nrow = S),
  "Weak Positive" = matrix(runif(S^2, 0.1, 0.4), nrow = S),
  "Strong Negative" = matrix(runif(S^2, -1, -0.6), nrow = S),
  "Weak Negative" = matrix(runif(S^2, -0.4, -0.1), nrow = S),
  "Mixed Weak" = generate_mixed_cov_matrix(S, c(-0.4, 0), c(0, 0.4)),
  "Mixed Strong" = generate_mixed_cov_matrix(S, c(-1, 0), c(0, 1)),
  "Null" = matrix(0, nrow = S, ncol = S)
)



prepare_plot_data <- function(cov_matrices) {
  plot_data <- lapply(names(cov_matrices), function(name) {
    matrix <- cov_matrices[[name]]
    df <- melt(matrix)
    df$Scenario <- name
    return(df)
  })
  combined_plot_data <- do.call(rbind, plot_data)
  return(combined_plot_data)
}


plot_data <- prepare_plot_data(cov_matrices)

p <- ggplot(plot_data, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white") +
  scale_fill_distiller(palette = "RdBu", direction = 1, na.value = "white") +
  labs(title = "Covariance Matrices",
       x = NULL,
       y = NULL,
       fill = "Covariance") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ Scenario, ncol = 3, scales = "free")+
  my_theme







###########################
# Simulation
###########################


num_simulations <- 100
set.seed(120)

# Function to simulate quantitative networks
sim_quantitative_network <- function(Net_type, S, C, aij_params, rho = 0) {
  A <- matrix(0, S, S)
  n_pairs <- S * (S - 1) / 2
  B <- runif(n_pairs) <= C
  if (Net_type == "random") {
    A[upper.tri(A)] <- B * rnorm(n_pairs, aij_params[1], aij_params[2])
    A <- t(A)
    A[upper.tri(A)] <- B * rnorm(n_pairs, aij_params[1], aij_params[2])
  } else if (Net_type == "predator-prey") {
    aij <- -abs(rnorm(n_pairs, aij_params[1], aij_params[2]))
    A[upper.tri(A)] <- B * aij
    A <- t(A)
    aij <- abs(rnorm(n_pairs, aij_params[1], aij_params[2]))
    A[upper.tri(A)] <- B * aij
  } else if (Net_type == "competition") {
    aij <- -abs(rnorm(n_pairs * 2, aij_params[1], aij_params[2]))
    A[upper.tri(A)] <- B * aij[1:n_pairs]
    A <- t(A)
    A[upper.tri(A)] <- B * aij[(n_pairs + 1):length(aij)]
  } else if (Net_type == "mutualistic") {
    aij <- abs(rnorm(n_pairs * 2, aij_params[1], aij_params[2]))
    A[upper.tri(A)] <- B * aij[1:n_pairs]
    A <- t(A)
    A[upper.tri(A)] <- B * aij[(n_pairs + 1):length(aij)]
  } else {
    stop("Incorrect network type")
  }
  diag(A) <- -(max(Re(eigen(A)$values)) + runif(S, 0.1))
  while (max(Re(eigen(A)$values)) > 0) {
    diag(A) <- -(max(Re(eigen(A)$values)) + runif(S, 0.1))
  }
  return(A)
}



calculate_metrics <- function(B, cov_delta_r) {
  n <- nrow(B)
  
  # Variance of Delta X
  var_delta_X <- numeric(n)
  for (i in 1:n) {
    var_sum <- 0
    for (k in 1:n) {
      for (l in 1:n) {
        var_sum <- var_sum + B[i, k] * B[i, l] * cov_delta_r[k, l]
      }
    }
    var_delta_X[i] <- var_sum
    
  }
  
  # Mean of the sum of each row of B
  mean_sum_B <- colMeans(B)
  
  # Mean of delta_r (assuming it is an average of a random distribution)
  mean_delta_r <- mean(cov_delta_r)
  
  # Covariance between sum of each row of B and delta_r
  cov_sum_B_delta_r <- cov(rowSums(B), cov_delta_r)
  
  # Calculate Abundance Change
  abundance_change <- mean_sum_B * mean_delta_r + cov_sum_B_delta_r
  
  return(list(Variance = var_delta_X, Abundance_Change = abundance_change))
}


################# Simulation

run_simulations_with_metrics <- function(cov_matrix, num_simulations, Net_type, S, C, aij_params) {
  results <- replicate(num_simulations, {
    interaction_matrix <- sim_quantitative_network(Net_type, S, C, aij_params)
    inverse_matrix <- solve(interaction_matrix)
    if (is.null(inverse_matrix)) {
      return(list(Variance = rep(NA, S), Abundance_Change = NA))
    }
    calculate_metrics(inverse_matrix, cov_matrix)
  }, simplify = FALSE)
  return(results)
}

# Generate results for each scenario
Net_type <- "predator-prey"  # Example network type
S <- 25  # Number of species
C <- 0.2  # Connectance
aij_params <- c(0, 0.1)  # Parameters for interaction strengths
num_simulations <- 100  # Number of simulations

simulation_results_with_metrics <- lapply(cov_matrices, function(cov_matrix) {
  run_simulations_with_metrics(cov_matrix, num_simulations, Net_type, S, C, aij_params)
})
############### Plot results


# Prepare the data for plotting
prepare_combined_plot_data <- function(simulation_results) {
  df_list_var <- lapply(names(simulation_results), function(scenario) {
    data <- simulation_results[[scenario]]
    df <- data.frame(
      Scenario = scenario,
      Variance = unlist(lapply(data, function(res) res$Variance)),
      Abundance_Change = unlist(lapply(data, function(res) res$Abundance_Change))
    )
    return(df)
  })
  combined_df_var <- bind_rows(df_list_var)
  combined_df_var <- combined_df_var %>% filter(!is.na(Variance) & !is.na(Abundance_Change))
  
  return(combined_df_var)
}

# Prepare the data
combined_plot_data <- prepare_combined_plot_data(simulation_results_with_metrics)

scenario_order <- c(
  "Strong Negative",
  "Weak Negative",
  "Null",
  "Weak Positive",
  "Strong Positive",
  "Mixed Weak",
  "Mixed Strong"
)

combined_plot_data$Scenario <- factor(combined_plot_data$Scenario, levels = scenario_order)


color_palette <- scales::brewer_pal(palette = "RdBu", direction = -1)(5)
mixed_color <- "violet"

scenario_colors <- setNames(c(color_palette, mixed_color, mixed_color), scenario_order)

# Plot variance in abundance change
p1 <- ggplot(combined_plot_data, aes(x = Scenario, y = Variance, fill = Scenario)) +
  geom_jitter(width = 0.2, alpha = 0.1, shape = 21, size = 2, aes(color = Scenario)) +
  scale_fill_manual(values = scenario_colors) +
  scale_color_manual(values = c(setNames(ifelse(scenario_order == "Null", "black", scenario_colors), scenario_order), "black")) +
  labs(x = " ", y = "Variance in Abundance Change") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  my_theme

# Plot abundance change
p2 <- ggplot(combined_plot_data, aes(x = Scenario, y = Abundance_Change, fill = Scenario)) +
  geom_jitter(width = 0.2, alpha = 0.1, shape = 21, size = 2, aes(color = Scenario)) +
  scale_fill_manual(values = scenario_colors) +
  scale_color_manual(values = c(setNames(ifelse(scenario_order == "Null", "black", scenario_colors), scenario_order), "black")) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(x = " ", y = "Abundance Change") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  my_theme

# Display plots
ggarrange(p1, p2, nrow = 2, ncol = 1, labels = LETTERS[1:2])
