# Load necessary libraries
library(clusterGeneration)
library(deSolve)
library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)

# Function to simulate quantitative interaction network
sim_quantitative_network <- function(Net_type, S, C, aij_params) {
  A <- matrix(0, S, S)
  n_pairs <- S * (S - 1) / 2
  B <- runif(n_pairs) <= C
  if (Net_type == "predator-prey") { 
    aij <- -abs(rnorm(n_pairs, aij_params[1], aij_params[2]))
    A[upper.tri(A)] <- B * aij
    A <- t(A)
    aij <- abs(rnorm(n_pairs, aij_params[1], aij_params[2]))
    A[upper.tri(A)] <- B * aij
  } 
  diag(A) <- -runif(S, min = 0, max = 1)
  while(max(Re(eigen(A)$values)) > 0){
    diag(A) <- -runif(S, min = 0, max = 1)
  }
  return(A)
}

# Function to simulate biomass dynamics with changing r based on the previous timestep
simulate_dynamic_r <- function(S, timesteps, cov_matrix, interval = 20, perturb_scale = 1) {
  A <- sim_quantitative_network("predator-prey", S = S, C = 0.2, aij_params = c(0, 0.5))
  B <- solve(A)
  init_biomass <- runif(S, min = 1, max = 10)
  r <- -A %*% init_biomass  # Initialize r based on initial biomass
  biomass <- init_biomass
  results <- data.frame(time = 1, t(biomass))
  colnames(results) <- c("time", paste0("species_", 1:S))
  
  # Iterate over time and update r using the previous timestep's r values as the mean
  for (t in seq(interval, timesteps, by = interval)) {
    r <- MASS::mvrnorm(n = 1, mu = r, Sigma = cov_matrix) * perturb_scale
    
    fw.model <- function(t, B, params) {
      with(as.list(c(B, params)), {
        B[B < 10^-8] <- 0
        dBdt <- (params$r + params$A %*% B) * B
        list(dBdt)
      })
    }
    
    out <- deSolve::ode(
      y = biomass, 
      times = seq(t - interval + 1, t), 
      func = fw.model, 
      parms = list(A = A, r = r),
      method = "ode45"
    )
    
    biomass <- as.numeric(out[nrow(out), -1])
    time_block_results <- data.frame(time = t, t(biomass))
    colnames(time_block_results) <- colnames(results)
    results <- rbind(results, time_block_results)
  }
  
  return(results)
}

# Generate covariance matrix
generate_cov_matrix <- function(S, alpha_d) {
  cor_matrix <- clusterGeneration::rcorrmatrix(d = S, alphad = alpha_d)
  sd_X <- rep(ifelse(alpha_d < 0.1, 1.3, 0.6), S)  
  cov_matrix <- diag(sd_X) %*% cor_matrix %*% diag(sd_X)
  return(cov_matrix)
}

# Simulation parameters
S <- 20  # Number of species
timesteps <- 1000
alpha_d_values <- c(0.001, 40)
nsim <- 200  # Number of simulations per alpha_d value

# Initialize lists to store results
variance_abundance <- data.frame()
cov_values <- data.frame()
biomass_data <- data.frame()

# Run multiple simulations for each alpha_d
for (alpha_d in alpha_d_values) {
  print(paste("Running simulations for alpha_d =", alpha_d))
  
  # Generate covariance matrix
  cov_matrix <- generate_cov_matrix(S, alpha_d)
  
  # Store covariance values from the generated covariance matrix
  cov_vals <- cov_matrix[upper.tri(cov_matrix)]
  cov_values <- rbind(cov_values, data.frame(covariance = cov_vals, alpha_d = as.factor(alpha_d)))
  
  for (sim in 1:nsim) {
    print(paste("sim", sim, "of", nsim, "- alpha_d =", alpha_d))
    # Run the simulation with dynamic r
    results <- simulate_dynamic_r(S = S, timesteps = timesteps, cov_matrix = cov_matrix, interval = 20)
    
    # Convert results to long format for biomass dynamics plot
    results_long <- results %>%
      pivot_longer(cols = starts_with("species_"), names_to = "species", values_to = "biomass") %>%
      mutate(sim_id = sim, alpha_d = as.factor(alpha_d))
    
    biomass_data <- rbind(biomass_data, results_long)
    
    # Summarize biomass changes for boxplot analysis
    summary_stats <- results_long %>%
      group_by(species) %>%
      arrange(time) %>%
      mutate(biomass_change = biomass - lag(biomass)) %>%  # Calculate biomass change between timesteps
      filter(!is.na(biomass_change)) %>%
      summarise(mean_change = mean(biomass_change, na.rm = TRUE),
                sd_change = sd(biomass_change, na.rm = TRUE)) %>%
      mutate(alpha_d = as.factor(alpha_d), sim_id = sim)
    
    variance_abundance <- rbind(variance_abundance, summary_stats)
  }
}

# Plot: Boxplot of mean change in biomass comparing two covariance structures
p1 <- ggplot(variance_abundance, aes(x = alpha_d, y = mean_change, fill = alpha_d)) +
  geom_boxplot() +
  labs(title = "Mean Biomass Change under Different Covariance Structures", x = "alpha_d", y = "Mean Change in Biomass") +
  theme_minimal()

# Plot: Boxplot of SD change in biomass comparing two covariance structures
p2 <- ggplot(variance_abundance, aes(x = alpha_d, y = sd_change, fill = alpha_d)) +
  geom_boxplot() +
  labs(title = "SD of Biomass Change under Different Covariance Structures", x = "alpha_d", y = "SD Change in Biomass") +
  theme_minimal()

# Plot: Distribution of covariance values for each alpha_d
p3 <- ggplot(cov_values, aes(x = covariance, fill = alpha_d)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~ alpha_d, scales = "free") +
  labs(title = "Distribution of Covariance Values for Different alpha_d", x = "Covariance", y = "Frequency") +
  theme_minimal()

# Plot: Biomass dynamics for selected simulations
p4 <- ggplot(biomass_data %>% filter(sim_id %in% c(1,2,3)), aes(x = time, y = biomass, color = species)) +
  geom_line() +
  facet_wrap(~ sim_id + alpha_d, ncol = 3) +
  labs(title = "Biomass Dynamics for Selected Simulations", x = "Time", y = "Biomass") +
  theme_minimal() +
  theme(legend.position = "none")  # Remove legend for clarity

# Display the plots
print(p1)
print(p2)
print(p3)
print(p4)


p2_violin <- ggplot(variance_abundance, aes(x = alpha_d, y = sd_change, fill = alpha_d)) +
  geom_violin(trim = FALSE) +
  labs(title = "SD of Biomass Change under Different Covariance Structures", x = "alpha_d", y = "SD Change in Biomass") +
  theme_minimal()

# Histogram: Distribution of SD change in biomass for each alpha_d
p2_histogram <- ggplot(variance_abundance, aes(x = sd_change, fill = alpha_d)) +
  geom_histogram(alpha = 0.6, position = "identity", bins = 30) +
  scale_fill_manual(values = c("salmon", "steelblue")) +  # Customize colors as needed
  labs(title = "Distribution of SD Change in Biomass by alpha_d", x = "SD Change in Biomass", y = "Frequency") +
  theme_minimal()

# Display the plots
print(p2_violin)
print(p2_histogram)
