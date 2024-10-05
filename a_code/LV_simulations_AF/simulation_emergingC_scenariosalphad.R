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

# Function to simulate dynamics with discrete perturbations and track delta_r values
simulate_with_discrete_perturbations <- function(S, timesteps, cov_matrix, perturb_interval = 20, perturb_scale = 1) {
  A <- sim_quantitative_network("predator-prey", S = S, C = 0.2, aij_params = c(0, 0.5))
  init_biomass <- runif(S, min = 1, max = 10)
  r <- -A %*% init_biomass
  
  r_new <- r
  biomass <- init_biomass
  
  results <- data.frame(time = 1, t(biomass))  
  colnames(results) <- c("time", paste0("species_", 1:S))
  
  delta_r_matrix <- matrix(0, nrow = timesteps, ncol = S)  
  
  for (t in seq(perturb_interval, timesteps, by = perturb_interval)) {
    delta_r <- MASS::mvrnorm(n = 1, mu = rep(0, S), Sigma = cov_matrix) * perturb_scale
    r_new <- r_new + delta_r
    delta_r_matrix[t, ] <- delta_r
    
    fw.model_block <- function(t, B, params) {
      with(as.list(c(B, params)), {
        B[B < 10^-8] <- 0  
        dBdt <- (params$r_new + params$A %*% B) * B  
        list(dBdt)
      })
    }
    
    out <- deSolve::ode(
      y = biomass, 
      times = seq(t - perturb_interval + 1, t), 
      func = fw.model_block, 
      parms = list(A = A, r_new = r_new),
      method = "ode45"
    )
    
    biomass <- as.numeric(out[nrow(out), -1])  
    time_block_results <- data.frame(time = t, t(biomass))  
    colnames(time_block_results) <- colnames(results)  
    
    results <- rbind(results, time_block_results)  
  }
  
  return(list(results = results, delta_r_matrix = delta_r_matrix))
}

# Function to generate covariance matrix based on alpha_d
generate_cov_matrix <- function(S, alpha_d) {
  cor_matrix <- clusterGeneration::rcorrmatrix(d = S, alphad = alpha_d)
  sd_X <- rep(1, S)
  cov_matrix <- diag(sd_X) %*% cor_matrix %*% diag(sd_X)
  return(cov_matrix)
}

# Simulation parameters
nsim <- 50
S <- 8
timesteps <- 1000

# List to store results for both alpha_d scenarios
scenarios <- list()

# Define the two alpha_d values
alpha_d_values <- c(0.001, 40)

# Run simulations for both alpha_d scenarios
for (alpha_d in alpha_d_values) {
  print(paste("Running simulations for alpha_d =", alpha_d))
  
  all_sim_results <- list()
  cov_values <- data.frame()  
  variance_covariance <- data.frame()  
  sd_covariance <- data.frame()
  mean_covariance <- data.frame()  
  
  for (sim in 1:nsim) {
    print(paste("sim", sim, "of", nsim))
    
    cov_matrix <- generate_cov_matrix(S, alpha_d)
    
    sim_result <- simulate_with_discrete_perturbations(S = S, timesteps = timesteps, cov_matrix = cov_matrix, perturb_scale = 1.2)
    result <- sim_result$results
    delta_r_matrix <- sim_result$delta_r_matrix
    
    result$sim_id <- sim
    result$alpha_d <- alpha_d
    
    all_sim_results[[sim]] <- result
    
    delta_r_cov_matrix <- cov(delta_r_matrix, use = "pairwise.complete.obs")
    
    cov_vals <- delta_r_cov_matrix[upper.tri(delta_r_cov_matrix)]
    cov_values <- rbind(cov_values, data.frame(sim_id = sim, covariance = cov_vals, alpha_d = alpha_d))
    
    variance_covariance <- rbind(variance_covariance, data.frame(sim_id = sim, var_cov = var(cov_vals), alpha_d = alpha_d))
    sd_covariance <- rbind(sd_covariance, data.frame(sim_id = sim, sd_cov = sd(cov_vals), alpha_d = alpha_d))
    mean_covariance <- rbind(mean_covariance, data.frame(sim_id = sim, mean_cov = mean(cov_vals), alpha_d = alpha_d))
  }
  
  combined_results <- do.call(rbind, all_sim_results)
  
  combined_long <- combined_results %>%
    pivot_longer(cols = starts_with("species_"), names_to = "species", values_to = "biomass")
  
  biomass_change_stats <- combined_long %>%
    group_by(species, sim_id) %>%
    arrange(time) %>%
    mutate(biomass_change = biomass - lag(biomass)) %>%
    filter(!is.na(biomass_change)) %>%
    summarise(mean_change = mean(biomass_change, na.rm = TRUE),
              variance_change = var(biomass_change, na.rm = TRUE)) %>%
    ungroup()
  
  biomass_var_cov <- biomass_change_stats %>%
    group_by(sim_id) %>%
    summarise(variance_abundance_change = mean(variance_change, na.rm = TRUE)) %>%
    inner_join(variance_covariance, by = "sim_id") %>%
    inner_join(sd_covariance, by = "sim_id") %>%
    inner_join(mean_covariance, by = "sim_id")
  
  scenarios[[as.character(alpha_d)]] <- list(
    combined_results = combined_results,
    combined_long = combined_long,
    cov_values = cov_values,
    biomass_var_cov = biomass_var_cov
  )
}

# Combine data from both alpha_d scenarios for comparisons
cov_values_combined <- do.call(rbind, lapply(scenarios, function(x) x$cov_values))
biomass_var_cov_combined <- do.call(rbind, lapply(scenarios, function(x) x$biomass_var_cov))

# Plot 1: Distribution of covariance values across alpha_d scenarios
p_cov <- ggplot(cov_values_combined, aes(x = covariance, fill = as.factor(alpha_d))) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~ alpha_d, scales = "free_x") +
  theme_minimal() +
  ggtitle("Distribution of Covariance Values for Different alpha_d Values") +
  labs(x = "Covariance", y = "Frequency")

# Plot 2: Variance of covariance vs. variance in abundance change across scenarios
p_var_cov <- ggplot(biomass_var_cov_combined, aes(x = var_cov, y = variance_abundance_change, color = as.factor(alpha_d))) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Variance of Covariance vs Variance in Abundance Change", x = "Variance of Covariance", y = "Variance in Abundance Change") +
  theme_minimal()

# Plot 3: Mean covariance vs. variance in abundance change across scenarios
p_mean_cov <- ggplot(biomass_var_cov_combined, aes(x = mean_cov, y = variance_abundance_change, color = as.factor(alpha_d))) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Mean Covariance vs Variance in Abundance Change", x = "Mean Covariance", y = "Variance in Abundance Change") +
  theme_minimal()

# Display plots
print(p_cov)
print(p_var_cov)
print(p_mean_cov)
