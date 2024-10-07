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


# Function to propagate covariances into biomass dynamics
compute_biomass_variance <- function(B, delta_r_cov_matrix) {
  S <- nrow(B)
  biomass_variance <- numeric(S)
  
  for (i in 1:S) {
    for (k in 1:S) {
      for (l in 1:S) {
        biomass_variance[i] <- biomass_variance[i] + B[i, k] * B[i, l] * delta_r_cov_matrix[k, l]
      }
    }
  }
  
  return(biomass_variance)
}

compute_biomass_covariance <- function(B, delta_r_cov_matrix) {
  S <- nrow(B)
  biomass_cov_matrix <- matrix(0, S, S)
  
  for (i in 1:S) {
    for (j in 1:S) {
      for (k in 1:S) {
        for (l in 1:S) {
          biomass_cov_matrix[i, j] <- biomass_cov_matrix[i, j] + B[i, k] * B[j, l] * delta_r_cov_matrix[k, l]
        }
      }
    }
  }
  
  return(biomass_cov_matrix)
}

# Function to simulate dynamics with discrete perturbations and track delta_r values
simulate_with_discrete_perturbations <- function(S, timesteps, cov_matrix, perturb_interval = 20, perturb_scale = 2) {
  A <- sim_quantitative_network("predator-prey", S = S, C = 0.2, aij_params = c(0, 0.5))
  B <- solve(A)
  init_biomass <- runif(S, min = 1, max = 10)
  r <- -A %*% init_biomass
  r_new <- r
  biomass <- init_biomass
  results <- data.frame(time = 1, t(biomass))
  colnames(results) <- c("time", paste0("species_", 1:S))
  delta_r_matrix <- matrix(0, nrow = timesteps, ncol = S)
  
  for (t in seq(perturb_interval, timesteps, by = perturb_interval)) {
    # Apply larger perturbations
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
  
  return(list(results = results, delta_r_matrix = delta_r_matrix, B = B))
}

# Function to generate covariance matrix based on emergent dynamics
generate_cov_matrix <- function(S, alpha_d) {
  cor_matrix <- clusterGeneration::rcorrmatrix(d = S, alphad = alpha_d)
  sd_X <- rep(1, S) * 1.5  # Increase standard deviation to scale up covariances
  cov_matrix <- diag(sd_X) %*% cor_matrix %*% diag(sd_X)
  return(cov_matrix)
}

# Simulation parameters
nsim <- 50
S <- 15  # Number of species (you can change back to 50)
timesteps <- 1000

alpha_d_values <- c(0.001, 40)
results_by_alpha <- list()

for (alpha_d in alpha_d_values) {
  print(paste("Running simulations for alpha_d =", alpha_d))
  
  all_sim_results <- list()
  all_cov_values <- data.frame()
  all_variance_cov <- data.frame()
  
  for (sim in 1:nsim) {
    print(paste("Sim", sim, "of", nsim))
    
    # Generate covariance matrix and simulate
    cov_matrix <- generate_cov_matrix(S, alpha_d)
    sim_result <- simulate_with_discrete_perturbations(S = S, timesteps = timesteps, cov_matrix = cov_matrix, perturb_scale = 4)  # Increased perturb_scale
    delta_r_matrix <- sim_result$delta_r_matrix
    B <- sim_result$B
    
    # Compute the covariance matrix of delta_r vectors
    delta_r_cov_matrix <- cov(delta_r_matrix, use = "pairwise.complete.obs")
    
    # Calculate variance and covariance in biomass dynamics
    biomass_variance <- compute_biomass_variance(B, delta_r_cov_matrix)
    biomass_covariance <- compute_biomass_covariance(B, delta_r_cov_matrix)
    
    # Store covariance values
    cov_vals <- delta_r_cov_matrix[upper.tri(delta_r_cov_matrix)]
    all_cov_values <- rbind(all_cov_values, data.frame(sim_id = sim, covariance = cov_vals))
    
    # Store biomass variance
    all_variance_cov <- rbind(all_variance_cov, data.frame(sim_id = sim, 
                                                           variance_abundance_change = mean(biomass_variance, na.rm = TRUE),
                                                           mean_cov = mean(cov_vals),
                                                           sd_cov = sd(cov_vals)))
  }
  
  results_by_alpha[[as.character(alpha_d)]] <- list(cov_values = all_cov_values, variance_cov = all_variance_cov)
}

# Combine results for plotting
combined_variance <- bind_rows(
  results_by_alpha[["0.001"]]$variance_cov %>% mutate(alpha_d = "0.001"),
  results_by_alpha[["40"]]$variance_cov %>% mutate(alpha_d = "40")
)

# Plot variance of abundance change vs. mean covariance for each alpha_d
p1 <- ggplot(combined_variance, aes(x = mean_cov, y = variance_abundance_change, color = alpha_d)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Variance of Abundance Change vs Mean Covariance", x = "Mean Covariance", y = "Variance of Abundance Change") +
  theme_minimal()

p1.2 <- ggplot(combined_variance, aes(x = sd_cov, y = variance_abundance_change, color = alpha_d)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(title = "Variance of Abundance Change vs Mean Covariance", x = "SD Covariance", y = "Variance of Abundance Change") +
  theme_minimal()

# Plot 2: Distribution of covariance values for each alpha_d
combined_cov_values <- bind_rows(
  results_by_alpha[["0.001"]]$cov_values %>% mutate(alpha_d = "0.001"),
  results_by_alpha[["40"]]$cov_values %>% mutate(alpha_d = "40")
)

p2 <- ggplot(combined_cov_values, aes(x = covariance, fill = alpha_d)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~ alpha_d, scales = "free") +
  theme_minimal() +
  labs(title = "Distribution of Covariance Values for Different alpha_d", x = "Covariance", y = "Frequency")

# Display plots
print(p1)
print(p1.2)
print(p2)
