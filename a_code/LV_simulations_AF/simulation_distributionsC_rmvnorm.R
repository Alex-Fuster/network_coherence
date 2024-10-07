# Load necessary libraries
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

# Generate symmetric correlation matrix for creating covariance matrices
createSymmetricCorMat <- function(S, distribution_type = "normal", mean = 0, sd = 1, beta_shape1 = NULL, beta_shape2 = NULL, ushape1 = NULL, ushape2 = NULL) {
  n_off_diag <- S * (S - 1) / 2
  
  if (distribution_type == "normal") {
    cor_values <- rnorm(n_off_diag, mean = mean, sd = sd)
  } else if (distribution_type == "beta") {
    cor_values <- rbeta(n_off_diag, shape1 = beta_shape1, shape2 = beta_shape2)
    cor_values <- 2 * (cor_values - 0.5)
  } else if (distribution_type == "u_shaped") {
    cor_values <- rbeta(n_off_diag, shape1 = ushape1, shape2 = ushape2)
    cor_values <- 2 * (cor_values - 0.5)
  } else {
    stop("Invalid distribution type specified.")
  }
  
  cor_values[cor_values < -1] <- -1
  cor_values[cor_values > 1] <- 1
  
  corMat <- diag(1, nrow = S, ncol = S)
  corMat[upper.tri(corMat)] <- cor_values
  corMat[lower.tri(corMat)] <- t(corMat)[lower.tri(corMat)]
  
  return(corMat)
}

# Define function to simulate biomass dynamics
simulate_dynamic_r <- function(S, timesteps, cov_matrix, interval = 20, perturb_scale = 1.5) {
  A <- sim_quantitative_network("predator-prey", S = S, C = 0.2, aij_params = c(0, 0.5))
  init_biomass <- runif(S, min = 1, max = 10)
  r <- as.numeric(-A %*% init_biomass)  # Initialize r based on initial biomass
  biomass <- init_biomass
  results <- data.frame(time = 1, t(biomass))
  colnames(results) <- c("time", paste0("species_", 1:S))
  
  for (t in seq(interval, timesteps, by = interval)) {
    r <- as.numeric(mvtnorm::rmvnorm(n = 1, mean = r, sigma = cov_matrix) * perturb_scale)
    
    fw.model <- function(t, B, params) {
      with(as.list(c(B, params)), {
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

# Parameters
S <- 20  # Number of species
timesteps <- 500
nsim <- 20
scenarios <- list(
  list("label" = "normal weak", "distribution_type" = "normal", "mean" = 0, "sd" = 0.1),
  list("label" = "normal strong", "distribution_type" = "normal", "mean" = 0, "sd" = 1.2),
  list("label" = "beta -", "distribution_type" = "beta", "beta_shape1" = 2, "beta_shape2" = 5),
  list("label" = "beta +", "distribution_type" = "beta", "beta_shape1" = 5, "beta_shape2" = 2),
  list("label" = "u-shaped", "distribution_type" = "u_shaped", "ushape1" = 0.5, "ushape2" = 0.5)
)

# Run simulations for each scenario
all_results <- data.frame()
cov_values <- data.frame()

for (scenario in scenarios) {
  print(paste("scenario", scenario))
  corMat <- do.call(createSymmetricCorMat, c(list(S = S), scenario[-1]))
  sd_X <- rep(1, S)
  cov_matrix <- diag(sd_X) %*% corMat %*% diag(sd_X)
  
  for (sim in 1:nsim) {
    
    print(paste("sim", sim, "of", nsim))
    results <- simulate_dynamic_r(S, timesteps, cov_matrix, interval = 20)
    
    results_long <- results %>%
      pivot_longer(cols = starts_with("species_"), names_to = "species", values_to = "biomass") %>%
      group_by(species) %>%
      arrange(time) %>%
      mutate(biomass_change = biomass - lag(biomass)) %>%
      filter(!is.na(biomass_change)) %>%
      summarise(mean_change = mean(biomass_change, na.rm = TRUE), sd_change = sd(biomass_change, na.rm = TRUE)) %>%
      mutate(label = scenario$label, sim_id = sim)
    all_results <- rbind(all_results, results_long)
  }
  
  cov_vals <- cov_matrix[upper.tri(cov_matrix)]
  cov_values <- rbind(cov_values, data.frame(covariance = cov_vals, label = scenario$label))
}

# Plotting
p1 <- ggplot(all_results, aes(x = label, y = sd_change, fill = label)) +
  geom_violin() +
  labs(title = "SD of Biomass Change Across Scenarios", x = "Scenario", y = "SD of Biomass Change") +
  theme_minimal()

p2 <- ggplot(all_results, aes(x = label, y = mean_change, fill = label)) +
  geom_boxplot() +
  labs(title = "Mean of Biomass Change Across Scenarios", x = "Scenario", y = "Mean of Biomass Change") +
  theme_minimal()

p3 <- ggplot(cov_values, aes(x = covariance, fill = label)) +
  geom_histogram(bins = 30, alpha = 0.7, position = "identity") +
  facet_wrap(~ label, scales = "free") +
  labs(title = "Covariance Distributions for Different Scenarios", x = "Covariance", y = "Frequency") +
  theme_minimal()

# Display plots
print(p1)
print(p2)
print(p3)
