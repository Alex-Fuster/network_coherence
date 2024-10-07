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

fw.model <- function(t, B, params) {
  with(as.list(c(B, params)), {
    B[B < 10^-8] <- 0
    dBdt <- (params$r[t+1,] + params$A %*% B) * B
    list(dBdt)
  })
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
  
  r <- MASS::mvrnorm(n = timesteps, mu = r, Sigma = cov_matrix)*perturb_scale
  out <- deSolve::ode(
    y = biomass, 
    times = seq(0, 100, length.out = timesteps), 
    func = fw.model, 
    parms = list(A = A, r = r)
  )
  
  return(as.data.frame(out))
}

# Generate covariance matrix
generate_cov_matrix <- function(S, alpha_d) {
  cor_matrix <- clusterGeneration::rcorrmatrix(d = S, alphad = alpha_d)
  sd_X <- rep(1, S)  
  cov_matrix <- diag(sd_X) %*% cor_matrix %*% diag(sd_X)
  return(cov_matrix)
}

generate_cov_matrix2 <- function(S, corr, sd_x) {
  cor_matrix <- matrix(1, nrow = S, ncol = S)
  cor_matrix[upper.tri(cor_matrix)] <- cor_matrix[lower.tri(cor_matrix)] <- corr
  sd_X <- rep(1, S)  
  cov_matrix <- diag(sd_X) %*% cor_matrix %*% diag(sd_X)
  return(cov_matrix)
}

# Simulation parameters
S <- 20  # Number of species
timesteps <- 1000
corr_values <- seq(0, 0.9, by = 0.1)
nsim <- 10  # Number of simulations per alpha_d value

# Initialize lists to store results
variance_abundance <- data.frame()
cov_values <- data.frame()
biomass_data <- data.frame()

# Run multiple simulations for each correlations
for (corr in corr_values) {
  print(paste("Running simulations for correlation =", corr))
  
  # Generate covariance matrix
  cov_matrix <- generate_cov_matrix2(S, corr = corr, sd_x = 1)
  
  # Store covariance values from the generated covariance matrix
 # cov_vals <- cov_matrix[upper.tri(cov_matrix)]
#  cov_values <- rbind(cov_values, data.frame(covariance = cov_vals, alpha_d = as.factor(alpha_d)))
  
  for (sim in 1:nsim) {
    print(paste("sim", sim, "of", nsim, "- correlation =", corr))
    # Run the simulation with dynamic r
    results <- simulate_dynamic_r(S = S, timesteps = timesteps, cov_matrix = cov_matrix, interval = 20)
    colnames(results) <- c("time", paste0("species", c(1:S)))
    
    # Convert results to long format for biomass dynamics plot
    results_long <- results %>%
      pivot_longer(cols = starts_with("species"), names_to = "species", values_to = "biomass") %>%
      mutate(sim_id = sim, correlation = as.factor(corr))
    
    biomass_data <- rbind(biomass_data, results_long)
  }
}

biomass_sd <- biomass_data %>%
  group_by(sim_id, species, correlation) %>%
  summarise(biomass_sd = sd(biomass)) %>%
  group_by(species, correlation) %>%
  summarise(mean_varX = mean(biomass_sd))


# Plot: Boxplot of mean change in biomass comparing two covariance structures
p1 <- ggplot(biomass_sd, aes(x = correlation, y = mean_varX, fill = correlation)) +
  geom_boxplot() +
  labs(title = "Mean Biomass SD under Different Covariance Structures", x = "corr", y = "Mean Var(X_i)", fill = "corr") +
  theme_minimal()

print(p1)