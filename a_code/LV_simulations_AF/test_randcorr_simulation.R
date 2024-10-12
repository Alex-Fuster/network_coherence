# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(deSolve)
library(mvtnorm)
library(randcorr)

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(deSolve)
library(mvtnorm)
library(randcorr)

# Define functions for dynamics and network generation
fw.model <- function(t, B, params) {
  with(as.list(c(B, params)), {
    B[B < 10^-8] <- 0 
    B[B < 0] <- 0 
    dBdt <- (params$r + params$A %*% B) * B
    dBdt[B > 100] <- 100  # Carrying capacity constraint
    list(dBdt)
  })
}

simulate_dynamics_c <- function(params, model, init_biomass = runif(params$S, min = 1, max = 10)) {
  init_biomass <- as.numeric(init_biomass)
  times <- seq(from = 1, to = params$maxt)
  out <- deSolve::ode(
    y = init_biomass,
    times = times,
    func = fw.model,
    parms = params
  ) %>% as.data.frame()
  return(out)
}

sim_quantitative_network <- function(S, C, aij_params) {
  A <- matrix(0, S, S)
  n_pairs <- S * (S - 1) / 2
  B <- runif(n_pairs) <= C
  
  aij <- -abs(rnorm(n_pairs, aij_params[1], aij_params[2]))
  A[upper.tri(A)] <- B * aij
  A <- t(A)
  aij <- abs(rnorm(n_pairs, aij_params[1], aij_params[2]))
  A[upper.tri(A)] <- B * aij
  diag(A) <- -runif(S, min = 0, max = 1)
  
  while (max(Re(eigen(A)$values)) > 0) {
    diag(A) <- -runif(S, min = 0, max = 1)
  }
  return(A)
}

simulate_dynamics_perturbed <- function(params, covMat, S, maxt, perturb_scale = 1) {
  init_biomass <- runif(S, min = 1, max = 10)
  r_pre <- -params$A %*% init_biomass
  
  dyn_params <- list(A = params$A, r = r_pre, S = S, maxt = maxt)
  pre_perturb <- simulate_dynamics_c(dyn_params, fw.model)
  equilibrium_pre <- as.numeric(pre_perturb[nrow(pre_perturb), -1])
  
  r_perturbation <- MASS::mvrnorm(n = 1, mu = r_pre, Sigma = covMat * perturb_scale)
  dyn_params$r <- matrix(r_perturbation, nrow = S, ncol = 1)
  
  post_perturb <- simulate_dynamics_c(dyn_params, fw.model, init_biomass = equilibrium_pre)
  equilibrium_post <- as.numeric(post_perturb[nrow(post_perturb), -1])
  
  delta_X <- equilibrium_post - equilibrium_pre
  return(list(delta_r = dyn_params$r - r_pre, delta_X = delta_X))
}

# Simulation parameters
S <- 8
C <- 0.2
aij_params <- c(0, 0.5)
sd_X <- 1
num_scenarios <- 100
num_simulations <- 20  # 20 simulations per correlation matrix scenario
maxt <- 1000

# Fixed interaction matrix A
A <- sim_quantitative_network(S, C, aij_params)

# Initialize results storage
results <- data.frame()

# Run simulations over multiple correlation scenarios
for (scenario in 1:num_scenarios) {
  print(paste("Scenario", scenario, "of", num_scenarios))
  
  # Generate random correlation matrix
  cor_matrix <- randcorr(S)
  covMat <- diag(sd_X, S) %*% cor_matrix %*% diag(sd_X, S)
  
  # Store summary stats of correlation matrix
  mean_cor <- mean(cor_matrix[upper.tri(cor_matrix)])
  var_cor <- var(cor_matrix[upper.tri(cor_matrix)])
  
  # Run multiple simulations per correlation matrix
  for (sim in 1:num_simulations) {
    sim_result <- simulate_dynamics_perturbed(list(A = A), covMat, S, maxt)
    
    # Store delta_X and correlation stats
    for (species in 1:S) {
      results <- rbind(results, data.frame(
        Scenario = scenario,
        Simulation = sim,
        Species = paste0("sp", species),
        Delta_X = sim_result$delta_X[species],
        Mean_Correlation = mean_cor,
        Var_Correlation = var_cor
      ))
    }
  }
}

# Compute SD of Delta_X for each species across simulations within each scenario
delta_summary <- results %>%
  group_by(Scenario, Species) %>%
  summarize(SD_Delta_X = sd(Delta_X), Mean_Correlation = first(Mean_Correlation), Var_Correlation = first(Var_Correlation), .groups = 'drop') %>%
  group_by(Scenario) %>%
  summarize(Mean_SD_Delta_X = mean(SD_Delta_X), Mean_Correlation = first(Mean_Correlation), Var_Correlation = first(Var_Correlation))

# Plot Mean Correlation vs. SD of abundance change
ggplot(delta_summary, aes(x = Mean_Correlation, y = Mean_SD_Delta_X)) +
  geom_point(alpha = 0.6, color = "blue") +
  labs(x = "Mean Correlation", y = "Mean SD in Abundance Change") +
  theme_minimal() +
  ggtitle("Mean SD in Abundance Change vs. Mean Correlation")

# Plot Variance of Correlation vs. SD of abundance change
ggplot(delta_summary, aes(x = Var_Correlation, y = Mean_SD_Delta_X)) +
  geom_point(alpha = 0.6, color = "red") +
  labs(x = "Variance of Correlation", y = "Mean SD in Abundance Change") +
  theme_minimal() +
  ggtitle("Mean SD in Abundance Change vs. Variance of Correlation")





# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(deSolve)
library(mvtnorm)
library(randcorr)

# Define functions for dynamics and network generation
fw.model <- function(t, B, params) {
  with(as.list(c(B, params)), {
    B[B < 10^-8] <- 0 
    B[B < 0] <- 0 
    dBdt <- (params$r + params$A %*% B) * B
    dBdt[B > 100] <- 100  # Carrying capacity constraint
    list(dBdt)
  })
}

simulate_dynamics_c <- function(params, model, init_biomass = runif(params$S, min = 1, max = 10)) {
  init_biomass <- as.numeric(init_biomass)
  times <- seq(from = 1, to = params$maxt)
  out <- deSolve::ode(
    y = init_biomass,
    times = times,
    func = fw.model,
    parms = params
  ) %>% as.data.frame()
  return(out)
}

sim_quantitative_network <- function(S, C, aij_params) {
  A <- matrix(0, S, S)
  n_pairs <- S * (S - 1) / 2
  B <- runif(n_pairs) <= C
  
  aij <- -abs(rnorm(n_pairs, aij_params[1], aij_params[2]))
  A[upper.tri(A)] <- B * aij
  A <- t(A)
  aij <- abs(rnorm(n_pairs, aij_params[1], aij_params[2]))
  A[upper.tri(A)] <- B * aij
  diag(A) <- -runif(S, min = 0, max = 1)
  
  while (max(Re(eigen(A)$values)) > 0) {
    diag(A) <- -runif(S, min = 0, max = 1)
  }
  return(A)
}

simulate_dynamics_perturbed <- function(params, covMat, S, maxt, perturb_scale = 1) {
  init_biomass <- runif(S, min = 1, max = 10)
  r_pre <- -params$A %*% init_biomass
  
  dyn_params <- list(A = params$A, r = r_pre, S = S, maxt = maxt)
  pre_perturb <- simulate_dynamics_c(dyn_params, fw.model)
  equilibrium_pre <- as.numeric(pre_perturb[nrow(pre_perturb), -1])
  
  r_perturbation <- MASS::mvrnorm(n = 1, mu = r_pre, Sigma = covMat * perturb_scale)
  dyn_params$r <- matrix(r_perturbation, nrow = S, ncol = 1)
  
  post_perturb <- simulate_dynamics_c(dyn_params, fw.model, init_biomass = equilibrium_pre)
  equilibrium_post <- as.numeric(post_perturb[nrow(post_perturb), -1])
  
  delta_X <- equilibrium_post - equilibrium_pre
  return(list(delta_r = dyn_params$r - r_pre, delta_X = delta_X))
}

# Simulation parameters
S <- 8
C <- 0.2
aij_params <- c(0, 0.5)
sd_X <- 1
num_scenarios <- 100
num_simulations <- 20  # 20 simulations per correlation matrix scenario
maxt <- 1000

# Fixed interaction matrix A
A <- sim_quantitative_network(S, C, aij_params)

# Initialize results storage
results <- data.frame()

# Run simulations over multiple correlation scenarios
for (scenario in 1:num_scenarios) {
  print(paste("Scenario", scenario, "of", num_scenarios))
  
  # Generate random correlation matrix
  cor_matrix <- randcorr(S)
  covMat <- diag(sd_X, S) %*% cor_matrix %*% diag(sd_X, S)
  
  # Store summary stats of correlation matrix
  mean_cor <- mean(cor_matrix[upper.tri(cor_matrix)])
  var_cor <- var(cor_matrix[upper.tri(cor_matrix)])
  
  # Run multiple simulations per correlation matrix
  for (sim in 1:num_simulations) {
    sim_result <- simulate_dynamics_perturbed(list(A = A), covMat, S, maxt)
    
    # Store delta_X and correlation stats
    for (species in 1:S) {
      results <- rbind(results, data.frame(
        Scenario = scenario,
        Simulation = sim,
        Species = paste0("sp", species),
        Delta_X = sim_result$delta_X[species],
        Mean_Correlation = mean_cor,
        Var_Correlation = var_cor
      ))
    }
  }
}

# Compute SD of Delta_X for each species across simulations within each scenario
delta_summary <- results %>%
  group_by(Scenario, Species) %>%
  summarize(SD_Delta_X = sd(Delta_X), Mean_Correlation = first(Mean_Correlation), Var_Correlation = first(Var_Correlation), .groups = 'drop') %>%
  group_by(Scenario) %>%
  summarize(Mean_SD_Delta_X = mean(SD_Delta_X), Mean_Correlation = first(Mean_Correlation), Var_Correlation = first(Var_Correlation))

# Plot Mean Correlation vs. SD of abundance change
ggplot(delta_summary, aes(x = Mean_Correlation, y = Mean_SD_Delta_X)) +
  geom_point(alpha = 0.6, color = "blue") +
  labs(x = "Mean Correlation", y = "Mean SD in Abundance Change") +
  theme_minimal() +
  ggtitle("Mean SD in Abundance Change vs. Mean Correlation")

# Plot Variance of Correlation vs. SD of abundance change
ggplot(delta_summary, aes(x = Var_Correlation, y = Mean_SD_Delta_X)) +
  geom_point(alpha = 0.6, color = "red") +
  labs(x = "Variance of Correlation", y = "Mean SD in Abundance Change") +
  theme_minimal() +
  ggtitle("Mean SD in Abundance Change vs. Variance of Correlation")



###################################


# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(deSolve)
library(mvtnorm)
library(randcorr)

# Define functions for dynamics and network generation
fw.model <- function(t, B, params) {
  with(as.list(c(B, params)), {
    B[B < 10^-8] <- 0 
    B[B < 0] <- 0 
    dBdt <- (params$r + params$A %*% B) * B
    dBdt[B > 100] <- 100  # Carrying capacity constraint
    list(dBdt)
  })
}

simulate_dynamics_c <- function(params, model, init_biomass = runif(params$S, min = 1, max = 10)) {
  init_biomass <- as.numeric(init_biomass)
  times <- seq(from = 1, to = params$maxt)
  out <- deSolve::ode(
    y = init_biomass,
    times = times,
    func = fw.model,
    parms = params
  ) %>% as.data.frame()
  return(out)
}

sim_quantitative_network <- function(S, C, aij_params) {
  A <- matrix(0, S, S)
  n_pairs <- S * (S - 1) / 2
  B <- runif(n_pairs) <= C
  
  aij <- -abs(rnorm(n_pairs, aij_params[1], aij_params[2]))
  A[upper.tri(A)] <- B * aij
  A <- t(A)
  aij <- abs(rnorm(n_pairs, aij_params[1], aij_params[2]))
  A[upper.tri(A)] <- B * aij
  diag(A) <- -runif(S, min = 0, max = 1)
  
  while (max(Re(eigen(A)$values)) > 0) {
    diag(A) <- -runif(S, min = 0, max = 1)
  }
  return(A)
}

simulate_dynamics_perturbed <- function(params, covMat, S, maxt, perturb_scale = 1) {
  init_biomass <- runif(S, min = 1, max = 10)
  r_pre <- -params$A %*% init_biomass
  
  dyn_params <- list(A = params$A, r = r_pre, S = S, maxt = maxt)
  pre_perturb <- simulate_dynamics_c(dyn_params, fw.model)
  equilibrium_pre <- as.numeric(pre_perturb[nrow(pre_perturb), -1])
  
  r_perturbation <- MASS::mvrnorm(n = 1, mu = r_pre, Sigma = covMat * perturb_scale)
  dyn_params$r <- matrix(r_perturbation, nrow = S, ncol = 1)
  
  post_perturb <- simulate_dynamics_c(dyn_params, fw.model, init_biomass = equilibrium_pre)
  equilibrium_post <- as.numeric(post_perturb[nrow(post_perturb), -1])
  
  delta_X <- equilibrium_post - equilibrium_pre
  return(list(delta_r = dyn_params$r - r_pre, delta_X = delta_X))
}

# Simulation parameters
S <- 8
C <- 0.2
aij_params <- c(0, 0.5)
sd_X <- 1
num_scenarios <- 100
num_simulations <- 20  # 20 simulations per correlation matrix scenario
maxt <- 1000

# Fixed interaction matrix A
A <- sim_quantitative_network(S, C, aij_params)

# Initialize results storage
results <- data.frame()

# Run simulations over multiple correlation scenarios
for (scenario in 1:num_scenarios) {
  print(paste("Scenario", scenario, "of", num_scenarios))
  
  # Generate random correlation matrix
  cor_matrix <- randcorr(S)
  covMat <- diag(sd_X, S) %*% cor_matrix %*% diag(sd_X, S)
  
  # Store summary stats of correlation matrix
  mean_cor <- mean(cor_matrix[upper.tri(cor_matrix)])
  var_cor <- var(cor_matrix[upper.tri(cor_matrix)])
  
  # Run multiple simulations per correlation matrix
  for (sim in 1:num_simulations) {
    sim_result <- simulate_dynamics_perturbed(list(A = A), covMat, S, maxt)
    
    # Store delta_X and correlation stats
    for (species in 1:S) {
      results <- rbind(results, data.frame(
        Scenario = scenario,
        Simulation = sim,
        Species = paste0("sp", species),
        Delta_X = sim_result$delta_X[species],
        Mean_Correlation = mean_cor,
        Var_Correlation = var_cor
      ))
    }
  }
}

# Compute SD of Delta_X for each species across simulations within each scenario
delta_summary <- results %>%
  group_by(Scenario, Species) %>%
  summarize(SD_Delta_X = sd(Delta_X), 
            Mean_Correlation = first(Mean_Correlation), 
            Var_Correlation = first(Var_Correlation), 
            .groups = 'drop')

# Reorder scenarios by mean correlation
delta_summary <- delta_summary %>%
  mutate(Scenario = factor(Scenario, levels = unique(Scenario[order(Mean_Correlation)])))

# Plot Mean Correlation vs. SD of abundance change with ordered x-axis
ggplot(delta_summary, aes(x = Scenario, y = SD_Delta_X)) +
  geom_boxplot(aes(fill = as.factor(round(Mean_Correlation, 2)))) +
  labs(x = "Scenario (Ordered by Mean Correlation)", y = "SD of Delta_X", fill = "Mean Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  ggtitle("SD of Abundance Change Across Scenarios Ordered by Mean Correlation")

# Reorder scenarios by variance of correlation
delta_summary <- delta_summary %>%
  mutate(Scenario = factor(Scenario, levels = unique(Scenario[order(Var_Correlation)])))

# Plot Variance of Correlation vs. SD of abundance change with ordered x-axis
ggplot(delta_summary, aes(x = Scenario, y = SD_Delta_X)) +
  geom_boxplot(aes(fill = as.factor(round(Var_Correlation, 2)))) +
  labs(x = "Scenario (Ordered by Variance of Correlation)", y = "SD of Delta_X", fill = "Variance of Correlation") +
  theme_minimal() +
  theme(axis.text.x = element_blank()) +
  ggtitle("SD of Abundance Change Across Scenarios Ordered by Variance of Correlation")

