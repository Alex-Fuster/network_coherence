library(ggplot2)
library(tidyr)
library(dplyr)
library(Matrix)
library(deSolve)
library(clusterGeneration)

# Function to simulate quantitative network
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

# Function to simulate dynamics with continuous press perturbations
# fw.model_with_perturb <- function(t, B, params) {
#   with(as.list(c(B, params)), {
#     B[B < 10^-8] <- 0 
#     B[B < 0] <- 0 # Prevent negative biomass
#     # Small perturbation at each step, influenced by covariance structure
#     delta_r <- MASS::mvrnorm(n = 1, mu = rep(0, length(B)), Sigma = params$covMat)
#     dBdt <- (params$r + params$A %*% B + delta_r) * B
#     list(dBdt)
#   })
# }

perturbation_counter <- 0

# Function to simulate dynamics with continuous press perturbations
fw.model_with_perturb <- function(t, B, params) {
  with(as.list(c(B, params)), {
    B[B < 10^-8] <- 0 
    B[B < 0] <- 0 
    
    # Increment the global counter
    perturbation_counter <<- perturbation_counter + 1
    
    # Apply perturbation every n function calls
    if (perturbation_counter %% 10 == 0) {  # Check if 50 function calls have passed
      delta_r <- MASS::mvrnorm(n = 1, mu = rep(0, length(B)), Sigma = params$covMat)
      # Scale delta_r relative to the current biomass
      delta_r <- delta_r * B * perturbation_scale # Scale delta_r by biomass
      print(paste("biomass values:", B))
      print(paste("delta_r values:",delta_r))
    } else {
      delta_r <- 0  # No perturbation in this timestep
    }
    
    dBdt <- (params$r + params$A %*% B + delta_r) * B
    list(dBdt)
  })
}

# Simulate dynamics with press perturbations over multiple timesteps
simulate_dynamics_press <- function(params, model, init_biomass = runif(params$S, min = 1, max = 10), timesteps = 1000) {
  init_biomass <- as.numeric(init_biomass)
  times <- seq(from = 1, to = timesteps)
  
  out <- deSolve::ode(
    y = init_biomass,
    times = times,
    func = model,
    parms = params,
    maxsteps = 50000, #50000, # Increase maxsteps to handle large amounts of work
    method = "vode", # Use a solver suitable for stiff systems
    atol = 1e-8,            # Set a tighter absolute tolerance 
    rtol = 1e-8             # Set a tighter relative tolerance
  ) %>% as.data.frame()
  
  return(out)
}

# Function to simulate response with continuous perturbations
simulate_response_press <- function(S, C, aij_params, mu_delta_r, sd_delta_r, sd_X, maxt, alpha_d) {
  
  perturbation_counter <- 0
  
  # Compute C matrix using rcorrmatrix
  cor_matrix <- clusterGeneration::rcorrmatrix(d = S, alphad = alpha_d)
  
  # Create covariance matrix
  covMat <- diag(sd_X) %*% cor_matrix %*% diag(sd_X)
  
  # Simulate interaction matrix A
  A <- sim_quantitative_network("predator-prey", S = S, C = C, aij_params = aij_params)
  
  equilibrium_pre <- runif(S, 1, 10)
  r_pre <- -A %*% equilibrium_pre
  dyn_params <- list(A = A, r = r_pre, S = S, maxt = maxt, covMat = covMat)
  
  # Simulate dynamics with press perturbations over time
  out <- simulate_dynamics_press(dyn_params, fw.model_with_perturb, init_biomass = equilibrium_pre, timesteps = 1000)
  
  return(out)
}

# Set up the alpha_d values and species number to test
alpha_d_values <- c(0.001, 40)
S <- 10  # Fixed number of species

# Run the simulation for a fixed number of species and varying alpha_d values
results <- data.frame()

# Initialize delta_r_values, cor_values, and cov_values
delta_r_values <- data.frame(species = character(0), alpha_d = numeric(0), delta_r = numeric(0))
cor_values <- data.frame(alpha_d = numeric(0), correlation = numeric(0))
cov_values <- data.frame(alpha_d = numeric(0), covariance = numeric(0))

nsim <- 1  # Number of simulations per scenario

perturbation_scale = 5

for (alpha_d in alpha_d_values) {
  print(paste("Running simulation for alpha_d =", alpha_d))
  
  for (j in 1:nsim) {
    result <- simulate_response_press(S = S, 
                                      C = 0.2, 
                                      aij_params = c(0, 0.5), 
                                      mu_delta_r = 0, 
                                      sd_delta_r = 0.5, 
                                      sd_X = rep(1, S), 
                                      maxt = 50000, #50000, 
                                      alpha_d = alpha_d)
    
    delta_r <- result$delta_r  # Capture delta_r values
    
    # Extract biomass dynamics for each species over time
    biomass_df <- result %>%
      pivot_longer(-time, names_to = "species", values_to = "biomass")
    
    # Calculate the total biomass change over time for each species
    biomass_df <- biomass_df %>%
      group_by(species) %>%
      summarize(sum_deltaX = sum(biomass - lag(biomass, default = biomass[1])),
                sd_deltaX = sd(biomass, na.rm = TRUE))  # Adding sd per simulation
    
    # Store results for sum_deltaX and sd_deltaX
    results <- rbind(
      results,
      data.frame(
        alpha_d = alpha_d,
        sum_deltaX = mean(biomass_df$sum_deltaX),
        sd_deltaX = mean(biomass_df$sd_deltaX)
      )
    )
    
    
  }
}

# Plot 1: Mean and Standard Deviation of Biomass Change for Each Scenario

# Sum of biomass changes across alpha_d
p1 <- ggplot(results, aes(x = as.factor(alpha_d), y = sum_deltaX)) +
  geom_boxplot(fill = "skyblue") +
  labs(title = "Mean of Sum Biomass Changes Across alpha_d", x = "alpha_d", y = "Mean Sum of Biomass Change") +
  theme_minimal()

# Standard deviation of biomass changes across alpha_d
p2 <- ggplot(results, aes(x = as.factor(alpha_d), y = sd_deltaX)) +
  geom_boxplot(fill = "orange") +
  labs(title = "Mean of Standard Deviation of Biomass Changes Across alpha_d", x = "alpha_d", y = "Mean Standard Deviation of Biomass Change") +
  theme_minimal()

print(p1)
print(p2)

###########################

#  Plot delta_r distributions per scenario
p3 <-  ggplot(delta_r_values, aes(x = delta_r)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  facet_wrap(~scenario, scales = "free") +  # Separate plots for each scenario
  theme_minimal() +
  ggtitle("Distribution of delta_r Across Scenarios") +
  labs(x = "delta_r", y = "Frequency")

# 2. Distribution of correlation values across alpha_d
p4 <- ggplot(cor_values, aes(x = correlation)) +
  geom_histogram(bins = 30, fill = "green", alpha = 0.7) +
  facet_wrap(~ alpha_d, scales = "free") +
  theme_minimal() +
  ggtitle("Distribution of Correlation Values Across alpha_d Values") +
  labs(x = "Correlation", y = "Frequency")

# 3. Distribution of covariance values across alpha_d
p5 <- ggplot(cov_values, aes(x = covariance)) +
  geom_histogram(bins = 30, fill = "red", alpha = 0.7) +
  facet_wrap(~ alpha_d, scales = "free") +
  theme_minimal() +
  ggtitle("Distribution of Covariance Values Across alpha_d Values") +
  labs(x = "Covariance", y = "Frequency")

print(p3)
print(p4)
print(p5)
