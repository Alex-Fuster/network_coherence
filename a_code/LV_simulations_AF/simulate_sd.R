library(ggplot2)


# Function to simulate quantitative network (interaction matrix)
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


# Function to simulate dynamics with a food web model
fw.model <- function(t, B, params) {
  with(as.list(c(B, params)), {
    # Ensure no negative biomass
    B[B < 10^-8] <- 0
    
    # Compute the rate of change in biomass (dB/dt)
    dBdt <- (params$r + params$A %*% B) * B
    list(dBdt)
  })
}

# Function to apply random perturbations
apply_random_perturbations <- function(biomass, perturbation_sd) {
  perturbation <- rnorm(length(biomass), mean = 0, sd = perturbation_sd)
  biomass <- biomass + perturbation
  biomass[biomass < 0] <- 0  # Ensure no negative biomass values
  return(biomass)
}

# Function to simulate dynamics over time with perturbations
simulate_dynamics_with_perturbations <- function(params, init_biomass = runif(params$S, min = 1, max = 10), perturbation_sd = 0.1, n_steps = 1000) {
  
  # Initialize biomass with random values
  init_biomass <- as.numeric(init_biomass)
  times <- seq(from = 1, to = n_steps)
  
  # Track biomass over time
  biomass_over_time <- matrix(NA, nrow = n_steps, ncol = params$S)
  biomass_over_time[1, ] <- init_biomass
  
  for (t in 2:n_steps) {
    
    #print(paste("step",t, "of", n_steps))
    # Apply dynamics
    out <- deSolve::ode(y = biomass_over_time[t-1, ], times = c(t-1, t), func = fw.model, parms = params)
    biomass <- out[2, -1]  # Extract the biomass at time t
    
    # Apply random perturbations
    biomass <- apply_random_perturbations(biomass, perturbation_sd)
    
    # Store biomass
    biomass_over_time[t, ] <- biomass
  }
  
  return(biomass_over_time)
}

# Set up parameters for simulation
params <- list(
  S = 20, # Number of species
  C = 0.2, # Connectance (fraction of possible links)
  aij_params = c(0, 0.5), # Range of interaction strengths
  maxt = 100 # Maximum time for the simulation
)

# Create the interaction matrix A (20 species with connectance 0.2)
params$A <- sim_quantitative_network("predator-prey",S = params$S, C = params$C, aij_params = params$aij_params)

# Initialize random growth rates for each species (r)
params$r <- runif(params$S, min = 0, max = 1)

# Simulate biomass dynamics with perturbations
biomass_over_time <- simulate_dynamics_with_perturbations(params, perturbation_sd = 0.1, n_steps = 1000)

# Plot biomass dynamics over time for each species

biomass_df <- as.data.frame(biomass_over_time)
biomass_df$time <- 1:1000
biomass_df_long <- tidyr::pivot_longer(biomass_df, cols = -time, names_to = "species", values_to = "biomass")

ggplot(biomass_df_long, aes(x = time, y = biomass, color = species)) +
  geom_line() +
  labs(title = "Biomass Dynamics Over Time", x = "Time", y = "Biomass") +
  theme_minimal()


### sd for each species

# Calculate the standard deviation of biomass for each species across the 1000 time steps
sd_X <- apply(biomass_over_time, 2, sd)

# View the standard deviation vector for each species
hist(sd_X)
