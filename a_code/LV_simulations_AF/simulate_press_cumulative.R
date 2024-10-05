simulate_with_covariance <- function(S, timesteps, cov_matrix, perturb_interval = 10, perturb_scale = 5) {
  # Initialize species parameters and biomasses
  A <- sim_quantitative_network("predator-prey", S = S, C = 0.2, aij_params = c(0, 0.5))
  init_biomass <- runif(S, min = 1, max = 10)
  r <- -A %*% init_biomass  # Initial r values based on biomass
  
  # Accumulating growth rates for each species (initialized as r)
  r_new <- r
  
  # Define model dynamics with periodic perturbations
  fw.model_with_predefined_cov <- function(t, B, params) {
    with(as.list(c(B, params)), {
      B[B < 10^-8] <- 0  # Avoid negative biomass
      step <- floor(t / perturb_interval) + 1  # Determine perturbation step
      
      # Apply delta_r every perturb_interval steps
      if (t %% perturb_interval == 0) {
        delta_r <- MASS::mvrnorm(n = 1, mu = rep(0, S), Sigma = params$cov_matrix)
        delta_r <- delta_r * perturb_scale  # Apply scaling to perturbation
        r_new <<- r_new + delta_r  # Accumulate the changes in r over time
        print(paste("Applying perturbation at step", t, ": delta_r =", paste(delta_r, collapse = ", ")))
      }
      
      # System dynamics evolve continuously based on r_new and A
      dBdt <- (r_new + params$A %*% B) * B  # Continuous dynamics from interactions
      list(dBdt)
    })
  }
  
  # Run the simulation
  out <- deSolve::ode(y = init_biomass, 
                      times = seq(1, timesteps), 
                      func = fw.model_with_predefined_cov, 
                      parms = list(A = A, r = r, cov_matrix = cov_matrix),
                      method = "ode45",  # Use a solver for non-stiff systems
                      atol = 1e-8, rtol = 1e-8)
  
  return(as.data.frame(out))
}

# Example covariance matrix for strong species interactions
covMat <- diag(S) * 0.5 + matrix(0.2, S, S)  # Adjust off-diagonal covariance

# Run the simulation with scaling applied
result <- simulate_with_covariance(S = 20, timesteps = 1000, cov_matrix = covMat, perturb_scale = 5)

# Plot the result to see changes in biomass over time
result_long <- result %>%
  pivot_longer(-time, names_to = "species", values_to = "biomass")

ggplot(result_long, aes(x = time, y = biomass, color = species)) +
  geom_line() +
  labs(title = "Biomass Dynamics with Continuous Interaction and Perturbations", x = "Time", y = "Biomass") +
  theme_minimal()


####################################################


# Function to simulate dynamics with discrete perturbations
simulate_with_discrete_perturbations <- function(S, timesteps, cov_matrix, perturb_interval = 10, perturb_scale = 5) {
  # Initialize species parameters and biomasses
  A <- sim_quantitative_network("predator-prey", S = S, C = 0.2, aij_params = c(0, 0.5))
  init_biomass <- runif(S, min = 1, max = 10)
  r <- -A %*% init_biomass  # Initial r values based on biomass
  
  # Accumulating growth rates for each species (initialized as r)
  r_new <- r
  biomass <- init_biomass
  
  # Create a data frame to store results, ensuring that column names match
  results <- data.frame(time = 1, t(biomass))  # Transpose biomass vector to match the shape
  
  # Set appropriate column names for the species (avoid name mismatch)
  colnames(results) <- c("time", paste0("species_", 1:S))
  
  # Iterate over time in discrete blocks, applying perturbations at specified intervals
  for (t in seq(perturb_interval, timesteps, by = perturb_interval)) {
    # Apply perturbation at each interval
    delta_r <- MASS::mvrnorm(n = 1, mu = rep(0, S), Sigma = cov_matrix) * perturb_scale
    r_new <- r_new + delta_r  # Update growth rates with perturbation
    
    # Define the ODE model using the new r_new
    fw.model_block <- function(t, B, params) {
      with(as.list(c(B, params)), {
        B[B < 10^-8] <- 0  # Avoid negative biomass
        dBdt <- (params$r_new + params$A %*% B) * B  # Use updated r_new
        list(dBdt)
      })
    }
    
    # Simulate the system for the next perturb_interval timesteps
    out <- deSolve::ode(
      y = biomass, 
      times = seq(t - perturb_interval + 1, t), 
      func = fw.model_block, 
      parms = list(A = A, r_new = r_new),
      method = "ode45"  # Use a solver for non-stiff systems
    )
    
    # Update biomass to the last step of the block and store results
    biomass <- as.numeric(out[nrow(out), -1])  # Extract the last row (new biomass)
    time_block_results <- data.frame(time = t, t(biomass))  # Transpose and create a data frame
    colnames(time_block_results) <- colnames(results)  # Ensure the column names match
    
    results <- rbind(results, time_block_results)  # Append the new results
    
    print(paste("Applied perturbation at time", t, ": delta_r =", paste(delta_r, collapse = ", ")))
  }
  
  return(results)
}

# Example covariance matrix for strong species interactions
covMat <- diag(S) * 0.5 + matrix(0.2, S, S)  # Adjust off-diagonal covariance

# Run the simulation with scaling applied and discrete perturbations
result <- simulate_with_discrete_perturbations(S = 20, timesteps = 1000, cov_matrix = covMat, perturb_scale = 1.2)

# Convert results to long format for plotting
result_long <- result %>%
  pivot_longer(-time, names_to = "species", values_to = "biomass")

# Plot the result to see changes in biomass over time
ggplot(result_long, aes(x = time, y = biomass, color = species)) +
  geom_line() +
  labs(title = "Biomass Dynamics with Discrete Perturbations", x = "Time", y = "Biomass") +
  theme_minimal()
