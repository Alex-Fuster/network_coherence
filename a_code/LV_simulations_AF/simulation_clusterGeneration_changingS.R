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

# Function to simulate dynamics
fw.model <- function(t, B, params) {
  with(as.list(c(B, params)), {
    B[B < 10^-8] <- 0 
    B[B < 0] <- 0 # Prevent negative biomass
    dBdt <- (params$r + params$A %*% B) * B
    list(dBdt)
  })
}


eventfun <- function(t, B, parms){
  with(as.list(B), {
    B[B < 1e-6] <- 0
    return(B)
  })
}

simulate_dynamics_c <- function(params, model, init_biomass = runif(params$S, min = 1, max = 10)) {
  init_biomass <- as.numeric(init_biomass)
  times <- seq(from = 1, to = params$maxt)
  out <- deSolve::ode(
    y = init_biomass,
    times = times,
    func = fw.model,
    parms = params,
    events = list(func = eventfun, time = times)
  ) %>% as.data.frame()
  return(out)
}


# Function to simulate response
simulate_response <- function(S, C, aij_params, mu_delta_r, sd_delta_r, sd_X, maxt, alpha_d){
  
  # Compute C matrix using rcorrmatrix
  cor_matrix <- clusterGeneration::rcorrmatrix(d = S, alphad = alpha_d)
  
  # Create covariance matrix
  covMat <- diag(sd_X) %*% cor_matrix %*% diag(sd_X)
  
  # Simulate interaction matrix A
  A <- sim_quantitative_network("predator-prey", S = S, C = C, aij_params = aij_params)
  
  equilibrium_pre <- runif(S, 1, 10)
  r_pre <- -A %*% equilibrium_pre
  dyn_params <- list(A = A, r = r_pre, S = S, maxt = maxt) 
  pre_perturb <- simulate_dynamics_c(dyn_params, fw.model)
  equilibrium_pre <- as.numeric(pre_perturb[nrow(pre_perturb), -1])
  
  # Generate delta_r
  delta_r <- MASS::mvrnorm(n = 1, mu =  rep(0, S), Sigma = covMat, empirical = FALSE)
  
  dyn_params$r <- dyn_params$r + delta_r 
  post_perturb <- simulate_dynamics_c(dyn_params, fw.model, init_biomass = equilibrium_pre)
  equilibrium_post <- as.numeric(post_perturb[nrow(post_perturb), -1])
  
  df <- data.frame(species = paste0("sp", c(1:S)), 
                   X_pre = equilibrium_pre, 
                   X_post = equilibrium_post)
  
  return(list(df = df, pre_perturb = pre_perturb, 
              post_perturb = post_perturb, 
              Sigma = covMat, 
              C_matrix = cor_matrix, 
              delta_r = delta_r))
}

# Set up the alpha_d values and species numbers to test
alpha_d_values <- c(0.1,20)
species_groups <- c(5, 20, 50)

# Run the simulation for different species and alpha_d values
results <- data.frame()

# Initialize delta_r_values, cor_values, and cov_values with the correct column names and types
delta_r_values <- data.frame(
  species_group = numeric(0), 
  alpha_d = numeric(0), 
  delta_r = numeric(0), 
  species = character(0)
)

cor_values <- data.frame(
  species_group = numeric(0),
  alpha_d = numeric(0),
  correlation = numeric(0)
)

cov_values <- data.frame(
  species_group = numeric(0),
  alpha_d = numeric(0),
  covariance = numeric(0)
)

nsim <- 20

for (S in species_groups) {
  for (alpha_d in alpha_d_values) {
    print(paste("Running simulation for S =", S, "alpha_d =", alpha_d))
    
    for (j in 1:nsim) {
      result <- simulate_response(S = S, 
                                  C = 0.2, 
                                  aij_params = c(0, 0.5), 
                                  mu_delta_r = 0, 
                                  sd_delta_r = 0.5, 
                                  sd_X = rep(1.5, S), 
                                  maxt = 100, 
                                  alpha_d = alpha_d)
      
      df <- result$df
      delta_r <- result$delta_r  # Capture delta_r values
      cor_matrix <- result$C_matrix  # Correlation matrix
      cov_matrix <- result$Sigma  # Covariance matrix
      
      # Store results for sum_deltaX and sd_deltaX
      results <- rbind(
        results,
        data.frame(
          species_group = S,
          alpha_d = alpha_d,
          sum_deltaX = sum(df$X_pre - df$X_post),
          sd_deltaX = sd(df$X_pre - df$X_post)
        )
      )
      
      # Store delta_r values for each species
      delta_r_df <- data.frame(
        species_group = rep(S, length(delta_r)),
        alpha_d = rep(alpha_d, length(delta_r)),
        delta_r = delta_r,
        species = paste0("sp", seq_len(S))
      )
      
      # Store the values from the correlation and covariance matrices
      n_upper_tri <- S * (S - 1) / 2  # Number of elements in the upper triangle of the matrix
      
      cor_df <- data.frame(
        species_group = rep(S, n_upper_tri),
        alpha_d = rep(alpha_d, n_upper_tri),
        correlation = cor_matrix[upper.tri(cor_matrix)]
      )
      
      cov_df <- data.frame(
        species_group = rep(S, n_upper_tri),
        alpha_d = rep(alpha_d, n_upper_tri),
        covariance = cov_matrix[upper.tri(cov_matrix)]
      )
      
      # Append to results data frames
      delta_r_values <- rbind(delta_r_values, delta_r_df)
      cor_values <- rbind(cor_values, cor_df)
      cov_values <- rbind(cov_values, cov_df)
    }
  }
}


# Plot results
# 1. Sum of biomass changes
p1 <- ggplot(results, aes(x = as.factor(alpha_d), y = sum_deltaX, fill = as.factor(species_group))) +
  geom_boxplot() +
  labs(title = "Sum of Biomass Changes Across alpha_d", x = "alpha_d", y = "Sum of Biomass Change", fill = "Species Group") +
  theme_minimal()

# 2. Standard deviation of biomass changes
p2 <- ggplot(results, aes(x = as.factor(alpha_d), y = sd_deltaX, fill = as.factor(species_group))) +
  geom_boxplot() +
  labs(title = "Standard Deviation of Biomass Changes Across alpha_d", x = "alpha_d", y = "Standard Deviation of Biomass Change", fill = "Species Group") +
  theme_minimal()

print(p1)
print(p2)

###########################

# Calculate global x-axis limits for delta_r, correlation, and covariance
delta_r_limits <- range(delta_r_values$delta_r, na.rm = TRUE)
correlation_limits <- range(cor_values$correlation, na.rm = TRUE)
covariance_limits <- range(cov_values$covariance, na.rm = TRUE)

# 1. Distribution of delta_r values across alpha_d
p3 <- ggplot(delta_r_values, aes(x = delta_r)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  facet_grid(species_group ~ alpha_d, scales = "free") +
  xlim(delta_r_limits) +  # Apply global limits
  theme_minimal() +
  ggtitle("Distribution of delta_r Across alpha_d Values") +
  labs(x = "delta_r", y = "Frequency")

# 2. Distribution of correlation values across alpha_d
p4 <- ggplot(cor_values, aes(x = correlation)) +
  geom_histogram(bins = 30, fill = "green", alpha = 0.7) +
  facet_grid(species_group ~ alpha_d, scales = "free") +
  xlim(correlation_limits) +  # Apply global limits
  theme_minimal() +
  ggtitle("Distribution of Correlation Values Across alpha_d Values") +
  labs(x = "Correlation", y = "Frequency")

# 3. Distribution of covariance values across alpha_d
p5 <- ggplot(cov_values, aes(x = covariance)) +
  geom_histogram(bins = 30, fill = "red", alpha = 0.7) +
  facet_grid(species_group ~ alpha_d, scales = "free") +
  xlim(covariance_limits) +  # Apply global limits
  theme_minimal() +
  ggtitle("Distribution of Covariance Values Across alpha_d Values") +
  labs(x = "Covariance", y = "Frequency")

# Print the plots
print(p3)
print(p4)
print(p5)


