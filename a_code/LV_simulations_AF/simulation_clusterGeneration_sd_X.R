# Load necessary libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(Matrix)
library(deSolve)
library(clusterGeneration)
library(mvtnorm)
library(gridExtra)

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

# Function to simulate dynamics
fw.model <- function(t, B, params) {
  with(as.list(c(B, params)), {
    B[B < 10^-8] <- 0 
    B[B < 0] <- 0 
    dBdt <- (params$r + params$A %*% B) * B
    list(dBdt)
  })
}

eventfun <- function(t, B, parms) {
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

# Function to simulate response and compute delta_r
simulate_response <- function(S, C, aij_params, sd_X, maxt, alpha_d) {
  
  cor_matrix <- clusterGeneration::rcorrmatrix(d = S, alphad = alpha_d)
  covMat <- diag(sd_X) %*% cor_matrix %*% diag(sd_X)
  
  A <- sim_quantitative_network("predator-prey", S = S, C = C, aij_params = aij_params)
  
  equilibrium_pre <- runif(S, 1, 10)
  r_pre <- -A %*% equilibrium_pre
  
  dyn_params <- list(A = A, r = r_pre, S = S, maxt = maxt) 
  pre_perturb <- simulate_dynamics_c(dyn_params, fw.model)
  equilibrium_pre <- as.numeric(pre_perturb[nrow(pre_perturb), -1])
  
  r_perturbation <- mvtnorm::rmvnorm(n = 1, mean = r_pre, sigma = covMat)
  delta_r <- as.vector(matrix(r_perturbation, nrow = S, ncol = 1) - r_pre)
  
  dyn_params$r <- matrix(r_perturbation, nrow = S, ncol = 1)  
  
  post_perturb <- simulate_dynamics_c(dyn_params, fw.model, init_biomass = equilibrium_pre)
  equilibrium_post <- as.numeric(post_perturb[nrow(post_perturb), -1])
  
  df <- data.frame(species = paste0("sp", c(1:S)), 
                   X_pre = equilibrium_pre, 
                   X_post = equilibrium_post)
  
  return(list(df = df, pre_perturb = pre_perturb, 
              post_perturb = post_perturb, 
              Sigma = covMat, 
              C_matrix = cor_matrix, 
              r_pre = r_pre, 
              delta_r = delta_r))
}

# Simulation parameters
alpha_d <- 40  # Fixed alpha_d for all scenarios
sd_X_values <- c(0.01, 0.1, 0.3, 0.5, 0.8, 1, 1.2, 1.5)  # Different sd_X scenarios
S <- 20  
nsim <- 20  
maxt <- 1000

# Initialize data frames to store results
results <- data.frame()
delta_r_values <- data.frame(species = character(0), sd_X = numeric(0), delta_r = numeric(0))
cor_values <- data.frame(sd_X = numeric(0), correlation = numeric(0))
cov_values <- data.frame(sd_X = numeric(0), covariance = numeric(0))

for (sd_X in sd_X_values) {
  print(paste("Running simulation for sd_X =", sd_X))
  
  for (j in 1:nsim) {
    result <- simulate_response(S = S, 
                                C = 0.2, 
                                aij_params = c(0, 0.5), 
                                sd_X = rep(sd_X, S), 
                                maxt = maxt, 
                                alpha_d = alpha_d)
    
    df <- result$df
    delta_r <- result$delta_r  
    cor_matrix <- result$C_matrix  
    cov_matrix <- result$Sigma  
    
    results <- rbind(
      results,
      data.frame(
        sd_X = sd_X,
        sum_deltaX = sum(df$X_pre - df$X_post),
        sd_deltaX = sd(df$X_pre - df$X_post),
        var_deltaX = var(df$X_pre - df$X_post)
      )
    )
    
    delta_r_df <- data.frame(
      sd_X = rep(sd_X, length(delta_r)),
      delta_r = delta_r,
      species = paste0("sp", seq_len(S))
    )
    
    n_upper_tri <- S * (S - 1) / 2  
    
    cor_df <- data.frame(
      sd_X = rep(sd_X, n_upper_tri),
      correlation = cor_matrix[upper.tri(cor_matrix)]
    )
    
    cov_df <- data.frame(
      sd_X = rep(sd_X, n_upper_tri),
      covariance = cov_matrix[upper.tri(cov_matrix)]
    )
    
    delta_r_values <- rbind(delta_r_values, delta_r_df)
    cor_values <- rbind(cor_values, cor_df)
    cov_values <- rbind(cov_values, cov_df)
  }
}

# Plot results
p1 <- ggplot(results, aes(x = as.factor(sd_X), y = sum_deltaX)) +
  geom_boxplot(fill = "skyblue") +
  labs(title = "Sum of Biomass Changes Across sd_X", x = "sd_X", y = "Sum of Biomass Change") +
  theme_minimal()

p2 <- ggplot(results, aes(x = as.factor(sd_X), y = sd_deltaX)) +
  geom_boxplot(fill = "orange") +
  labs(title = "Standard Deviation of Biomass Changes Across sd_X", x = "sd_X", y = "Standard Deviation of Biomass Change") +
  theme_minimal()

delta_r_limits <- range(delta_r_values$delta_r, na.rm = TRUE)
correlation_limits <- range(cor_values$correlation, na.rm = TRUE)
covariance_limits <- range(cov_values$covariance, na.rm = TRUE)

p3 <- ggplot(delta_r_values, aes(x = delta_r)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  facet_wrap(~ sd_X, scales = "free") +
  xlim(delta_r_limits) +
  theme_minimal() +
  ggtitle("Distribution of delta_r Across sd_X Values") +
  labs(x = "delta_r", y = "Frequency")

p4 <- ggplot(cor_values, aes(x = correlation)) +
  geom_histogram(bins = 30, fill = "green", alpha = 0.7) +
  facet_wrap(~ sd_X, scales = "free") +
  xlim(correlation_limits) +
  theme_minimal() +
  ggtitle("Distribution of Correlation Values Across sd_X Values") +
  labs(x = "Correlation", y = "Frequency")

p5 <- ggplot(cov_values, aes(x = covariance)) +
  geom_histogram(bins = 30, fill = "red", alpha = 0.7) +
  facet_wrap(~ sd_X, scales = "free") +
  xlim(covariance_limits) +
  theme_minimal() +
  ggtitle("Distribution of Covariance Values Across sd_X Values") +
  labs(x = "Covariance", y = "Frequency")

# Display all plots
gridExtra::grid.arrange(p1, p2, p3, p4, p5, ncol = 2)
