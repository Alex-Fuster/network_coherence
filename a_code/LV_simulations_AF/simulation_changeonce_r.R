
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






# Simulation with constant variances and varying covariances

simulate_response <- function(S, C, aij_params, sd_X, maxt, alpha_d, perturb_scale) {
  # Generate a correlation matrix with desired structure
  cor_matrix <- clusterGeneration::rcorrmatrix(d = S, alphad = alpha_d)
  covMat <- diag(sd_X) %*% cor_matrix %*% diag(sd_X)
  
  
  A <- sim_quantitative_network("predator-prey", S = S, C = C, aij_params = aij_params)
  equilibrium_pre <- runif(S, 1, 10)
  r_pre <- -A %*% equilibrium_pre
  
  dyn_params <- list(A = A, r = r_pre, S = S, maxt = maxt) 
  pre_perturb <- simulate_dynamics_c(dyn_params, fw.model)
  equilibrium_pre <- as.numeric(pre_perturb[nrow(pre_perturb), -1])
  
  r_perturbation <-  MASS::mvrnorm(n = 1, mu = r_pre, Sigma = covMat*perturb_scale)
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





# Parameters for two contrasting scenarios
alpha_d_values <- c(0.001, 40)  # Two different alpha_d values
cov_scale <- 1  # Fixed scale since we want to focus on different alpha_d scenarios
S <- 8  # Number of species set to 8
nsim <- 50  # Number of simulations
maxt <- 1000

# Initialize data frames to store results
results <- data.frame()
delta_r_values <- data.frame(species = character(0), alpha_d = numeric(0), delta_r = numeric(0))
cor_values <- data.frame(alpha_d = numeric(0), correlation = numeric(0))
cov_values <- data.frame(alpha_d = numeric(0), covariance = numeric(0))

for (alpha_d in alpha_d_values) {
  print(paste("Running simulation for alpha_d =", alpha_d))
  
  for (j in 1:nsim) {
    result <- simulate_response(S = S, 
                                C = 0.2, 
                                aij_params = c(0, 0.5), 
                                sd_X = rep(1, times = S), 
                                maxt = maxt, 
                                alpha_d = alpha_d,
                                perturb_scale = 1)
    
    df <- result$df
    delta_r <- result$delta_r  
    cor_matrix <- result$C_matrix  
    cov_matrix <- result$Sigma  
    
    results <- rbind(
      results,
      data.frame(
        alpha_d = alpha_d,
        sum_deltaX = sum(df$X_pre - df$X_post),
        sd_deltaX = sd(df$X_pre - df$X_post),
        var_deltaX = var(df$X_pre - df$X_post)
      )
    )
    
    delta_r_df <- data.frame(
      alpha_d = rep(alpha_d, length(delta_r)),
      delta_r = delta_r,
      species = paste0("sp", seq_len(S))
    )
    
    n_upper_tri <- S * (S - 1) / 2  
    
    cor_df <- data.frame(
      alpha_d = rep(alpha_d, n_upper_tri),
      correlation = cor_matrix[upper.tri(cor_matrix)]
    )
    
    cov_df <- data.frame(
      alpha_d = rep(alpha_d, n_upper_tri),
      covariance = cov_matrix[upper.tri(cov_matrix)]
    )
    
    delta_r_values <- rbind(delta_r_values, delta_r_df)
    cor_values <- rbind(cor_values, cor_df)
    cov_values <- rbind(cov_values, cov_df)
  }
}

# Plot results with separate panels for each alpha_d
p1 <- ggplot(results, aes(x = as.factor(alpha_d), y = sum_deltaX, fill = as.factor(alpha_d))) +
  geom_boxplot() +
  scale_fill_manual(values = c("#6A0DAD", "#EAD7F5")) +  # Colors for the two alpha_d values
  labs(title = "Sum of Biomass Changes Across alpha_d Scenarios", x = "alpha_d", y = "Sum of Biomass Change") +
  theme_minimal() +
  theme(legend.position = "none")

p2 <- ggplot(results, aes(x = as.factor(alpha_d), y = sd_deltaX, fill = as.factor(alpha_d))) +
  geom_boxplot() +
  scale_fill_manual(values = c("#6A0DAD", "#EAD7F5")) +
  labs(title = "Standard Deviation of Biomass Changes Across alpha_d Scenarios", x = "alpha_d", y = "Standard Deviation of Biomass Change") +
  theme_minimal() +
  theme(legend.position = "none")

delta_r_limits <- range(delta_r_values$delta_r, na.rm = TRUE)
correlation_limits <- range(cor_values$correlation, na.rm = TRUE)
covariance_limits <- range(cov_values$covariance, na.rm = TRUE)

p3 <- ggplot(delta_r_values, aes(x = delta_r)) +
  geom_histogram(bins = 30, fill = "blue", alpha = 0.7) +
  facet_wrap(~ alpha_d, scales = "free") +
  xlim(delta_r_limits) +
  theme_minimal() +
  ggtitle("Distribution of delta_r Across alpha_d Scenarios") +
  labs(x = "delta_r", y = "Frequency")

p4 <- ggplot(cor_values, aes(x = correlation)) +
  geom_histogram(bins = 30, fill = "green", alpha = 0.7) +
  facet_wrap(~ alpha_d, scales = "free") +
  xlim(correlation_limits) +
  theme_minimal() +
  ggtitle("Distribution of Correlation Values Across alpha_d Scenarios") +
  labs(x = "Correlation", y = "Frequency")

p5 <- ggplot(cov_values, aes(x = covariance)) +
  geom_histogram(bins = 30, fill = "red", alpha = 0.7) +
  facet_wrap(~ alpha_d, scales = "free") +
  xlim(covariance_limits) +
  theme_minimal() +
  ggtitle("Distribution of Covariance Values Across alpha_d Scenarios") +
  labs(x = "Covariance", y = "Frequency")

# Arrange all plots
gridExtra::grid.arrange(p1, p2, p3, p4, p5, ncol = 2)

