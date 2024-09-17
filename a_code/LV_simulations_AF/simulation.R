# Load necessary libraries
library(ggplot2)
library(deSolve)
library(Matrix)
library(dplyr)
library(tidyr)

# Define parameters for dynamic simulation
params <- list(
  S = 10, # number of species
  C = 0.2, # connectance of the food web
  aij_params = c(0, 0.5), # minimum and maximum interspecific effect of predation
  mu_delta_r = 0, # mean responses to perturbation
  sd_delta_r = 0.5, # sd of response to perturbation
  covMatrix_type = "positive", # type of covariance (mixed or all positive)
  sd_X = rep(1, 10), # standard deviations of responses for each species
  maxt = 100 # maximum time for dynamics
)

# Function to create a correlation matrix
createCorMat <- function(S, type){
  if (type == "positive"){
    cor <- runif((S*(S-1)/2), 0, 1)
  } else if (type == "mixed"){
    cor <- runif((S*(S-1)/2), -1, 1)
  } else {
    stop("type not defined")
  }
  corMat <- matrix(1, nrow = S, ncol = S)
  corMat[upper.tri(corMat)] <- cor
  corMat <- corMat * t(corMat)
  return(corMat)
}

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

# Function to simulate response
simulate_response <- function(S, C, aij_params, mu_delta_r, sd_delta_r, covMatrix_type, sd_X, maxt){
  corMat <- createCorMat(S, covMatrix_type)
  covMat <- as.matrix(Matrix::nearPD(sd_X %*% t(sd_X) * corMat)$mat)
  A <- sim_quantitative_network("predator-prey", S = S, C = C, aij_params = aij_params)
  
  # Simulate pre-perturbation equilibrium
  dyn_params <- list(A = A, r = rep(1, S), S = S, maxt = maxt) 
  pre_perturb <- simulate_dynamics_c(dyn_params, fw.model)
  equilibrium_pre <- as.numeric(pre_perturb[nrow(pre_perturb), -1])
  
  # Calculate response to perturbation
  delta_r <- MASS::mvrnorm(n = 1, mu = rnorm(S, mu_delta_r, sd_delta_r), Sigma = covMat, empirical = FALSE)
  dyn_params$r <- dyn_params$r + delta_r 
  
  # Simulate post-perturbation equilibrium
  post_perturb <- simulate_dynamics_c(dyn_params, fw.model, init_biomass = equilibrium_pre)
  equilibrium_post <- as.numeric(post_perturb[nrow(post_perturb), -1])
  
  df <- data.frame(species = paste0("sp", c(1:S)), X_pre = equilibrium_pre, X_post = equilibrium_post)
  return(list(df = df, pre_perturb = pre_perturb, post_perturb = post_perturb, Sigma = covMat))
}




########################################################## DEMONSTRATION

# Run the simulation for demonstration
result <- simulate_response(params$S, params$C, params$aij_params, params$mu_delta_r, params$sd_delta_r, params$covMatrix_type, params$sd_X, params$maxt)





# Plot species biomass dynamics
pre_perturb_long <- result$pre_perturb %>%
  pivot_longer(-time, names_to = "species", values_to = "biomass")

post_perturb_long <- result$post_perturb %>%
  pivot_longer(-time, names_to = "species", values_to = "biomass")

p1 <- ggplot(pre_perturb_long, aes(x = time, y = biomass, color = species)) +
  geom_line() +
  ggtitle("Pre-Perturbation Biomass Dynamics")

p2 <- ggplot(post_perturb_long, aes(x = time, y = biomass, color = species)) +
  geom_line() +
  ggtitle("Post-Perturbation Biomass Dynamics")

print(p1)
print(p2)

# Plot the covariance matrix
p3 <- result$Sigma %>%
  as_tibble() %>%
  rownames_to_column("Var1") %>%
  pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
  ggplot(aes(Var1, Var2)) +
  geom_tile(aes(fill = value)) +
  scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-1,1)*max(abs(result$Sigma))) +
  ggtitle("Covariance Matrix")

print(p3)


########################################################## SIMULATION

# Save necessary data to make these plots
out <- data.frame()

params$mu_delta_r <- 0
params$sd_delta_r <- 0.1
for (i in 1:100){
  X <- with(params, simulate_response(S, C, aij_params, mu_delta_r, sd_delta_r, covMatrix_type, sd_X, maxt)$df)
  out <- rbind(
    out,
    data.frame(covtype = "mixed", perturbation = "weak", sum_deltaX = sum(X$X_pre - X$X_post), sd_deltaX = sd(X$X_pre - X$X_post), X_pre = X$X_pre, X_post = X$X_post)
  )
}

params$mu_delta_r <- 0
params$sd_delta_r <- 0.5
for (i in 1:100){
  X <- with(params, simulate_response(S, C, aij_params, mu_delta_r, sd_delta_r, covMatrix_type, sd_X, maxt)$df)
  out <- rbind(
    out,
    data.frame(covtype = "mixed", perturbation = "strong", sum_deltaX = sum(X$X_pre - X$X_post), sd_deltaX = sd(X$X_pre - X$X_post), X_pre = X$X_pre, X_post = X$X_post)
  )
}

# Plot results of perturbation scenarios
p4 <- ggplot(out) +
  geom_boxplot(aes(x = covtype, y = sum_deltaX, fill = perturbation)) +
  lims(y = c(-100, 100)) +
  ggtitle("Sum of Biomass Changes Across Scenarios")

p5 <- ggplot(out) +
  geom_boxplot(aes(x = covtype, y = sd_deltaX, fill = perturbation)) +
  lims(y = c(0, 25)) +
  ggtitle("Standard Deviation of Biomass Changes Across Scenarios")

# Display scenario plots
print(p4)
print(p5)


########################################

