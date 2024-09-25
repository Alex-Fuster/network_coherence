

# Set up parameters for simulation
params <- list(
  S = 20, # Number of species
  C = 0.2, # Connectance (fraction of possible links)
  aij_params = c(0, 0.5), # Range of interaction strengths
  maxt = 100 # Maximum time for the simulation
)


# Function to create a correlation matrix with different distributions
createCorMat <- function(S, distribution_type = "normal", 
                         mean = 0, sd = 1, 
                         beta_shape1 = NULL, beta_shape2 = NULL, 
                         ushape1 = NULL, ushape2 = NULL) {
  n_off_diag <- S * (S - 1) / 2
  
  # Generate correlation values based on the chosen distribution
  if (distribution_type == "normal") {
    cor_values <- rnorm(n_off_diag, mean = mean, sd = sd)
  } else if (distribution_type == "beta") {
    cor_values <- rbeta(n_off_diag, shape1 = beta_shape1, shape2 = beta_shape2)
    cor_values <- 2 * (cor_values - 0.5)  # Map Beta values to [-1, 1]
  } else if (distribution_type == "u_shaped") {
    cor_values <- rbeta(n_off_diag, shape1 = ushape1, shape2 = ushape2)
    cor_values <- 2 * (cor_values - 0.5)  # Map Beta values to [-1, 1]
  } else {
    stop("Invalid distribution type specified.")
  }
  
  if (distribution_type == "normal") {
    cor_values <- scale(cor_values, center = mean(cor_values), scale = sd(cor_values))
    cor_values <- cor_values * sd + mean
  }
  
  cor_values[cor_values < -1] <- -1
  cor_values[cor_values > 1] <- 1
  
  corMat <- matrix(1, nrow = S, ncol = S)
  corMat[upper.tri(corMat)] <- cor_values
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
simulate_response <- function(S, C, aij_params, mu_delta_r, sd_delta_r, sd_X, maxt, distribution_type, mean, sd, beta_shape1, beta_shape2, ushape1, ushape2){
  
  # Compute C matrix
  cor_matrix <- createCorMat(S = S, 
                             distribution_type = distribution_type, 
                             mean = mean, sd = sd, 
                             beta_shape1 = beta_shape1, beta_shape2 = beta_shape2, 
                             ushape1 = ushape1, ushape2 = ushape2)
  
  # Adjust C matrix
  cor_matrix <- as.data.frame(cor_matrix)
  cor_matrix$Var1 <- paste0("V", rownames(cor_matrix))
  cor_matrix <- cor_matrix %>%
    pivot_longer(-Var1, names_to = "Var2", values_to = "value")
  cor_matrix$Var1 <- factor(cor_matrix$Var1, levels = unique(cor_matrix$Var1))
  cor_matrix$Var2 <- factor(cor_matrix$Var2, levels = rev(unique(cor_matrix$Var2)))
  
  # Pivot the long-format dataframe back to wide format
  wide_cor_matrix <- cor_matrix %>%
    pivot_wider(names_from = Var2, values_from = value) %>%
    select(-Var1)
  
  # Convert back to a matrix
  recovered_matrix <- as.matrix(wide_cor_matrix)
  
  #covMat <- as.matrix(Matrix::nearPD(sd_X %*% t(sd_X) * recovered_matrix)$mat)
  covMat <- diag(sd_X)%*%recovered_matrix%*%diag(sd_X)
  covMat <- as.matrix(Matrix::nearPD(covMat)$mat)
  A <- sim_quantitative_network("predator-prey", S = S, C = C, aij_params = aij_params)
  
  equilibrium_pre <- runif(S, 1, 10)
  r_pre <- -A %*% equilibrium_pre
  dyn_params <- list(A = A, r = r_pre, S = S, maxt = maxt) 
  pre_perturb <- simulate_dynamics_c(dyn_params, fw.model)
  equilibrium_pre <- as.numeric(pre_perturb[nrow(pre_perturb), -1])
  
  delta_r <- MASS::mvrnorm(n = 1, mu = rnorm(S, mu_delta_r, sd_delta_r), Sigma = covMat, empirical = FALSE)
  # delta_r <- MASS::mvrnorm(n = 1, mu = rep(0, S), Sigma = covMat, empirical = FALSE)
  dyn_params$r <- dyn_params$r + delta_r 
  
  post_perturb <- simulate_dynamics_c(dyn_params, fw.model, init_biomass = equilibrium_pre)
  equilibrium_post <- as.numeric(post_perturb[nrow(post_perturb), -1])
  
  df <- data.frame(species = paste0("sp", c(1:S)), X_pre = equilibrium_pre, X_post = equilibrium_post)
  return(list(df = df, pre_perturb = pre_perturb, post_perturb = post_perturb, Sigma = covMat, C_matrix = recovered_matrix, delta_r = delta_r))
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



# Function to simulate dynamics over time and compute species-specific standard deviations
compute_species_sd <- function(params, n_steps = 1000, perturbation_sd = 0.1) {
  biomass_over_time <- simulate_dynamics_with_perturbations(params, perturbation_sd = perturbation_sd, n_steps = n_steps)
  species_sd <- apply(biomass_over_time, 2, sd)  # Compute the standard deviation for each species
  return(species_sd)
}



# Step 1: Create a dataframe for scenarios with the correct classification
scenarios <- data.frame(
  scenario = c("Normal weak", "Normal strong", 
               "Beta + weak", "Beta + strong", 
               "Beta - weak", "Beta - strong", 
               "U-shaped weak", "U-shaped strong"),
  distribution = c("normal", "normal", 
                   "beta", "beta", 
                   "beta", "beta", 
                   "u_shaped", "u_shaped"),
  strength = c("weak", "strong", 
               "weak_positive", "strong_positive", 
               "weak_negative", "strong_negative", 
               "weak", "strong"),
  mean = c(0, 0, NA, NA, NA, NA, NA, NA),
  sd = c(0.1, 0.4, NA, NA, NA, NA, NA, NA),
  beta_shape1 = c(NA, NA, 3, 5, 2, 2, NA, NA),
  beta_shape2 = c(NA, NA, 2, 2, 3, 5, NA, NA),
  ushape1 = c(NA, NA, NA, NA, NA, NA, 0.5, 0.05),
  ushape2 = c(NA, NA, NA, NA, NA, NA, 0.5, 0.05)
)


# Step 2: Run simulation and store the results for each scenario
results <- data.frame()
delta_r_values <- data.frame()  # To store delta_r values
correlation_plots <- list()
covariance_plots <- list()
biomass_plots <- list()

nsim <- 15
n_steps <- 1000

for (i in 1:nrow(scenarios)) {
  scenario <- scenarios[i, ]
  print(paste("Running scenario", scenario$scenario, ":", scenario$distribution, scenario$strength))
  
  for (j in 1:nsim) {
    # Step 1: Simulate biomass dynamics and compute sd_X (species-specific standard deviations)
    sd_X <- compute_species_sd(params, n_steps = n_steps, perturbation_sd = 0.1)
    
    # Step 2: Run response simulation with computed sd_X
    result <- simulate_response(S = 20, 
                                C = 0.2, 
                                aij_params = c(0, 0.5), 
                                mu_delta_r = 0, 
                                sd_delta_r = 0.5, 
                                sd_X = sd_X,  # Use species-specific standard deviations
                                maxt = 100,
                                distribution_type = scenario$distribution,
                                mean = scenario$mean,
                                sd = scenario$sd,
                                beta_shape1 = scenario$beta_shape1,
                                beta_shape2 = scenario$beta_shape2,
                                ushape1 = scenario$ushape1,
                                ushape2 = scenario$ushape2)
    
    df <- result$df
    delta_r <- result$delta_r  # Capture delta_r values
    
    # Store results for each simulation
    results <- rbind(
      results,
      data.frame(
        scenario = scenario$scenario,
        distribution = scenario$distribution,
        strength = scenario$strength,
        sum_deltaX = sum(df$X_pre - df$X_post),
        sd_deltaX = sd(df$X_pre - df$X_post)
      )
    )
    
    # Store delta_r values for later analysis
    delta_r_values <- rbind(
      delta_r_values,
      data.frame(scenario = scenario$scenario, delta_r = delta_r)
    )
    
    # Randomly pick one simulation to visualize correlation, covariance matrices, and biomass dynamics
    if (j == 1) {  # Only visualize the first run
      # 1. Plot Correlation Matrix
      cor_matrix <- result$C_matrix
      cor_matrix_long <- cor_matrix %>%
        as.data.frame() %>%
        rownames_to_column("Var1") %>%
        pivot_longer(-Var1, names_to = "Var2", values_to = "value")
      
      cor_plot <- ggplot(cor_matrix_long, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile() +
        scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-1, 1)) +
        ggtitle(paste("Correlation Matrix -", scenario$scenario)) +
        theme_minimal()
      correlation_plots[[scenario$scenario]] <- cor_plot
      
      # 2. Plot Covariance Matrix
      cov_matrix <- result$Sigma
      cov_matrix_long <- cov_matrix %>%
        as.data.frame() %>%
        rownames_to_column("Var1") %>%
        pivot_longer(-Var1, names_to = "Var2", values_to = "value")
      
      cov_plot <- ggplot(cov_matrix_long, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile() +
        scale_fill_distiller(palette = "RdBu", direction = -1) +
        ggtitle(paste("Covariance Matrix -", scenario$scenario)) +
        theme_minimal()
      covariance_plots[[scenario$scenario]] <- cov_plot
      
      # 3. Plot Biomass Dynamics (Pre and Post Perturbation)
      pre_perturb <- result$pre_perturb
      post_perturb <- result$post_perturb
      
      # Combine pre- and post-perturbation dynamics
      pre_perturb$phase <- "Pre-perturbation"
      post_perturb$phase <- "Post-perturbation"
      pre_perturb$time <- 1:nrow(pre_perturb)
      post_perturb$time <- 1:nrow(post_perturb)
      
      biomass_df <- rbind(pre_perturb, post_perturb) %>%
        tidyr::pivot_longer(-c(phase, time), names_to = "species", values_to = "biomass")
      
      biomass_plot <- ggplot(biomass_df, aes(x = time, y = biomass, color = species)) +
        geom_line() +
        facet_wrap(~phase) +
        ggtitle(paste("Biomass Dynamics -", scenario$scenario)) +
        theme_minimal()
      biomass_plots[[scenario$scenario]] <- biomass_plot
    }
  }
}

# Step 3: Plot the results organized by distribution and strength

# Plot for sum of biomass changes
p4 <- ggplot(results, aes(x = scenario, y = sum_deltaX, fill = distribution)) +
  geom_boxplot() +
  labs(title = "Sum of Biomass Changes Across Scenarios", x = "Scenario", y = "Sum of Biomass Change") +
  scale_x_discrete(limits = c("Normal weak", "Normal strong", 
                              "Beta + weak", "Beta + strong", 
                              "Beta - weak", "Beta - strong", 
                              "U-shaped weak", "U-shaped strong")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot for standard deviation of biomass changes
p5 <- ggplot(results, aes(x = scenario, y = sd_deltaX, fill = distribution)) +
  geom_boxplot() +
  labs(title = "Standard Deviation of Biomass Changes Across Scenarios", x = "Scenario", y = "Standard Deviation of Biomass Change") +
  scale_x_discrete(limits = c("Normal weak", "Normal strong", 
                              "Beta + weak", "Beta + strong", 
                              "Beta - weak", "Beta - strong", 
                              "U-shaped weak", "U-shaped strong")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(0,20)

print(p4)
print(p5)

# Step 4: Visualize one simulation's correlation, covariance, and biomass dynamics
print(correlation_plots[["Normal weak"]])  # Example for "Normal weak"
print(covariance_plots[["Normal weak"]])
print(biomass_plots[["Normal weak"]])


############################


# Set up parameters for simulation
params <- list(
  S = 20, # Number of species
  C = 0.2, # Connectance (fraction of possible links)
  aij_params = c(0, 0.5), # Range of interaction strengths
  maxt = 100 # Maximum time for the simulation
)

# Step 2: Run simulation and store the results for each scenario
results <- data.frame()
delta_r_values <- data.frame()  # To store delta_r values
correlation_plots <- list()
covariance_plots <- list()
biomass_plots <- list()

nsim <- 15
n_steps <- 1000

for (i in 1:nrow(scenarios)) {
  scenario <- scenarios[i, ]
  print(paste("Running scenario", scenario$scenario, ":", scenario$distribution, scenario$strength))
  
  for (j in 1:nsim) {
    # Step 1: Create the interaction matrix A before each simulation
    params$A <- sim_quantitative_network("predator-prey", S = params$S, C = params$C, aij_params = params$aij_params)
    
    # Simulate biomass dynamics and compute sd_X (species-specific standard deviations)
    sd_X <- compute_species_sd(params, n_steps = n_steps, perturbation_sd = 0.1)
    
    # Step 2: Run response simulation with computed sd_X
    result <- simulate_response(S = 20, 
                                C = 0.2, 
                                aij_params = c(0, 0.5), 
                                mu_delta_r = 0, 
                                sd_delta_r = 0.5, 
                                sd_X = sd_X,  # Use species-specific standard deviations
                                maxt = 100,
                                distribution_type = scenario$distribution,
                                mean = scenario$mean,
                                sd = scenario$sd,
                                beta_shape1 = scenario$beta_shape1,
                                beta_shape2 = scenario$beta_shape2,
                                ushape1 = scenario$ushape1,
                                ushape2 = scenario$ushape2)
    
    df <- result$df
    delta_r <- result$delta_r  # Capture delta_r values
    
    # Store results for each simulation
    results <- rbind(
      results,
      data.frame(
        scenario = scenario$scenario,
        distribution = scenario$distribution,
        strength = scenario$strength,
        sum_deltaX = sum(df$X_pre - df$X_post),
        sd_deltaX = sd(df$X_pre - df$X_post)
      )
    )
    
    # Store delta_r values for later analysis
    delta_r_values <- rbind(
      delta_r_values,
      data.frame(scenario = scenario$scenario, delta_r = delta_r)
    )
    
    # Randomly pick one simulation to visualize correlation, covariance matrices, and biomass dynamics
    if (j == 1) {  # Only visualize the first run
      # 1. Plot Correlation Matrix
      cor_matrix <- result$C_matrix
      cor_matrix_long <- cor_matrix %>%
        as.data.frame() %>%
        rownames_to_column("Var1") %>%
        pivot_longer(-Var1, names_to = "Var2", values_to = "value")
      
      cor_plot <- ggplot(cor_matrix_long, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile() +
        scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-1, 1)) +
        ggtitle(paste("Correlation Matrix -", scenario$scenario)) +
        theme_minimal()
      correlation_plots[[scenario$scenario]] <- cor_plot
      
      # 2. Plot Covariance Matrix
      cov_matrix <- result$Sigma
      cov_matrix_long <- cov_matrix %>%
        as.data.frame() %>%
        rownames_to_column("Var1") %>%
        pivot_longer(-Var1, names_to = "Var2", values_to = "value")
      
      cov_plot <- ggplot(cov_matrix_long, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile() +
        scale_fill_distiller(palette = "RdBu", direction = -1) +
        ggtitle(paste("Covariance Matrix -", scenario$scenario)) +
        theme_minimal()
      covariance_plots[[scenario$scenario]] <- cov_plot
      
      # 3. Plot Biomass Dynamics (Pre and Post Perturbation)
      pre_perturb <- result$pre_perturb
      post_perturb <- result$post_perturb
      
      # Combine pre- and post-perturbation dynamics
      pre_perturb$phase <- "Pre-perturbation"
      post_perturb$phase <- "Post-perturbation"
      pre_perturb$time <- 1:nrow(pre_perturb)
      post_perturb$time <- 1:nrow(post_perturb)
      
      biomass_df <- rbind(pre_perturb, post_perturb) %>%
        tidyr::pivot_longer(-c(phase, time), names_to = "species", values_to = "biomass")
      
      biomass_plot <- ggplot(biomass_df, aes(x = time, y = biomass, color = species)) +
        geom_line() +
        facet_wrap(~phase) +
        ggtitle(paste("Biomass Dynamics -", scenario$scenario)) +
        theme_minimal()
      biomass_plots[[scenario$scenario]] <- biomass_plot
    }
  }
}

# Step 3: Plot the results organized by distribution and strength

# Plot for sum of biomass changes
p4 <- ggplot(results, aes(x = scenario, y = sum_deltaX, fill = distribution)) +
  geom_boxplot() +
  labs(title = "Sum of Biomass Changes Across Scenarios", x = "Scenario", y = "Sum of Biomass Change") +
  scale_x_discrete(limits = c("Normal weak", "Normal strong", 
                              "Beta + weak", "Beta + strong", 
                              "Beta - weak", "Beta - strong", 
                              "U-shaped weak", "U-shaped strong")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Plot for standard deviation of biomass changes
p5 <- ggplot(results, aes(x = scenario, y = sd_deltaX, fill = distribution)) +
  geom_boxplot() +
  labs(title = "Standard Deviation of Biomass Changes Across Scenarios", x = "Scenario", y = "Standard Deviation of Biomass
