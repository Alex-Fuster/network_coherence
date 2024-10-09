# Load necessary libraries
library(ggplot2)
library(deSolve)
library(Matrix)
library(dplyr)
library(tidyr)
library(ggpubr)


# pars for plotting


my_theme<-theme(axis.text=element_text(size=12),
                axis.title = element_text(size = 14),
                legend.text=element_text(size=10),
                legend.title = element_text(size=12),
                plot.title = element_text(face="bold",size=14,margin=margin(0,0,20,0),hjust = 0.5),
                axis.title.y = element_text(hjust = 0.5),
                axis.title.x = element_text(hjust = 0.5))

# Define parameters for dynamic simulation
params <- list(
  S = 20, # number of species
  C = 0.2, # connectance of the food web
  aij_params = c(0, 0.5), # minimum and maximum interspecific effect of predation
  mu_delta_r = 0, # mean responses to perturbation
  sd_delta_r = 0.5, # sd of response to perturbation
  sd_X = rep(0.1, 20), # standard deviations of responses for each species
  maxt = 100 # maximum time for dynamics
)

# Function to create a correlation matrix with specified mean and sd
createCorMat <- function(S, target_mean = C_mean, target_sd = C_sd){
  # Define the number of off-diagonal elements
  n_off_diag <- S * (S - 1) / 2
  
  cor_values <- runif(n_off_diag, -1, 1)
  
  # Adjust mean and sd to match the target
  cor_values <- scale(cor_values, center = mean(cor_values), scale = sd(cor_values))  # Standardize values
  cor_values <- cor_values * C_sd + C_mean  # Rescale to desired mean and sd
  
  # Ensure values stay between -1 and 1
  cor_values[cor_values < -1] <- -1
  cor_values[cor_values > 1] <- 1
  
  # Create a symmetric matrix with the correlation values
  corMat <- matrix(1, nrow = S, ncol = S)
  corMat[upper.tri(corMat)] <- cor_values
  corMat <- corMat * t(corMat)  # Make the matrix symmetric
  
  
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
simulate_response <- function(S, C, aij_params, mu_delta_r, sd_delta_r, sd_X, maxt, C_mean, C_sd){
  
  # compute C matrix
  cor_matrix <- createCorMat(S, target_mean = C_mean, target_sd = C_sd)
  
  # adjust C matrix
  
  cor_matrix <- as.data.frame(cor_matrix)
  cor_matrix$Var1 = paste0("V", rownames(cor_matrix))
  cor_matrix = cor_matrix %>%
    pivot_longer(-Var1, names_to = "Var2", values_to = "value") #|> na.omit()
  cor_matrix$Var1 = factor(cor_matrix$Var1, levels = unique(cor_matrix$Var1))
  cor_matrix$Var2 = factor(cor_matrix$Var2, levels = rev(unique(cor_matrix$Var2)))
  
  ## Pivot the long-format dataframe back to wide format
  wide_cor_matrix <- cor_matrix %>%
    pivot_wider(names_from = Var2, values_from = value)
  
  ## Remove the row identifier column (Var1) to convert it back to a matrix
  wide_cor_matrix <- wide_cor_matrix %>%
    select(-Var1)
  
  ## Convert the dataframe back to a matrix
  recovered_matrix <- as.matrix(wide_cor_matrix)
  
  
  
  covMat <- as.matrix(Matrix::nearPD(sd_X %*% t(sd_X) * recovered_matrix)$mat)
  A <- sim_quantitative_network("predator-prey", S = S, C = C, aij_params = aij_params)
  
  # Simulate pre-perturbation equilibrium
  equilibrium_pre <- runif(S, 1, 10)
  r_pre <- -A %*% equilibrium_pre
  dyn_params <- list(A = A, r = r_pre, S = S, maxt = maxt) 
  pre_perturb <- simulate_dynamics_c(dyn_params, fw.model)
  equilibrium_pre <- as.numeric(pre_perturb[nrow(pre_perturb), -1])
  
  # Calculate response to perturbation
  delta_r <- MASS::mvrnorm(n = 1, mu = rnorm(S, mu_delta_r, sd_delta_r), Sigma = covMat, empirical = FALSE)
  dyn_params$r <- dyn_params$r + delta_r 
  
  # Simulate post-perturbation equilibrium
  post_perturb <- simulate_dynamics_c(dyn_params, fw.model, init_biomass = equilibrium_pre)
  equilibrium_post <- as.numeric(post_perturb[nrow(post_perturb), -1])
  
  df <- data.frame(species = paste0("sp", c(1:S)), X_pre = equilibrium_pre, X_post = equilibrium_post)
  return(list(df = df, pre_perturb = pre_perturb, post_perturb = post_perturb, Sigma = covMat, C_matrix = recovered_matrix))
}




########################################################## DEMONSTRATION

# Run the simulation for demonstration
C_mean = 0
C_sd = 0.3


result <- simulate_response(params$S, 
                            params$C, 
                            params$aij_params, 
                            params$mu_delta_r, 
                            params$sd_delta_r, 
                            params$sd_X, 
                            params$maxt,
                            C_mean = C_mean,
                            C_sd = C_sd)

while (all(result$post_perturb[nrow(result$post_perturb),]>0)) {
  result <- simulate_response(params$S, 
                              params$C, 
                              params$aij_params, 
                              params$mu_delta_r, 
                              params$sd_delta_r, 
                              params$sd_X, 
                              params$maxt,
                              C_mean = C_mean,
                              C_sd = C_sd)
}


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
# p3 <- result$Sigma %>%
#   as_tibble() %>%
#   rownames_to_column("Var1") %>%
#   pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
#   ggplot(aes(Var1, Var2)) +
#   geom_tile(aes(fill = value)) +
#   scale_fill_distiller(palette = "RdBu", direction = -1, limits = c(-1,1)*max(abs(result$Sigma))) +
#   ggtitle("Covariance Matrix")

#print(p3)

########################################################## SIMULATION

# Save necessary data to make these plots
out <- data.frame()



params$mu_delta_r <- 0
params$sd_delta_r <- 0.1

C_mean = 0
C_sd = 0.3

for (i in 1:100){
  X <- with(params, simulate_response(params$S, 
                                      params$C, 
                                      params$aij_params, 
                                      params$mu_delta_r, 
                                      params$sd_delta_r, 
                                      params$sd_X, 
                                      params$maxt,
                                      C_mean = C_mean,
                                      C_sd = C_sd)$df)
  out <- rbind(
    out,
    data.frame(covtype = "mixed", perturbation = "weak", sum_deltaX = sum(X$X_pre - X$X_post), sd_deltaX = sd(X$X_pre - X$X_post), X_pre = X$X_pre, X_post = X$X_post)
  )
}

params$mu_delta_r <- 0
params$sd_delta_r <- 1
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


################################################################
# SIMULATING EFFECTS OF MEAN AND SD FROM NORMALLY DISTRIBUTED ENC
################################################################

# Prepare to store results
out <- data.frame()

# Define the range of values to test
mu_delta_r_values <- seq(-3, 3, by = 1)    # mu_delta_r from -15 to 15 by 1
sd_delta_r_values <- seq(0, 3, by = 0.5)     # sd_delta_r from 0 to 5 by 0.2

# Loop through each combination of mu_delta_r and sd_delta_r
for (mu_delta_r in mu_delta_r_values) {
  
  print(paste(mu_delta_r, "to", length(mu_delta_r_values), "mu values"))
  
  for (sd_delta_r in sd_delta_r_values) {
    # Update parameters for the current scenario
    params$covMatrix_type <- "mixed"
    params$mu_delta_r <- mu_delta_r
    params$sd_delta_r <- sd_delta_r
    
    print(paste(sd_delta_r, "to", length(sd_delta_r_values), "sd values"))
    
    # Run 100 simulations for each scenario
    for (i in 1:100) {
      X <- with(params, simulate_response(S, C, aij_params, mu_delta_r, sd_delta_r, covMatrix_type, sd_X, maxt)$df)
      
      # Store the results
      out <- rbind(
        out,
        data.frame(
          mu_delta_r = mu_delta_r,
          sd_delta_r = sd_delta_r,
          sum_deltaX = sum(X$X_pre - X$X_post), 
          sd_deltaX = sd(X$X_pre - X$X_post)
        )
      )
    }
  }
}

# Calculate the mean sum_deltaX and sd_deltaX for each scenario
results_summary <- out %>%
  group_by(mu_delta_r, sd_delta_r) %>%
  summarise(mean_sum_deltaX = mean(sum_deltaX), mean_sd_deltaX = mean(sd_deltaX))




# Find the maximum absolute value in mean_sd_deltaX
max_abs_value_sum <- max(abs(results_summary$mean_sum_deltaX), na.rm = TRUE)
max_abs_value_sd <- max(abs(results_summary$mean_sd_deltaX), na.rm = TRUE)

# Define the RdBu palette with continuous colors
colors <- colorRampPalette(brewer.pal(11, "RdBu"))(100)  # Create a continuous color gradient

# Plot using scale_fill_gradientn with a symmetric color scale around 0

p1 <- ggplot(results_summary, aes(x = sd_delta_r, y = mu_delta_r, fill = mean_sum_deltaX)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = colors,            # Use the continuous RdBu colors
    limits = c(-max_abs_value_sum, max_abs_value_sum),  # Set symmetric limits around 0
    values = scales::rescale(c(-max_abs_value_sum, 0, max_abs_value_sum)),  # Rescale colors to center 0
    na.value = "white"          # Set NA values to white
  ) +
  labs(title = "Mean sum_deltaX Across Scenarios",
       x = "sd_delta_r",
       y = "mu_delta_r") +
  theme_minimal()


p2 <- ggplot(results_summary, aes(x = sd_delta_r, y = mu_delta_r, fill = mean_sd_deltaX)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = colors,            # Use the continuous RdBu colors
    limits = c(-max_abs_value_sd, max_abs_value_sd),  # Set symmetric limits around 0
    values = scales::rescale(c(-max_abs_value_sd, 0, max_abs_value_sd)),  # Rescale colors to center 0
    na.value = "white"          # Set NA values to white
  ) +
  labs(title = "Mean sd_deltaX Across Scenarios",
       x = "sd_delta_r",
       y = "mu_delta_r") +
  theme_minimal()

# Display the plot

ggarrange(p1,
          p2,
          ncol = 1,
          nrow = 2)


########################################################################
# SAME BUT CREATING ENC PLOTS
########################################################################

# Function to compute and plot ENC pattern for each scenario
plot_enc_pattern <- function(a_matrix, c_matrix, scenario_name) {
  # Multiply the A matrix by the C matrix to get the ENC matrix
  enc_matrix <- a_matrix * c_matrix
  
  # Convert the matrices to data frames
  enc_matrix_df <- as.data.frame(enc_matrix)
  a_matrix_df <- as.data.frame(a_matrix)
  
  # Set ENC values to NA where A matrix has value 0
  enc_matrix_df[a_matrix_df == 0] <- NA
  
  # Melt the dataframe to long format
  enc_matrix_melt <- enc_matrix_df %>%
    rownames_to_column(var = "Var1") %>%
    gather(key = "Var2", value = "value", -Var1) %>%
    drop_na(value)
  
  # Determine the limits for the x-axis
  min_enc <- min(enc_matrix_melt$value, na.rm = TRUE)
  max_enc <- max(enc_matrix_melt$value, na.rm = TRUE)
  
  # Create the plot
  p <- ggplot(enc_matrix_melt, aes(x = value)) +
    geom_histogram(aes(fill = ..x..), bins = 30, alpha = 0.7) +
    scale_fill_distiller(palette = "RdBu", direction = 1) +
    scale_x_continuous(limits = c(min_enc, max_enc)) +
    labs(x = "Co-response", y = "Frequency", title = paste("ENC for", scenario_name)) +
    theme_classic() +
    theme(legend.position = "bottom", legend.title = element_blank(),
          plot.title = element_text(hjust = 0.5)) +
    my_theme
  
  return(p)
}

# Initialize a list to store plots
plot_list <- list()

# Define the range of values to test
mu_delta_r_values <- seq(-3, 3, by = 1)    # mu_delta_r from -3 to 3 by 1
sd_delta_r_values <- seq(0, 3, by = 0.5)   # sd_delta_r from 0 to 3 by 0.5

# Loop through each combination of mu_delta_r and sd_delta_r
for (mu_delta_r in mu_delta_r_values) {
  
  print(paste(mu_delta_r, "to", length(mu_delta_r_values), "mu values"))
  
  for (sd_delta_r in sd_delta_r_values) {
    # Update parameters for the current scenario
    params$covMatrix_type <- "mixed"
    params$mu_delta_r <- mu_delta_r
    params$sd_delta_r <- sd_delta_r
    
    print(paste(sd_delta_r, "to", length(sd_delta_r_values), "sd values"))
    
    # Run 100 simulations for each scenario
    for (i in 1:100) {
      # Simulate response for the current parameters
      response <- with(params, simulate_response(S, C, aij_params, mu_delta_r, sd_delta_r, covMatrix_type, sd_X, maxt))
      X <- response$df
      a_matrix <- response$pre_perturb  # Example: Use pre-perturbation A matrix
      c_matrix <- response$Sigma  # Example: Use covariance matrix
      
      # Store the results
      out <- rbind(
        out,
        data.frame(
          mu_delta_r = mu_delta_r,
          sd_delta_r = sd_delta_r,
          sum_deltaX = sum(X$X_pre - X$X_post), 
          sd_deltaX = sd(X$X_pre - X$X_post)
        )
      )
    }
    
    # After simulations, generate ENC plot for the current scenario
    scenario_name <- paste("mu_delta_r:", mu_delta_r, "sd_delta_r:", sd_delta_r)
    enc_plot <- plot_enc_pattern(a_matrix, c_matrix, scenario_name)
    
    # Save the plot in the list
    plot_list[[scenario_name]] <- enc_plot
  }
}

# Example to print one plot from the list
print(plot_list[["mu_delta_r: -3 sd_delta_r: 1"]])

print(plot_list[["mu_delta_r: 0 sd_delta_r: 1"]])

print(plot_list[["mu_delta_r: 3 sd_delta_r: 1"]])

print(plot_list[["mu_delta_r: -3 sd_delta_r: 3"]])
