
###############################################
# testing effects of different sd values
###############################################

# Prepare to store results
out <- data.frame()

# Define the range of sd values to test
sd_values <- seq(0, 4, by = 0.2)  # Explore sd from 0 to 4 in steps of 0.2

# Loop through each sd value and run the simulations
for (sd in sd_values) {
  # Set the parameters for mixed covariance scenario
  params$covMatrix_type <- "mixed"
  params$mu_delta_r <- 0
  params$sd_delta_r <- sd
  
  # Run simulations for a range of sd values
  for (i in 1:100) {
    X <- with(params, simulate_response(S, C, aij_params, mu_delta_r, sd_delta_r, covMatrix_type, sd_X, maxt)$df)
    
    # Store the results
    out <- rbind(
      out,
      data.frame(sd = sd, 
                 covtype = "mixed", 
                 perturbation = ifelse(sd <= 0.5, "weak", "strong"), 
                 sum_deltaX = sum(X$X_pre - X$X_post), 
                 sd_deltaX = sd(X$X_pre - X$X_post), 
                 X_pre = X$X_pre, 
                 X_post = X$X_post)
    )
  }
}

# Plot the results to show the effect of different sd values
p1 <- ggplot(out, aes(x = sd, y = sum_deltaX, color = perturbation)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Effect of Different SD Values on Biomass Changes", 
       x = "Standard Deviation (sd)", 
       y = "Sum of Biomass Changes (sum_deltaX)")

p2 <- ggplot(out, aes(x = sd, y = X_post, color = perturbation)) +
  geom_point() +
  geom_smooth(method = "loess") +
  labs(title = "Final Biomass (X_post) Across Different SD Values", 
       x = "Standard Deviation (sd)", 
       y = "Final Biomass (X_post)")
