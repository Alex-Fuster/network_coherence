# Initialize output dataframe
out <- data.frame()

# Define the range of C_sd values to test
C_sd_values <- seq(0.1, 3, by = 0.2)  # C_sd from 0.1 to 1, in steps of 0.1

C_mean = 0

# Run the simulations
for (C_sd in C_sd_values) {
  for (i in 1:50) {
    # Run the simulation for each combination of C_sd and iteration
    result <- simulate_response(params$S, 
                                params$C, 
                                params$aij_params, 
                                params$mu_delta_r, 
                                params$sd_delta_r, 
                                params$sd_X, 
                                params$maxt,
                                C_mean = 0,
                                C_sd = C_sd)
    
    # Store the results
    out <- rbind(
      out,
      data.frame(C_sd = C_sd, 
                 sum_deltaX = sum(result$df$X_pre - result$df$X_post), 
                 sd_deltaX = sd(result$df$X_pre - result$df$X_post))
    )
  }
}

# Plot the results for biomass change (sum_deltaX)
p4 <- ggplot(out, aes(x = as.factor(C_sd), y = sum_deltaX, fill = as.factor(C_sd))) +
  geom_boxplot() +
  scale_fill_brewer(palette = "RdBu", direction = -1) +
  labs(title = "Sum of Biomass Changes Across Different C_sd Values", 
       x = "C_sd", 
       y = "Sum of Biomass Changes") +
  theme_minimal() +
  my_theme

# Plot the results for standard deviation in biomass change (sd_deltaX)
p5 <- ggplot(out, aes(x = as.factor(C_sd), y = sd_deltaX, fill = as.factor(C_sd))) +
  geom_boxplot() +
  scale_fill_brewer(palette = "RdBu", direction = -1) +
  labs(title = "Standard Deviation of Biomass Changes Across Different C_sd Values", 
       x = "C_sd", 
       y = "Standard Deviation of Biomass Changes") +
  theme_minimal() +
  my_theme

# Display the plots
print(p4)
print(p5)