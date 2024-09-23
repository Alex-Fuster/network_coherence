library(ggpubr)

# Plot biomass dynamics per scenario
nsim <- 3  # number of simulations per scenario
plots <- list()  # list to store plots

for (i in 1:nrow(scenarios)) {
  scenario <- scenarios[i, ]
  
  for (j in 1:nsim) {
    result <- simulate_response(S = 20, 
                                C = 0.2, 
                                aij_params = c(0, 0.5), 
                                mu_delta_r = 0, 
                                sd_delta_r = 0.5, 
                                sd_X = rep(1, 20), 
                                maxt = 100,
                                distribution_type = scenario$distribution,
                                mean = scenario$mean,
                                sd = scenario$sd,
                                beta_shape1 = scenario$beta_shape1,
                                beta_shape2 = scenario$beta_shape2,
                                ushape1 = scenario$ushape1,
                                ushape2 = scenario$ushape2)
    
    # Extract pre- and post-perturbation dynamics
    pre_perturb <- result$pre_perturb
    post_perturb <- result$post_perturb
    
    # Combine pre- and post-perturbation dynamics into a single dataframe
    pre_perturb$phase <- "Pre-perturbation"
    post_perturb$phase <- "Post-perturbation"
    
    # Add time column
    pre_perturb$time <- seq(1, nrow(pre_perturb))
    post_perturb$time <- seq(1, nrow(post_perturb))
    
    # Combine both phases
    dynamics_df <- rbind(pre_perturb, post_perturb)
    
    # Set the phase order explicitly: pre first, post second
    dynamics_df$phase <- factor(dynamics_df$phase, levels = c("Pre-perturbation", "Post-perturbation"))
    
    # Melt the dataframe to long format for ggplot
    dynamics_df_long <- tidyr::pivot_longer(dynamics_df, cols = -c(phase, time), names_to = "species", values_to = "biomass")
    
    # Create plot for biomass dynamics
    p <- ggplot(dynamics_df_long, aes(x = time, y = biomass, color = species)) +
      geom_line() +
      facet_wrap(~phase) +
      labs(title = paste("Biomass Dynamics -", scenario$scenario),
           x = "Time", y = "Biomass") +
      theme_minimal() +
      theme(legend.position = "none")
    
    # Store the plot in the list with correct scenario names
    plots[[paste(scenario$scenario, "Sim", j)]] <- p
  }
}

# Example: Print or save a plot from the list
print(plots[["Normal weak Sim 1"]])

# Arrange the plots for all scenarios
ggarrange(
  plots[["Normal weak Sim 1"]],
  plots[["Normal strong Sim 1"]],
  plots[["Beta + weak Sim 1"]],
  plots[["Beta + strong Sim 1"]],
  plots[["Beta - weak Sim 1"]],
  plots[["Beta - strong Sim 1"]],
  plots[["U-shaped weak Sim 1"]],
  plots[["U-shaped strong Sim 1"]],
  
  ncol = 2,
  nrow = 4
)
