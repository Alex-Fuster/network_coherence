
# Plot ENC frequencies ----------------------------------------------------
library(ggpubr)
library(patchwork)

# Set theme
theme_set(theme_pubr() +
            theme(panel.grid.major.x = element_line()))

# Source code for ENC frequencies
source("a_code/example_bees/06-correlations.R")


# Descriptive stats for interacting and non-int. species ------------------
stats_interactions <- 
interaction_matrix |> 
  filter(weighted_interaction > 0) |>             # Filter only interacting species
  drop_na(correlation) |>                         # Drop NAs
  summarise(mean_correlation = mean(correlation), # Calculate mean and sd of the correlation
            sd_correlation = sd(correlation))


stats_noninteractions <- 
  interaction_matrix |> 
  filter(weighted_interaction == 0) |>            # Filter only non-interacting species
  drop_na(correlation) |>                         # Drop NAs
  summarise(mean_correlation = mean(correlation), # Calculate mean and sd of the correlation
            sd_correlation = sd(correlation))



# Plot for interacting species
correlations_histogram_interacting <- 
  interaction_matrix |> 
  filter(weighted_interaction > 0) |>          # Filter only interacting species
  ggplot(aes(x = correlation)) + 
  geom_histogram(aes(fill = after_stat(x)),bins = 15) +
  labs(title = "Distribution of ENC frequencies \nfor interacting species",
       x = "Correlation between species (R)",
       y = "Frequency",
       fill = "correlation") +
  theme_minimal() +
  scale_fill_distiller(palette = "RdBu", limits = c(-1.1,1.1), direction = 1)+
  geom_vline(xintercept = stats_interactions$mean_correlation) + # Add mean and sd lines
  geom_vline(xintercept = stats_interactions$mean_correlation - stats_interactions$sd_correlation, lty = 2) +
  geom_vline(xintercept = stats_interactions$mean_correlation + stats_interactions$sd_correlation, lty = 2) +
  theme(legend.position = "none")


# Plot for non-interacting species

correlations_histogram_noninteracting <- 
  interaction_matrix |> 
  filter(weighted_interaction == 0) |>
 ggplot(aes(x = correlation)) +
  geom_histogram(aes(fill = after_stat(x)),bins = 15) +
  labs(title = "Distribution of ENC frequencies \nfor non-interacting species",
       x = "Correlation between species (R)",
       y = "Frequency",
       fill = "correlation") +
  theme_minimal() +
  scale_fill_distiller(palette = "RdBu", limits = c(-1.1,1.1), direction = 1)+
  geom_vline(xintercept = stats_noninteractions$mean_correlation) +
  geom_vline(xintercept = stats_noninteractions$mean_correlation - stats_noninteractions$sd_correlation, lty = 2) + 
  geom_vline(xintercept = stats_noninteractions$mean_correlation + stats_noninteractions$sd_correlation, lty = 2) +
  theme(legend.position = "none")


# add them together in one plot and save it as a png file in the c_outputs/figures/bees folder
correlations_histogram_interacting / correlations_histogram_noninteracting + 
  plot_annotation(tag_levels = "a")
ggsave("c_outputs/figures/bees/correlations_histogram_interactingVSnoninteracting.png", 
       width = 5, height = 8)
