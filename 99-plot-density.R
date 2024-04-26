combined_df |> 
  ggplot(aes(x = environ_values,
             color = group,
             group = species)) +
  geom_density() +
  labs(title = NULL,
       x = "Mean Temperature (°C)",
       y = "Frequency") +
  theme_minimal() +
  theme(panel.background = element_rect(fill = "white", color = "white"),
        plot.background = element_rect(fill = "white"))

ggsave(filename = "c_outputs/figures/bees/density_plot.png" , 
       width = 6, 
       height = 4)
