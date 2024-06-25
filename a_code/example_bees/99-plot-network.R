# Packages ----------------------------------------------------------------

source("./a_code/functions/ipak.R")

ipak(c("bipartite","gridGraphics","ggplot2","ggplotify"))

# Rearranging the interaction matrix --------------------------------------

# Source last script
source("a_code/example_bees/06-correlations.R")

# Rearranging the interaction matrix --------------------------------------


interaction_matrix_fig <- 
  interaction_matrix |> 
  select(-correlation) |> # exclude the correlation column
  pivot_wider(names_from = spp_pollinators, # make the matrix wider
              values_from = weighted_interaction, # fill the matrix with the weighted interactions
              values_fill = 0) |> # fill the NA values with 0
  column_to_rownames(var = "spp_plants") |> # make the plants the row names
  as.matrix() # make it a matrix

# Plotting the weighted interactions -------------------------------------------

visweb_interaction_matrix <- grid.grabExpr(grid.echo(function() visweb(interaction_matrix_fig, 
                                                      circles=TRUE,  
                                                      boxes=TRUE,  
                                                      labsize=1, 
                                                      circle.max=3, 
                                                      text="no")))

# shift the left axis labels to the right

visweb_interaction_matrix[["children"]][["graphics-plot-1-left-axis-labels-1"]][["x"]] <- unit(0.8, units = "in")

grid.newpage(); grid.draw(visweb_interaction_matrix)

ggsave(filename = "c_outputs/figures/bees/interaction_matrix_visweb.png", 
       plot = visweb_interaction_matrix,
       width = 7.5,
       height = 8)


# Plot binary network -----------------------------------------------------

interaction_matrix_fig <-
  interaction_matrix |>
  select(-correlation) |>   # Exclude the correlation column
  mutate(binary_interaction = case_when(weighted_interaction > 0 ~ 1,
                                          .default = weighted_interaction)) |>  # Binarize the matrix
  select(-weighted_interaction)                   # Exclude the weighted_interaction column

binary_network <- ggplot(data = interaction_matrix_fig) +
  geom_tile(aes(x = spp_plants, y = spp_pollinators, 
                fill = as.character(binary_interaction)), col = "lightgrey") +
  theme_minimal() +
  ggtitle("A matrix\n\nInteraction network") +
  theme(axis.text.x = element_text(angle = 90, 
                                   vjust = 0.5, hjust = 1, 
                                   face = "italic"),
        axis.text.y = element_text(face = "italic"),
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = c("white", "black")) +
  labs(fill = "Interaction", x = "", y = "")


ggsave(filename = "c_outputs/figures/bees/interaction_matrix_binary.png", 
       width = 8.23,
       height = 6.82)
