# Packages ----------------------------------------------------------------

source("./a_code/functions/ipak.R")

ipak(c("bipartite","gridGraphics","ggplot2","ggplotify"))

# Rearranging the interaction matrix --------------------------------------

interaction_matrix_fig <- interaction_matrix[,-4]

interaction_matrix_fig <- interaction_matrix_fig |> 
  pivot_wider(names_from = spp_pollinators, 
              values_from = weighted_interaction, 
              values_fill = 0) |> 
  column_to_rownames(var = "spp_plants") |> 
  as.matrix()

saveRDS(interaction_matrix_fig, file = "meu_objeto.rda")
interaction_matrix_fig <- readRDS("meu_objeto.rda")

# Viewing the network -----------------------------------------------------

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

