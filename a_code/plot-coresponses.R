# script to make species' co-response curves from the time series model


# load libraries ---------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
theme_set(theme_pubr() +
            theme(panel.grid.major.x = element_line()))

# set the data sizes and seed --------------------------------------------------

set.seed(12)
npops = 29
tsl = 33

load("~/Documents/GitHub/groundfish-data-analysis/data/year_geom_means.Rdata")
rm(Year_Geom_Means, Year_Geom_Means_rare, Year_Geom_Means_SE)
# remove Notacanthidae, which was 0 in year 1
Year_Geom_Means_all = Year_Geom_Means_all[,-which(colnames(Year_Geom_Means_all) %in% "NOTACANTHIDAE")]

# Population growth rates ######################################################

# take the derivative of each population trend ----

# load the model
mod1 = readRDS(paste0("outputs/gam_hierarchical_gp.rds"))

source("~/Documents/GitHub/hierarchical-lpi/scripts/plot_mvgam_trend_custom.R")
trend_vals = list()
derivs_ls = list()
for(i in 1:npops){
  trend_vals[[i]] = plot_mvgam_trend_custom(mod1, derivatives = TRUE, series = i)
  derivs_ls[[i]] = trend_vals[[i]]$derivs
}
derivs = do.call(rbind, derivs_ls)

# each species' median derivative 
derivs_pops = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, median, na.rm = TRUE)))
derivs_pops_lower = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, quantile, probs = .05, na.rm = TRUE)))
derivs_pops_upper = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, quantile, probs = .95, na.rm = TRUE)))

# Interactions between species #################################################

# read diet matrix
A <- read.csv("~/Documents/GitHub/network_coherence/b_data/clean/globi_adjacencymatrix.csv", row.names = 1)
which(rownames(A) == "Notacanthidae")
A = A[-21,-21]

# reformat the matrix to get a list of species
Along = A |> mutate(sp1 = rownames(A)) |>
  pivot_longer(cols = -sp1, names_to = "sp2")
Along$sp2 = gsub("\\.", " ", Along$sp2)

# identify interacting partners vs. non-interacting species
non_interacting = Along |> dplyr::filter(value == 0) |> distinct()
interacting = Along |> dplyr::filter(value == 1) |> distinct()

# plot co-response curves for interacting species ------------------------------

# convert to matrix to more easily feed into the function
interacting = as.matrix(interacting[,-3])

# assign names to the deriv list
names(derivs_ls) <- colnames(Year_Geom_Means_all) |> 
  stringr::str_to_sentence()
names(derivs_ls) <- gsub("_", " ", names(derivs_ls))

# change the deriv names to match the interaction matrix
names(derivs_ls) <- gsub("Lycodes esmarki", "Lycodes esmarkii", names(derivs_ls))
names(derivs_ls) <- gsub("Nezumia bairdi", "Nezumia bairdii", names(derivs_ls))
names(derivs_ls) <- gsub("Notacanthus nasus", "Notacanthus chemnitzii", names(derivs_ls))
names(derivs_ls) <- gsub("Raja spinicauda", "Bathyraja spinicauda", names(derivs_ls))
names(derivs_ls) <- gsub("Raja radiata", "Amblyraja radiata", names(derivs_ls))
names(derivs_ls) <- gsub("Sebastes marinus", "Sebastes norvegicus", names(derivs_ls))
names(derivs_ls) <- gsub("Urophycis chesteri", "Phycis chesteri", names(derivs_ls))


# function to plot coresponses -------------------------------------------------

plot_coresponses = function(interaction_partners_vector){
  
  x = unname(interaction_partners_vector)
  
  df = data.frame(
    "species" = unlist(c(rep(x[1], tsl), rep(x[2], tsl))),
    # get mean rate per year
    "rate" = c(apply(derivs_ls[[x[1]]], 2, mean),
               apply(derivs_ls[[x[2]]], 2, mean))
  )
  
  p = ggpubr::ggdensity(data = na.omit(df), x = "rate", 
                        fill = "species", col = "species", 
                        add = "mean", rug = TRUE) +
    coord_cartesian(xlim = c(-.5,.5)) +
    labs(x = "Population growth rate", y = "Density",
         fill = "Species", col = "Species")
  return(p)
}

apply(interacting, 1, plot_coresponses)

