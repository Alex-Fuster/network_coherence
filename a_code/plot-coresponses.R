# script to make species' co-response curves from the time series model


# load libraries ---------------------------------------------------------------

library(dplyr)
library(tidyr)
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
# this model output is in the shared GDrive at:
#https://drive.google.com/drive/folders/1MaZpuz5NS2urYeKnSVhxJRUFPdv6Zy4y?usp=drive_link
mod1 = readRDS(paste0("~/Documents/GitHub/hierarchical-lpi/outputs/gam_hierarchical_gp.rds"))

source("a_code/plot_mvgam_trend_custom.R")
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


## Covariance matrix ###########################################################

time <- rownames(Year_Geom_Means_all) %>% as.numeric()
time <- time-min(time)
time_m <- as.matrix(time)

YData <- Year_Geom_Means_all 
# center on the baseline biomass and standardise by the mean
for(i in 1:ncol(YData)){
  YData[,i] = (YData[,i] - YData[1,i])/(mean(YData[,i], na.rm = TRUE))
}
biomass <- YData |> as.data.frame()

# format into long
biomass$time = time
dat = pivot_longer(biomass, cols = -c(time), names_to = "series", values_to = "y")
dat$series <- as.factor(dat$series)
dat$time <- as.integer(dat$time)
data_train = dat

# get species correlations
sp_correlations = lv_correlations(mod1)
covmatrix = sp_correlations$mean_correlations

# Interactions between species #################################################

# read diet matrix
A <- read.csv("b_data/clean/globi_adjacencymatrix.csv", row.names = 1)
which(rownames(A) == "Notacanthidae")
A = A[-21,-21]
A = as.matrix(A)
colnames(A) = gsub("\\.", " ", colnames(A))
rownames(A) = gsub("\\.", " ", rownames(A))

# match up the names to order the matrices
colnames(covmatrix) <- colnames(covmatrix) |> stringr::str_to_sentence()
colnames(covmatrix) <- gsub("_", " ", colnames(covmatrix))
# change the covmatrix names to match the interaction matrix
colnames(covmatrix) <- gsub("Lycodes esmarki", "Lycodes esmarkii", colnames(covmatrix))
colnames(covmatrix) <- gsub("Nezumia bairdi", "Nezumia bairdii", colnames(covmatrix))
colnames(covmatrix) <- gsub("Notacanthus nasus", "Notacanthus chemnitzii", colnames(covmatrix))
colnames(covmatrix) <- gsub("Raja spinicauda", "Bathyraja spinicauda", colnames(covmatrix))
colnames(covmatrix) <- gsub("Raja radiata", "Amblyraja radiata", colnames(covmatrix))
colnames(covmatrix) <- gsub("Sebastes marinus", "Sebastes norvegicus", colnames(covmatrix))
colnames(covmatrix) <- gsub("Urophycis chesteri", "Phycis chesteri", colnames(covmatrix))

rownames(covmatrix) <- rownames(covmatrix) |> stringr::str_to_sentence()
rownames(covmatrix) <- gsub("_", " ", rownames(covmatrix))
# change the covmatrix names to match the interaction matrix
rownames(covmatrix) <- gsub("Lycodes esmarki", "Lycodes esmarkii", rownames(covmatrix))
rownames(covmatrix) <- gsub("Nezumia bairdi", "Nezumia bairdii", rownames(covmatrix))
rownames(covmatrix) <- gsub("Notacanthus nasus", "Notacanthus chemnitzii", rownames(covmatrix))
rownames(covmatrix) <- gsub("Raja spinicauda", "Bathyraja spinicauda", rownames(covmatrix))
rownames(covmatrix) <- gsub("Raja radiata", "Amblyraja radiata", rownames(covmatrix))
rownames(covmatrix) <- gsub("Sebastes marinus", "Sebastes norvegicus", rownames(covmatrix))
rownames(covmatrix) <- gsub("Urophycis chesteri", "Phycis chesteri", rownames(covmatrix))

# order the interaction matrix by the names in covmatrix so both matrices match
A = A[rownames(covmatrix), colnames(covmatrix)]

## mask the correlation matrix by the interaction matrix

# histogram of interacting species 
A[which(A == 0)] <- NA
covmatrix_interacting = covmatrix*A

# histogram of non-interacting species
A[which(is.na(A))] <- 0
A[which(A == 1)] <- NA
A[which(A == 0)] <- 1
covmatrix_noninteracting = covmatrix*A

# histogram of interacting species' correlations
cov_interacting <- as.vector(covmatrix_interacting) |> na.omit() |> as.vector()
corr_df = data.frame("value" = cov_interacting,
                     "group" = rep("A", length(cov_interacting)))
(correlations_histogram_interacting = ggplot(data = corr_df, aes(x = value)) +
    geom_histogram(aes(fill = after_stat(x)),bins = 15) +
    scale_fill_distiller(palette = "RdBu", limits = c(-1.1,1.1), direction = 1) +
    labs(y = "Frequency", x = "Correlation between species (R)", fill = "R",
         title = "Interacting species") +
    geom_vline(xintercept = mean(corr_df$value)) +
    geom_vline(xintercept = mean(corr_df$value) - sd(corr_df$value), lty = 2) + 
    geom_vline(xintercept = mean(corr_df$value) + sd(corr_df$value), lty = 2)
)

# histogram of non-interacting species' correlations
cov_noninteracting <- as.vector(covmatrix_noninteracting) |> na.omit() |> as.vector()
corr_df2 = data.frame("value" = cov_noninteracting,
                     "group" = rep("A", length(cov_noninteracting)))
(correlations_histogram_noninteracting = ggplot(data = corr_df2, aes(x = value)) +
    geom_histogram(aes(fill = after_stat(x)),bins = 15) +
    scale_fill_distiller(palette = "RdBu", limits = c(-1.1,1.1), direction = 1) +
    labs(y = "Frequency", x = "Correlation between species (R)", fill = "R",
         title = "Non-interacting species") +
    geom_vline(xintercept = mean(corr_df2$value)) +
    geom_vline(xintercept = mean(corr_df2$value) - sd(corr_df2$value), lty = 2) + 
    geom_vline(xintercept = mean(corr_df2$value) + sd(corr_df2$value), lty = 2) +
    theme(legend.position = "none")
)

# add them together in one plot
correlations_histogram_interacting / correlations_histogram_noninteracting + 
  plot_annotation(tag_levels = "a")
ggsave("c_outputs/figures/correlations_histogram_interactingVSnoninteracting.png", 
       width = 5, height = 8)


ggpubr::ggdensity(data = corr_df_all, x = "value", 
                  fill = "interacting", col = "interacting",
                                           add = "mean") +
                      labs(x = "Correlation between species (R)", y = "Density",
                           fill = "Species", col = "Species") +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2")
ggsave("c_outputs/figures/correlations_density_interactingVSnoninteracting.png", 
       width = 5, height = 5)
