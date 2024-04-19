# Script to evaluate species' contributions to network coherence
# for the Newfoundland groundfish community example based on a multivariate 
# hierarchical generalized additive model 

# Author: Katherine Hébert

# libraries ----

library(here)
library(dplyr)
library(tidyr)
library(mvgam)
library(tidybayes)
library(ggplot2)
library(ggpubr)
library(corrplot)
library(patchwork)

# set ggplot theme
theme_set(theme_pubr() +
            theme(panel.grid.major.x = element_line()))

set.seed(12)

# import the model
mod1 = readRDS(paste0("c_outputs/gam_hierarchical_gp.rds"))

# prepare the data =============================================================

# groundfish-data-analysis contains the data folder that is copied from a clone of:
# https://github.com/eric-pedersen/groundfish-data-analysis
load("~/Documents/GitHub/groundfish-data-analysis/data/year_geom_means.Rdata")
rm(Year_Geom_Means, Year_Geom_Means_rare, Year_Geom_Means_SE)

time <- rownames(Year_Geom_Means_all) %>% as.numeric()
time <- time-min(time)
time_m <- as.matrix(time)

# remove Notacanthidae, which was 0 in year 1
Year_Geom_Means_all = Year_Geom_Means_all[,-which(colnames(Year_Geom_Means_all) %in% "NOTACANTHIDAE")]

npops <- ncol(Year_Geom_Means_all)
tsl <- nrow(Year_Geom_Means_all)

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

## Community coherence =========================================================

# take the derivative of each population trend ----

# load this mvgam function I customised to access the derivatives more easily
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
colnames(derivs_pops) = colnames(Year_Geom_Means_all)
matplot(derivs_pops, x = time+1981, type = "l")


## Covariance between populations' trends ======================================

# get species correlations
data_train = dat
sp_correlations = lv_correlations(mod1)

# edit the species names to be prettier
colnames(sp_correlations$mean_correlations) = gsub("_", " ", colnames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()
rownames(sp_correlations$mean_correlations) = gsub("_", " ", rownames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()

# save mean correlation in an object for easier use
corrs = sp_correlations$mean_correlations

# make a df to store results
df = data.frame("species" = rownames(corrs),
                "mu_corr" = NA,
                "sd_corr" = NA)
df$mu_corr = apply(corrs, 2, median, na.rm = TRUE) # median! not mean. 
df$sd_corr = apply(corrs, 2, sd, na.rm = TRUE)

# medians, to see + and - contributions
df$species = factor(df$species, 
                    levels = df$species[order(df$mu_corr)])
ggplot(data = df,
       aes(x = mu_corr, y = species)) +
  geom_vline(xintercept = mean(mu_corr, na.rm = T), lwd = .2) +
  geom_point(aes(fill = mu_corr), size = 3, pch = 21) +
  scale_fill_distiller(palette = "RdBu")

# absolute value version
df$abs_mu = abs(df$mu_corr)
df$species = factor(df$species, 
                    levels = df$species[order(df$abs_mu)])
ggplot(data = df,
       aes(x = abs_mu, y = species)) +
  geom_point(aes(fill = abs_mu), size = 3, pch = 21) +
  scale_fill_viridis_c(option = "plasma") +
  labs(x = "Contribution to community correlation\n(Absolute median correlation with all other species)",
       fill = "Median correlation", y = "") +
  theme(axis.text.y = element_text(face = "italic"))

## Contributions as the variance of each species to overall community variance ----

# make a df to store results
df = data.frame("species" = rownames(corrs),
                "var_sp" = NA)
df$var_sp = derivs_pops |> apply(2, var, na.rm = TRUE)

# order by variance
df$species = factor(df$species, 
                    levels = df$species[order(df$var_sp)])
# plot!
ggplot(data = df,
       aes(x = var_sp, y = species)) +
  geom_point(aes(fill = var_sp), size = 3, pch = 21) +
  scale_fill_viridis_c(option = "plasma") +
  labs(x = "Contribution to community coherence\n(Variance of each species' growth rate trend)",
       fill = "Variance", y = "") +
  theme(axis.text.y = element_text(face = "italic"))
