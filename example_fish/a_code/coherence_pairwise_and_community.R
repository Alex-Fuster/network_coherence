# Script to calculate the covariance between population abundance trends and
# community-level coherence across all populations, based on a multivariate 
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

# prepare the data =============================================================

# groundfish-data-analysis contains the data folder that is copied from a clone of:
# https://github.com/eric-pedersen/groundfish-data-analysis
load("groundfish-data-analysis/data/year_geom_means.Rdata")
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

################################################################################
# check number of knots to use in the GAM 

knots = ceiling(tsl/4)

kcheck_gam = mgcv::gam(y ~ s(time, bs = "tp", k = 9) + s(series, bs = 're', k = npops),
                       data = data_train)
summary(kcheck_gam)

kcheck_gam = mgcv::gam(y ~ s(time, bs = "tp", k = 18) + s(series, bs = 're', k = npops),
                       data = data_train)
summary(kcheck_gam)
mgcv::k.check(kcheck_gam)
mgcv::gam.check(kcheck_gam)

################################################################################

# hierarchical gam on all populations ==========================================

# prepare the priors 

mvgam_prior <- mvgam(data = data_train,
                     formula = y ~ 
                       # global smoother for all pops over time
                       s(time, bs = "tp", k = knots) + 
                       # random intercept per group
                       s(series, bs = 're', k = npops),
                     family = "gaussian",
                     trend_model = 'GP',
                     chains = 3,
                     use_stan = TRUE,
                     prior_simulation = TRUE)

# record the priors
test_priors <- get_mvgam_priors(y ~ 
                                  # global smoother for all pops over time
                                  s(time, bs = "tp", k = knots) + 
                                  # random intercept per group
                                  s(series, bs = 're', k = npops),
                                family = "gaussian",
                                data = data_train,
                                trend_model = 'GP',
                                use_stan = TRUE)
# look at the priors
plot(mvgam_prior, type = 'smooths')

# train the model on data ======================================================

# only rerun the model if something changes, because it takes a long time...
# mod1 <- mvgam(data = data_train,
#               formula =  y ~ s(time, bs = "tp", k = knots) + 
#                 s(series, bs = "re"),
#               use_lv = TRUE,
#               family = "gaussian",
#               trend_model = 'GP',
#               use_stan = TRUE,
#               chains = 3, 
#               burnin = 5000,
#               samples = 10000
# )
# saveRDS(mod1, paste0("c_outputs/gam_hierarchical_gp.rds")) 

mod1 = readRDS(paste0("c_outputs/gam_hierarchical_gp.rds"))

# to view the Stan model file:
m <- mod1
# save the model file
cmdstanr::write_stan_file(m$model_file, "stan")
summary(mod1)

# check convergence
# rstan::stan_trace(mod1$model_output, 'b')

## Covariance between populations' trends ======================================
                  
# get species correlations
data_train = dat
sp_correlations = lv_correlations(mod1)
saveRDS(sp_correlations, "c_outputs/gam_hierarchical_species_correlations.rds")

# edit the species names to be prettier
colnames(sp_correlations$mean_correlations) = gsub("_", " ", colnames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()
rownames(sp_correlations$mean_correlations) = gsub("_", " ", rownames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()

# Plot as histogram
sp_correlations$mean_correlations[which(lower.tri(sp_correlations$mean_correlations))] |>
  hist(col = "grey20", border = "white", lwd = .2,
       xlab = "Corrélation entre espèces",
       ylab = "Fréquence",
       main = "", cex = 3)

corrs = sp_correlations$mean_correlations[which(lower.tri(sp_correlations$mean_correlations))]                  

# Plot as a heatmap of pairwise correlations
png(height=1800, width=1800, file="c_outputs/figures/heatmap_species_associations.png", type = "cairo")
corrplot::corrplot(sp_correlations$mean_correlations, 
                   type = "lower",
                   method = "color", 
                   tl.cex = 2.5, cl.cex = 3, tl.col = "black", font = 3)
dev.off()

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
derivs_pops_lower = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, quantile, probs = .05, na.rm = TRUE)))
derivs_pops_upper = do.call(cbind, lapply(derivs_ls, FUN = function(x) apply(x, 2, quantile, probs = .95, na.rm = TRUE)))
matplot(derivs_pops, x = time+1981, type = "l")
abline(h = 0)
abline(v = 1990, lty = 2) # start of the cod collapse according to Pedersen 2017
abline(v = 2000, lty = 2) # start of the recovery phase (beginning of the 2000s)

# get each species' mean derivative over the whole time series ----

(temporal_trend = data.frame(
  species = colnames(Year_Geom_Means_all),
  mu_deriv = derivs_pops |> apply(2, median, na.rm = TRUE),
  lower = derivs_pops |> apply(2, quantile, probs = .05, na.rm = TRUE),
  upper = derivs_pops |> apply(2, quantile, probs = .95, na.rm = TRUE),
  sd = derivs_pops |> apply(2, sd, na.rm = TRUE),
  mu = derivs_pops |> apply(2, mean, na.rm = TRUE)
))
temporal_trend$species = gsub("_", " ", temporal_trend$species) |> stringr::str_to_sentence()
temporal_trend$species =  factor(temporal_trend$species,
                                 levels = temporal_trend$species[order(temporal_trend$mu_deriv)])
colnames(derivs_pops) = colnames(Year_Geom_Means_all)
colnames(derivs_pops_lower) = colnames(Year_Geom_Means_all)
colnames(derivs_pops_upper) = colnames(Year_Geom_Means_all)

derivs_pops_df = derivs_pops |> 
  as.data.frame() |> 
  mutate(year = time + 1981) |>
  pivot_longer(cols = -year)
derivs_pops_df_lower = derivs_pops_lower |> 
  as.data.frame() |> 
  mutate(year = time + 1981) |>
  pivot_longer(cols = -year, values_to = "lower")
derivs_pops_df_upper = derivs_pops_upper |> 
  as.data.frame() |> 
  mutate(year = time + 1981) |>
  pivot_longer(cols = -year, values_to = "upper")
derivs_pops_df = full_join(derivs_pops_df, derivs_pops_df_lower) |>
  full_join(derivs_pops_df_upper)
# set first time step's derivative to 0
derivs_pops_df$value[which(derivs_pops_df$year == 1981)] <- 0
derivs_pops_df$lower[which(derivs_pops_df$year == 1981)] <- 0
derivs_pops_df$upper[which(derivs_pops_df$year == 1981)] <- 0

derivs_pops_df$name = gsub("_", " ", derivs_pops_df$name) |> 
  stringr::str_to_sentence()
derivs_pops_df$name =  factor(derivs_pops_df$name,
                              levels = temporal_trend$species[order(temporal_trend$mu_deriv)])

derivs_without1981 = dplyr::filter(derivs_pops_df, year != 1981)

#### average and variance of assemblage per year ----

(temporal_trend_yearly = data.frame(
  year = rownames(Year_Geom_Means_all),
  mu_deriv = derivs_pops |> apply(1, median, na.rm = TRUE),
  lower = derivs_pops |> apply(1, quantile, probs = .05, na.rm = TRUE),
  upper = derivs_pops |> apply(1, quantile, probs = .95, na.rm = TRUE),
  sd = derivs_pops |> apply(1, sd, na.rm = TRUE),
  mu = derivs_pops |> apply(1, mean, na.rm = TRUE)
))
temporal_trend_yearly$year = as.character(temporal_trend_yearly$year)
saveRDS(temporal_trend_yearly, "c_outputs/temporal_trend_yearly.rds")

# coherence ----

(plot_coherence = ggplot(data = temporal_trend_yearly) +
   geom_line(aes(y = sd, x = as.numeric(year), col = sd), lwd = 1.5) +
   geom_point(aes(y = sd, x = as.numeric(year), col = sd), 
              size = .9) +
   labs(y = "Community variability (σ)", 
        x = "",
        col = "σ") +
   scale_color_distiller(palette = "YlGnBu", direction = 1, limits = c(0, .17)) +
   coord_cartesian(ylim = c(0, .18)) +
   theme_pubr() +
   theme(panel.grid.major = element_line())
)
