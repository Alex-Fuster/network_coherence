# rerun the model for the different time chunks to get correlation matrices
# for each stage of the collapse... right?

# phases to do ----
# until 1990: pre-collapse
# 1990 to 2000: collapse
# 2000 to end: recovery/post-collapse

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

data_precollapse = dat[which(dat$time %in% 0:9),] # pre-1990
data_duringcollapse = dat[which(dat$time %in% 10:19),] # 1990-1999
data_postcollapse = dat[which(dat$time %in% 20:33),] # 2000-2013

## PRE-COLLAPSE ################################################################
################################################################################
# check number of knots to use in the GAM 

knots = ceiling(length(unique(data_precollapse$time))/4)

kcheck_gam = mgcv::gam(y ~ s(time, bs = "tp", k = 6) + s(series, bs = 're', k = npops),
                       data = data_precollapse)
mgcv::gam.check(kcheck_gam)
summary(kcheck_gam)

knots = 6

################################################################################

# hierarchical gam on all populations ==========================================

# prepare the priors 

mvgam_prior <- mvgam(data = data_precollapse,
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
                                data = data_precollapse,
                                trend_model = 'GP',
                                use_stan = TRUE)
# look at the priors
plot(mvgam_prior, type = 'smooths')

# train the model on data ======================================================

# only rerun the model if something changes, because it takes a long time...
mod1 <- mvgam(data = data_precollapse,
              formula =  y ~ s(time, bs = "tp", k = knots) +
                s(series, bs = "re"),
              use_lv = TRUE,
              family = "gaussian",
              trend_model = 'GP',
              use_stan = TRUE,
              chains = 3,
              burnin = 500,#0,
              samples = 1000#0
)
saveRDS(mod1, paste0("c_outputs/gam_hierarchical_gp_precollapse.rds")) 

mod1 = readRDS(paste0("c_outputs/gam_hierarchical_gp_precollapse.rds"))

# to view the Stan model file:
m <- mod1
# save the model file
cmdstanr::write_stan_file(m$model_file, "stan")
summary(mod1)

## Covariance between populations' trends ======================================

# get species correlations
data_train = data_precollapse
sp_correlations = lv_correlations(mod1)

# edit the species names to be prettier
colnames(sp_correlations$mean_correlations) = gsub("_", " ", colnames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()
rownames(sp_correlations$mean_correlations) = gsub("_", " ", rownames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()

# Plot as a heatmap of pairwise correlations
png(height=1800, width=1800, file="c_outputs/fish-example/figures/heatmap_species_associations-precollapse.png", type = "cairo")
corrplot::corrplot(sp_correlations$mean_correlations, 
                   type = "lower",
                   method = "color", 
                   tl.cex = 2.5, cl.cex = 3, tl.col = "black", font = 3)
dev.off()

## DURING-COLLAPSE #############################################################
################################################################################
# check number of knots to use in the GAM 

knots = ceiling(length(unique(data_duringcollapse$time))/4)

kcheck_gam = mgcv::gam(y ~ s(time, bs = "tp", k = 5) + s(series, bs = 're', k = npops),
                       data = data_duringcollapse)
mgcv::gam.check(kcheck_gam)
summary(kcheck_gam)

knots = 5

################################################################################

# hierarchical gam on all populations ==========================================

# prepare the priors 

mvgam_prior <- mvgam(data = data_duringcollapse,
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
                                data = data_duringcollapse,
                                trend_model = 'GP',
                                use_stan = TRUE)
# look at the priors
plot(mvgam_prior, type = 'smooths')

# train the model on data ======================================================

# only rerun the model if something changes, because it takes a long time...
mod1 <- mvgam(data = data_duringcollapse,
              formula =  y ~ s(time, bs = "tp", k = knots) +
                s(series, bs = "re"),
              use_lv = TRUE,
              family = "gaussian",
              trend_model = 'GP',
              use_stan = TRUE,
              chains = 3,
              burnin = 500,#0,
              samples = 1000#0
)
saveRDS(mod1, paste0("c_outputs/gam_hierarchical_gp_duringcollapse.rds")) 

mod1 = readRDS(paste0("c_outputs/gam_hierarchical_gp_duringcollapse.rds"))

# to view the Stan model file:
m <- mod1
# save the model file
cmdstanr::write_stan_file(m$model_file, "stan")
summary(mod1)

## Covariance between populations' trends ======================================

# get species correlations
data_train = data_duringcollapse
sp_correlations = lv_correlations(mod1)

# edit the species names to be prettier
colnames(sp_correlations$mean_correlations) = gsub("_", " ", colnames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()
rownames(sp_correlations$mean_correlations) = gsub("_", " ", rownames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()

# Plot as a heatmap of pairwise correlations
png(height=1800, width=1800, file="c_outputs/fish-example/figures/heatmap_species_associations-duringcollapse.png", type = "cairo")
corrplot::corrplot(sp_correlations$mean_correlations, 
                   type = "lower",
                   method = "color", 
                   tl.cex = 2.5, cl.cex = 3, tl.col = "black", font = 3)
dev.off()


## POST-COLLAPSE #############################################################
################################################################################
# check number of knots to use in the GAM 

knots = ceiling(length(unique(data_postcollapse$time))/4)

kcheck_gam = mgcv::gam(y ~ s(time, bs = "tp", k = 4) + s(series, bs = 're', k = npops),
                       data = data_postcollapse)
mgcv::gam.check(kcheck_gam)
summary(kcheck_gam)

knots = 4

################################################################################

# hierarchical gam on all populations ==========================================

# prepare the priors 

mvgam_prior <- mvgam(data = data_postcollapse,
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
                                data = data_postcollapse,
                                trend_model = 'GP',
                                use_stan = TRUE)
# look at the priors
plot(mvgam_prior, type = 'smooths')

# train the model on data ======================================================

# only rerun the model if something changes, because it takes a long time...
mod1 <- mvgam(data = data_postcollapse,
              formula =  y ~ s(time, bs = "tp", k = knots) +
                s(series, bs = "re"),
              use_lv = TRUE,
              family = "gaussian",
              trend_model = 'GP',
              use_stan = TRUE,
              chains = 3,
              burnin = 500,#0,
              samples = 1000#0
)
saveRDS(mod1, paste0("c_outputs/gam_hierarchical_gp_postcollapse.rds")) 

mod1 = readRDS(paste0("c_outputs/gam_hierarchical_gp_postcollapse.rds"))

# to view the Stan model file:
m <- mod1
# save the model file
cmdstanr::write_stan_file(m$model_file, "stan")
summary(mod1)

## Covariance between populations' trends ======================================

# get species correlations
data_train = data_postcollapse
sp_correlations = lv_correlations(mod1)

# edit the species names to be prettier
colnames(sp_correlations$mean_correlations) = gsub("_", " ", colnames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()
rownames(sp_correlations$mean_correlations) = gsub("_", " ", rownames(sp_correlations$mean_correlations)) |> 
  stringr::str_to_sentence()

# Plot as a heatmap of pairwise correlations
png(height=1800, width=1800, file="c_outputs/fish-example/figures/heatmap_species_associations-postcollapse.png", type = "cairo")
corrplot::corrplot(sp_correlations$mean_correlations, 
                   type = "lower",
                   method = "color", 
                   tl.cex = 2.5, cl.cex = 3, tl.col = "black", font = 3)
dev.off()


### triplot of coherence before, during, and after collapse

pre = readRDS(paste0("c_outputs/gam_hierarchical_gp_precollapse.rds"))
during = readRDS(paste0("c_outputs/gam_hierarchical_gp_duringcollapse.rds"))
post = readRDS(paste0("c_outputs/gam_hierarchical_gp_postcollapse.rds"))

pre_correlations = lv_correlations(pre)
during_correlations = lv_correlations(during)
post_correlations = lv_correlations(post)


# Plot as histogram
pre_correlations$mean_correlations[which(lower.tri(pre_correlations$mean_correlations))] |>
  hist(col = "grey20", border = "white", lwd = .2,
       xlab = "Corrélation entre espèces",
       ylab = "Fréquence",
       main = "", cex = 3)

pre_corrs = data.frame("group" = "pre", "value" = pre_correlations$mean_correlations[which(lower.tri(pre_correlations$mean_correlations))])
during_corrs = data.frame("group" = "during", 
                          "value" = during_correlations$mean_correlations[which(lower.tri(during_correlations$mean_correlations))])
post_corrs = data.frame("group" = "post", "value" = post_correlations$mean_correlations[which(lower.tri(post_correlations$mean_correlations))])
corrs = rbind(pre_corrs, during_corrs, post_corrs)
corrs$group = factor(corrs$group, levels = c("pre", "during", "post"))

(A = ggplot(data = pre_corrs, aes(x = value)) +
    geom_histogram(aes(fill = after_stat(x)),bins = 15, col = "black", lwd = .1) +
    geom_vline(xintercept =  mean(pre_corrs$value)) +
    geom_vline(xintercept =  mean(pre_corrs$value) - sd(pre_corrs$value), lty = 2) +
    geom_vline(xintercept =  mean(pre_corrs$value) + sd(pre_corrs$value), lty = 2) +
    scale_fill_distiller(palette = "RdBu", limits = c(-1.1,1.1), direction = 1) +
    labs(y = "Frequency", 
         #x = "Correlation between species (R)", 
         title = "Pre-collapse (1980-1989)",
         fill = "R") +
    theme(legend.position = "none") +
    coord_cartesian(xlim = c(-1.1,1.1))
)
(B = ggplot(data = during_corrs, aes(x = value)) +
    geom_histogram(aes(fill = after_stat(x)),bins = 15, col = "black", lwd = .1) +
    geom_vline(xintercept =  mean(during_corrs$value)) +
    geom_vline(xintercept =  mean(during_corrs$value) - sd(during_corrs$value), lty = 2) +
    geom_vline(xintercept =  mean(during_corrs$value) + sd(during_corrs$value), lty = 2) +
    scale_fill_distiller(palette = "RdBu", limits = c(-1.1,1.1), direction = 1) +
    labs(y = "Frequency", 
         #x = "Correlation between species (R)", 
         title = "Collapse (1990-1999)",
         fill = "R") +
    theme(legend.position = "none")+
    coord_cartesian(xlim = c(-1.1,1.1))
)

(C = ggplot(data = post_corrs, aes(x = value)) +
    geom_histogram(aes(fill = after_stat(x)), bins = 15, col = "black", lwd = .1) +
    geom_vline(xintercept =  mean(post_corrs$value)) +
    geom_vline(xintercept =  mean(post_corrs$value) - sd(post_corrs$value), lty = 2) +
    geom_vline(xintercept =  mean(post_corrs$value) + sd(post_corrs$value), lty = 2) +
    scale_fill_distiller(palette = "RdBu", limits = c(-1.1,1.1), direction = 1) +
    labs(y = "Frequency", 
         x = "Correlation between species (R)", 
         title = "Post-collapse (2000-2013)",
         fill = "R") +
    theme(legend.position = "none")+
    coord_cartesian(xlim = c(-1.1,1.1))
)

A / B / C
ggsave("c_outputs/fish-example/figures/coherence_preduringpostcollapse.png", width = 4.25, height = 7.64)


