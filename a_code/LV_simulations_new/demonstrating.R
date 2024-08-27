# 4 scenarios:
# 1: positive correlations (cor between 0 and 1) among species responses, highly variable responses (sd_delta_r = 0.5)
# 2: positive correlations (cor between 0 and 1) among species responses, weakly variable responses (sd_delta_r = 0.1)
# 3: mixed correlations (cor between -1 and 1) among species responses, highly variable responses (sd_delta_r = 0.5)
# 4: mixed correlations (cor between -1 and 1) among species responses, weakly variable responses (sd_delta_r = 0.1)

# load functions
library(Matrix)
library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
source("generateCovMat.R")

S <- 25

# Step 1: build cov Matrix
# build positive correlation and covariance matrices.
positive_corMat <- createCorMat(S, "positive")
plot_matrix(positive_corMat, "Positive correlations")
# to get the covariances, we need to set the standard deviations of responses for each species (Cov(X,Y) = sd_X*sd_Y*Cor(X,Y))
# let's set it to 1 for all species
sd <- rep(1, S) 
positive_covMat <- as.matrix(Matrix::nearPD(sd %*% t(sd) * positive_corMat)$mat)
plot_matrix(positive_covMat, "Positive covariances")

# build a mixed correlation and covariance matrices
mixed_corMat <- createCorMat(S, "mixed")
plot_matrix(mixed_corMat, "Positive correlations")
# to get the covariances, we need to set the standard deviations of responses for each species (Cov(X,Y) = sd_X*sd_Y*Cor(X,Y))
# let's set it to 1 for all species
sd <- rep(1, S) 
mixed_covMat <- as.matrix(Matrix::nearPD(sd %*% t(sd) * mixed_corMat)$mat)
plot_matrix(mixed_covMat, "Positive covariances")

# Step 2: Generate food web
C = 0.2 # connectance
A <- sim_quantitative_network("predator-prey", S = S, C = C, aij_params = c(0, 0.5)) # minimum effect is 0 and max effect is 0.1
B <- ginv(A)
plot_matrix(A, "Food web (A)")
plot_matrix(B, "Total effect (B)")

# Step 3: Generate r_pre and r_post
mu_delta_r <- 0 # average response
sd_delta_r <- 0.5 # variation among responses
equilibrium_pre <- runif(S, 1, 10)
r_pre <- -A %*% Biomass_at_equilibrium
delta_r <- mvrnorm(n = 1, mu = rnorm(S, mu_delta_r, sd_delta_r), Sigma = mixed_covMat, empirical = FALSE)
r_post <- r_pre + delta_r

# Step 4: get new equilbrium
equilibrium_post <- -B%*%r_post
delta_X <- equilibrium_post - equilibrium_pre
mean(delta_X)
sd(delta_X)
