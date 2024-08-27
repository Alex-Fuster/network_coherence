# load functions
library(Matrix)
library(MASS)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
source("LV_simulations_new/functions.R")

# parameters
params <- list(
  S = 10, # number of species
  C = 0.2, # connectance of the food web
  aij_params = c(0, 0.5), # minimum and maximum interspecific effect of predation
  mu_delta_r = 0, # mean responses to perturbation
  sd_delta_r = 0.5, # sd of response to perturbation
  covMatrix_type = "positive", # type of covariance (mixed or all positive)
  sd_X = rep(1, 10) # standard deviations of responses for each species
)

simulate_response <- function(S, C, aij_params, mu_delta_r, sd_delta_r, covMatrix_type, sd_X){
  corMat <- createCorMat(S, covMatrix_type)
  covMat <- as.matrix(Matrix::nearPD(sd_X %*% t(sd_X) * corMat)$mat)
  A <- sim_quantitative_network("predator-prey", S = S, C = C, aij_params = aij_params) # minimum effect is 0 and max effect is 0.1
  B <- ginv(A)
  equilibrium_pre <- runif(S, 1, 10)
  r_pre <- -A %*% equilibrium_pre
  delta_r <- mvrnorm(n = 1, mu = rnorm(S, mu_delta_r, sd_delta_r), Sigma = covMat, empirical = FALSE)
  r_post <- r_pre + delta_r
  equilibrium_post <- -B%*%r_post
  df <- data.frame(species = paste0("sp", c(1:S)), X_pre = equilibrium_pre, X_post = equilibrium_post)
  return(df)
}

out <- data.frame()
# Scenario 1: positive correlations (cor between 0 and 1) among species responses, highly variable responses (sd_delta_r = 0.5)
for (i in 1:100){
  X <- with(params, simulate_response(S, C, aij_params, mu_delta_r, sd_delta_r, covMatrix_type, sd_X))
  out <- rbind(
    out,
    data.frame(covtype = "positive", perturbation = "strong", sum_deltaX = sum(X$X_pre - X$X_post), sd_deltaX = sd(X$X_pre - X$X_post))
  )
}
# Scenario 2: positive correlations (cor between 0 and 1) among species responses, weakly variable responses (sd_delta_r = 0.1)
params$sd_delta_r <- 0.1
for (i in 1:100){
  X <- with(params, simulate_response(S, C, aij_params, mu_delta_r, sd_delta_r, covMatrix_type, sd_X))
  out <- rbind(
    out,
    data.frame(covtype = "positive", perturbation = "weak", sum_deltaX = sum(X$X_pre - X$X_post), sd_deltaX = sd(X$X_pre - X$X_post))
  )
}
# Scenario 3: mixed correlations (cor between 0 and 1) among species responses, weakly variable responses (sd_delta_r = 0.1)
params$covMatrix_type <- "mixed"
for (i in 1:100){
  X <- with(params, simulate_response(S, C, aij_params, mu_delta_r, sd_delta_r, covMatrix_type, sd_X))
  out <- rbind(
    out,
    data.frame(covtype = "mixed", perturbation = "weak", sum_deltaX = sum(X$X_pre - X$X_post), sd_deltaX = sd(X$X_pre - X$X_post))
  )
}
# Scenario 4: mixed correlations (cor between 0 and 1) among species responses, highly variable responses (sd_delta_r = 0.5)
params$sd_delta_r <- 0.5
for (i in 1:100){
  X <- with(params, simulate_response(S, C, aij_params, mu_delta_r, sd_delta_r, covMatrix_type, sd_X))
  out <- rbind(
    out,
    data.frame(covtype = "mixed", perturbation = "strong", sum_deltaX = sum(X$X_pre - X$X_post), sd_deltaX = sd(X$X_pre - X$X_post))
  )
}

# Scenario 5: mixed correlations (cor between 0 and 1) among species responses, highly variable responses (sd_delta_r = 0.5), and positive on average (mu_delta_r = 0.5)
params$mu_delta_r <- 0.5
for (i in 1:100){
  X <- with(params, simulate_response(S, C, aij_params, mu_delta_r, sd_delta_r, covMatrix_type, sd_X))
  out <- rbind(
    out,
    data.frame(covtype = "mixed", perturbation = "strong positive", sum_deltaX = sum(X$X_pre - X$X_post), sd_deltaX = sd(X$X_pre - X$X_post))
  )
}

# Scenario 6: positive correlations (cor between 0 and 1) among species responses, highly variable responses (sd_delta_r = 0.5), and positive on average (mu_delta_r = 0.5)
params$covMatrix_type <- "positive"
for (i in 1:100){
  X <- with(params, simulate_response(S, C, aij_params, mu_delta_r, sd_delta_r, covMatrix_type, sd_X))
  out <- rbind(
    out,
    data.frame(covtype = "positive", perturbation = "strong positive", sum_deltaX = sum(X$X_pre - X$X_post), sd_deltaX = sd(X$X_pre - X$X_post))
  )
}

ggplot(out) +
  geom_boxplot(aes(x = covtype, y = sum_deltaX, fill = perturbation)) +
  lims(y = c(-100, 100))

ggplot(out) +
  geom_boxplot(aes(x = covtype, y = sd_deltaX, fill = perturbation)) +
  lims(y = c(0, 25))
