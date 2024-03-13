### simulate species dynamics using predefined parameters
### uses the function simulate_dynamics.R

### TODO: calculate network coherence using procrustes 

set.seed(4321)

# packages
library(deSolve)
library(faux)
library(tidyr)
library(MASS)

# list of parameters for Lotka-Volterra simulations
params <- list(
  
  # number of species
  S = 10, 
  
  # simulating interaction networks 
  Net_type = "predator-prey", # type of network: predator-prey, mutualistic, competition
  C = 0.1, # connectance 
  aij_params = c(0, 1), # mean and standard deviation of interaction strength
  efficiency = 0.5, # efficiency of predators to transform prey biomass 
  rho = 0, 
  
  # simulating species responses to an environmental perturbation
  NC = 0, # network coherence 
  delta_r_params = c(1, 0.2), # mean and standard deviation of the changes in species growth rates after the perturbation (squared)
  prop_neg = 0.5, # proportion of delta r dans are negative
  
  # simulating species dynamics 
  maxt = 100 # number of time steps before and after the perturbation
)


# simulate species dynamics using the above parameters and a home-made function 
source("a_code/simulate_dynamics.R")

out <- simulate_dynamics(params)


# plot species dynamics
plot_dynamic <- function(out){
  plot(y = out[,2], x = 1:nrow(out), type = "l", ylim = c(0, max(out[,c(2:ncol(out))])), xlab = "time", ylab = "Biomass")
  cols <- sample(1:ncol(out))
  for (i in 3:ncol(out)){
    lines(y = out[,i], x = 1:nrow(out), col = cols[i])
  }
}
plot_dynamic(out)


# calculate differences in biomass before and after perturbation 
delta_biomass <- sum(out[params$maxt,-1]) - sum(out[nrow(out),-1]) 
delta_biomass




