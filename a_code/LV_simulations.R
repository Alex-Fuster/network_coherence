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
  delta_r_params = c(0, 1), # mean and standard deviation of the changes in species growth rates after the perturbation (squared)
  prop_neg = 0.5, # proportion of delta r dans are negative
  
  # simulating species dynamics 
  maxt = 100 # number of time steps before and after the perturbation
)


# simulate species dynamics using the above parameters and a home-made function 
source("a_code/simulate_dynamics.R")

out <- simulate_dynamics(params)


# plot species dynamics
plot_dynamic <- function(dyn){
  plot(y = dyn[,2], x = 1:nrow(dyn), type = "l", ylim = c(0, max(dyn[,c(2:ncol(dyn))])), xlab = "time", ylab = "Biomass")
  cols <- sample(1:ncol(dyn))
  for (i in 3:ncol(dyn)){
    lines(y = dyn[,i], x = 1:nrow(dyn), col = cols[i])
  }
}
plot_dynamic(out$dyn)
NC(out$delta_r, out$A)

# calculate differences in biomass before and after perturbation 
delta_biomass <- sum(out[params$maxt,-1]) - sum(out[nrow(out),-1]) 
delta_biomass

# repeat simulation 100 times
df <- c()
for (i in 1:1000) {
  out <- simulate_dynamics(params)
  NC <- net_coherence(out$delta_r, out$A)
  delta_biomass <- out$dyn[params$maxt,-1] - out$dyn[nrow(out$dyn),-1]
  extinctions <- sum(out$dyn[nrow(out$dyn),-1] < 0.01)
  df <- rbind(df,
              c(NC = NC, delta_biom_sum = sum(delta_biomass), delta_biom_sd = sd(delta_biomass), extinctions = extinctions))
}


