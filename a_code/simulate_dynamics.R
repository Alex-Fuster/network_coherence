### simulate species dynamics before and after an environmental perturbation using the concept of network coherence 
### function currently implemented for unipartite networks only 

### TODO: relax the assumption of symmetry in the effect of one species on another 
### see Dom's code which simulate coefficients from a multivariate normal distribution with a fixed covariance between species effects 

### TODO: find a better way to simulate delta r whose correlation with the inverse of the adjacency matrix is network coherence

### TODO: double-check the lotka-volterra model for potential mistakes 

####################
  

### step 1: simulate a binary interaction matrix using the following arguments:
### S: number of species 
### C: connectance (proportion of links that are realized)

sim_binary_network <- function(S, C) {
  
  # empty binary matrix (square matrix because we simulate an unipartite network)
  A <- matrix(0, nrow = S, ncol = S)

  # number of interactions that are realized (C = L/S^2)
  L <- round(C * S^2)

  # sample interactions (all on the upper triangle to facilitate calculations)
  L_ij <- sample(c(1:floor((S^2-S)/2)), L)

  # add binary interactions in the adjacency matrix 
  A[upper.tri(A)][L_ij] <- 1

  return(A)
}



### step 2: simulate a quantitative interaction matrix using the following arguments:
### Net_type: network type (predator-prey, competition, mutualistic) 
### S: number of species 
### C: connectance 
### rho correlation between aij and aji
### parametrization inspired from Allessina & Tang 2012 (https://doi.org/10.1038/nature10832)
### efficiency: energy efficiency of predators (for predator-prey networks)

sim_quantitative_network <- function(Net_type, S, C, aij_params, efficiency = 1, rho = 0) {
  
  # simulate binary network
  B <- sim_binary_network(S, C)
  
  # simulate interspecific coefficients in pairs depending on the network type
  if (Net_type == "random") {
    
    # for random network, interactions are drawn in pairs, but the sign does not matter
    interactions <- MASS::mvrnorm(n = S * (S-1) / 2,
                                  mu = c(aij_params[1], aij_params[1]),
                                  Sigma = aij_params[2]^2 * matrix(c(1, rho, rho, 1), 2, 2))
    
  } else if (Net_type == "predator-prey") { 
    
    # for predator-prey network, the effect of prey on predator is positive and the effect of
    # predator on prey is negative.
    # there is a efficiency in transforming prey biomass into predator biomass
    # if i eats j: aji = -efficiency * aij 
    
    # effects of predator on prey:
    interactions <- -abs(rnorm(n = S * (S-1) / 2 , mean = aij_params[1], sd = aij_params[2]))
    # effects of prey on predator:
    interactions <- cbind(interactions, interactions * -efficiency)
    
  } else if (Net_type == "competition") {
    
    # interactions are drawn in pairs and we put everything as negative
    
    # for random network, interactions are drawn in pairs, but the sign does not matter
    interactions <- MASS::mvrnorm(n = S * (S-1) / 2,
                                  mu = c(aij_params[1], aij_params[1]),
                                  Sigma = aij_params[2]^2 * matrix(c(1, rho, rho, 1), 2, 2))
    interactions <- -abs(interactions)
    
    
  } else if (Net_type == "mutualistic") {
    
    # interactions are drawn in pairs and we put everything as positive
    
    # for random network, interactions are drawn in pairs, but the sign does not matter
    interactions <- MASS::mvrnorm(n = S * (S-1) / 2,
                                  mu = c(aij_params[1], aij_params[1]),
                                  Sigma = aij_params[2]^2 * matrix(c(1, rho, rho, 1), 2, 2))
    interactions <- abs(interactions)
    
  } else { stop("Incorrect network type") }
  
  # build a completely filled matrix
  A <- matrix(0, S, S)
  A[upper.tri(A)] <- interactions[,1]
  A <- t(A)
  A[upper.tri(A)] <- interactions[,2]
  
  # remove non-interactins
  A <- A * (B + t(B))
  
  # make sure that the equilibrium is stable
  diag(A) <- sqrt(C*S*sd(A)) * (1 + rho) - 1
  
  return(A)
}
    


### step 3: simulate the effect of the perturbation on species intrinsic growth rates using the following arguments:
### A: quantitative adjacency matrix 
### NC: network coherence (NC = cor(A-1, C matrix))
### delta_r_params: mean and standard deviation of the differences in intrinsic growth rate (squared)
### prop_neg: proportion of delta r dans are negative

sim_delta_r <- function(A, NC, delta_r_params, prop_neg) {
  
  # inverse the matrix of quantitative interactions
  A_inv <- solve(A)

  # simulate the vector of squared co-responses to the environment (diagonal of the C matrix) 
  # the elements of the C matrix are the pairwise products of delta r
  # the diagonals elements are the squared delta r for each species
  # simulating the diagonal of the C matrix that correlates with the one of the A matrix allows us to respect the structure of the C matrix 
  C_diag <- faux::rnorm_pre(x = diag(A_inv), mu = delta_r_params[1], sd = delta_r_params[2], r = NC)

  # calculate species delta r and assign negative values randomly
  delta_r <- sqrt(C_diag)
  delta_r <- delta_r * ifelse(rbinom(length(delta_r), size = 1, p = prop_neg) == 1, -1, 1)
  
  return(delta_r)
}


### step 4: define lotka-volterra model using the following arguments:
### t: time step
### A: quantitative adjacency matrix
### params: list of parameters for all of the subfunctions

fw.model <- function(t, B, params) {
  with(as.list(c(B, params)), {
    B[B < 10^-8] <- 0 # prevent numerical problems
    dBdt <- t(r + A %*% B)*B
    list(dBdt)
  })
}


### step 5: define dynamical function of the system with the following arguments
### params: list of parameters
### model: lotka-volterra model (defined at step 4)
### init_biomass: initial biomasses of each species

simulate_dynamics_c <- function(parms, model, init_biomass = runif(length(parms$r), min = 1, max = 10)){
 
  # define matrix of initial biomasses
  out <- matrix(c(0, init_biomass), ncol = params$S + 1)
  colnames(out) <- c("time", as.character(1:params$S))
  
  # define time steps 
  times <- seq(from = 1, to = params$maxt)
  
  # dynamical model 
  out <- rbind(out[-nrow(out),],
               ode(
                 func = fw.model,
                 y = init_biomass,
                 times = times,
                 parms = parms
               ) %>%
                 as.data.frame())
  return(out)
}


### step 6: wrapping function with the following arguments:
### params: list of parameters for all of the subfunctions
### model: food-web model (defined at step 4)

simulate_dynamics <- function(params, model = fw.model) {

    with(params, {
      
      # simulate quantitative adjacency matrix 
      A <- sim_quantitative_network(Net_type, S, C, aij_params, efficiency, rho) 
      
      # define initial intrinsic growth rates of species
      # to make sure all species have a positive biomass at equilibrium
      Biomass_at_equilibrium <- runif(S, 1, 10)
      r <- -A%*%Biomass_at_equilibrium
                 
      # simulate parameters before the perturbation 
      dyn_params <- list(A = A, r = r, S = S)
      pre_perturb <- simulate_dynamics_c(dyn_params, fw.model)
      
      # get delta r
      delta_r <- sim_delta_r(A, NC, delta_r_params, prop_neg) 
        
      # calculate new r
      new_r <- r + delta_r
      dyn_params$r <- new_r
      
      # find values at equilibrium 
      equil <- as.numeric(pre_perturb[nrow(pre_perturb), -1])
      
      # simulate dynamics after the perturbation 
      post_perturb <- simulate_dynamics_c(dyn_params, fw.model, init_biomass = equil)
      
      out <- rbind(pre_perturb, post_perturb)
      out[out < 0] <- 0
      
    return(out)
    
    })
}

