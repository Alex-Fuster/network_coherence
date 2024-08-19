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
  # we divide by two because we only fill the upper triangle
  L <- round(C/2 * S^2)

  # sample interactions (all on the upper triangle to facilitate calculations)
  L_ij <- sample(c(1:floor((S^2-S)/2)), L)

  # add binary interactions in the adjacency matrix 
  A[upper.tri(A)][L_ij] <- 1

  return(A)
}



### step 2: simulate a quantitative interaction matrix using the following arguments:
### Net_type: network type (random, predator-prey, competition, mutualistic) 
### S: number of species 
### C: connectance 
### rho correlation between aij and aji
### parametrization inspired from Allesina & Tang 2012 (https://doi.org/10.1038/nature10832)

sim_quantitative_network <- function(Net_type, S, C, aij_params) {
  
  # initialize network
  A <- matrix(0, S, S)
  
  # number of species pairs
  n_pairs <- S*(S-1)/2
  # which ones are interacting
  B <- runif(n_pairs) <= C
  
  # simulate interspecific coefficients in pairs depending on the network type
  if (Net_type == "random") {

    # add independent interaction strengths
    A[upper.tri(A)] <- B * rnorm(n_pairs, aij_params[1], aij_params[2])
    A <- t(A)
    A[upper.tri(A)] <- B * rnorm(n_pairs, aij_params[1], aij_params[2])
    
  } else if (Net_type == "predator-prey") { 
    
    # for predator-prey network, the effect of prey on predator is positive and the effect of
    # predator on prey is negative.
    
    # effects of predator on prey:
    aij <- -abs(rnorm(n_pairs, aij_params[1], aij_params[2]))
    A[upper.tri(A)] <- B * aij
    A <- t(A)
    # effects of prey on predator:
    aij <- abs(rnorm(n_pairs, aij_params[1], aij_params[2]))
    A[upper.tri(A)] <- B * aij
    
  } else if (Net_type == "competition") {
    
    # interactions are drawn in a negative half normal distribution
    aij <- -abs(rnorm(n_pairs*2, aij_params[1], aij_params[2]))
    A[upper.tri(A)] <- B * aij[1:n_pairs]
    A <- t(A)
    A[upper.tri(A)] <- B * aij[(n_pairs+1):length(aij)]
    
  } else if (Net_type == "mutualistic") {
    
    # interactions are drawn in a negative half normal distribution
    aij <- abs(rnorm(n_pairs*2, aij_params[1], aij_params[2]))
    A[upper.tri(A)] <- B * aij[1:n_pairs]
    A <- t(A)
    A[upper.tri(A)] <- B * aij[(n_pairs+1):length(aij)]
    
  } else { stop("Incorrect network type") }
  diag(A) <- -(max(Re(eigen(A)$values)) + runif(S, 0.1))
  
  # make sure that the equilibrium is stable
  while(max(Re(eigen(A)$values)) > 0){
    diag(A) <- -(max(Re(eigen(A)$values)) + runif(S, 0.1))
  }
  
  return(A)
}
    


### step 3: simulate the effect of the perturbation on species intrinsic growth rates using the following arguments:
### A: quantitative adjacency matrix 
### NC: network coherence (NC = cor(A-1, C matrix))
### delta_r_params: mean and standard deviation of the differences in intrinsic growth rate (squared)
### prop_neg: proportion of delta r dans are negative


# sim_delta_r <- function(A, NC_parms, delta_r_params) {
# 
#   # inverse the matrix of quantitative interactions
#   A_inv <- solve(A)
#   n <- nrow(A)
#   # simulate NCs
#   NC <- rnorm(n, NC_parms[1], NC_parms[2])
#   delta_r <- rep(0,n)
#   for (i in 1:n){
#     delta_r[i] <- sum(NC[-i] * A_inv[-i,i])
#   }
#   
#   # center and scale
#   delta_r <- delta_r - mean(delta_r) + delta_r_params[1]
#   delta_r <- delta_r/sd(delta_r) * delta_r_params[2]
#   
#   return(delta_r)
#}

sim_delta_r <- function(A, delta_r_params, cov_type = "weak") {
  n <- nrow(A)
  delta_r <- rep(0, n)
  
  if (cov_type == "weak_negative") {
    # Small perturbations, mostly around zero
    delta_r <- -abs(rnorm(n, delta_r_params[1], 0.1))
    
  }else if (cov_type == "weak_positive") {
    # Larger perturbations
    delta_r <- abs(rnorm(n, delta_r_params[1], 0.1))
    
  } else if (cov_type == "strong_negative") {
    # Larger perturbations
    delta_r <- -abs(rnorm(n, delta_r_params[1], 1))
    
  } else if (cov_type == "strong_positive") {
    # Strong positive perturbations
    delta_r <- abs(rnorm(n, delta_r_params[1], 1))
    
  } else if (cov_type == "strong_mixed") {
    # Strong positive and negative perturbations
    delta_r <- sample(c(-1, 1), n, replace = TRUE) * rnorm(n, delta_r_params[1], 1)
    
  } else if (cov_type == "weak_mixed") {
    # Small mixed perturbations
    delta_r <- sample(c(-0.1, 0.1), n, replace = TRUE) * rnorm(n, delta_r_params[1], 0.1)
  }
  
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

# simulate_dynamics <- function(params, model = fw.model) {
# 
#     with(params, {
#       
#       # simulate quantitative adjacency matrix 
#       A <- sim_quantitative_network(Net_type, S, C, aij_params) 
#       
#       # define initial intrinsic growth rates of species
#       # to make sure all species have a positive biomass at equilibrium
#       Biomass_at_equilibrium <- runif(S, 1, 10)
#       r <- -A%*%Biomass_at_equilibrium
#                  
#       # simulate parameters before the perturbation 
#       dyn_params <- list(A = A, r = r, S = S)
#       pre_perturb <- simulate_dynamics_c(dyn_params, fw.model)
#       
#       # get delta r
#       delta_r <- sim_delta_r(A, NC_parms, delta_r_params) 
#       #delta_r <- rnorm(S, delta_r_params[1], delta_r_params[2]) 
# 
#       # calculate new r
#       new_r <- r + delta_r
#       dyn_params$r <- new_r
#       
#       #########################################
#       
#       delta_r_vectors <- cbind(r, new_r)
#       C <- compute_cov_matrix(delta_r_vectors) # compute C matrix from delta_r vectors
#       
#       ##########################################
#       
#       # find values at equilibrium 
#       equil <- as.numeric(pre_perturb[nrow(pre_perturb), -1])
#       
#       # simulate dynamics after the perturbation 
#       post_perturb <- simulate_dynamics_c(dyn_params, fw.model, init_biomass = equil)
#       
#       dyn <- rbind(pre_perturb, post_perturb)
#       dyn[dyn < 0] <- 0
#       
#       out <- list(dyn = dyn,
#                   A = A,
#                   r = r,
#                   new_r = new_r,
#                   delta_r = delta_r,
#                   C = C)
#       
#     return(out)
#     
#     })
# }


simulate_dynamics <- function(params, cov_type = "weak", model = fw.model) {
  with(params, {
    
    # Simulate quantitative interaction matrix
    A <- sim_quantitative_network(Net_type, S, C, aij_params)
    
    # Define initial intrinsic growth rates (r) of species
    Biomass_at_equilibrium <- runif(S, 1, 10)
    r <- -A %*% Biomass_at_equilibrium
    
    # Simulate dynamics before the perturbation
    dyn_params <- list(A = A, r = r, S = S)
    pre_perturb <- simulate_dynamics_c(dyn_params, fw.model)
    
    # Simulate perturbation (delta_r) based on desired covariance structure
    delta_r <- sim_delta_r(A, delta_r_params, cov_type = cov_type)
    
    # Calculate new growth rates (r_new)
    new_r <- r + delta_r
    dyn_params$r <- new_r
    
    # Compute covariance matrix based on r and r_new
    delta_r_vectors <- cbind(r, new_r)
    C <- compute_cov_matrix(delta_r_vectors)  # This will capture the covariance based on delta_r
    
    # Find equilibrium values and simulate dynamics after perturbation
    equil <- as.numeric(pre_perturb[nrow(pre_perturb), -1])
    post_perturb <- simulate_dynamics_c(dyn_params, fw.model, init_biomass = equil)
    
    dyn <- rbind(pre_perturb, post_perturb)
    dyn[dyn < 0] <- 0
    
    # Return simulation output
    out <- list(dyn = dyn, A = A, r = r, new_r = new_r, delta_r = delta_r, C = C)
    
    return(out)
  })
}



######################################

compute_cov_matrix <- function(delta_r_vectors) {
  S <- nrow(delta_r_vectors)
  cov_matrix <- matrix(NA, nrow = S, ncol = S)
  for (i in 1:S) {
    for (j in 1:S) {
      cov_matrix[i, j] <- cov(delta_r_vectors[i, ], delta_r_vectors[j, ])
    }
  }
  return(cov_matrix)
}


# Convert A matrix to binary matrix
convert_to_binary_matrix <- function(A_matrix) {
  binary_A_matrix <- A_matrix
  binary_A_matrix[binary_A_matrix > 0] <- 1
  binary_A_matrix[binary_A_matrix < 0] <- 1
  return(binary_A_matrix)
}

# Calculate ENC pattern by element-wise multiplication of binary A and C matrices
calculate_enc_pattern <- function(A_matrix, C_matrix) {
  binary_A_matrix <- convert_to_binary_matrix(A_matrix)
  ENC_matrix <- binary_A_matrix * C_matrix
  return(ENC_matrix)
}

######################################
