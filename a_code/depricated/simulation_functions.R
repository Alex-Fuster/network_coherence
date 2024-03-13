##TODO
# Different type of network: random, mutualistic, competition, niche model
# Clean a little: turn into function + documentation
# Make sure the parametrization of LV model is adequate
# Implement the control of NC using faux::rnorm_pre()

# function simulating the dynamic
simulate_dynamic <- function(params){
  with(params, {
    if (Net_type == "predator-prey") {
      metaweb <- make_metaweb(nsp, C) # create a binary matrix with only the predator->prey interaction
      nonbasal_sp <- which(colSums(metaweb) != 0) # get which species are basal (no prey)
      # growth rate positive for basal species
      r <- runif(nsp, r_basal[1], r_basal[2])
      # negative growth rate for non basal
      r[nonbasal_sp] <- runif(length(nonbasal_sp), r_nonbasal[1], r_nonbasal[2])
    } else {
      stop("Incorrect network type")
    }
    # add interactions strength to matrix
    Net <- add_weights(metaweb, Net_type, interaction[1], interaction[2], efficiency)
    dyn_params <- list(fw = Net, r = r)
    pre_perturb <- simulate_dynamics_c(dyn_params, maxt, fw.model)
    delta_r <- get_delta_r(Net, NC, delta_r[1], delta_r[2])
    new_r <- r + delta_r
    dyn_params$r <- new_r
    equil <- as.numeric(pre_perturb[nrow(pre_perturb),-1])
    post_perturb <- simulate_dynamics_c(dyn_params, maxt, fw.mode, init_biomass = equil)
    out <- rbind(pre_perturb, post_perturb)
  })
  return(out)
}

add_weights <- function(web, Net_type, min_weight, max_weight, efficiency = NULL) {
  if (Net_type == "random"){
    
  } else if (Net_type == "predator-prey") {
    Net <- metaweb * runif(nrow(metaweb)^2, min_weight, max_weight)
    Net <- Net - t(Net) * efficiency
    diag(Net) <- runif(nsp, -0.01, 0)
  } else {
    stop("Incorrect network type")
  }
  return(Net)
}

# function creating the binary matrix of interaction
make_metaweb <- function(S, C){
  A <- matrix(0, nrow = S, ncol = S)
  n_possible_ints <- (S*(S-1))/2
  n_ints <- round(n_possible_ints* 2 * C)
  ints <- sample(c(1:n_possible_ints), n_ints)
  A[upper.tri(A)][ints] <- 1
  return(A)
}

# LV model
fw.model <- function (t, B, params) {
  with(as.list(c(B, params)), {
    B[B < 10^-8] <- 0 # prevent numerical problems
    dBdt <- t(r + fw %*% B)*B
    list(dBdt)
  })
}

# get delta_r based on NC
get_delta_r <- function(Net, NC, mean, sd) {
  rnorm(nrow(nsp), mean, sd)
}

# build interaction matrix
build_A <- function(S, C, d, sigma, rho, Net_type = "random"){
  # sample coefficients in pairs
  interactions <- MASS::mvrnorm(n = S * (S-1) / 2,
                                mu = c(0, 0),
                                Sigma = sigma^2 * matrix(c(1, rho, rho, 1), 2, 2))
  # build a completely filled matrix
  A <- matrix(0, S, S)
  A[upper.tri(A)] <- interactions[,1]
  A <- t(A)
  A[upper.tri(A)] <- interactions[,2]
  if (Net_type == "mutualistic"){
    A[A<0] <- -A[A<0]
    } else if (Net_type == "competition"){
      A[A>0] <- -A[A>0]
      }
  # determine which connections to retain
  Connections <- (matrix(runif(S * S), S, S) <= C) * 1 
  Connections[lower.tri(Connections)] <- 0
  diag(Connections) <- 0
  Connections <- Connections + t(Connections)
  A <- A * Connections
  diag(A) <- -d
  return(A)
}
