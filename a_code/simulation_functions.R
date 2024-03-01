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
      
      # add interactions strength to matrix
      efficiency <- 0.5
      Net <- metaweb * runif(nsp^2, -0.5, 0)
      Net <- Net - t(Net) * efficiency
      diag(Net) <- runif(nsp, -0.01, 0)
    } else {
      stop("Incorrect network type")
    }
    dyn_params <- list(fw = Net, r = r)
    pre_perturb <- simulate_dynamics_c(dyn_params, maxt, fw.model)
    delta_r <- rnorm(nsp, delta_r[1], delta_r[2])
    new_r <- r + delta_r
    parms$r <- new_r
    equil <- as.numeric(pre_perturb[nrow(pre_perturb),-1])
    post_perturb <- simulate_dynamics_c(parms, maxt, fw.mode, init_biomass = equil)
    out <- rbind(pre_perturb, post_perturb)
  })
  return(out)
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
