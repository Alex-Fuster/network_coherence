set.seed(4321)

library(deSolve)
library(tidyr)
library(MASS)

# What to return:
# The dynamic of the system
params <- list(
  Net_type = "predator-prey", # type of network
  nsp = 10, # number of species
  C = 0.1, # connectance
  NC = 0, # Network coherence
  r_basal = c(0.25, 0.75), # Minimum and maximum growth rate of basal species
  r_nonbasal = c(-0.5, 0), # Minimum and maximum growth rate of non-basal species
  efficiency = 0.5, # Efficiency of predator to transform prey biomass
  interaction = c(-0.5, 0), # Minimum and maximum interaction strength
  delta_r = c(0, 1), # Mean and standard deviation of the perturbation
  maxt = 1000 # Number of timestep before and after perturbation
)



# Make it easy to change throughout
# TODO make it functions
nsp <- 10
C <- 0.1
maxt <- 1000

# params: food web, b
fw.model <- function (B, params) {
  with(as.list(c(B, params)), {
    B[B < 10^-8] <- 0 # prevent numerical problems
    dBdt <- t(r + fw %*% B)*B
    list(dBdt)
  })
}

# need a function that creates realistic interaction matrix
make_metaweb <- function(S, C){
  A <- matrix(0, nrow = S, ncol = S)
  n_possible_ints <- (S*(S-1))/2
  n_ints <- round(n_possible_ints* 2 * C)
  ints <- sample(c(1:n_possible_ints), n_ints)
  A[upper.tri(A)][ints] <- 1
  return(A)
}

metaweb <- make_metaweb(nsp, C) # web of static and binary interactions

nonbasal_sp <- which(colSums(metaweb) != 0)

# growth rate positive for basal species
r <- runif(nsp, 0.25, 0.75)
# negative growth rate for non basal
r[nonbasal_sp] <- runif(length(nonbasal_sp), -0.5, 0)

efficiency <- 0.5
FW <- metaweb * runif(nsp^2, -0.5, 0)
FW <- FW - t(FW) * efficiency
diag(FW) <- runif(nsp, -0.01, 0)

parms <- list(fw = FW, r = r)

simulate_dynamics_c <- function(parms, maxt, model, init_biomass = runif(length(parms$r))){
  t <- 0 
  nsp <- length(parms$r)
  out <- matrix(c(0, init_biomass), ncol=nsp+1)
  colnames(out) <- c("time", as.character(1:nsp))
  times <- seq(from=1, to=maxt)
  out <- rbind(out[-nrow(out),],
                ode(
                func=fw.model,
                y=init_biomass,
                times=times,
                parms=parms
                ) %>%
                 as.data.frame())
    return(out)
}
pre_perturb <- simulate_dynamics_c(parms, maxt, fw.model)

# C matrix: correlation in how b changes
# coherence: correlation between C and FW
delta_r <- rnorm(nsp, 0.2, 0.001)
# C <- delta_r %*% t(delta_r)
new_r <- r + delta_r
parms$r <- new_r

equil <- as.numeric(pre_perturb[nrow(pre_perturb),-1])
post_perturb <- simulate_dynamics_c(parms, maxt, fw.mode, init_biomass = equil)

out <- rbind(pre_perturb, post_perturb)


plot_dynamic <- function(out){
  plot(y = out[,2], x = 1:nrow(out), type = "l", ylim = c(0, max(out[,c(2:ncol(out))])), xlab = "time", ylab = "Biomass")
  cols <- sample(1:ncol(out))
  for (i in 3:ncol(out)){
    lines(y = out[,i], x = 1:nrow(out), col = cols[i])
  }
}
plot_dynamic(out)

compute_NC <- function(A, delta_r){
  A_inv <- ginv(A)
  NC <- c()
  for (i in c(1:length(delta_r))){
    NC <- c(NC, cor(A_inv[i,-i], delta_r[-i]))
  }
  return(NC)
}

delta_prod <- sum(out[nrow(out),-1]) - sum(equil)

# # To see
# compute_NC <- function(A, C) {
#   pc_A <- dudi.pco(as.dist(A), scannf = FALSE, full = TRUE)
#   pc_C <- dudi.pco(as.dist(C), scannf = FALSE, full = TRUE)
#   procrustes_result <- protest(pc_A, pc_C)
#   return(procrustes_result$ss)
# }
# 
# C <- FW + rnorm(nrow(FW)*ncol(FW), mean = 0, sd = 0.1)
# 
# A <- t(as.matrix(as.dist(C)))
# diag(A) <- diag(C)
# compute_NC(FW, A)
# 
# cor(as.vector(C), as.vector(FW))
# cor(as.vector(A), as.vector(FW))
