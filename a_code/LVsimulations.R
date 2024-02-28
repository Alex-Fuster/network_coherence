library(deSolve)
library(vegan)
library(ade4)

# Make it easy to change throughout
# TODO make it functions
nsp <- 10
nbasal <- 4

# params: food web, b
fw.model <- function (t, B, params) {
  with(as.list(c(B, params)), {
    B[B < 10^-8] <- 0 # prevent numerical problems
    dBdt <- t(r + fw %*% B)*B
    list(dBdt)
  })
}

r <- c(runif(2, 0.25, 0.75), runif(3, -0.5, 0))

metaweb <- matrix(c(rep(0,10), 1, 1, rep(0,3), 1, 1, rep(0,5), 1, 1,0), nrow =5 ) # web of static and binary interactions
efficiency <- 0.5
FW <- metaweb * runif(25, -0.5, 0)
FW <- FW - t(FW) * efficiency
diag(FW) <- runif(5, -0.01, 0)

parms <- list(fw = FW, r = r)

simulate_dynamics_c <- function(parms, maxt, model, init_biomass = runif(length(parms$r))){
  t <- 0 
  nsp <- length(parms$r)
  out <- matrix(c(0, init_biomass), ncol=6)
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
pre_perturb <- simulate_dynamics_c(parms, 100, fw.model)

# C matrix: correlation in how b changes
# coherence: correlation between C and FW
delta_r <- rnorm(5, 0, 0.1)
C <- delta_r %*% t(delta_r)
new_r <- r + delta_r
parms$r <- new_r

equil <- as.numeric(pre_perturb[nrow(pre_perturb),-1])
post_perturb <- simulate_dynamics_c(parms, 100, fw.mode, init_biomass = equil)

out <- rbind(pre_perturb, post_perturb)

plot(y = out[,2], x = 1:200, type = "l", ylim = c(0, max(out[,c(2:6)])))
lines(y = out[,3], x = 1:200, col = "red")
lines(y = out[,4], x = 1:200, col = "blue")
lines(y = out[,5], x = 1:200, col = "violet")
lines(y = out[,6], x = 1:200, col = "green")

# To see
compute_NC <- function(A, C) {
  pc_A <- dudi.pco(as.dist(A), scannf = FALSE, full = TRUE)
  pc_C <- dudi.pco(as.dist(C), scannf = FALSE, full = TRUE)
  procrustes_result <- protest(pc_A, pc_C)
  return(procrustes_result$ss)
}

C <- FW + rnorm(nrow(FW)*ncol(FW), mean = 0, sd = 0.1)

A <- t(as.matrix(as.dist(C)))
diag(A) <- diag(C)
compute_NC(FW, A)

cor(as.vector(C), as.vector(FW))
cor(as.vector(A), as.vector(FW))
