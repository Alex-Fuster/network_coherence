library(MASS)
library(tibble)
library(ggplot2)
library(tidyr)

# There are 3 entities: mu, sd, the correlation matrix. What do they do?
# I think mu would control the average strength of the perturbance and sd among species
# The correlation matrix controls the correlation between responses of species.

createCorMat <- function(S, type){
  if (type == "positive"){
    cor <- runif(S, 0, 1)
  }
  else if (type == "mixed"){
    cor <- runif(S, -1, 1)
  } else {
    stop("type not defined")
  }
  corMat <- matrix(1, nrow = S, ncol = S)
  corMat[upper.tri(corMat)] <- cor
  corMat <- corMat * t(corMat)
}

plot_matrix <- function(CovMatrix, title=""){
  colnames(CovMatrix) <- rownames(CovMatrix) <- paste0("sp", c(1:nrow(CovMatrix)))
  CovMatrix %>%
    as_tibble(rownames  = NA) %>%
    rownames_to_column("Var1") %>%
    pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
    ggplot(aes(Var1, Var2)) +
    geom_tile(aes(fill = value)) +
    scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(-1,1)*max(abs(CovMatrix))) +
    labs(title = title, x= "", y = "")
}

# positive mu, low sd, positive corMat
S <- 5
mu <- rnorm(S, 0, 1)
sd <- rexp(S, 1)

corMat <- createCorMat(S, type = "positive")

covMat <- as.matrix(Matrix::nearPD(sd %*% t(sd) * corMat)$mat)

delta_r <- mvrnorm(n = 1, mu = mu, Sigma = covMat, empirical = FALSE)
mean(delta_r)
sd(delta_r)
plot_covmatrix(covMat)

# playing with sd, it seems like it change the covariances. The larger sd, the larger the covariance.
# Makes sense since cov(X,Y) = sd_X*sd_Y*cor(X,Y)
# The smaller sd, the closer to the mean responses (mu) each species will have
sd <- rexp(S, 50)
covMat <- sd %*% t(sd) * corMat
delta_r <- mvrnorm(n = 1, mu = mu, Sigma = covMat, empirical = FALSE)

plot_covmatrix(covMat)


