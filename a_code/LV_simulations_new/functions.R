createCorMat <- function(S, type){
  if (type == "positive"){
    cor <- runif((S*(S-1)/2), 0, 1)
  }
  else if (type == "mixed"){
    cor <- runif((S*(S-1)/2), -1, 1)
  } else {
    stop("type not defined")
  }
  corMat <- matrix(1, nrow = S, ncol = S)
  corMat[upper.tri(corMat)] <- cor
  corMat <- corMat * t(corMat)
  return(corMat)
}

plot_matrix <- function(CovMatrix, title=""){
  species <- paste0("sp", c(1:nrow(CovMatrix)))
  colnames(CovMatrix) <- rownames(CovMatrix) <- species
  CovMatrix %>%
    as_tibble(rownames  = NA) %>%
    rownames_to_column("Var1") %>%
    pivot_longer(-Var1, names_to = "Var2", values_to = "value") %>%
    mutate(Var1 = factor(Var1,levels= species), Var2 = factor(Var2, levels = species)) %>%
    ggplot(aes(Var1, Var2)) +
    geom_tile(aes(fill = value)) +
    scale_fill_distiller(palette = "RdBu", direction = 1, limits = c(-1,1)*max(abs(CovMatrix))) +
    labs(title = title, x= "", y = "")
}

### step 2: simulate a quantitative interaction matrix using the following arguments:
### Net_type: network type (random, predator-prey, competition, mutualistic) 
### S: number of species 
### C: connectance 
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
  diag(A) <- -runif(S, min = 0, max = 1)
  
  # make sure that the equilibrium is stable
  while(max(Re(eigen(A)$values)) > 0){
    diag(A) <- -runif(S, min = 0, max = 1)
  }
  
  return(A)
}
