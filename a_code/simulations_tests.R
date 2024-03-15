# testing function to see where it fails

for (i in c(1:500)) {
  sim_quantitative_network("random", 10, 0.2, c(0,1)) 
}

while (max(pre_perturb) == 100) {
  # simulate quantitative adjacency matrix 
  A <- sim_quantitative_network(Net_type, S, C, c(0,0.5), efficiency, rho) 
  
  # define initial intrinsic growth rates of species
  # to make sure all species have a positive biomass at equilibrium
  Biomass_at_equilibrium <- runif(S, 1, 10)
  r <- -A%*%Biomass_at_equilibrium
  
  # simulate parameters before the perturbation 
  dyn_params <- list(A = A, r = r, S = S)
  pre_perturb <- simulate_dynamics_c(dyn_params, fw.model)
  
}

df<-data.frame()
for (i in c(1:1000)) {
  Net_type <- sample(c("random", "predator-prey", "mutualistic", "competition"),1)
  # simulate quantitative adjacency matrix 
  A <- sim_quantitative_network(Net_type, 25, 0.25, c(0,0.5)) 
  NC_parms <- c(runif(1,-5,5), runif(1,0,1))
  delta_r_params<- c(0,1)
  
  delta_r <- sim_delta_r(A, NC_parms, delta_r_params)
 
  NC <- net_coherence2(delta_r, A)
  
  df <- rbind(df,
              data.frame(Net_type, expected_mean = NC_parms[1], expected_sd = NC_parms[2], observed_mean = mean(NC, na.rm = T), observed_sd = sd(NC, na.rm = T))
  ) 
}
ggplot(df) +
  geom_point(aes(x = expected_mean, y = observed_mean, colour = Net_type))

ggplot(df) +
  geom_point(aes(x = expected_sd, y = observed_sd, colour = Net_type))
