### simulate species dynamics using predefined parameters
### uses the function simulate_dynamics.R

### TODO: calculate network coherence using procrustes 

set.seed(4321)

# packages
library(deSolve)
library(faux)
library(tidyr)
library(MASS)
library(ggplot2)
library(ggpubr)

# list of parameters for Lotka-Volterra simulations
params <- list(
  
  # number of species
  S = 100, 
  
  # simulating interaction networks 
  Net_type = "predator-prey", # type of network: predator-prey, mutualistic, competition
  C = 0.2, # connectance 
  aij_params = c(0, 0.1), # mean and standard deviation of interaction strength
  efficiency = 0.5, # efficiency of predators to transform prey biomass 
  rho = 0, 
  
  # simulating species responses to an environmental perturbation
  NC_parms = 0, # network coherence 
  delta_r_params = c(0, 1), # mean and standard deviation of the changes in species growth rates after the perturbation (squared)
  prop_neg = 0.5, # proportion of delta r dans are negative
  
  # simulating species dynamics 
  maxt = 100 # number of time steps before and after the perturbation
)


# simulate species dynamics using the above parameters and a home-made function 
source("a_code/simulate_dynamics.R")

df <- data.frame()
for (Net_type in c("random", "predator-prey", "mutualistic", "competition")) {
  params$Net_type <- Net_type
  params$NC_parms <- c(runif(1,-1,1),runif(1,0,1))
  for (i in 1:300) {
    out <- simulate_dynamics(params)
    NC <- net_coherence3(out$delta_r, out$A)
    delta_biomass <- out$dyn[nrow(out$dyn),] - out$dyn[params$maxt,]
    extinctions <- sum(out$dyn[nrow(out$dyn),-1] < 0.01)
    df <- rbind(df,
                data.frame(Net_type = Net_type, NC_mean = mean(NC, na.rm = T), NC_sd = sd(NC, na.rm = T), delta_biom_sum = sum(delta_biomass), delta_biom_sd = sd(delta_biomass), extinctions = extinctions))
  }
}



# pars for plotting:

my_theme = theme(axis.text=element_text(size=12),
                axis.title = element_text(size = 14),
                legend.text=element_text(size=14),
                legend.title = element_text(size=14),
                plot.title = element_text(face="bold",size=14,margin=margin(0,0,20,0),hjust = 0.5),
                axis.title.y = element_text(hjust = 0.5),
                axis.title.x = element_text(hjust = 0.5))


ggplot(df) +
  geom_point(aes(NC_mean, NC_sd, colour = Net_type), alpha = 0.25) + 
  labs(x = "Mean ENC", y = "SD ENC", colour = "Network type")+
  theme_classic()+
  my_theme

ggsave("c_outputs/figures/100_sd_mean.png", height = 5, width = 6)

p1 <- ggplot(df) +
  geom_point(aes(NC_mean, delta_biom_sum, colour = Net_type), alpha = 0.25) +
  geom_smooth(aes(NC_mean, delta_biom_sum, colour = Net_type), se = F) +
  labs(x = "Mean ENC", y = expression(paste(italic(Delta)," total biomass")),
         colour = "Network type") +
  theme_classic()+
  my_theme


p2 <- ggplot(df) +
  geom_point(aes(NC_mean, delta_biom_sd, colour = Net_type), alpha = 0.25) +
  geom_smooth(aes(NC_mean, delta_biom_sd, colour = Net_type), se = F) +
  scale_y_continuous(trans = "log", breaks = c(0.2, 1, 5), limits = c(0.2, 10)) +
  labs(x = "Mean ENC", y = "SD in biomass", colour = "Network type") +
  theme_classic()+
  my_theme


p3 <- ggplot(df) +
  geom_point(aes(NC_mean, extinctions, colour = Net_type), alpha = 0.25) +
  geom_smooth(aes(NC_mean, extinctions, colour = Net_type), se = F) +
  labs(x = "Mean ENC", y = "Number of extinctions", colour = "Network type") +
  theme_classic()+
  my_theme


p4 <- ggplot(df) +
  geom_point(aes(NC_sd, delta_biom_sum, colour = Net_type), alpha = 0.25) +
  geom_smooth(aes(NC_sd, delta_biom_sum, colour = Net_type), se = F) +
  labs(x = "SD ENC", y = expression(paste(italic(Delta)," total biomass")), colour = "Network type") +
  theme_classic()+
  my_theme


p5 <- ggplot(df) +
  geom_point(aes(NC_sd, delta_biom_sd, colour = Net_type), alpha = 0.25) +
  geom_smooth(aes(NC_sd, delta_biom_sd, colour = Net_type), se = F) +
  scale_y_continuous(trans = "log", breaks = c(0.2, 1, 5), limits = c(0.2, 10)) +
  labs(x = "SD ENC", y = "SD in biomass", colour = "Network type") +
  theme_classic()+
  my_theme


p6 <- ggplot(df) +
  geom_point(aes(NC_sd, extinctions, colour = Net_type), alpha = 0.25) +
  geom_smooth(aes(NC_sd, extinctions, colour = Net_type), se = F) +
  labs(x = "SD ENC", y = "Number of extinctions", colour = "Network type") +
  theme_classic()+
  my_theme

#(p1 + p2 + p3) / (p4 + p5 + p6) + plot_layout(guides = "collect")


ggarrange(
  p1, p2, p3,
  p4, p5, p6,
  
  nrow = 2,
  ncol = 3,
  
  labels = LETTERS[1:6],
  
  common.legend = TRUE
)


ggsave("c_outputs/figures/LVsimulations_cov_100sp.png", height = 9, width = 10)



params$Net_type <- "random"
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


