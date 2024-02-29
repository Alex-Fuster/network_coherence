# Demonstrate Dom's equations with the fish case study
# (using the grand banks ecosim model)

# read the data ----

# interaction effect matrix
A_grand_banks <- read.table("data/A_grand_banks.txt")
# growth rates
delta_r <- read.table("data/B_grand_banks.txt")

# take the inverse of the interaction effect matrix ----
inv_A = solve(as.matrix(A_grand_banks)) # mid-90s food web

# plot the matrices
heatmap(as.matrix(A_grand_banks), Rowv = NA, Colv = NA)
heatmap(inv_A, Rowv = NA, Colv = NA)

# number of species
n = length(delta_r[,1])

# calculate change in abundance based on interaction effect and abundances
species_effect = c()
cov_effect = c()
delta_X = c()
self_effect = c()
for(i in 1:n){
  
  # take the alphas (interaction effects) of all species except species i on itself
  alphas_ij = inv_A[i,-i]
  
  # calculate the covariance between interaction effects and species' responses (r)
  cov_effect[i] = cov(alphas_ij, delta_r[-i,1])
  
  # calculate the effect of all other species on the abundance of species i
  species_effect[i] = (n-1)*(mean(alphas_ij) + mean(delta_r[-i,1]) + cov_effect[i])

  # calculate the change in abundance of species i, depending on the effect of itself
  # and all other species in the community
  self_effect[i] = inv_A[i,i]*delta_r[i,1]
  
  delta_X[i] = self_effect[i] + species_effect[i]
}


# plot the values (X axis is species)
par(mfrow = c(1,3))

plot(cov_effect, pch = 16, 
     main = "Covariance effect per species", 
     xlab = "species")
abline(h = 0)

plot(species_effect, pch = 16, 
     main = "Other species' net effect on each species' abundance",
     xlab = "species")
abline(h = 0)

plot(delta_X, pch = 16, 
     main = "Change in abundance per species (delta X)",
     xlab = "species")
abline(h = 0)

par(mfrow=c(1,1))
plot(self_effect ~ species_effect, ylim = c(-2000, 500))
abline(v = 0)
abline(h = 0)


# interaction effects
temp = inv_A
diag(temp) = 0
interaction_effet_mean = apply(temp, 1, mean)
interaction_effet_var = apply(temp, 1, var)

temp = A_grand_banks
diag(temp) = 0
interaction_effet_mean = apply(temp, 1, mean)
interaction_effet_var = apply(temp, 1, var)

bnA = A_grand_banks*0
bnA[A_grand_banks>0] = 1
deg = apply(bnA,2,sum)
deg

plot(interaction_effet_mean ~ deg, ylim = c(-1, 1))
