# draft scripts for species contributions that i will not push to GH

mean_lowertri = function(x) {mean(x[which(lower.tri(x))])}
mu_all = mean_lowertri(corrs)
mu_without = corrs[-10,-10] |> mean_lowertri()

var_lowertri = function(x) {var(x[which(lower.tri(x))])}
var_all = var_lowertri(corrs)
var_without = corrs[-10,-10] |> var_lowertri()

sd_lowertri = function(x) {sd(x[which(lower.tri(x))])}
sd_all = sd_lowertri(corrs)

df = data.frame("species" = rownames(corrs),
                "mean_corr_withoutthissp" = NA,
                "var_corr_withoutthissp" = NA,
                "sd_corr_withoutthissp" = NA)
for(i in 1:npops){
  df[i, "mean_corr_withoutthissp"] = corrs[-i,-i] |> mean_lowertri()
  df[i, "var_corr_withoutthissp"] = corrs[-i,-i] |> var_lowertri()
  df[i, "sd_corr_withoutthissp"] = corrs[-i,-i] |> sd_lowertri()
}

df$species = factor(df$species, 
                    levels = df$species[order(df$mean_corr_withoutthissp)])
ggplot(data = df,
       aes(x = mean_lowertri(corrs) - mean_corr_withoutthissp, y = species)) +
  geom_vline(xintercept = 0, lwd = .2) +
  geom_point(aes(fill = mean_corr_withoutthissp), size = 3, pch = 21) +
  #geom_vline(xintercept = mean_lowertri(corrs)) +
  scale_fill_distiller(palette = "Spectral") +
  labs(fill = "Δ Community correlation", x = "Contribution to Community correlation", y = "")


df$species = factor(df$species, 
                    levels = df$species[order(df$var_corr_withoutthissp)])
ggplot(data = df,
       aes(x = var_lowertri(corrs) - var_corr_withoutthissp, y = species)) +
  geom_vline(xintercept = 0, lwd = .2) +
  geom_point(aes(fill = var_corr_withoutthissp), size = 3, pch = 21) +
  scale_fill_distiller(palette = "Spectral")


df$species = factor(df$species, 
                    levels = df$species[order(df$sd_corr_withoutthissp)])
ggplot(data = df,
       aes(x = sd_lowertri(corrs) - sd_corr_withoutthissp, y = species)) +
  geom_vline(xintercept = 0, lwd = .2) +
  geom_point(aes(fill = sd_corr_withoutthissp), size = 3, pch = 21) +
  scale_fill_distiller(palette = "Spectral")

