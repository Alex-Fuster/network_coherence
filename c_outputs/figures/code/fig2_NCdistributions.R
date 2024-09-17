library(ggplot2)
y1 <- rnorm(1000, -1, 0.4)
y2 <- rnorm(1000, 0, 0.4)
y3 <- rnorm(1000, 1, 0.4)
y4 <- rnorm(1000, -1, 1)
y5 <- rnorm(1000, 0, 1)
y6 <- rnorm(1000, 1, 1)
df <- data.frame(
  NC = as.factor(rep(c(1:6), each = 1000)),
  meanNC = as.factor(rep(c(1,0,1,0,1,0), each = 1000)),
  sdNC = as.factor(rep(c(0,0,0,1,1,1), each = 1000)),
  x = c(y1, y2, y3, y4, y5, y6)
)

ggplot(df, aes(x, group = NC, color = meanNC, linetype = sdNC)) +
  geom_density(alpha = 0.4, size = 1) +
  xlim(-3, 3) +
  theme_void() +
  scale_color_manual(values = c("#69d1c5ff", "#ca005bff")) +
  theme(legend.position = "none")

ggsave("~/Documents/Projects/network_coherence/c_outputs/figures/svg/NCdistributions.png", width = 3, height = 1)
