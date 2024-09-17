library(ggplot2)
y1 <- rnorm(1000, 5, 4)
y2 <- rnorm(1000, 6, 2.5)
y3 <- rnorm(1000, 7, 2.5)
y4 <- rnorm(1000, 2, 2)
df <- data.frame(
  sp = as.factor(rep(c(1:4), each = 1000)),
  x = c(y1, y2, y3, y4)
)


ggplot(df, aes(x, color = sp, fill = sp)) +
  geom_density(alpha = 0.4) +
  xlim(-3, 12) +
  theme_void() +
  scale_color_manual(values = c("#8d18e6ff", "#ca005bff", "#69d1c5ff", "#e0ca3cff")) +
  scale_fill_manual(values = c("#8d18e6ff", "#ca005bff", "#69d1c5ff", "#e0ca3cff")) +
  theme(legend.position = "none")

ggsave("~/Documents/Projects/network_coherence/c_outputs/figures/svg/response_curves_conceptual.png", width = 12, height = 5)
