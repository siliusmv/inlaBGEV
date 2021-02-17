library(inlaBGEV)
library(dplyr)
library(tidyr)
library(ggplot2)

# In this script we compare the different Kullback-Leibler distances that have been discussed.
# We also use these to compare the PC priors for different extreme value distributions

# Compute KLD for different distributions
df = data.frame(ξ = seq(.001, .999, by = .001))
df$gp = kld_gp(df$ξ)
df$gev = kld_gev(df$ξ, approx = FALSE)
df$gev_approx = kld_gev(df$ξ, approx = TRUE)
df$bgev = kld_bgev(df$ξ, approx = FALSE)
df$bgev_approx = kld_bgev(df$ξ, approx = TRUE)

# Plot KLD
tidyr::pivot_longer(df, starts_with(c("g", "b"))) %>%
  dplyr::filter(ξ < .4) %>%
  ggplot() +
  geom_line(aes(x = ξ, y = value, col = name, group = name))
# There is almost no difference between the approximated KLD and the numerically computed KLD

# Compute PC prior densities
df2 = data.frame(ξ = seq(.001, .999, by = .001))
df2$gev = pc_gev(df2$ξ, λ = 7, approx = FALSE)
df2$gev_approx = pc_gev(df2$ξ, λ = 7, approx = TRUE)
df2$gp = pc_gp(df2$ξ, λ = 7)
df2$bgev = pc_bgev(df2$ξ, λ = 7, approx = FALSE)
df2$bgev_approx = pc_bgev(df2$ξ, λ = 7, approx = TRUE)

# Plot the densities
tidyr::pivot_longer(df2, starts_with(c("g", "b"))) %>%
  ggplot() +
  geom_line(aes(x = ξ, y = value, col = name, group = name))
# With λ = 7, all the PC priors are practically identical. However, there are
# some considerable differences between the PC prior and its approximation
# for the GEV distribution

plot = df2 %>%
  dplyr::select(-gev_approx, -bgev_approx) %>%
  pivot_longer(starts_with(c("g", "b"))) %>%
  dplyr::mutate(name = factor(name, levels = c("gev", "bgev", "gp"),
                              labels = c("GEV", "BGEV", "GP"))) %>%
  ggplot() +
  geom_line(aes(x = ξ, y = value, linetype = name, group = name)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  theme_bw() +
  labs(x = "$\\xi$", y = "Density", linetype = "Distribution")
plot

#tikz_plot(file.path(here::here(), "inst", "extdata", "pc-priors.pdf"),
#          print(plot), width = 10, height = 7, view = TRUE)
