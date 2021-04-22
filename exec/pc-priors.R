library(inlaBGEV)
library(dplyr)
library(tidyr)
library(ggplot2)

# In this script we compare the different Kullback-Leibler distances that have been discussed.
# We also use these to compare the PC priors for different extreme value distributions

# Compute KLD for different distributions
df = data.frame(ξ = seq(.001, .999, by = .001))
df$gp = kld_gp(df$ξ)
df$gev = kld_gev(df$ξ)
df$bgev = kld_bgev(df$ξ)

# Plot KLD
plot = tidyr::pivot_longer(df, starts_with(c("g", "b"))) %>%
  dplyr::filter(ξ < .4) %>%
  ggplot() +
  geom_line(aes(x = ξ, y = value, col = name, group = name))
if (interactive()) print(plot)
# There is almost no difference between the approximated KLD and the numerically computed KLD

# Compute PC prior densities
λ = 1
df2 = data.frame(ξ = seq(.001, .999, by = .001))
df2$gev = pc_gev(df2$ξ, λ = λ)
df2$gp = pc_gp(df2$ξ, λ = λ)
df2$bgev = pc_bgev(df2$ξ, λ = λ)

# Plot the densities
plot = tidyr::pivot_longer(df2, starts_with(c("g", "b"))) %>%
  ggplot() +
  geom_line(aes(x = ξ, y = value, col = name, group = name))
if (interactive()) print(plot)
# With λ = 7, all the PC priors are practically identical. However, there are
# some considerable differences between the PC prior and its approximation
# for the GEV distribution

# Create a plot for the article
plot = df2 %>%
  pivot_longer(starts_with(c("g", "b"))) %>%
  dplyr::mutate(name = factor(name, levels = c("gev", "bgev", "gp"),
                              labels = c("GEV", "BGEV", "GP"))) %>%
  ggplot() +
  geom_line(aes(x = ξ, y = value, linetype = name, group = name)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  theme_bw() +
  labs(x = "$\\xi$", y = "Density", linetype = "Distribution")

tikz_plot(file.path(here::here(), "results", "pc-priors.pdf"),
          plot, width = 10, height = 7, view = TRUE)


df3 = list()
λ_vec = c(.1, .5, 1, 2.5, 5, 10)
for (i in seq_along(λ_vec)) {
  df3[[i]] = data.frame(ξ = seq(.005, .999, by = .001))
  df3[[i]]$gev = pc_gev(df3[[i]]$ξ, λ = λ_vec[i])
  df3[[i]]$gp = pc_gp(df3[[i]]$ξ, λ = λ_vec[i])
  df3[[i]]$bgev = pc_bgev(df3[[i]]$ξ, λ = λ_vec[i])
  df3[[i]]$λ = λ_vec[i]
  message(i)
}
df3 = do.call(rbind, df3)

# Create a plot for the article
plot = df3 %>%
  pivot_longer(starts_with(c("g", "b"))) %>%
  dplyr::mutate(name = factor(name, levels = c("gev", "bgev", "gp"),
                              labels = c("GEV", "BGEV", "GP")),
                λ = factor(λ, levels = λ_vec, labels = paste("$\\lambda =", λ_vec, "$"))) %>%
  ggplot() +
  geom_line(aes(x = ξ, y = value, linetype = name, group = name)) +
  scale_linetype_manual(values = c("solid", "dashed", "dotted")) +
  theme_bw() +
  labs(x = "$\\xi$", y = "Density", linetype = "Distribution") +
  facet_wrap(~λ, scales = "free_y") +
  theme(text = element_text(size = 14))

tikz_plot(file.path(here::here(), "results", "pc-priors.pdf"),
          plot, width = 10, height = 7, view = TRUE)
