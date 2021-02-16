library(dplyr)
library(ggplot2)
library(tidyr)
library(inlaBGEV)
library(sf)

# In this script we estimate GEV parameters locally at all weather stations with
# long enough time series. We then visualise the estimates to see if a spatially
# varying model is needed for parameter estimation. Additionally we compare the
# parameter estimates with explanatory variables to look for trends.
# Finally, we standardise the observations using the standard deviation of
# heavy precipitation at each location. We then estimate parameters and look
# for trends once more.


# functions ==================================================================
bootstrap_gev_fit = function(x, n, q = c(.025, .25, .5, .75, .975), ...) {
  pars = matrix(NA, nrow = n, ncol = length(x$estimate))
  for (b in 1:n) {
    y = sample(x$data, length(x$data), replace = TRUE)
    pars[b, ] = gev_fit(y, ...)$estimate
  }
  stats = data_stats(t(pars), q)
  stats$par = names(x$estimate)
  stats
}

# Estimate parameters for different time periods with PWM ============================
min_years = 4
α = .5
β = .8
hour_vec = c(1, 3, 6)
covariate_names = c("x", "y", "height", "dist_sea", "precipitation", "wetdays")

params = list()
for (i in seq_along(hour_vec)) {
  # Prepare the data
  data = dplyr::filter(observations, n_hours == hour_vec[i], n_years >= min_years)
  standardisation_stats = get_stats_for_standardisation(data, covariate_names)
  data = standardise(data, standardisation_stats)

  # Estimate parameters
  ids = unique(data$id)
  fits = list()
  set.seed(1)
  for (j in seq_along(ids)) {
      obs = dplyr::filter(data, id == ids[j])$value
      fits[[j]] = gev_fit(obs, method = "pwm", α = α, β = β)
  }

  # Estimate confidence intervals using bootstrap
  bootstraps = list()
  set.seed(1)
  for (j in seq_along(fits)) {
    bootstraps[[j]] = bootstrap_gev_fit(fits[[j]], n = 1000, method = "pwm", α = α, β = β)
    bootstraps[[j]]$estimator = fits[[j]]$estimate
    bootstraps[[j]]$id = ids[j]
  }
  bootstraps = do.call(rbind, bootstraps)

  # Add explanatory variables to the estimated variables
  explanatory_variables = data %>%
    dplyr::select(-value, -year) %>%
    dplyr::distinct(id, .keep_all = TRUE)
  bootstraps = dplyr::left_join(bootstraps, explanatory_variables, by = "id")

  params[[i]] = bootstraps
  message("Done with ", hour_vec[i], " hours")
}
params = do.call(rbind, params)

# Plot parameter estimates for all stations ====================
plot = params %>%
  dplyr::filter(!is.na(mean)) %>%
  # Remove one stations for purely aesthetic reasons. The station with id SN49085
  # has four observations, and for n_hours = 6, one of them is really large.
  # This causes a very high estimate and a really wide confidence interval for the spread,
  # making it harder to see the patterns in the parameter estimates from other stations.
  # You are welcome to uncomment the line below if you wish
  dplyr::filter(id != "SN49085") %>%
  dplyr::mutate(
    par = factor(
      par, levels = c("q", "s", "ξ"),
      labels = c("$\\hat{q}_\\alpha$", "$\\hat{s}_\\beta$", "$\\hat{\\xi}$")),
    n_hours = factor(n_hours, levels = hour_vec,
                     labels = paste(hour_vec, "hours"))) %>%
  ggplot(aes(x = as.numeric(factor(id)))) +
  geom_point(aes(y = estimator), size = .2) +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = .3) +
  facet_grid(par ~ n_hours, scales = "free_y") +
  labs(x = "Location nr.", y = "Value") +
  theme_bw() +
  theme(text = element_text(size = 20),
        strip.text.y = element_text(angle = 0))
plot

# tikz_plot(file.path(here::here(), "inst", "extdata", "pwm_gev_parameters.pdf"),
#           print(plot),
#           width = 10, height = 7, view = TRUE)

# Plot parameter estimates with covariates ==============================
q_plot = params %>%
  dplyr::filter(par %in% c("q"), !is.na(mean)) %>%
  dplyr::mutate(n_hours = factor(n_hours, levels = hour_vec,
                                 labels = paste(hour_vec, "hours"))) %>%
  tidyr::pivot_longer(all_of(covariate_names)) %>%
  dplyr::mutate(name = factor(
    name, levels = covariate_names,
    labels = c("Mean daily\nprecipitation", "Mean yearly\nwetdays", "Distance\nfrom sea",
               "Altitude", "Easting", "Norting"))) %>%
  ggplot() +
  geom_point(aes(y = estimator, x = value), size = .2) +
  geom_smooth(aes(y = estimator, x = value), formula = y ~ x, method = "lm",
              col = "black", size = .5) +
  facet_grid(n_hours ~ name, scales = "free_y") +
  labs(y = "$\\hat{q}_\\alpha$", x = "Standardised value") +
  theme_bw() +
  theme(text = element_text(size = 20),
        strip.text.y = element_text(angle = 0),
        axis.title.y = element_text(angle = 0, vjust = .5))
q_plot

s_plot = params %>%
  dplyr::filter(par %in% c("s"), !is.na(mean)) %>%
  dplyr::mutate(n_hours = factor(n_hours, levels = hour_vec,
                                 labels = paste(hour_vec, "hours"))) %>%
  tidyr::pivot_longer(all_of(covariate_names)) %>%
  dplyr::mutate(name = factor(
    name, levels = covariate_names,
    labels = c("Mean daily\nprecipitation", "Mean yearly\nwetdays", "Distance\nfrom sea",
               "Altitude", "Easting", "Norting"))) %>%
  ggplot() +
  geom_point(aes(y = estimator, x = value), size = .2) +
  geom_smooth(aes(y = estimator, x = value), formula = y ~ x, method = "lm",
              col = "black", size = .5) +
  labs(y = "$\\hat{s}_\\beta$", x = "Standardised value") +
  facet_grid(n_hours ~ name, scales = "free_y") +
  theme_bw() +
  theme(text = element_text(size = 20),
        strip.text.y = element_text(angle = 0),
        axis.title.y = element_text(angle = 0, vjust = .5))
s_plot

# tikz_plot(file.path(here::here(), "inst", "extdata", "pwm_gev_parameters_with_covariates.pdf"),
#           expression = {print(q_plot); print(s_plot)},
#           width = 10, height = 7, view = TRUE)





# Standardise the response using the standard deviation of heavy precipitation =========
min_years = 4
α = .5
β = .8
hour_vec = c(1, 3, 6)
covariate_names = c("x", "y", "height", "dist_sea", "precipitation", "wetdays")

params = list()
for (i in seq_along(hour_vec)) {
  # Prepare the data
  data = dplyr::filter(observations, n_hours == hour_vec[i], n_years >= min_years)
  standardisation_stats = get_stats_for_standardisation(data, covariate_names)
  data = standardise(data, standardisation_stats)
  sd_vals = dplyr::filter(estimated_sd, n_hours == hour_vec[i], n_years >= min_years)
  data = dplyr::left_join(data, sd_vals, by = c("id", "n_hours", "n_years"))
  data$value = data$value / data$sd

  # Estimate parameters
  ids = unique(data$id)
  fits = list()
  set.seed(1)
  for (j in seq_along(ids)) {
      obs = dplyr::filter(data, id == ids[j])$value
      fits[[j]] = gev_fit(obs, method = "pwm", α = α, β = β)
  }

  # Estimate confidence intervals using bootstrap
  bootstraps = list()
  set.seed(1)
  for (j in seq_along(fits)) {
    bootstraps[[j]] = bootstrap_gev_fit(fits[[j]], n = 1000, method = "pwm", α = α, β = β)
    bootstraps[[j]]$estimator = fits[[j]]$estimate
    bootstraps[[j]]$id = ids[j]
  }
  bootstraps = do.call(rbind, bootstraps)

  # Add explanatory variables to the estimated variables
  explanatory_variables = data %>%
    dplyr::select(-value, -year) %>%
    dplyr::distinct(id, .keep_all = TRUE)
  bootstraps = dplyr::left_join(bootstraps, explanatory_variables, by = "id")

  params[[i]] = bootstraps
  message("Done with ", hour_vec[i], " hours")
}
params = do.call(rbind, params)

# Plot parameter estimates for all stations ====================

plot = params %>%
  dplyr::filter(!is.na(mean)) %>%
  # Remove one stations for purely aesthetic reasons. The station with id SN49085
  # has four observations, and for n_hours = 6, one of them is really large.
  # This causes a very high estimate and a really wide confidence interval for the spread,
  # making it harder to see the patterns in the parameter estimates from other stations.
  # You are welcome to uncomment the line below if you wish
  #dplyr::filter(id != "SN49085") %>%
  dplyr::mutate(
    par = factor(
      par, levels = c("q", "s", "ξ"),
      labels = c("$\\hat{q}_\\alpha$", "$\\hat{s}_\\beta$", "$\\hat{\\xi}$")),
    n_hours = factor(n_hours, levels = hour_vec,
                     labels = paste(hour_vec, "hours"))) %>%
  ggplot(aes(x = as.numeric(factor(id)))) +
  geom_point(aes(y = estimator), size = .2) +
  geom_linerange(aes(ymin = `2.5%`, ymax = `97.5%`), alpha = .3) +
  facet_grid(par ~ n_hours, scales = "free_y") +
  labs(x = "Location nr.", y = "Value") +
  theme_bw() +
  theme(text = element_text(size = 20),
        strip.text.y = element_text(angle = 0))
plot

#tikz_plot(file.path(here::here(), "inst", "extdata", "pwm_gev_parameters_twostep.pdf"),
#          print(plot),
#          width = 10, height = 7, view = TRUE)

# Plot parameter estimates with covariates ==============================
q_plot = params %>%
  dplyr::filter(par %in% c("q"), !is.na(mean)) %>%
  dplyr::mutate(n_hours = factor(n_hours, levels = hour_vec,
                                 labels = paste(hour_vec, "hours"))) %>%
  tidyr::pivot_longer(all_of(covariate_names)) %>%
  dplyr::mutate(name = factor(
    name, levels = covariate_names,
    labels = c("Mean daily\nprecipitation", "Mean yearly\nwetdays", "Distance\nfrom sea",
               "Altitude", "Easting", "Norting"))) %>%
  ggplot() +
  geom_point(aes(y = estimator, x = value), size = .2) +
  geom_smooth(aes(y = estimator, x = value), formula = y ~ x, method = "lm",
              col = "black", size = .5) +
  facet_grid(n_hours ~ name, scales = "free_y") +
  labs(y = "$\\hat{q}_\\alpha$", x = "Standardised value") +
  theme_bw() +
  theme(text = element_text(size = 20),
        strip.text.y = element_text(angle = 0),
        axis.title.y = element_text(angle = 0, vjust = .5))
q_plot

s_plot = params %>%
  dplyr::filter(par %in% c("s"), !is.na(mean)) %>%
  dplyr::mutate(n_hours = factor(n_hours, levels = hour_vec,
                                 labels = paste(hour_vec, "hours"))) %>%
  tidyr::pivot_longer(all_of(covariate_names)) %>%
  dplyr::mutate(name = factor(
    name, levels = covariate_names,
    labels = c("Mean daily\nprecipitation", "Mean yearly\nwetdays", "Distance\nfrom sea",
               "Altitude", "Easting", "Norting"))) %>%
  ggplot() +
  geom_point(aes(y = estimator, x = value), size = .2) +
  geom_smooth(aes(y = estimator, x = value), formula = y ~ x, method = "lm",
              col = "black", size = .5) +
  labs(y = "$\\hat{s}_\\beta$", x = "Standardised value") +
  facet_grid(n_hours ~ name, scales = "free_y") +
  theme_bw() +
  theme(text = element_text(size = 20),
        strip.text.y = element_text(angle = 0),
        axis.title.y = element_text(angle = 0, vjust = .5))
s_plot

#tikz_plot(file.path(here::here(), "inst", "extdata",
#                    "pwm_gev_parameters_twostep_with_covariates.pdf"),
#          expression = {print(q_plot); print(s_plot)},
#          width = 10, height = 7, view = TRUE)
