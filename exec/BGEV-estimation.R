library(dplyr)
library(INLA)
library(parallel)
library(sf)
library(inlaBGEV)
library(ggplot2)
library(patchwork)

# In this script we estimate return level maps for 1, 3, and 6 hour precipitation,
# using the two-step model

hour_vec = c(1, 3, 6, 12, 24) # Which aggregation lengths are we examining?
α = .5; β = .8 # Probabilities used in the location and spread parameters
min_sd_years = 4L # Minimum number of years before we use the computed SD values
return_level_period = 20 # Period we are computing return levels for
n_sd_samples = 1 # Number of samples drawn from the distribution of the SD
num_cores = 1 # Number of cores used for parallel computations

# A list containing covariate_names for location, spread and tail parameter
covariate_names = list(c("precipitation", "height", "x", "y", "dist_sea"),
                       c("x", "y", "dist_sea"), NULL)

stats = list()
for (i in seq_along(hour_vec)) {
  n_hours = hour_vec[i]

  # Filter out the data of interest and standardise the covariates in the observations data
  # and in the prediction data
  data = dplyr::left_join(observations, estimated_sd, by = c("id", "n_hours", "n_years")) %>%
    dplyr::filter(n_hours == !!n_hours)
  standardisation_stats = get_stats_for_standardisation(data, unique(unlist(covariate_names)))
  data = standardise(data, standardisation_stats)
  prediction_data = standardise(prediction_grid, standardisation_stats)

  # Create the mesh used for modelling, and define the prior for the
  # spatial Gaussian field
  mesh = inla_mesh(data)
  spde = inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(75, .05),
    prior.sigma = c(.5, .05))

  # Data used for modelling the SD at all observation locations
  sd_df = data %>%
    dplyr::distinct(id, .keep_all = TRUE) %>%
    dplyr::filter(n_years >= min_sd_years) %>%
    dplyr::mutate(log_sd = log(sd))

  # Create a prior for the spatial Gaussian field used for modelling the SD
  sd_spde = inla.spde2.pcmatern(
    mesh = mesh,
    prior.sigma = c(.5, .05),
    prior.range = c(75, .05))

  # Estimate s^*
  sd_stack = inla_stack(sd_df, covariate_names[[2]], response_name = "log_sd", spde = sd_spde)
  sd_inla_formula = as.formula(paste("log_sd ~ -1 + intercept +",
                                     paste(covariate_names[[2]], collapse = " + "),
                                     "+ f(matern_field, model = sd_spde)"))
  sd_inla_args = inla_default_args("gaussian")
  sd_inla_args$formula = sd_inla_formula
  sd_inla_args$control.predictor$A = INLA::inla.stack.A(sd_stack)
  sd_inla_args$data = INLA::inla.stack.data(sd_stack)
  sd_inla_args$data$sd_spde = sd_spde
  sd_res = do.call(inla, sd_inla_args)

  # Run the two-step model
  samples = twostep_modelling(
    data = data,
    sd_model = sd_res,
    covariate_names = covariate_names,
    response_name = "value",
    n_sd_samples = n_sd_samples,
    prediction_data = prediction_data,
    spde = spde,
    num_cores = num_cores,
    α = α,
    β = β)

  # Compute parameter stats and return level stats at all locations
  stats[[i]] = inla_stats(
    sample_list = lapply(samples, `[[`, "samples"),
    data = prediction_data,
    covariate_names = list(covariate_names[[1]], NULL, NULL),
    s_list = lapply(samples, `[[`, "s_est"),
    mesh = mesh,
    family = "bgev",
    n_batches = 20,
    fun = function(pars) {
      locscale_pars = locspread_to_locscale(pars$q, pars$s, pars$ξ, α, β)
      return_level_gev(return_level_period, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ)
    })

  message("Done with ", n_hours, " hours")
  message("Number of succesful runs: ", length(samples), " of ", n_sd_samples)
}

saveRDS(stats, file.path(here::here(), "inst", "extdata", "return-level-stats.rds"))

# Plot the results ===================================================

# Return levels =======================
my_breaks = c(6, 10, 14, 18, 22, 26)
p1 = stats[[1]]$fun %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(breaks = my_breaks, use_tex = TRUE, size = .3)
p1[[1]] = p1[[1]] + labs(title = "1 hour precipitation")

my_breaks = seq(10, by = 6, length = 6)
p2 = stats[[2]]$f %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(breaks = my_breaks, use_tex = TRUE, size = .3)
p2[[1]] = p2[[1]] + labs(title = "3 hour precipitation")

my_breaks = seq(16, by = 10, length = 6)
p3 = stats[[3]]$f %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(breaks = my_breaks, use_tex = TRUE, size = .3)
p3[[1]] = p3[[1]] + labs(title = "6 hour precipitation")

my_breaks = seq(20, by = 15, length = 6)
p4 = stats[[4]]$f %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(breaks = my_breaks, use_tex = TRUE, size = .3)
p4[[1]] = p4[[1]] + labs(title = "12 hour precipitatio")

my_breaks = seq(30, by = 20, length = 6)
p5 = stats[[5]]$f %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(breaks = my_breaks, use_tex = TRUE, size = .3)
p5[[1]] = p5[[1]] + labs(title = "24 hour precipitatio")


text_size = 8
myplot = patchwork::wrap_plots(p1, p2, p3, nrow = 3) *
  theme(text = element_text(size = text_size))

tikz_plot(file.path(here::here(), "inst", "extdata", "return-level-maps.pdf"),
          myplot, width = 7, height = 10, view = TRUE)

# BGEV parameters ===========================
pq = stats[[1]]$q %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(use_tex = TRUE, size = .3, breaks = c(8, 10, 12, 14, 16))
pq[[1]] = pq[[1]] + labs(title = "1 hour precipitation", fill = "Posterior mean $q_\\alpha$")

ps = stats[[1]]$s %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(use_tex = TRUE, size = .3, breaks = c(1, 1.5, 2, 2.5, 3))
ps[[1]] = ps[[1]] + labs(fill = "Posterior mean $s_\\beta$")
#ps[[1]] = ps[[1]] + labs(title = "1 hour precipitation, $s_\\beta$")

p1 = pq[[1]] + ps[[1]]

pq = stats[[2]]$q %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(use_tex = TRUE, size = .3, breaks = c(15, 20, 25, 30, 35))
pq[[1]] = pq[[1]] + labs(title = "3 hour precipitation", fill = "Posterior mean $q_\\alpha$")

ps = stats[[2]]$s %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(use_tex = TRUE, size = .3, breaks = c(2, 2.5, 3, 3.5, 4))
ps[[1]] = ps[[1]] + labs(fill = "Posterior mean $s_\\beta$")

p2 = pq[[1]] + ps[[1]]

pq = stats[[3]]$q %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(use_tex = TRUE, size = .3, breaks = c(24, 30, 36, 42, 48))
pq[[1]] = pq[[1]] + labs(title = "6 hour precipitation", fill = "Posterior mean $q_\\alpha$")

ps = stats[[3]]$s %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(use_tex = TRUE, size = .3, breaks = c(3, 3.5, 4, 4.5, 5))
ps[[1]] = ps[[1]] + labs(fill = "Posterior mean $s_\\beta$")

p3 = pq[[1]] + ps[[1]]

text_size = 8
myplot = patchwork::wrap_plots(p1, p2, p3, nrow = 3) *
  theme(text = element_text(size = text_size))

tikz_plot(file.path(here::here(), "inst", "extdata", "BGEV-parameter-maps.pdf"),
          myplot, width = 7, height = 10, view = TRUE)


# Tail parameter summary =========================

for (i in seq_along(hour_vec)) {
  message(paste(hour_vec[i], "hour ξ:"))
  print(stats[[i]]$ξ)
}
