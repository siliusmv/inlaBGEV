library(dplyr)
library(INLA)
library(parallel)
library(sf)
library(inlaBGEV)
library(ggplot2)
library(patchwork)

# In this script we estimate return level maps for 1, 3, and 6 hour precipitation,
# using the two-step procedure of sparately modelling the spread in space


hour_vec = c(1, 3, 6)
α = .5
β = .8
min_sd_years = 6L # Minimum number of years before we use the computed SD values
return_level_period = 20 # Period we are computing return levels for
n_sd_samples = 20 # Number of samples drawn from the distribution of the SD
num_cores = 6 # Number of cores used for parallel computations

# A list containing covariate_names for location, spread and tail parameter
covariate_names = list(c("precipitation", "height", "x", "y", "dist_sea", "wetdays"),
                       c("x", "y", "height", "dist_sea"), NULL)

stats = list()
for (i in seq_along(hour_vec)) {
  n_hours = hour_vec[i]

  # Filter out the data of interest and standardise the covariates in the observations data
  # and in the prediction data
  data = dplyr::left_join(observations, estimated_sd, by = c("id", "n_hours", "n_years")) %>%
    dplyr::filter(n_hours == !!n_hours)
  standardisation_stats = get_stats_for_standardisation(data, covariate_names[[1]])
  data = standardise(data, standardisation_stats)
  prediction_data = standardise(prediction_grid, standardisation_stats)

  # Create the mesh used for modelling, and define the prior for the
  # spatial Gaussian field
  mesh = inla_mesh(data)
  spde = inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(75, .05),
    prior.sigma = c(.5, .05),
    constr = TRUE)

  # Data used for modelling the SD at all observation locations
  sd_df = data %>%
    dplyr::distinct(id, .keep_all = TRUE) %>%
    dplyr::filter(n_years >= min_sd_years) %>%
    dplyr::mutate(log_sd = log(sd))

  # Create a prior for the spatial Gaussian field used for modelling the SD
  sd_spde = inla.spde2.pcmatern(
    mesh = mesh,
    prior.sigma = c(1, .05),
    prior.range = c(75, .05),
    constr = TRUE)

  # Prepare to use R-INLA for modelling the SD
  sd_stack = inla_stack(sd_df, covariate_names[[2]], response_name = "log_sd", spde = sd_spde)
  sd_inla_formula = as.formula(paste("log_sd ~ -1 + intercept +",
                                     paste(covariate_names[[2]], collapse = " + "),
                                     "+ f(matern_field, model = sd_spde)"))
  sd_inla_args = inla_default_args("gaussian")
  sd_inla_args$formula = sd_inla_formula
  sd_inla_args$control.predictor$A = INLA::inla.stack.A(sd_stack)
  sd_inla_args$data = INLA::inla.stack.data(sd_stack)
  sd_inla_args$data$sd_spde = sd_spde

  # Run R-INLA
  s_res = do.call(inla, sd_inla_args)

  # Sample from the distribution of the SD at all observation locations and prediction locations
  sd_prediction_data = sd_df %>%
    st_transform(st_crs(prediction_data)) %>%
    dplyr::select(all_of(names(prediction_data))) %>%
    rbind(prediction_data)
  set.seed(1)
  log_sd_samples = inla.posterior.sample(s_res, n = n_sd_samples, seed = 1)
  log_sd_pars = inla_gaussian_pars(
    samples = log_sd_samples,
    data = sd_prediction_data,
    covariate_names = covariate_names[[2]],
    mesh = mesh,
    coords = st_geometry(sd_prediction_data))
  sd_samples = rnorm(length(log_sd_pars$μ), log_sd_pars$μ, 1 / sqrt(log_sd_pars$τ)) %>%
    matrix(nrow = nrow(log_sd_pars$μ)) %>%
    exp()

  # Run R-inla to estimate the BGEV-parameters once for each of the SD samples
  samples = parallel::mclapply(
    X = seq_len(n_sd_samples),
    mc.cores = num_cores,
    mc.preschedule = FALSE,
    FUN = function(i) {
      res = tryCatch({
        inla_bgev(
          data = data,
          s_est = sd_samples[, i][1:nrow(sd_df)],
          covariate_names = list(covariate_names[[1]], NULL, NULL),
          response_name = "value",
          spde = spde,
          α = α,
          β = β)},
        error = function(e) NULL)
      message("Done with iter nr. ", i)
      if (is.null(res)) return(NULL)
      set.seed(1)
      samples = inla.posterior.sample(100, res, seed = 1)
      list(const = res$standardising_const, samples = samples)
    })

  # Sometimes, INLA might have some numerical problems. Remove the bad
  # models
  bad_samples = which(sapply(samples, is.null))
  if (any(bad_samples)) {
    samples = samples[-bad_samples]
    sd_samples = sd_samples[, -bad_samples]
  }

  # Compute the SD multiplied with the standardising const at all prediction locations
  s_est = lapply(
    seq_along(samples),
    function(i) samples[[i]]$const * sd_samples[-(1:nrow(sd_df)), i])

  # Compute parameter stats and return level stats at all locations
  stats[[i]] = inla_bgev_stats(
    sample_list = lapply(samples, `[[`, "samples"),
    data = prediction_data,
    covariate_names = list(covariate_names[[1]], NULL, NULL),
    s_list = s_est,
    mesh = mesh,
    n_batches = 50,
    fun = function(pars) {
      locscale_pars = locspread_to_locscale(pars$q, pars$s, pars$ξ, α, β)
      return_level_gev(return_level_period, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ)
    })

  message("Done with ", n_hours, " hours")
  message("Number of succesful runs: ", length(samples), " of ", n_sd_samples)
}

#saveRDS(stats, file.path(here::here(), "inst", "extdata", "return-level-stats.rds"))

# Plot the results ===================================================

my_breaks = c(8, 12, 16, 20, 24, 28)
p1 = stats[[1]]$fun %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(breaks = my_breaks, CI_breaks = my_breaks, use_tex = TRUE, size = .3)
p1[[1]] = p1[[1]] + labs(title = "1 hour precipitation\n20 year return level")

my_breaks = c(18, 24, 30, 36, 42, 48)
p2 = stats[[2]]$f %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(breaks = my_breaks, CI_breaks = my_breaks, use_tex = TRUE, size = .3)
p2[[1]] = p2[[1]] + labs(title = "3 hour precipitation\n20 year return level")

my_breaks = c(20, 30, 40, 50, 60, 70)
p3 = stats[[3]]$f %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(breaks = my_breaks, CI_breaks = my_breaks, use_tex = TRUE, size = .3)
p3[[1]] = p3[[1]] + labs(title = "6 hour precipitation\n20 year return level")

text_size = 8
myplot = patchwork::wrap_plots(p1, p2, p3, nrow = 3) *
  theme(text = element_text(size = text_size))

#tikz_plot(file.path(here::here(), "inst", "extdata", "1-3-6-hour-return-level-maps.pdf"),
#          print(myplot), width = 7, height = 10)
