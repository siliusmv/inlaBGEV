library(sf)
library(dplyr)
library(inlaBGEV)
library(INLA)

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
  data = dplyr::left_join(observations, estimated_sd, by = c("id", "n_hours", "n_years")) %>%
    dplyr::filter(n_hours == !!n_hours)
  standardisation_stats = get_stats_for_standardisation(data, covariate_names[[1]])
  data = standardise(data, standardisation_stats)

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

  # Sample from the distribution of the SD at all observation locations
  set.seed(1)
  log_sd_samples = inla.posterior.sample(s_res, n = n_sd_samples, seed = 1)
  log_sd_pars = inla_gaussian_pars(
    samples = log_sd_samples,
    data = dplyr::distinct(data, id, .keep_all = TRUE),
    covariate_names = covariate_names[[2]],
    mesh = mesh,
    coords = st_geometry(dplyr::distinct(data, id, .keep_all = TRUE)))
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
          s_est = sd_samples[, i],
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
    s_samples = s_samples[-bad_samples]
  }

  # Compute the SD multiplied with the standardising const at all observation locations
  s_est = lapply(
    seq_along(samples),
    function(i) samples[[i]]$const * sd_samples[, i])

  # Compute sampled parameters at all observation locations
  est_pars = list()
  for (j in seq_along(samples)) {
    est_pars[[j]] = inla_bgev_pars(
      samples = samples[[j]]$samples,
      data = dplyr::distinct(data, id, .keep_all = TRUE),
      covariate_names = list(covariate_names[[1]], NULL, NULL),
      s_est = s_est[[j]],
      mesh = mesh,
      coords = st_geometry(dplyr::distinct(data, id, .keep_all = TRUE)))
  }
  est_pars = purrr::transpose(est_pars)
  for (j in seq_along(est_pars)) {
    if (is.null(dim(est_pars[[j]][[1]]))) {
      est_pars[[j]] = do.call(c, est_pars[[j]])
    } else {
      est_pars[[j]] = do.call(cbind, est_pars[[j]])
    }
  }

  stats[[i]] = list()

  twcrps = list()
  for (j in seq_along(unique(data$id))) {
    id = unique(data$id)[j]
    obs = dplyr::filter(data, id == !!id)$value
    locscale_pars = locspread_to_locscale(est_pars$q[j, ], est_pars$s[j, ], est_pars$ξ[j, ], α, β)
    twcrps[[j]] = twcrps_gev(obs, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ, .9)
  }
  stats[[i]]$two_step = data_stats(as.numeric(unlist(twcrps)))

  res = inla_bgev(
    data = data,
    covariate_names = covariate_names,
    response_name = "value",
    spde = spde,
    α = α,
    β = β)

  samples = inla.posterior.sample(100 * length(samples), res)

  est_pars = inla_bgev_pars(
    samples = samples,
    data = dplyr::distinct(data, id, .keep_all = TRUE),
    covariate_names = covariate_names,
    s_est = rep(res$standardising_const, length(unique(data$id))),
    mesh = mesh,
    coords = st_geometry(dplyr::distinct(data, id, .keep_all = TRUE)))

  twcrps = list()
  for (j in seq_along(unique(data$id))) {
    id = unique(data$id)[j]
    obs = dplyr::filter(data, id == !!id)$value
    locscale_pars = locspread_to_locscale(est_pars$q[j, ], est_pars$s[j, ], est_pars$ξ[j, ], α, β)
    twcrps[[j]] = twcrps_gev(obs, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ, .9)
  }
  stats[[i]]$joint_modelling = data_stats(as.numeric(unlist(twcrps)))

  message("Done with ", n_hours, " hours")
  message("Number of succesful runs: ", length(samples), " of ", mynum)
}

saveRDS(stats, file.path(here::here(), "inst", "extdata", "twcrps-comparison.rds"))

# CRPS
for (i in seq_along(all_stats)) {
  message("=========================================\n",
          hour_vec[i], " hour(s)\n",
          "=========================================")
  message("With standardisation:")
  print(all_stats[[i]]$standardised_crps)
  message("With direct modelling:")
  print(all_stats[[i]]$direct_crps)
}
# This is great. The standardisation approach is better above the 90% quantile!

# Tail parameter estimation
for (i in seq_along(all_stats)) {
  message("=========================================\n",
          hour_vec[i], " hour(s)\n",
          "=========================================")
  message("With standardisation:")
  print(summary(all_stats[[i]]$standardised_tail))
  message("With direct modelling:")
  print(summary(all_stats[[i]]$direct_tail))
}
# The directly estimated ξ seems to be too small
