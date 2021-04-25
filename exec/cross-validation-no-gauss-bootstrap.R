library(sf)
library(dplyr)
library(inlaBGEV)
library(INLA)

# In this script we perform k-fold cross-validation on the joint model and on the
# two-step model. The performance of the foreacsts are evaluated using StwCRPS, where
# the weight function is an indicator function that is one for larger quantiles than
# the p_0 quantile

hour_vec = c(1, 3, 6, 12, 24) # Which aggregation lengths are we examining?
α = .5; β = .8 # Probabilities used in the location and spread parameters
min_sd_years = 4L # Minimum number of years before we use the computed SD values
n_sd_samples = 100 # Number of samples drawn from the distribution of the SD
num_cores = 20 # Number of cores used for parallel computations
n_folds = 5 # number of folds for cross-validation
p0 = .9 # Threshold used in the stwCRPS

# A list containing covariate_names for location, spread and tail parameter
covariate_names = list(c("precipitation", "height", "x", "y", "dist_sea"),
                       c("x", "y", "dist_sea"), NULL)

stats = list()
for (i in seq_along(hour_vec)) {

  # Create the empty list that will contain all the result
  stats[[i]] = list(
    in_sample = list(),
    out_of_sample = list())
  n_hours = hour_vec[i]

  # Filter out the data of interest and standardise the covariates in the observations data
  data = dplyr::left_join(observations, estimated_sd, by = c("id", "n_hours", "n_years")) %>%
    dplyr::filter(n_hours == !!n_hours)
  standardisation_stats = get_mean_and_sd_stats(data, unique(unlist(covariate_names)))
  data = standardise(data, standardisation_stats)

  # Create the folds for cross-validation
  set.seed(1)
  ids = unique(data$id)
  folds = rep(seq_len(n_folds), length(ids) / n_folds)
  folds = c(folds, sample.int(n_folds, length(ids) - length(folds)))
  folds = sample(folds, length(folds))

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
    prior.sigma = c(1, .05),
    prior.range = c(75, .05))

  # Estimate s^*
  message("Start in-sample evaluation")
  sd_inla_args2 = inla_default_args("gaussian")
  sd_inla_args2$data = sd_df
  sd_inla_args2$data$intercept = 1
  sd_inla_args2$formula = as.formula(paste("log_sd ~ -1 + intercept +",
                                           paste(covariate_names[[2]], collapse = " + ")))
  sd_res2 = do.call(inla, sd_inla_args2)

  # Perform in-sample estimation using all the data, twostep no gaussian
  set.seed(1)
  twostep_nogaussian = twostep_modelling(
    data = data,
    sd_model = sd_res2,
    covariate_names = covariate_names,
    response_name = "value",
    diagonal = .05,
    n_sd_samples = 1,
    spde = spde,
    num_cores = 1,
    α = α,
    β = β)
  message("Done with in-sample two-step model, no gaussian")

  # Estimate parameters at all locations
  params = list()
  for (k in seq_along(twostep_nogaussian)) {
    params[[k]] = inla_bgev_pars(
      samples = twostep_nogaussian[[k]]$samples,
      data = data,
      covariate_names = list(covariate_names[[1]], NULL, NULL),
      s_est = twostep_nogaussian[[k]]$s_est,
      mesh = mesh,
      coords = st_geometry(dplyr::distinct(data, id)))
  }
  params = purrr::transpose(params)
  for (k in seq_along(params)) params[[k]] = do.call(cbind, params[[k]])

  # Compute stwCRPS
  for (k in seq_along(unique(data$id))) {
    id = unique(data$id)[k]
    obs = dplyr::filter(data, id == !!id)$value
    locscale_pars = locspread_to_locscale(params$q[k, ], params$s[k, ], params$ξ[k, ], α, β)
    stats[[i]]$in_sample[[k]] = stwcrps_bgev(
      obs, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ, p0)
  }

  message("Start with the out-of-sample modelling")
  for (j in seq_len(n_folds)) {

    # Split the data up into in-sample and out-of-sample data for the current folds
    in_fold_sd_df = dplyr::filter(sd_df, id %in% ids[folds != j])
    in_fold_data = dplyr::filter(data, id %in% ids[folds != j])
    out_of_fold_data = dplyr::filter(data, id %in% ids[folds == j])

    # Estimate s^*
    sd_inla_args2 = inla_default_args("gaussian")
    sd_inla_args2$data = in_fold_sd_df
    sd_inla_args2$data$intercept = 1
    sd_inla_args2$formula = as.formula(paste("log_sd ~ -1 + intercept +",
                                             paste(covariate_names[[2]], collapse = " + ")))
    sd_res2 = do.call(inla, sd_inla_args2)

    # Perform out-of-sample estimation with the two-step model, with no gaussian
    set.seed(1)
    twostep_nogaussian = twostep_modelling(
      data = in_fold_data,
      prediction_data = dplyr::distinct(out_of_fold_data, id, .keep_all = TRUE),
      sd_model = sd_res2,
      covariate_names = covariate_names,
      response_name = "value",
      diagonal = .05,
      n_sd_samples = 1,
      verbose = FALSE,
      spde = spde,
      num_cores = 1,
      α = α,
      β = β)
    message("Done with out-of-sample two-step model no gaussian for fold ", j)

    # Estimate parameters at all out-of-fold locations
    params = list()
    for (k in seq_along(twostep_nogaussian)) {
      params[[k]] = inla_bgev_pars(
        samples = twostep_nogaussian[[k]]$samples,
        data = dplyr::distinct(out_of_fold_data, id, .keep_all = TRUE),
        covariate_names = list(covariate_names[[1]], NULL, NULL),
        s_est = twostep_nogaussian[[k]]$s_est,
        mesh = mesh,
        coords = st_geometry(dplyr::distinct(out_of_fold_data, id)))
    }
    params = purrr::transpose(params)
    for (k in seq_along(params)) params[[k]] = do.call(cbind, params[[k]])

    # Compute stwCRPS
    stats[[i]]$out_of_sample[[j]] = list()
    for (k in seq_along(unique(out_of_fold_data$id))) {
      id = unique(out_of_fold_data$id)[k]
      obs = dplyr::filter(out_of_fold_data, id == !!id)$value
      locscale_pars = locspread_to_locscale(params$q[k, ], params$s[k, ], params$ξ[k, ], α, β)
      stats[[i]]$out_of_sample[[j]][[k]] = stwcrps_bgev(
        obs, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ, p0)
    }

    message("Done with fold nr. ", j)
  }

  # Print results for this specific aggregation period
  message("=========================================\n",
          hour_vec[i], " hour(s)\n",
          "=========================================")

}

saveRDS(stats, file.path(here::here(), "results", "cross-validation-no-gauss-bootstrap.rds"))
