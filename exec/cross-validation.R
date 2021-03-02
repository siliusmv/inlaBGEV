library(sf)
library(dplyr)
library(inlaBGEV)
library(INLA)

# In this script we perform k-fold cross-validation on the joint model and on the
# two-step model. The performance of the foreacsts are evaluated using twCRPS, where
# the weight function is an indicator function that is one for larger quantiles than
# the p_0 quantile

hour_vec = c(1, 3, 6) # Which aggregation lengths are we examining?
α = .5; β = .8 # Probabilities used in the location and spread parameters
min_sd_years = 4L # Minimum number of years before we use the computed SD values
n_sd_samples = 20 # Number of samples drawn from the distribution of the SD
num_cores = 20 # Number of cores used for parallel computations
n_folds = 5 # number of folds for cross-validation
p0 = .9 # Threshold used in the twCRPS

# A list containing covariate_names for location, spread and tail parameter
covariate_names = list(c("precipitation", "height", "x", "y", "dist_sea"),
                       c("x", "y", "dist_sea"), NULL)
stats = list()
for (i in seq_along(hour_vec)) {
  # Create the empty list that will contain all the result
  stats[[i]] = list(
    in_sample_twostep = list(),
    out_of_sample_twostep = list(),
    in_sample_joint = list(),
    out_of_sample_joint = list())
  n_hours = hour_vec[i]

  # Filter out the data of interest and standardise the covariates in the observations data
  data = dplyr::left_join(observations, estimated_sd, by = c("id", "n_hours", "n_years")) %>%
    dplyr::filter(n_hours == !!n_hours)
  standardisation_stats = get_stats_for_standardisation(data, unique(unlist(covariate_names)))
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

  # Perform in-sample estimation using all the data
  message("Start in-sample evaluation")
  sd_stack = inla_stack(sd_df, covariate_names[[2]],
                        response_name = "log_sd", spde = sd_spde)
  sd_inla_formula = as.formula(paste("log_sd ~ -1 + intercept +",
                                     paste(covariate_names[[2]], collapse = " + "),
                                     "+ f(matern_field, model = sd_spde)"))
  sd_inla_args = inla_default_args("gaussian")
  sd_inla_args$formula = sd_inla_formula
  sd_inla_args$control.predictor$A = INLA::inla.stack.A(sd_stack)
  sd_inla_args$data = INLA::inla.stack.data(sd_stack)
  sd_inla_args$data$sd_spde = sd_spde
  sd_inla_args$control.family = list(hyper = list(prec = list(init = 1e10, fixed = TRUE)))
  s_res = do.call(inla, sd_inla_args)

  # Sample from the posterior of log(s^*) and transform back to s^*
  set.seed(1)
  log_sd_samples = inla.posterior.sample(s_res, n = n_sd_samples, seed = 1)
  log_sd_pars = inla_gaussian_pars(
    samples = log_sd_samples,
    data = dplyr::distinct(data, id, .keep_all = TRUE),
    covariate_names = covariate_names[[2]],
    mesh = mesh,
    coords = st_geometry(dplyr::distinct(data, id, .keep_all = TRUE)))
  in_sample_sd_samples = rnorm(length(log_sd_pars$μ), log_sd_pars$μ, 1 / sqrt(1e10)) %>%
    matrix(nrow = nrow(log_sd_pars$μ)) %>%
    exp()

  # Run R-INLA once for each sample of s^*, and merge the posterior samples
  in_sample_samples = parallel::mclapply(
    X = seq_len(n_sd_samples),
    mc.cores = num_cores,
    mc.preschedule = FALSE,
    FUN = function(i) {
      res = tryCatch({
        inla_bgev(
          data = data,
          s_est = in_sample_sd_samples[, i],
          covariate_names = list(covariate_names[[1]], NULL, NULL),
          response_name = "value",
          spde = spde,
          α = α,
          β = β)},
        error = function(e) NULL)
      if (is.null(res) || !res$convergence) return(NULL)
      set.seed(1)
      samples = inla.posterior.sample(100, res, seed = 1)
      list(const = res$standardising_const, samples = samples)
    })
  message("Done with in-sample two-step model")

  # Sometimes, R-INLA might have some numerical problems. Remove the bad models
  bad_samples = which(sapply(in_sample_samples, is.null))
  if (any(bad_samples)) {
    in_sample_samples = in_sample_samples[-bad_samples]
    in_sample_sd_samples = in_sample_sd_samples[, -bad_samples]
  }

  # Run the joint model on the data
  in_sample_res2 = tryCatch(inla_bgev(
    data = data,
    covariate_names = covariate_names,
    response_name = "value",
    spde = spde,
    α = α,
    β = β), error = function(e) NULL)
  message("Done with in-sample joint model")

  # Sample from the posterior of the joint model
  if (!is.null(in_sample_res2)) {
    in_sample_samples2 = inla.posterior.sample(100 * length(in_sample_samples), in_sample_res2)
  }

  message("Start with the out-of-sample modelling")
  for (j in seq_len(n_folds)) {

    # Split the data up into in-sample and out-of-sample data for the current folds
    in_sample_sd_df = dplyr::filter(sd_df, id %in% ids[folds != j])
    in_sample_data = dplyr::filter(data, id %in% ids[folds != j])
    leave_out_data = dplyr::distinct(data, id, .keep_all = TRUE) %>%
      dplyr::filter(id %in% ids[folds == j])

    # Prepare to use R-INLA for modelling the SD based on data from the in-sample folds
    sd_stack = inla_stack(in_sample_sd_df, covariate_names[[2]],
                          response_name = "log_sd", spde = sd_spde)
    sd_inla_formula = as.formula(paste("log_sd ~ -1 + intercept +",
                                       paste(covariate_names[[2]], collapse = " + "),
                                       "+ f(matern_field, model = sd_spde)"))
    sd_inla_args = inla_default_args("gaussian")
    sd_inla_args$formula = sd_inla_formula
    sd_inla_args$control.predictor$A = INLA::inla.stack.A(sd_stack)
    sd_inla_args$data = INLA::inla.stack.data(sd_stack)
    sd_inla_args$data$sd_spde = sd_spde
    sd_inla_args$control.family = list(hyper = list(prec = list(init = 1e10, fixed = TRUE)))
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
    sd_samples = rnorm(length(log_sd_pars$μ), log_sd_pars$μ, 1 / sqrt(1e10)) %>%
      #log_sd_pars$μ %>%
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
            data = in_sample_data,
            s_est = sd_samples[, i][folds != j],
            covariate_names = list(covariate_names[[1]], NULL, NULL),
            response_name = "value",
            spde = spde,
            α = α,
            β = β)},
          error = function(e) NULL)
        if (is.null(res) || !res$convergence) return(NULL)
        set.seed(1)
        samples = inla.posterior.sample(100, res, seed = 1)
        list(const = res$standardising_const, samples = samples)
      })
    message("Done with out-of-sample two-step model for fold ", j)

    # Sometimes, R-INLA might have some numerical problems. Remove the bad models
    bad_samples = which(sapply(samples, is.null))
    if (any(bad_samples)) {
      samples = samples[-bad_samples]
      sd_samples = sd_samples[, -bad_samples]
    }

    # Compute the SD multiplied with the standardising const at all leave_out observation locations
    s_est = lapply(
      seq_along(samples),
      function(i) samples[[i]]$const * sd_samples[, i][folds == j])

    # Compute sampled parameters at all leave-out locations
    est_pars = list()
    for (k in seq_along(samples)) {
      est_pars[[k]] = inla_bgev_pars(
        samples = samples[[k]]$samples,
        data = leave_out_data,
        covariate_names = list(covariate_names[[1]], NULL, NULL),
        s_est = s_est[[k]],
        mesh = mesh,
        coords = st_geometry(leave_out_data))
    }
    est_pars = purrr::transpose(est_pars)
    for (k in seq_along(est_pars)) {
      if (is.null(dim(est_pars[[k]][[1]]))) {
        est_pars[[k]] = do.call(c, est_pars[[k]])
      } else {
        est_pars[[k]] = do.call(cbind, est_pars[[k]])
      }
    }

    # Compute twCRPS at all leave-out locations
    twcrps = list()
    for (k in seq_along(unique(leave_out_data$id))) {
      id = unique(leave_out_data$id)[k]
      obs = dplyr::filter(data, id == !!id)$value
      locscale_pars = locspread_to_locscale(est_pars$q[k, ], est_pars$s[k, ],
                                            est_pars$ξ[k, ], α, β)
      twcrps[[k]] = twcrps_gev(obs, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ, p0)
    }
    stats[[i]]$out_of_sample_twostep[[j]] = twcrps

    # Use joint modelling, out-of-sample
    res2 = tryCatch(inla_bgev(
      data = in_sample_data,
      covariate_names = covariate_names,
      response_name = "value",
      spde = spde,
      α = α,
      β = β), error = function(e) NULL)
    message("Done with out-of-sample joint model for fold ", j)

    if (!is.null(res2)) {
      samples2 = inla.posterior.sample(100 * length(samples), res2)

      # Compute sampled parameters at all leave-out locations
      est_pars2 = inla_bgev_pars(
        samples = samples2,
        data = leave_out_data,
        covariate_names = covariate_names,
        s_est = rep(res2$standardising_const, nrow(leave_out_data)),
        mesh = mesh,
        coords = st_geometry(leave_out_data))

      # Compute twCRPS at all leave-out locations
      twcrps2 = list()
      for (k in seq_along(unique(leave_out_data$id))) {
        id = unique(data$id)[k]
        obs = dplyr::filter(data, id == !!id)$value
        locscale_pars = locspread_to_locscale(est_pars2$q[k, ], est_pars2$s[k, ],
                                              est_pars2$ξ[k, ], α, β)
        twcrps2[[k]] = twcrps_gev(obs, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ, p0)
      }
      stats[[i]]$out_of_sample_joint[[j]] = twcrps2
    }

    message("Compute in-sample twCRPS for fold ", j)
    # Here we just do all the same things once more, but now for the models that are
    # trained on all the data

    s_est = lapply(
      seq_along(in_sample_samples),
      function(i) in_sample_samples[[i]]$const * in_sample_sd_samples[, i][folds == j])

    est_pars = list()
    for (k in seq_along(in_sample_samples)) {
      est_pars[[k]] = inla_bgev_pars(
        samples = in_sample_samples[[k]]$samples,
        data = leave_out_data,
        covariate_names = list(covariate_names[[1]], NULL, NULL),
        s_est = s_est[[k]],
        mesh = mesh,
        coords = st_geometry(leave_out_data))
    }
    est_pars = purrr::transpose(est_pars)
    for (k in seq_along(est_pars)) {
      if (is.null(dim(est_pars[[k]][[1]]))) {
        est_pars[[k]] = do.call(c, est_pars[[k]])
      } else {
        est_pars[[k]] = do.call(cbind, est_pars[[k]])
      }
    }

    twcrps = list()
    for (k in seq_along(unique(leave_out_data$id))) {
      id = unique(leave_out_data$id)[k]
      obs = dplyr::filter(data, id == !!id)$value
      locscale_pars = locspread_to_locscale(est_pars$q[k, ], est_pars$s[k, ],
                                            est_pars$ξ[k, ], α, β)
      twcrps[[k]] = twcrps_gev(obs, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ, p0)
    }
    stats[[i]]$in_sample_twostep[[j]] = twcrps

    if (!is.null(in_sample_res2)) {
      est_pars2 = inla_bgev_pars(
        samples = in_sample_samples2,
        data = leave_out_data,
        covariate_names = covariate_names,
        s_est = rep(in_sample_res2$standardising_const, nrow(leave_out_data)),
        mesh = mesh,
        coords = st_geometry(leave_out_data))

      twcrps2 = list()
      for (k in seq_along(unique(leave_out_data$id))) {
        id = unique(data$id)[k]
        obs = dplyr::filter(data, id == !!id)$value
        locscale_pars = locspread_to_locscale(est_pars2$q[k, ], est_pars2$s[k, ],
                                              est_pars2$ξ[k, ], α, β)
        twcrps2[[k]] = twcrps_gev(obs, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ, p0)
      }
      stats[[i]]$in_sample_joint[[j]] = twcrps2
    }

    message("Done with fold nr. ", j)
  }

  # Print results for this specific aggregation period
  message("Number of succesful runs: ", length(samples), " of ", n_sd_samples)
  message("=========================================\n",
          hour_vec[i], " hour(s)\n",
          "=========================================")
  message("In sample, twostep:")
  print(data_stats(unlist(stats[[i]]$in_sample_twostep)))
  message("In sample, joint:")
  print(data_stats(unlist(stats[[i]]$in_sample_joint)))
  message("Out of sample, twostep:")
  print(data_stats(unlist(stats[[i]]$out_of_sample_twostep)))
  message("Out of sample, joint:")
  print(data_stats(unlist(stats[[i]]$out_of_sample_joint)))
}

saveRDS(stats, file.path(here::here(), "inst", "extdata", "cross-validation.rds"))

# Print the final results
for (i in seq_along(stats)) {
  for (j in seq_along(stats[[i]])) {
    stats[[i]][[j]] = data_stats(unlist(stats[[i]][[j]]))
  }
}
for (i in seq_along(stats)) {
  message("=========================================\n",
          hour_vec[i], " hour(s)\n",
          "=========================================")
  message("In sample, twostep:")
  print(stats[[i]]$in_sample_twostep)
  message("In sample, joint:")
  print(stats[[i]]$in_sample_joint)
  message("Out of sample, twostep:")
  print(stats[[i]]$out_of_sample_twostep)
  message("Out of sample, joint:")
  print(stats[[i]]$out_of_sample_joint)
}
