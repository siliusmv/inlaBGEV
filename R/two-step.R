
# This script contains a very large function used for applying the two-step model
# for modelling with the bGEV distribution

#' @export
twostep_modelling = function(data,
                             sd_model,
                             covariate_names,
                             response_name,
                             n_sd_samples = 1,
                             prediction_data = NULL,
                             spde = NULL,
                             num_cores = 1,
                             num_samples = 500,
                             verbose = FALSE,
                             α = .5,
                             β = .8,
                             ...) {

  # Sample from the distribution of the SD at all observation locations and prediction locations
  if (is.null(prediction_data)) {
    all_data = dplyr::distinct(data, id, .keep_all = TRUE)
    prediction_indices = seq_len(nrow(dplyr::distinct(data, id, .keep_all = TRUE)))
  } else {
    all_data = dplyr::distinct(data, id, .keep_all = TRUE) %>%
      sf::st_transform(sf::st_crs(prediction_data)) %>%
      dplyr::select(dplyr::all_of(names(prediction_data))) %>%
      rbind(prediction_data)
    prediction_indices = seq_len(nrow(all_data))[-seq_len(nrow(dplyr::distinct(data, id, .keep_all = TRUE)))]
  }
  if (n_sd_samples == 1) {
    log_sd_samples = INLA::inla.posterior.sample(1000, sd_model, seed = 1, parallel.configs = FALSE)
    log_sd_stats = inla_stats(
      sample_list = list(log_sd_samples),
      data = all_data,
      covariate_names = covariate_names[[2]],
      mesh = spde$mesh,
      verbose = verbose,
      fun = function(pars) exp(rnorm(length(pars$μ), pars$μ, 1 / sqrt(pars$τ))),
      n_batches = 20,
      family = "gaussian")
    sd_samples = log_sd_stats$fun$mean
  } else {
    log_sd_samples = INLA::inla.posterior.sample(n_sd_samples, sd_model, seed = 1, parallel.configs = FALSE)
    log_sd_pars = inla_gaussian_pars(
      samples = log_sd_samples,
      data = all_data,
      covariate_names = covariate_names[[2]],
      mesh = spde$mesh,
      coords = sf::st_geometry(all_data))
    sd_samples = matrix(exp(log_sd_pars$μ), nrow = nrow(log_sd_pars$μ))
  }
  location_indices = as.numeric(factor(data$id))
  # Run R-inla to estimate the BGEV-parameters once for each of the SD samples
  samples = parallel::mclapply(
    X = seq_len(n_sd_samples),
    mc.cores = num_cores,
    mc.preschedule = FALSE,
    FUN = function(i) {
      if (n_sd_samples == 1) {
        s_est = sd_samples
      } else {
        s_est = sd_samples[, i]
      }
      res = tryCatch({
        inla_bgev(
          data = data,
          s_est = s_est[location_indices],
          covariate_names = list(covariate_names[[1]], NULL, NULL),
          response_name = response_name,
          spde = spde,
          α = α,
          β = β,
          ...)},
        error = function(e) NULL)
      if (verbose) message("Done with iter nr. ", i)
      if (is.null(res)) return(NULL)
      samples = INLA::inla.posterior.sample(ceiling(num_samples / n_sd_samples), res, seed = 1, parallel.configs = FALSE)
      list(s_est = s_est[prediction_indices] * res$standardising_const,
           samples = samples,
           standardising_const = res$standardising_const,
           s_coeffs = log_sd_pars$coeffs[, i])
    })

  # Sometimes, R-INLA might have some numerical problems. Remove the bad models
  bad_samples = which(sapply(samples, is.null))
  if (any(bad_samples)) {
    samples = samples[-bad_samples]
  }
  if (length(samples) == 0) stop("All runs were faulty")

  samples
}

#' @export
twostep_gev = function(data,
                       sd_model,
                       covariate_names,
                       response_name,
                       n_sd_samples = 1,
                       prediction_data = NULL,
                       spde = NULL,
                       num_cores = 1,
                       num_samples = 500,
                       verbose = FALSE,
                       α = .5,
                       β = .8,
                       ...) {

  # Sample from the distribution of the SD at all observation locations and prediction locations
  if (is.null(prediction_data)) {
    all_data = dplyr::distinct(data, id, .keep_all = TRUE)
    prediction_indices = seq_len(nrow(dplyr::distinct(data, id, .keep_all = TRUE)))
  } else {
    all_data = dplyr::distinct(data, id, .keep_all = TRUE) %>%
      sf::st_transform(sf::st_crs(prediction_data)) %>%
      dplyr::select(dplyr::all_of(names(prediction_data))) %>%
      rbind(prediction_data)
    prediction_indices = seq_len(nrow(all_data))[-seq_len(nrow(dplyr::distinct(data, id, .keep_all = TRUE)))]
  }
  if (n_sd_samples == 1) {
    log_sd_samples = INLA::inla.posterior.sample(1000, sd_model, seed = 1, parallel.configs = FALSE)
    log_sd_stats = inla_stats(
      sample_list = list(log_sd_samples),
      data = all_data,
      covariate_names = covariate_names[[2]],
      mesh = spde$mesh,
      verbose = verbose,
      fun = function(pars) exp(rnorm(length(pars$μ), pars$μ, 1 / sqrt(pars$τ))),
      n_batches = 20,
      family = "gaussian")
    sd_samples = log_sd_stats$fun$mean
  } else {
    log_sd_samples = INLA::inla.posterior.sample(n_sd_samples, sd_model, seed = 1, parallel.configs = FALSE)
    log_sd_pars = inla_gaussian_pars(
      samples = log_sd_samples,
      data = all_data,
      covariate_names = covariate_names[[2]],
      mesh = spde$mesh,
      coords = sf::st_geometry(all_data))
    sd_samples = matrix(exp(log_sd_pars$μ), nrow = nrow(log_sd_pars$μ))
  }
  location_indices = as.numeric(factor(data$id))
  # Run R-inla to estimate the BGEV-parameters once for each of the SD samples
  samples = parallel::mclapply(
    X = seq_len(n_sd_samples),
    mc.cores = num_cores,
    mc.preschedule = FALSE,
    FUN = function(i) {
      if (n_sd_samples == 1) {
        s_est = sd_samples
      } else {
        s_est = sd_samples[, i]
      }
      res = tryCatch({
        inla_gev(
          data = data,
          s_est = s_est[location_indices],
          covariate_names = covariate_names[[1]],
          response_name = response_name,
          spde = spde,
          ...)},
        error = function(e) NULL)
      if (verbose) message("Done with iter nr. ", i)
      if (is.null(res)) return(NULL)
      samples = INLA::inla.posterior.sample(ceiling(num_samples / n_sd_samples), res, seed = 1, parallel.configs = FALSE)
      gev.scale.xi = res$.args$control.family[[1]]$gev.scale.xi
      list(s_est = s_est[prediction_indices] * res$standardising_const, samples = samples, gev.scale.xi = gev.scale.xi)
    })

  # Sometimes, R-INLA might have some numerical problems. Remove the bad models
  bad_samples = which(sapply(samples, is.null))
  if (any(bad_samples)) {
    samples = samples[-bad_samples]
  }
  if (length(samples) == 0) stop("All runs were faulty")

  samples
}
