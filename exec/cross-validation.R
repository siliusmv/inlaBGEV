library(sf)
library(dplyr)
library(inlaBGEV)
library(INLA)

# In this script we perform k-fold cross-validation on the joint model and on the
# two-step model. The performance of the foreacsts are evaluated using stwCRPS, where
# the weight function is an indicator function that is one for larger quantiles than
# the p_0 quantile

hour_vec = c(1, 3, 6, 12, 24) # Which aggregation lengths are we examining?
α = .5; β = .8 # Probabilities used in the location and spread parameters
min_sd_years = 4L # Minimum number of years before we use the computed SD values
n_sd_samples = 100 # Number of samples drawn from the distribution of the SD
num_cores = 100 # Number of cores used for parallel computations
n_folds = 10 # number of folds for cross-validation
p0 = .9 # Threshold used in the stwCRPS

# A list containing covariate_names for location, spread and tail parameter
covariate_names = list(c("precipitation", "height", "x", "y", "dist_sea"),
                       c("x", "y", "dist_sea"), NULL)

stats = list()
for (i in seq_along(hour_vec)) {

  # Create the empty list that will contain all the result
  stats[[i]] = list(
    in_sample_twostep = list(),
    out_of_sample_twostep = list(),
    in_sample_twostep_one = list(),
    out_of_sample_twostep_one = list(),
    in_sample_joint = list(),
    out_of_sample_joint = list())
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
  sd_res = do.call(inla, sd_inla_args)

  # Perform in-sample estimation using all the data
  set.seed(1)
  twostep = twostep_modelling(
    data = data,
    sd_model = sd_res,
    covariate_names = covariate_names,
    response_name = "value",
    diagonal = .05,
    n_sd_samples = n_sd_samples,
    spde = spde,
    num_cores = num_cores,
    α = α,
    β = β)
  message("Done with in-sample two-step model")

  # Estimate parameters at all locations
  params = list()
  for (k in seq_along(twostep)) {
    params[[k]] = inla_bgev_pars(
      samples = twostep[[k]]$samples,
      data = data,
      covariate_names = list(covariate_names[[1]], NULL, NULL),
      s_est = twostep[[k]]$s_est,
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
    stats[[i]]$in_sample_twostep[[k]] = stwcrps_bgev(
      obs, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ, p0)
  }

  # Perform in-sample estimation using all the data
  set.seed(1)
  twostep_one = twostep_modelling(
    data = data,
    sd_model = sd_res,
    covariate_names = covariate_names,
    response_name = "value",
    diagonal = .05,
    n_sd_samples = 1,
    spde = spde,
    num_cores = 1,
    α = α,
    β = β)
  message("Done with in-sample two-step model")

  # Estimate parameters at all locations
  params = list()
  for (k in seq_along(twostep_one)) {
    params[[k]] = inla_bgev_pars(
      samples = twostep_one[[k]]$samples,
      data = data,
      covariate_names = list(covariate_names[[1]], NULL, NULL),
      s_est = twostep_one[[k]]$s_est,
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
    stats[[i]]$in_sample_twostep_one[[k]] = stwcrps_bgev(
      obs, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ, p0)
  }

  # Run the joint model on the data
  joint = tryCatch(inla_bgev(
    data = data,
    covariate_names = covariate_names,
    response_name = "value",
    diagonal = .1,
    spde = spde,
    α = α,
    β = β), error = function(e) NULL)
  message("Done with in-sample joint model")

  if (!is.null(joint)) {
    # Sample from the posterior of the joint model
    set.seed(1)
    samples = inla.posterior.sample(
      sum(sapply(twostep, function(x) length(x$samples))), joint, seed = 1)

    # Estimate parameters at all locations
    params = inla_bgev_pars(
      samples = samples,
      data = data,
      covariate_names = covariate_names,
      s_est = rep(joint$standardising_const, length(unique(data$id))),
      mesh = mesh,
      coords = st_geometry(dplyr::distinct(data)))

    # Compute stwCRPS
    for (k in seq_along(unique(data$id))) {
      id = unique(data$id)[k]
      obs = dplyr::filter(data, id == !!id)$value
      locscale_pars = locspread_to_locscale(params$q[k, ], params$s[k, ], params$ξ[k, ], α, β)
      stats[[i]]$in_sample_joint[[k]] = stwcrps_bgev(
        obs, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ, p0)
    }
  }

  message("Start with the out-of-sample modelling")
  for (j in seq_len(n_folds)) {

    # Split the data up into in-sample and out-of-sample data for the current folds
    in_fold_sd_df = dplyr::filter(sd_df, id %in% ids[folds != j])
    in_fold_data = dplyr::filter(data, id %in% ids[folds != j])
    out_of_fold_data = dplyr::filter(data, id %in% ids[folds == j])

    # Estimate s^*
    sd_stack = inla_stack(in_fold_sd_df, covariate_names[[2]],
                          response_name = "log_sd", spde = sd_spde)
    sd_inla_formula = as.formula(paste("log_sd ~ -1 + intercept +",
                                       paste(covariate_names[[2]], collapse = " + "),
                                       "+ f(matern_field, model = sd_spde)"))
    sd_inla_args = inla_default_args("gaussian")
    sd_inla_args$formula = sd_inla_formula
    sd_inla_args$control.predictor$A = INLA::inla.stack.A(sd_stack)
    sd_inla_args$data = INLA::inla.stack.data(sd_stack)
    sd_inla_args$data$sd_spde = sd_spde
    sd_res = do.call(inla, sd_inla_args)

    # Perform out-of-sample estimation with the two-step model, with bootsrtapping
    set.seed(1)
    twostep = twostep_modelling(
      data = in_fold_data,
      prediction_data = dplyr::distinct(out_of_fold_data, id, .keep_all = TRUE),
      sd_model = sd_res,
      covariate_names = covariate_names,
      response_name = "value",
      diagonal = .05,
      n_sd_samples = n_sd_samples,
      verbose = FALSE,
      spde = spde,
      num_cores = num_cores,
      α = α,
      β = β)
    message("Done with out-of-sample two-step model for fold ", j)

    # Estimate parameters at all out-of-fold locations
    params = list()
    for (k in seq_along(twostep)) {
      params[[k]] = inla_bgev_pars(
        samples = twostep[[k]]$samples,
        data = dplyr::distinct(out_of_fold_data, id, .keep_all = TRUE),
        covariate_names = list(covariate_names[[1]], NULL, NULL),
        s_est = twostep[[k]]$s_est,
        mesh = mesh,
        coords = st_geometry(dplyr::distinct(out_of_fold_data, id)))
    }
    params = purrr::transpose(params)
    for (k in seq_along(params)) params[[k]] = do.call(cbind, params[[k]])

    # Compute stwCRPS
    stats[[i]]$out_of_sample_twostep[[j]] = list()
    for (k in seq_along(unique(out_of_fold_data$id))) {
      id = unique(out_of_fold_data$id)[k]
      obs = dplyr::filter(out_of_fold_data, id == !!id)$value
      locscale_pars = locspread_to_locscale(params$q[k, ], params$s[k, ], params$ξ[k, ], α, β)
      stats[[i]]$out_of_sample_twostep[[j]][[k]] = stwcrps_bgev(
        obs, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ, p0)
    }

    # Perform out-of-sample estimation with the two-step model, without bootstrapping
    set.seed(1)
    twostep_one = twostep_modelling(
      data = in_fold_data,
      prediction_data = dplyr::distinct(out_of_fold_data, id, .keep_all = TRUE),
      sd_model = sd_res,
      covariate_names = covariate_names,
      response_name = "value",
      diagonal = .05,
      n_sd_samples = 1,
      verbose = FALSE,
      spde = spde,
      num_cores = 1,
      α = α,
      β = β)
    message("Done with out-of-sample two-step model for fold ", j)

    # Estimate parameters at all out-of-fold locations
    params = list()
    for (k in seq_along(twostep_one)) {
      params[[k]] = inla_bgev_pars(
        samples = twostep_one[[k]]$samples,
        data = dplyr::distinct(out_of_fold_data, id, .keep_all = TRUE),
        covariate_names = list(covariate_names[[1]], NULL, NULL),
        s_est = twostep_one[[k]]$s_est,
        mesh = mesh,
        coords = st_geometry(dplyr::distinct(out_of_fold_data, id)))
    }
    params = purrr::transpose(params)
    for (k in seq_along(params)) params[[k]] = do.call(cbind, params[[k]])

    # Compute stwCRPS
    stats[[i]]$out_of_sample_twostep_one[[j]] = list()
    for (k in seq_along(unique(out_of_fold_data$id))) {
      id = unique(out_of_fold_data$id)[k]
      obs = dplyr::filter(out_of_fold_data, id == !!id)$value
      locscale_pars = locspread_to_locscale(params$q[k, ], params$s[k, ], params$ξ[k, ], α, β)
      stats[[i]]$out_of_sample_twostep_one[[j]][[k]] = stwcrps_bgev(
        obs, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ, p0)
    }

    # Use joint modelling, out-of-sample
    joint = tryCatch(inla_bgev(
      data = in_fold_data,
      covariate_names = covariate_names,
      response_name = "value",
      diagonal = .1,
      spde = spde,
      α = α,
      β = β), error = function(e) NULL)
    message("Done with out-of-sample joint model for fold ", j)

    if (!is.null(joint)) {
      set.seed(1)
      samples = inla.posterior.sample(
        sum(sapply(twostep, function(x) length(x$samples))), joint, seed = 1)

      # Compute sampled parameters at all leave-out locations
      params = inla_bgev_pars(
        samples = samples,
        data = dplyr::distinct(out_of_fold_data, id, .keep_all = TRUE),
        covariate_names = covariate_names,
        s_est = rep(joint$standardising_const, length(unique(out_of_fold_data$id))),
        mesh = mesh,
        coords = st_geometry(dplyr::distinct(out_of_fold_data, id)))

      # Compute stwCRPS at all leave-out locations
      stats[[i]]$out_of_sample_joint[[j]] = list()
      for (k in seq_along(unique(out_of_fold_data$id))) {
        id = unique(out_of_fold_data$id)[k]
        obs = dplyr::filter(out_of_fold_data, id == !!id)$value
        locscale_pars = locspread_to_locscale(params$q[k, ], params$s[k, ], params$ξ[k, ], α, β)
        stats[[i]]$out_of_sample_joint[[j]][[k]] = stwcrps_bgev(
          obs, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ, p0)
      }
    }

    message("Done with fold nr. ", j)
  }

  # Print results for this specific aggregation period
  message("=========================================\n",
          hour_vec[i], " hour(s)\n",
          "=========================================")
  message("In sample, twostep with bootstrapping:")
  print(data_stats(unlist(stats[[i]]$in_sample_twostep)))
  message("In sample, twostep without bootstrapping:")
  print(data_stats(unlist(stats[[i]]$in_sample_twostep_one)))
  message("In sample, joint:")
  print(data_stats(unlist(stats[[i]]$in_sample_joint)))
  message("Out of sample, twostep with bootstrapping:")
  print(data_stats(unlist(stats[[i]]$out_of_sample_twostep)))
  message("Out of sample, twostep without bootstrapping:")
  print(data_stats(unlist(stats[[i]]$out_of_sample_twostep_one)))
  message("Out of sample, joint:")
  print(data_stats(unlist(stats[[i]]$out_of_sample_joint)))
}

saveRDS(stats, file.path(here::here(), "results", "cross-validation.rds"))

# Print the final results
q = c(.025, .25, .5, .75, .975, .99, 1)
for (i in seq_along(stats)) {
  message("=========================================\n",
          hour_vec[i], " hour(s)\n",
          "=========================================")
  message("In sample, twostep with bootstrapping:")
  print(data_stats(unlist(stats[[i]]$in_sample_twostep), q))
  message("In sample, twostep without bootstrapping:")
  print(data_stats(unlist(stats[[i]]$in_sample_twostep_one), q))
  message("In sample, joint:")
  print(data_stats(unlist(stats[[i]]$in_sample_joint), q))
  message("Out of sample, twostep with bootstrapping:")
  print(data_stats(unlist(stats[[i]]$out_of_sample_twostep), q))
  message("Out of sample, twostep without bootstrapping:")
  print(data_stats(unlist(stats[[i]]$out_of_sample_twostep_one), q))
  message("Out of sample, joint:")
  print(data_stats(unlist(stats[[i]]$out_of_sample_joint), q))
}

# Print the results in a latex-friendly tabular format
df = lapply(
  seq_along(stats),
  function(i) {
    res = c(
      base::mean(unlist(stats[[i]]$out_of_sample_joint)),
      base::mean(unlist(stats[[i]]$out_of_sample_twostep)),
      base::mean(unlist(stats[[i]]$out_of_sample_twostep_one)),
      base::mean(unlist(stats[[i]]$in_sample_joint)),
      base::mean(unlist(stats[[i]]$in_sample_twostep),
      base::mean(unlist(stats[[i]]$in_sample_twostep_one))))
    best_indx = c(1, 3) + as.numeric(res[c(2, 4)] < res[c(1, 3)])
    res = format(res, digits = 3)
    res[best_indx] = paste0("\\bm{", res[best_indx], "}")
    res = paste0("\\(", res, "\\)")
    res = data.frame(res)
    names(res) = paste(hour_vec[i], ifelse(hour_vec[i] == 1, "hour", "hours"))
    res
  }
) %>%
  do.call(cbind, .)
table = list(c("\n& & ", paste(names(df), collapse = " & "), "\\\\\n"))
for (i in 1:nrow(df)) {
  table[[i + 1]] = c(case_when(i == 1 ~ "\\midrule\nOut-of-sample",
                               i == 4 ~ "\\midrule\nIn-sample",
                               TRUE ~ ""), " &",
                     case_when(i %in% c(1, 4) ~ "Joint",
                               i %in% c(2, 5) ~ "Two-step with bootstrap",
                               i %in% c(3, 6) ~ "Two-step"), " model &",
                     paste(df[i, ], collapse = " & "), "\\\\\n")
}
cat(unlist(table))
