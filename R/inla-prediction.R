
# This script contains functions for estimating parameters at various locations
# using output from the R-INLA package

#' @export
inla_stats = function(sample_list,
                      data,
                      covariate_names,
                      s_list = NULL,
                      mesh = NULL,
                      n_batches = 1L,
                      verbose = TRUE,
                      fun = NULL,
                      gev.scale.xi = .1,
                      quantiles = c(.025, .25, .5, .75, .975),
                      family = c("bgev", "gaussian", "gev")) {
  # Divide the data into batches. This is necessary if you use a laptop with "finite" RAM
  batch_index = floor(seq(0, nrow(data), length.out = n_batches + 1))
  stats = vector("list", n_batches)
  pb = get_progress_bar(n_batches)
  for (b in seq_len(n_batches)) {
    rows = (batch_index[b] + 1):batch_index[b + 1]
    pars = inla_multiple_sample_pars(
      sample_list = sample_list,
      data = data,
      covariate_names = covariate_names,
      rows = rows,
      s_list = s_list,
      gev.scale.xi = gev.scale.xi,
      mesh = mesh,
      fun = fun,
      family = family)
    stats[[b]] = compute_data_stats(pars, quantiles)
    if (verbose) pb$tick()
  }
  pb$terminate()
  # Some tidying of the data
  stats = purrr::transpose(stats)
  for (i in seq_along(stats)) stats[[i]] = do.call(rbind, stats[[i]])
  stats$ξ = stats$ξ[1, ]
  stats
}

#' @export
inla_multiple_sample_pars = function(sample_list,
                                     data,
                                     covariate_names,
                                     s_list = NULL,
                                     mesh = NULL,
                                     rows = seq_len(nrow(data)),
                                     fun = NULL,
                                     gev.scale.xi = .1,
                                     family = c("bgev", "gaussian", "gev")) {
  pars = list()
  if (any(class(data) %in% c("sf", "sfc"))) {
    mydata = sf::st_drop_geometry(data)[rows, ]
    mycoords = sf::st_geometry(data)[rows]
  } else {
    mydata = data[rows, ]
    mycoords = NULL
  }
  for (i in seq_along(sample_list)) {
    # Compute posterior samples for the parameters at all locations
    if (family[1] == "bgev") {
      pars[[i]] = inla_bgev_pars(
        samples = sample_list[[i]],
        data = mydata,
        covariate_names = covariate_names,
        fun = fun,
        s_est = s_list[[i]][rows],
        mesh = mesh,
        coords = mycoords)
    } else if (family[1] == "gaussian") {
      pars[[i]] = inla_gaussian_pars(
        samples = sample_list[[i]],
        data = mydata,
        covariate_names = covariate_names,
        fun = fun,
        mesh = mesh,
        coords = mycoords)
    } else {
      stop("unknown family type")
    }
  }
  # Compute stats for the parameters and possibly the function values at all locations
  pars = purrr::transpose(pars)
  for (i in seq_along(pars)) {
    if (is.null(dim(pars[[i]][[1]]))) {
      pars[[i]] = do.call(c, pars[[i]])
    } else {
      pars[[i]] = do.call(cbind, pars[[i]])
    }
  }
  pars
}


#' @export
inla_gaussian_pars = function(samples,
                              data,
                              covariate_names,
                              mesh = NULL,
                              fun = NULL,
                              coords = NULL) {

  # Extract design matrix
  X = dplyr::select(data, tidyselect::all_of(unique(unlist(covariate_names)))) %>%
    dplyr::distinct(.keep_all = TRUE)
  if (is(X, "sf") || is(X, "sfc")) X = sf::st_drop_geometry(X)
  X = cbind(intercept = 1, as.matrix(X))

  # Extract all coefficients from the samples
  coeffs = inla_gaussian_coeffs(samples, covariate_names)

  # Compute parameters at all locations and for all samples
  if (is.null(dim(coeffs$μ))) {
    μ = coeffs$μ
  } else {
    μ = X[, c("intercept", covariate_names)] %*% coeffs$μ
    #μ_intercept = matrix(X[, "intercept"], nrow = nrow(X)) %*% coeffs$μ[1, ]
    #μ_no_intercept = X[, covariate_names] %*% coeffs$μ[-1, ]
  }
  τ = matrix(rep(coeffs$τ, each = nrow(X)), ncol = length(samples))

  # If there is a spatial Gaussian field in the model, sample from it
  # and add it to the location parameter
  is_matern_field = any(attributes(samples)$.contents$tag == "matern_field")
  #if (is_matern_field) matern = inla_sample_matern_field(samples, mesh, coords)
  if (is_matern_field) μ = μ + inla_sample_matern_field(samples, mesh, coords)

  pars = list(μ = μ, τ = τ, coeffs = coeffs$μ)
  # Compute some function of the sampled parameters, e.g. return level
  if (!is.null(fun)) {
    if (is.function(fun)) {
      pars$fun = matrix(fun(pars), nrow = nrow(X))
    } else {
      for (j in seq_along(fun)) {
        pars[[names(fun)[j]]] = matrix(fun[[j]](pars), nrow = nrow(X))
      }
    }
  }
  pars
}

#' @export
inla_bgev_pars = function(samples,
                          data,
                          covariate_names,
                          s_est = NULL,
                          fun = NULL,
                          mesh = NULL,
                          coords = NULL) {

  # Extract design matrix
  X = dplyr::select(data, tidyselect::all_of(unique(unlist(covariate_names)))) %>%
    dplyr::distinct(.keep_all = TRUE)
  if (is(X, "sf") || is(X, "sfc")) X = sf::st_drop_geometry(X)
  X = cbind(intercept = 1, as.matrix(X))

  # Extract all coefficients from the samples
  coeffs = inla_bgev_coeffs(samples, covariate_names)

  # Compute parameters at all locations and for all samples
  q = matrix(X[, c("intercept", covariate_names[[1]])], nrow = nrow(X)) %*% coeffs$q
  if (is.null(dim(coeffs$s))) {
    s = matrix(rep(coeffs$s, each = nrow(X)), ncol = length(samples))
  } else {
    s = matrix(X[, "intercept"], nrow = nrow(X)) %*% coeffs$s[1, ] *
      exp(matrix(X[, covariate_names[[2]]], nrow = nrow(X)) %*% coeffs$s[-1, ])
  }
  ξ = matrix(rep(coeffs$ξ, each = nrow(X)), ncol = length(samples))

  # If there is a spatial Gaussian field in the model, sample from it
  # and add it to the location parameter
  is_matern_field = any(attributes(samples)$.contents$tag == "matern_field")
  #if (is_matern_field) q = q + inla_sample_matern_field(samples, mesh, coords)
  if (is_matern_field) {
    matern = inla_sample_matern_field(samples, mesh, coords)
    q = q + matern
  }

  # If the response had been standardised, compute un-standardised parameters
  if (!is.null(s_est)) {
    s = s * matrix(rep(s_est, length(samples)), ncol = length(samples))
    q = q * matrix(rep(s_est, length(samples)), ncol = length(samples))
    if (is_matern_field) {
      matern = matern * matrix(rep(s_est, length(samples)), ncol = length(samples))
    }
  }

  if (is_matern_field) {
    pars = list(q = q, s = s, ξ = ξ, matern = matern)
  } else {
    pars = list(q = q, s = s, ξ = ξ)
  }
  # Compute some function of the sampled parameters, e.g. return level
  if (!is.null(fun)) {
    if (is.function(fun)) {
      pars$fun = matrix(fun(pars), nrow = nrow(X))
    } else {
      for (j in seq_along(fun)) {
        pars[[names(fun)[j]]] = matrix(fun[[j]](pars), nrow = nrow(X))
      }
    }
  }
  pars
}

inla_bgev_coeffs = function(samples, covariate_names) {
  contents = attributes(samples)$.contents
  q_covariates = c("intercept", covariate_names[[1]])
  q_indx = sapply(q_covariates, function(cov) which(contents$tag == cov))
  q_start = contents$start[q_indx]
  q_samples = sapply(samples, function(s) s$latent[q_start])
  ξ_indx = grep("tail", names(samples[[1]]$hyperpar))
  ξ_samples = sapply(samples, function(s) s$hyperpar[ξ_indx])
  s_indx = grep("spread", names(samples[[1]]$hyperpar))
  s_samples = sapply(samples, function(s) s$hyperpar[s_indx])
  list(q = q_samples, s = s_samples, ξ = ξ_samples)
}

inla_gaussian_coeffs = function(samples, covariate_names) {
  contents = attributes(samples)$.contents
  covariate_names = unique(c("intercept", covariate_names))
  μ_index = sapply(covariate_names, function(cov) which(contents$tag == cov))
  μ_start = contents$start[μ_index]
  μ_samples = sapply(samples, function(s) s$latent[μ_start])
  τ_indx = grep("Precision", names(samples[[1]]$hyperpar))
  τ_samples = sapply(samples, function(s) s$hyperpar[τ_indx])
  list(μ = μ_samples, τ = τ_samples)
}

inla_sample_matern_field = function(samples, mesh, coords) {
  contents = attributes(samples)$.contents
  indx = which(contents$tag == "matern_field")
  if (length(indx) == 0) return(0)
  A = INLA::inla.spde.make.A(mesh, sf::st_coordinates(sf::st_transform(coords, get_proj_xy())))
  l = contents$length[indx]
  start = contents$start[indx]
  matern_samples = sapply(samples, function(s) s$latent[start:(start + l - 1)])
  matrix((A %*% matern_samples)@x, nrow = nrow(A), ncol = ncol(matern_samples))
}
