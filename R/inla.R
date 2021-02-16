
#' @export
inla_default_args = function(family = "bgev", α = .5, β = .8) {
  res = list(
    family = family,
    control.compute = list(
      # This allows us to always do predictions via sampling
      config = TRUE,
      # And this allows us for easy cross-validation
      cpo = TRUE, dic = TRUE, waic = TRUE),
    # use identity link and compute marginals
    control.predictor = list(compute = TRUE, link = 1),
    # Make inla more stable
    control.inla = list(tolerance = 1e-12, h = 1e-3),
    # There are some race-conditions happening when num.threads > 1, making INLA less stable
    num.threads = 1)
  if (family == "bgev") {
    # Tell INLA what we are using for α and β
    res$control.family = list(control.bgev = list(q.location = α, q.spread = β))
  }
  res
}

#' @export
inla_location_prior = function(coeffs, precision = 10) {
  fix_lengths(coeffs, precision)
  mean = list(); prec = list()
  if ("(Intercept)" %in% names(coeffs)) {
    names(coeffs) = sub("\\(Intercept\\)", "intercept", names(coeffs))
  }
  for (i in seq_along(coeffs)) {
    mean[[names(coeffs)[i]]] = coeffs[i]
    prec[[names(coeffs)[i]]] = precision[i]
  }
  list(mean = mean, prec = prec)
}

#' @export
inla_tail_prior = function(lims = c(0, .5), init = .1, λ = 7, fixed = FALSE) {
  list(prior = "pc.gevtail",
       fixed = fixed,
       initial = get_internal_bgev_tail(init, lims),
       param = c(λ, lims))
}

get_internal_bgev_tail = function(x, lims = c(0, .5)) log((x - lims[1]) / (lims[2] - x))
get_outer_bgev_tail = function(x, lims = c(0, .5)) lims[1] + diff(lims) * exp(x) / (1 + exp(x))

#' @export
inla_spread_prior = function(coeffs, precision = 10) {
  fix_lengths(coeffs, precision)
  prior = list()
  prior$spread = list(
    prior = "loggamma",
    # These parameters result in a gamma prior for the spread
    # with mean `coeffs` and precision `precision`
    param = precision[1] * coeffs[1]^c(2, 1))
  if (length(coeffs) > 1) {
    for (i in 2:length(coeffs)) {
      prior[[paste0("beta", i - 1)]] = list(
        prior = "gaussian",
        param = as.numeric(c(coeffs[i], precision[i])))
    }
  }
  prior
}

#' @export
inla_stack = function(df, covariate_names, response_name, spde = NULL, tag = "tag") {
  stack_response = list(df[[response_name]])
  names(stack_response) = response_name
  df$intercept = 1
  effects_df = df[, c("intercept", unique(unlist(covariate_names)))]
  if (is(effects_df, "sf") || is(effects_df, "sfc")) effects_df = sf::st_drop_geometry(effects_df)
  effects = list(as.data.frame(effects_df))
  names(effects[[1]]) = c("intercept", unique(unlist(covariate_names)))
  A = list(1)
  if (!is.null(spde)) {
    coords = sf::st_coordinates(sf::st_transform(sf::st_as_sf(df), get_proj_xy()))
    A = c(list(INLA::inla.spde.make.A(spde$mesh, coords)), A)
    effects = c(list(matern_field = seq_len(spde$mesh$n)), effects)
  }
  INLA::inla.stack(data = stack_response, A = A, effects = effects, tag = tag)
}

#' @export
inla_bgev = function(data,
                     covariate_names,
                     response_name,
                     s_est = rep(1, nrow(data)),
                     spde = NULL,
                     α = .5,
                     β = .8,
                     ...) {
  # Standardise the observations by dividing on the estimated spread, and then
  # dividing by a constant so the difference between 5% quantile and 95% quantile is 1
  data[[response_name]] = data[[response_name]] / s_est
  standardising_const = diff(quantile(data[[response_name]], c(.05, .95)))
  data[[response_name]] = data[[response_name]] / standardising_const

  # Create the design matrix
  X = dplyr::select(data, tidyselect::all_of(unique(unlist(covariate_names))))
  if (is(X, c("sf", "sfc"))) X = sf::st_drop_geometry(X)
  X = as.matrix(X)

  # Create data structures used by R-INLA
  stack = inla_stack(data, covariate_names, response_name = response_name, spde = spde)
  mdata = INLA::inla.mdata(
    data[[response_name]], X[, covariate_names[[2]]], X[, covariate_names[[3]]])

  # Create prior(s) for the location coefficient(s)
  if (length(covariate_names[[1]]) > 0) {
    q_formula = as.formula(paste(response_name, "~", paste(covariate_names[[1]], collapse = " + ")))
    est_q_coeffs = quantreg::rq(q_formula, α, data)$coefficients
  } else {
    est_q_coeffs = c("intercept" = median(data[[response_name]]))
  }
  control.fixed = inla_location_prior(est_q_coeffs, precision = 10)

  # Create prior(s) for the spread coefficient(s)
  spread_mean = diff(quantile(data[[response_name]], c(β / 2, 1 - β / 2)))
  n_s = length(covariate_names[[2]])
  spread_par = spread_mean
  spread_precision = 10
  if (n_s > 0) {
    spread_par = c(spread_par, rep(0, n_s))
    spread_precision = c(spread_precision, rep(1, n_s))
  }
  hyper = inla_spread_prior(spread_par, spread_precision)

  # Create a prior for the tail parameter
  hyper$tail = inla_tail_prior(lims = c(0, .5), λ = 7)

  # Create the formula needed by R-INLA
  if (is.null(spde)) {
    inla_formula = as.formula(paste("mdata ~ -1 + ",
                                    paste(c("intercept", covariate_names[[1]]), collapse = " + ")))
  } else {
    inla_formula = as.formula(paste("mdata ~ -1 + ",
                                    paste(c("intercept", covariate_names[[1]]), collapse = " + "),
                                    "+ f(matern_field, model = spde)"))
  }

  # Create a list of all the arguments needed for running R-INLA
  inla_args = inla_default_args("bgev", α = α, β = β)
  inla_args$formula = inla_formula
  inla_args$control.predictor$A = INLA::inla.stack.A(stack)
  inla_args$data = INLA::inla.stack.data(stack)
  inla_args$data$mdata = mdata
  inla_args$data$spde = spde
  inla_args$control.fixed = control.fixed
  inla_args$control.family$hyper = hyper

  # Add more arguments to INLA if they were provided in the ellipsis (...)
  more_args = list(...)
  for (name in names(more_args)) inla_args[[name]] = more_args[[name]]

  # Run R-INLA
  res = do.call(INLA::inla, inla_args)

  res$standardising_const = standardising_const
  res
}

#' @export
inla_mesh = function(coords,
                     convex = -.15, # The non-convex boundary
                     offset = c(10, 100), # lengths from points to boundary?
                     cutoff = 5, # Minimum edge length (when multiple points are close)
                     max.edge = c(25, 50)) {# Max edge length for inner/outer area
  coords = sf::st_transform(coords, get_proj_xy())
  boundary = inla.nonconvex.hull(points = sf::st_coordinates(coords), convex = convex)
  mesh_crs = sp::CRS(as.character(get_proj_xy())[1])
  INLA::inla.mesh.2d(
    sf::st_coordinates(coords),
    boundary = boundary,
    max.edge = max.edge,
    cutoff = cutoff,
    offset = offset,
    crs = mesh_crs)
}
