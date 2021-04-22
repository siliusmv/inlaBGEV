
# This script contains functions for computing summary statistics on large data-sets that are
# outputs from the functions in `inla-prediction.R`


#' @export
compute_data_stats = function(pars, q = c(.025, .25, .5, .75, .975)) {
  varname = as.character(substitute(pars))
  if (!is.list(pars)) {
    pars = list(pars)
    names(pars) = varname
  }
  res = list()
  for (name in names(pars)) {
    stats = data_stats(pars[[name]], q)
    stats$par = name
    res[[name]] = stats
  }
  res
}

#' @export
data_stats = function(x, q = c(.025, .25, .5, .75, .975)) {
  if (all(is.null(dim(x)))) {
    quantiles = quantile(x, probs = q, na.rm = TRUE)
    μ = mean(x, na.rm = TRUE)
    σ = sd(x, na.rm = TRUE)
    res = c(mean = μ, sd = σ, quantiles)
    res = as.data.frame(t(res))
  } else {
    quantiles = matrixStats::rowQuantiles(x, probs = q, na.rm = TRUE)
    μ = matrixStats::rowMeans2(x, na.rm = TRUE)
    σ = matrixStats::rowSds(x, na.rm = TRUE)
    res = cbind(mean = μ, sd = σ, quantiles)
    res = as.data.frame(res)
  }
  res
}
