
# This script provides the necessary functions for working with the GEV distribution


#' @export
pgev = function(x, μ, σ, ξ) {
  fix_lengths(x, μ, σ, ξ)
  ifelse(ξ == 0,
         exp(-exp(- (x - μ) / σ)),
         exp(-pmax(0, 1 + ξ * (x - μ) / σ) ^ (-1 / ξ)))
}

#' @export
qgev = function(p, μ, σ, ξ) {
  fix_lengths(p, μ, σ, ξ)
  ifelse(ξ == 0,
         μ - σ * log(-log(p)),
         μ - σ * (1 / ξ) * (1 - (- log(p)) ^ (-ξ)))
}

#' @export
rgev = function(n, μ, σ, ξ) {
  lengths = sapply(list(μ, σ, ξ), length)
  if (any(lengths > n)) stop("Bad input lengths")
  qgev(runif(n), μ, σ, ξ)
}

#' @export
dgev = function(x, μ, σ, ξ, log = FALSE) {
  fix_lengths(x, μ, σ, ξ)
  res = ifelse(ξ == 0,
               -exp(- (x - μ) / σ),
               -pmax(0, 1 + ξ * (x - μ) / σ) ^ (-1 / ξ))
  res = res - log(σ) +
    ifelse(ξ == 0,
           - (x - μ) / σ,
           ifelse(1 + ξ * (x - μ) / σ > 0,
                  - (1 / ξ + 1) * log(1 + ξ * (x - μ) / σ),
                  -Inf))
  if (!log) res = exp(res)
  res
}

#' @export
return_level_gev = function(period, μ, σ, ξ) {
  if (any(period <= 1)) warning("invalid period")
  p = ifelse(period > 1, 1 - 1 / period, NA)
  qgev(p, μ, σ, ξ)
}

#' @export
locscale_to_locspread = function(μ, σ, ξ, α = .5, β = .8) {
  fix_lengths(μ, σ, ξ)
  s = σ * (ℓ(1 - β / 2, ξ) - ℓ(β / 2, ξ)) / ifelse(ξ == 0, 1, ξ)
  q = μ + s * ℓ(α, ξ) / (ℓ(1 - β / 2, ξ) - ℓ(β / 2, ξ))
  list(q = q, s = s, ξ = ξ)
}

#' @export
locspread_to_locscale = function(q, s, ξ, α = .5, β = .8) {
  fix_lengths(q, s, ξ)
  μ = q - s * ℓ(α, ξ) / (ℓ(1 - β / 2, ξ) - ℓ(β / 2, ξ))
  σ = s / (ℓ(1 - β / 2, ξ) - ℓ(β / 2, ξ)) * ifelse(ξ == 0, 1, ξ)
  list(μ = μ, σ = σ, ξ = ξ)
}

ℓ = function(x, ξ) {
  ifelse(ξ == 0, -log(-log(x)), (-log(x))^-ξ - 1)
}


#' @export
stwcrps_gev = function(y, μ, σ, ξ, p) {
  S = abs(expected_twcrps_gev(μ, σ, ξ, p))
  twcrps = twcrps_gev(y, μ, σ, ξ, p)
  twcrps / S + log(S)
}

#' @export
twcrps_gev = function(y, μ, σ, ξ, p) {
  F = function(x) sapply(x, function(z) mean(pgev(z, μ, σ, ξ)))
  quantiles = qgev(p, μ, σ, ξ)
  if (length(quantiles) == 1) {
    y_min = quantiles
  } else {
    y_min = uniroot(function(x) F(x) - p, lower = min(quantiles), upper = max(quantiles))$root
  }
  p_max = .999
  y_max = max(qgev(p_max, μ, σ, ξ), max(y) + 1)
  res = rep(0, length(y))
  for (i in seq_along(y)) {
    if (y[i] < y_min) {
      res[i] = integrate(function(x) (1 - F(x))^2, lower = y_min, upper = y_max)$value
    } else if (y[i] < y_max) {
      res[i] = integrate(function(x) F(x)^2, lower = y_min, upper = y[i])$value
      res[i] = res[i] + integrate(function(x) (1 - F(x))^2, lower = y[i], upper = y_max)$value
    } else {
      res[i] = integrate(function(x) F(x)^2, lower = y_min, upper = y_max)$value
    }
  }
  res = res + (y_min - y) * (ifelse(y <= y_min, 1, 0) - p)^2
  res
}

#' @export
expected_twcrps_gev = function(μ, σ, ξ, p,
                                μ_true = μ, σ_true = σ, ξ_true = ξ) {
  p_min = .00001
  y_min = min(qgev(p_min, μ_true, σ_true, ξ_true))
  p_max = .99999
  y_max = max(qgev(p_max, μ_true, σ_true, ξ_true))
  if (length(c(μ_true, σ_true, ξ_true)) == 3) {
    density = function(x) dgev(x, μ_true, σ_true, ξ_true)
  } else {
    density = function(x) sapply(x, function(z) mean(dgev(z, μ_true, σ_true, ξ_true)))
  }
  integrate(function(y) density(y) * twcrps_gev(y, μ, σ, ξ, p),
            lower = y_min, upper = y_max)$value
}


