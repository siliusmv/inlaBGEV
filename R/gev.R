
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
twcrps_gev = function(y, μ, σ, ξ, p) {
  if (p < 0 || p > 1) stop("p must be a probability between 0 and 1")
  res = mapply(twcrps_gev_one_par, μ = μ, σ = σ, ξ = ξ, MoreArgs = list(y = y, p = p))
  if (!is.null(dim(res)) && ncol(res) == 1) res = as.numeric(res)
  res
}

#' @export
stwcrps_gev = function(y, μ, σ, ξ, p, num_cores = 6) {
  f = function(x, μ, σ, ξ, p) twcrps_gev(x, μ, σ, ξ, p) * dgev(x, μ, σ, ξ)
  g_lower = function(u, ...) f(log(u), ...) / u
  g_upper = function(u, ...) f(-log(u), ...) / u
  one_stwcrps = function(y, μ, σ, ξ, p, p_a, p_b) {
    expected_score = integrate(
      function(u, ...) g_lower(u, ...) + g_upper(u, ...), lower = 0, upper = 1,
      μ = μ, σ = σ, ξ = ξ, p = p)$value
    twcrps_gev(y, μ, σ, ξ, p) / abs(expected_score) + log(abs(expected_score))
  }
  fix_lengths(μ, σ, ξ)
  parallel::mcmapply(
    FUN = function(y, ...) one_stwcrps(y, ...),
    mc.cores = num_cores,
    μ = μ, σ = σ, ξ = ξ,
    MoreArgs = list(y = y, p = p))
}

twcrps_gev_one_par = function(y, μ, σ, ξ, p) {
  Fy = pgev(y, μ, σ, ξ)
  Ei = function(x) {
    res = rep(0, length(x))
    res[x > -700] = gsl::expint_Ei(x[x > -700])
    res
  }
  p1 = pmax(p, Fy)
  if (ξ == 0) {
    res = (y - μ) * (2 * p1 - p^2 - 1) + 2 * σ * p1 * log(-log(p1)) +
      σ * (Ei(2 * log(p)) - 2 * Ei(log(p1)) - digamma(1) - log(2))
    if (p > 0) res = res - σ * p^2 * log(-log(p))
  } else {
    res = (y - μ + σ / ξ) * (2 * p1 - p^2 - 1) - (σ / ξ) * gamma(1 - ξ) * (
      2^ξ * pgamma(-2 * log(p), 1 - ξ) - 2 * pgamma(-log(p1), 1 - ξ))
  }
  res
}

#' @export
locscale_to_locspread = function(μ, σ, ξ, α = .5, β = .25) {
  fix_lengths(μ, σ, ξ)
  s = σ * (ℓ(1 - β / 2, ξ) - ℓ(β / 2, ξ)) / ifelse(ξ == 0, 1, ξ)
  q = μ + s * ℓ(α, ξ) / (ℓ(1 - β / 2, ξ) - ℓ(β / 2, ξ))
  list(q = q, s = s, ξ = ξ)
}

#' @export
locspread_to_locscale = function(q, s, ξ, α = .5, β = .25) {
  fix_lengths(q, s, ξ)
  μ = q - s * ℓ(α, ξ) / (ℓ(1 - β / 2, ξ) - ℓ(β / 2, ξ))
  σ = s / (ℓ(1 - β / 2, ξ) - ℓ(β / 2, ξ)) * ifelse(ξ == 0, 1, ξ)
  list(μ = μ, σ = σ, ξ = ξ)
}

ℓ = function(x, ξ) {
  ifelse(ξ == 0, -log(-log(x)), (-log(x))^-ξ - 1)
}
