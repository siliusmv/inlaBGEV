
#' @export
pbgev = function(x, μ, σ, ξ, p_a = .1, p_b = .2, s = 5) {
  fix_lengths(x, μ, σ, ξ, p_a, p_b, s)
  g = get_gumbel_par(μ, σ, ξ, p_a, p_b)
  a = qgev(p_a, μ, σ, ξ)
  b = qgev(p_b, μ, σ, ξ)
  p = pbeta((x - a) / (b - a), s, s)
  pgev(x, μ, σ, ξ) ^ p * pgev(x, g$μ, g$σ, 0) ^ (1 - p)
}

#' @export
qbgev = function(p, μ, σ, ξ, p_a = .1, p_b = .2, s = 5) {
  fix_lengths(p, μ, σ, ξ, p_a, p_b, s)
  res = vector("numeric", length(p))
  gumbel = which(p <= p_a)
  frechet = which(p >= p_b)
  mixing = setdiff(seq_along(p), c(gumbel, frechet))
  if (any(gumbel)) {
    g = get_gumbel_par(μ[gumbel], σ[gumbel], ξ[gumbel], p_a[gumbel], p_b[gumbel])
    res[gumbel] = qgev(p[gumbel], g$μ, g$σ, 0)
  }
  if (any(frechet)) res[frechet] = qgev(p[frechet], μ[frechet], σ[frechet], ξ[frechet])
  if (any(mixing)) {
    res[mixing] = qbgev_mixing(p[mixing], μ[mixing], σ[mixing], ξ[mixing],
                               p_a[mixing], p_b[mixing], s[mixing])
  }
  res
}

#' @export
rbgev = function(n, μ, σ, ξ, p_a = .1, p_b = .2, s = 5) {
  lengths = sapply(list(μ, σ, ξ), length)
  if (any(lengths > n)) stop("Bad input lengths")
  qbgev(runif(n), μ, σ, ξ, p_a, p_b, s)
}

#' @export
dbgev = function(x, μ, σ, ξ, p_a = .1, p_b = .2, s = 5, log = FALSE) {
  fix_lengths(x, μ, σ, ξ, p_a, p_b, s)
  a = qgev(p_a, μ, σ, ξ)
  b = qgev(p_b, μ, σ, ξ)
  res = vector("numeric", length(x))
  gumbel = which(x <= a)
  frechet = which(x >= b)
  mixing = setdiff(seq_along(x), c(gumbel, frechet))
  if (any(gumbel)) {
    g = get_gumbel_par(μ[gumbel], σ[gumbel], ξ[gumbel], p_a[gumbel], p_b[gumbel])
    res[gumbel] = dgev(x[gumbel], g$μ, g$σ, 0, log = log)
  }
  if (any(frechet)) res[frechet] = dgev(x[frechet], μ[frechet], σ[frechet], ξ[frechet], log = log)
  if (any(mixing)) {
    res[mixing] = dbgev_mixing(x[mixing], μ[mixing], σ[mixing], ξ[mixing],
                               p_a[mixing], p_b[mixing], s[mixing], log = log)
  }
  res
}

#' @export
twcrps_bgev = function(y, μ, σ, ξ, p, p_b = .2) {
  if (p < p_b) stop("twcrps_bgev is not implemented for p < p_b")
  twcrps_gev(y, μ, σ, ξ, p)
}

#' @export
return_level_bgev = function(period, μ, σ, ξ, p_a = .1, p_b = .2, s = 5) {
  if (any(period <= 1)) warning("invalid period")
  p = ifelse(period > 1, 1 - 1 / period, NA)
  qbgev(p, μ, σ, ξ, p_a, p_b, s)
}

dbgev_mixing = function(x, μ, σ, ξ, p_a = .1, p_b = .2, s = 5, log = FALSE) {
  g = get_gumbel_par(μ, σ, ξ, p_a, p_b)
  a = qgev(p_a, μ, σ, ξ)
  b = qgev(p_b, μ, σ, ξ)
  if (any(x <= a | x >= b)) stop("x is outside the domain for mixing")
  p = pbeta((x - a) / (b - a), s, s)
  p_der = dbeta((x - a) / (b - a), s, s) / (b - a)
  term1 = - p_der * (1 + ξ * (x - μ) / σ) ^ (-1 / ξ)
  term2 = p / σ * (1 + ξ * (x - μ) / σ) ^ (-1 / ξ - 1)
  term3 = p_der * exp(- (x - g$μ) / g$σ)
  term4 = (1 - p) / g$σ * exp(- (x - g$μ) / g$σ)
  term0 = p * log(pgev(x, μ, σ, ξ)) + (1 - p) * log(pgev(x, g$μ, g$σ, 0))
  res = term0 + log(term1 + term2 + term3 + term4)
  if (!log) res = exp(res)
  res
}

qbgev_mixing = function(p, μ, σ, ξ, p_a = .1, p_b = .2, s = 5, lower = 0, upper = 100) {
  if (any(p <= p_a | p >= p_b)) stop("p is outside the domain for mixing")
  res = vector("numeric", length(p))
  for (i in seq_along(p)) {
    f = function(x) (pbgev(x, μ, σ, ξ, p_a, p_b, s) - p)[i]
    sol = uniroot(f, lower = lower, upper = upper, extendInt = "upX")
    res[i] = sol$root
  }
  res
}

get_gumbel_par = function(μ, σ, ξ, p_a = .1, p_b = .2) {
  if (any(ξ < 0)) stop("ξ must be nonnegative")
  a = qgev(p_a, μ, σ, ξ)
  b = qgev(p_b, μ, σ, ξ)
  σ2 = (b - a) / log(log(p_a) / log(p_b))
  μ2 = a + σ2 * log(-log(p_a))
  list(μ = μ2, σ = σ2)
}
