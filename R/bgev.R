
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
stwcrps_bgev = function(y, μ, σ, ξ, p, p_a = .1, p_b = .2, ...) {
  S = abs(expected_twcrps_bgev(μ, σ, ξ, p, p_a, p_b, ...))
  twcrps = twcrps_bgev(y, μ, σ, ξ, p, p_b)
  # We need to repeat the constants S if we have multiple observations to test against
  if (length(y) > 1) {
    # If we also have multiple parameter combinations then twcrps is a matrix
    if (!is.null(dim(twcrps))) {
      S = matrix(rep(S, each = length(y)), nrow = length(y), ncol = length(S))
    }
  }
  twcrps / S + log(S)
}

expected_twcrps_bgev = function(μ, σ, ξ, p, p_a = .1, p_b = .2, n_steps = 3) {
  #Fy = pgev(y, μ, σ, ξ)
  Ei = function(x) {
    res = rep(0, length(x))
    res[x > -700] = gsl::expint_Ei(x[x > -700])
    res
  }
  #p1 = pmax(p, Fy)
  a = qgev(p_a, μ, σ, ξ); b = qgev(p_b, μ, σ, ξ)
  gumbel_par = get_gumbel_par(μ, σ, ξ, p_a, p_b)
  μ2 = gumbel_par$μ; σ2 = gumbel_par$σ

  # Compute expected value of y * (-p^2 - 1) ===========================
  # integral over the Frechet and Gumbel domains
  gumbel_expected_value = μ2 * p_a - σ2 * p_a * log(-log(p_a)) + σ2 * Ei(log(p_a))
  frechet_expected_value = (μ - σ / ξ) * (1 - p_b) - σ / ξ * gamma(1 - ξ) *
    (pgamma(0, 1 - ξ) - pgamma(-log(p_b), 1 - ξ))
  # Approximate integral over the blended domain using linear interpolation of log(F(x))
  blended_expected_value = 0
  for (i in seq_len(n_steps)) {
    x_i = a + (i - 1) * (b - a) / n_steps
    x_iplus = a + i * (b - a) / n_steps
    h_i = log(pbgev(x_i, μ, σ, ξ, p_a, p_b))
    h_iplus = log(pbgev(x_iplus, μ, σ, ξ, p_a, p_b))
    blended_expected_value = blended_expected_value + (x_iplus - x_i) / (h_iplus - h_i) *
      (exp(h_iplus) * (h_iplus - 1) - exp(h_i) * (h_i - 1)) +
      (x_i * h_iplus - x_iplus * h_i) / (h_iplus - h_i) * (exp(h_iplus) - exp(h_i))
  }
  # Compute the expected value, and multiply with the constant
  tmp1 = (gumbel_expected_value + frechet_expected_value + blended_expected_value) * (-p^2 - 1)

  # Compute the expected value of 2y * max(p0, F(y)) ===========================
  # integral over the Frechet domain, ending at p
  frechet_integral_to_p = (μ - σ / ξ) * (p - p_b) - σ / ξ * gamma(1 - ξ) *
    (pgamma(-log(p), 1 - ξ) - pgamma(-log(p_b), 1 - ξ)) 
  # Compute the expected value, and multiply with the constant
  tmp2 = (gumbel_expected_value + blended_expected_value + frechet_integral_to_p) * 2 * p +
    (μ - σ / ξ) * (1 - p^2) +
    2^ξ * σ / ξ * gamma(1 - ξ) * (pgamma(-2 * log(p), 1 - ξ) - pgamma(0, 1 - ξ))

  # Compute the expected value of 2(-μ + σ/ξ) * max(p0, F(y)) ====================
  tmp3 = (-μ + σ / ξ) * (1 + p^2)

  # Compute the expected value of 2σ/ξ * Γ_l(-log(max(p0, F(y))), 1 - ξ) ===========
  approx_gamma_int = get_integral_of_gamma_l_with_xi(p, max(.5, ξ))
  tmp4 = 2 * (σ / ξ) * (p * gamma(1 - ξ) * pgamma(-log(p), 1 - ξ) + approx_gamma_int(ξ))

  # Compute the expected value of some constants ===========
  tmp5 = (-μ + σ / ξ) * (-p^2 - 1) - σ / ξ * 2^ξ * gamma(1 - ξ) * pgamma(-2 * log(p), 1 - ξ)

  # Sum everything and return it
  tmp1 + tmp2 + tmp3 + tmp4 + tmp5
}

get_integral_of_gamma_l_with_xi = function(p, ξ_max = .5) {
  if (ξ_max > .5) warning("I don't know if this gives a good approximation for ξ > 0.5")
  ξ_vals = seq(0, ξ_max, length = 5)
  int_vals = vector("numeric", length(ξ_vals))
  f = function(x, ξ) gamma(1 - ξ) * pgamma(-log(x), 1 - ξ)
  for (i in seq_along(ξ_vals)) {
    int_vals[i] = integrate(f, p, 1, ξ = ξ_vals[i])$value
  }
  coeffs = lm(y~poly(ξ, 4, raw = TRUE), data = data.frame(ξ = ξ_vals, y = int_vals))$coefficients
  function(ξ) {
    res = coeffs[1] + coeffs[2] * ξ + coeffs[3] * ξ^2 + coeffs[4] * ξ^3 + coeffs[5] * ξ^4
    unname(res)
  }
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
