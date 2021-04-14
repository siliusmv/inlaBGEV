
#' @export
pc_gp = function(ξ, λ = 7) {
  sqrt(2) * λ * exp(-sqrt(2) * λ * ξ / sqrt(1 - ξ)) * (1 - ξ / 2) / (1 - ξ)^1.5
}

#' @export
pc_gev = function(ξ, λ = 7) {
  # abs() is because of possible numerical unstability
  d = sqrt(2 * abs(kld_gev(ξ)))
  λ * exp(-λ * d) / d * abs(kld_gev_derivative(ξ))
}

#' @export
pc_bgev = function(ξ, λ = 7, p_a = .1, p_b = .2) {
  # abs() is because of possible numerical unstability
  d = sqrt(2 * abs(kld_bgev(ξ, p_a = p_a, p_b = p_b)))
  λ * exp(-λ * d) / d * abs(kld_bgev_derivative(ξ, p_a = p_a, p_b = p_b))
}

#' @export
kld_gp = function(ξ) {
  if (any(ξ < 0 | ξ >= 1)) stop("Invalid ξ")
  ξ^2 / (1 - ξ)
}

#' @export
kld_gev = function(ξ) {
  if (any(ξ < 0 | ξ >= 1)) stop("Invalid ξ")
  int = sapply(ξ, function(ξ) {
    if (ξ == 0) {
      1
    } else {
      integrate(function(x) exp((1 - (-log(x))^-ξ) / ξ), lower = 0, upper = 1)$value
    }
  })
  -1 + (1 + ξ) * digamma(1) + (gamma(1 - ξ) - 1) / ξ + int
}

#' @export
kld_bgev = function(ξ, p_a = .1, p_b = .2, ...) {
  if (any(ξ < 0 | ξ >= 1)) stop("Invalid ξ")
  res = kld_bgev_gumbel(ξ, p_a, p_b) + kld_bgev_frechet(ξ, p_a, p_b)
  res = res + sapply(seq_along(ξ), function(i) kld_bgev_mixing(ξ[i], p_a, p_b, ...))
  res
}

kld_bgev_frechet = function(ξ, p_a = .1, p_b = .2) {
  int = sapply(seq_along(ξ), function(i) {
    if (ξ[i] == 0) {
      0
    } else {
      integrate(function(x) exp((1 - (-log(x))^-ξ[i]) / ξ[i]), p_b, 1)$value
    }
  })
  - Γ_l(-log(p_b), 2) +
    (ξ + 1) * (-p_b * log(-log(p_b)) + gsl::expint_Ei(log(p_b)) + digamma(1)) +
    (Γ_l(-log(p_b), 1 - ξ) - (1 - p_b)) / ξ + int
}

kld_bgev_gumbel = function(ξ, p_a = .1, p_b = .2) {
  c1 = (1 / ξ) * ((-log(p_b))^-ξ - (-log(p_a))^-ξ) / log(log(p_a) / log(p_b))
  c2 = log(-log(p_a)) / ξ * ((-log(p_b))^-ξ - (-log(p_a))^-ξ) / log(log(p_a) / log(p_b)) -
    1 / ξ + (-log(p_a))^-ξ / ξ
  -Γ_u(-log(p_a), 2) + p_a * log(-log(p_a)) - gsl::expint_Ei(log(p_a)) +
    Γ_u(-log(p_a), c1 + 1) * exp(-c2) - c1 * (p_a * log(-log(p_a)) - gsl::expint_Ei(log(p_a))) +
      p_a * (c2 + log(ξ) + log(log(log(p_a) / log(p_b))) - log((-log(p_b))^-ξ - (-log(p_a))^-ξ))
}

Γ_l = function(x, α) pgamma(x, α) * gamma(α)
Γ_u = function(x, α) gamma(α) * (1 - pgamma(x, α))

kld_bgev_mixing = function(ξ, p_a = .1, p_b = .2, μ = 0, σ = 1) {
  if (length(ξ) > 1) stop("This can only take one ξ")
  f = function(x) {
    d = dbgev(x, μ, σ, ξ)
    l = (dbgev(x, μ, σ, ξ, log = TRUE) - dgev(x, μ, σ, 0, log = TRUE))
    ifelse(d == 0, 0, d * l)
  }
  if (ξ == 0) {
    0
  } else {
    lower = qgev(p_a, μ, σ, ξ)
    upper = qgev(p_b, μ, σ, ξ)
    integrate(f, lower = lower, upper = upper)$value
  }
}

kld_bgev_derivative = function(ξ, ...) {
  numDeriv::grad(kld_bgev, ξ, ...)
}

kld_gev_derivative = function(ξ, ...) {
  f = function(u, ξ) {
    exp((1 - (-log(u))^-ξ) / ξ) *
      (ξ^-2 * (-1 + (-log(u))^-ξ) + ξ^-1 * (-log(u))^-ξ * log(-log(u)))
  }
  int_derivative = sapply(ξ, function(ξ) integrate(f, lower = 0, upper = 1, ξ = ξ)$value)
  digamma(1) - gamma(1 - ξ) * digamma(1 - ξ) / ξ - (gamma(1 - ξ) - 1) / ξ^2 + int_derivative
}
