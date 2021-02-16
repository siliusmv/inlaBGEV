
#' @export
pc_gp = function(ξ, λ = 7) {
  sqrt(2) * λ * exp(-sqrt(2) * λ * ξ / sqrt(1 - ξ)) * (1 - ξ / 2) / (1 - ξ)^1.5
}

#' @export
pc_gev = function(ξ, λ = 7, approx = TRUE) {
  # abs() is because of possible numerical unstability if approx = FALSE
  d = sqrt(2 * abs(kld_gev(ξ, approx = approx)))
  λ * exp(-λ * d) / d * abs(kld_gev_derivative(ξ, approx = approx))
}

#' @export
pc_bgev = function(ξ, λ = 7, p_a = .1, p_b = .2, approx = TRUE) {
  # abs() is because of possible numerical unstability if approx = FALSE
  d = sqrt(2 * abs(kld_bgev(ξ, p_a = p_a, p_b = p_b, approx = approx)))
  λ * exp(-λ * d) / d * abs(kld_bgev_derivative(ξ, p_a = p_a, p_b = p_b, approx = approx))
}

#' @export
kld_gp = function(ξ) {
  if (any(ξ < 0 | ξ >= 1)) stop("Invalid ξ")
  ξ^2 / (1 - ξ)
}

#' @export
kld_gev = function(ξ, approx = TRUE) {
  if (any(ξ < 0 | ξ >= 1)) stop("Invalid ξ")
  if (approx) {
    # These uncommented lines show how we estimate the coefficients below
    # fff = function(x, ξ) exp((1 - (-log(x))^-ξ) / ξ)
    # ff = function(x) sapply(x, function(x) integrate(fff, 0, 1, ξ = x)$value)
    # #x = seq(.001, 1, by = .001)
    # x = c(seq(.0001, .1, by = .0001), seq(.15, 1, by = .05))
    # df = data.frame(x = x, f = ff(x) - 1)
    # coeffs = lm(f ~ -1 + poly(x, 3, raw = TRUE), data = df)$coefficients
    # coeffs
    a = -.0870732
    b = .2553301
    c = -.4088320
    d = 1
    int = a * ξ^3 + b * ξ^2 + c * ξ + d
  } else {
    int = sapply(ξ, function(ξ) {
      if (ξ == 0) {
        1
      } else {
        integrate(function(x) exp((1 - (-log(x))^-ξ) / ξ), lower = 0, upper = 1)$value
      }
    })
  }
  -1 + (1 + ξ) * digamma(1) + (gamma(1 - ξ) - 1) / ξ + int
}

#' @export
kld_bgev = function(ξ, p_a = .1, p_b = .2, approx = TRUE, ...) {
  if (any(ξ < 0 | ξ >= 1)) stop("Invalid ξ")
  res = kld_bgev_gumbel(ξ, p_a, p_b) + kld_bgev_frechet(ξ, p_a, p_b, approx = approx)
  if (approx) {
    res = res + kld_bgev_mixing_approx(ξ, p_a, p_b)
  } else {
    res = res + sapply(seq_along(ξ), function(i) kld_bgev_mixing(ξ[i], p_a, p_b, ...))
  }
  res
}

kld_bgev_frechet = function(ξ, p_a = .1, p_b = .2, approx = TRUE) {
  if (approx) {
    d = 1 + p_b * (log(p_b) - 1)
    if (p_b == .2) {
      a = -.007052692
      b = .035104600
      c = -.103604470
    } else {
      fff = function(x, ξ) exp((1 - (-log(x))^-ξ) / ξ)
      ff = function(x) sapply(x, function(x) integrate(fff, p_b, 1, ξ = x)$value)
      x = seq(.001, 1, by = .001)
      df = data.frame(x = x, f = ff(x) - d)
      coeffs = lm(f ~ -1 + poly(x, 3, raw = TRUE), data = df)$coefficients
      a = coeffs[3]
      b = coeffs[2]
      c = coeffs[1]
    }
    int = a * ξ^3 + b * ξ^2 + c * ξ + d
  } else {
    int = sapply(seq_along(ξ), function(i) {
      if (ξ[i] == 0) {
        0
      } else {
        integrate(function(x) exp((1 - (-log(x))^-ξ[i]) / ξ[i]), p_b, 1)$value
      }
    })
  }
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

kld_bgev_mixing_approx = function(ξ, p_a = .1, p_b = .2) {
  f1 = kld_bgev_mixing(.5, p_a, p_b)
  f2 = kld_bgev_mixing(1, p_a, p_b)
  b = (f2 - (1 / .5)^2 * f1) * (.5 / 1) / (.5 - 1)
  a = f1 / .5^2 - b / .5
  a * ξ^2 + b * ξ
}

kld_bgev_derivative = function(ξ, approx = TRUE, ...) {
  numDeriv::grad(kld_bgev, ξ, approx = approx, ...)
}

kld_gev_derivative = function(ξ, approx = TRUE, ...) {
  if (approx) {
    a = -.0870732
    b = .2553301
    c = -.4088320
    int_derivative = 3 * a * ξ^2 + 2 * b * ξ + c
  } else {
    f = function(u, ξ) {
      exp((1 - (-log(u))^-ξ) / ξ) *
        (ξ^-2 * (-1 + (-log(u))^-ξ) + ξ^-1 * (-log(u))^-ξ * log(-log(u)))
    }
    int_derivative = sapply(ξ, function(ξ) integrate(f, lower = 0, upper = 1, ξ = ξ)$value)
  }
  digamma(1) - gamma(1 - ξ) * digamma(1 - ξ) / ξ - (gamma(1 - ξ) - 1) / ξ^2 + int_derivative
}
