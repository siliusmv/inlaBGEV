test_that("pbgev is the inverse of qbgev", {
  n1 = 20
  n2 = 20
  μ = rnorm(n1)
  σ = exp(rnorm(n1))
  ξ = runif(n1)
  p = matrix(runif(n1 * n2), n1, n2)
  p2 = p
  for (i in seq_len(n2)) {
    p2[, i] = qbgev(p[, i], μ, σ, ξ) %>%
      pbgev(μ, σ, ξ)
  }
  expect_equal(p, p2, tolerance = 1e-4)
})

test_that("dbgev integrates to pbgev", {
  n = 5
  μ = rnorm(n)
  σ = exp(rnorm(n))
  ξ = runif(n)
  integrals = rep(0, n)
  p_low = runif(n, 0, .1)
  p_high = runif(n, .9, 1)
  for (i in seq_len(n)) {
    lower = qbgev(p_low[i], μ[i], σ[i], ξ[i])
    upper = qbgev(p_high[i], μ[i], σ[i], ξ[i])
    integrals[i] = integrate(function(x) dbgev(x, μ[i], σ[i], ξ[i]), lower, upper)$value
  }
  expect_equal(integrals, p_high - p_low, tolerance = 1e-3)
})

test_that("get_gumbel_par() gives a = a and b = b for frechét and gumbel", {
  n = 100
  μ = rnorm(n)
  σ = exp(rnorm(n))
  ξ = runif(n)
  p_a = runif(1, 0, .2)
  p_b = runif(1, p_a, .4)
  par = get_gumbel_par(μ, σ, ξ, p_a, p_b)
  expect_equal(qgev(p_a, μ, σ, ξ), qgev(p_a, par$μ, par$σ, 0))
  expect_equal(qgev(p_b, μ, σ, ξ), qgev(p_b, par$μ, par$σ, 0))
})

test_that("twcrps_bgev is correctly computed", {
  n = 1
  μ = rnorm(n)
  σ = exp(rnorm(n))
  ξ = runif(n)
  p = .9
  y = rgev(5, μ, σ, ξ)
  ρ_p = function(u, p) ifelse(u >= 0, p * u, (p - 1) * u)
  F_inv = function(p) qbgev(p, μ, σ, ξ)
  quantile_integrand = function(p, y) 2 * ρ_p(y - F_inv(p), p)
  score1 = sapply(y, function(y) integrate(quantile_integrand, p, 1, y = y)$value)
  score2 = twcrps_bgev(y, μ, σ, ξ, p)
  expect_equal(score1, score2, tolerance = .1 * max(score1))
})

test_that("expected_twcrps_bgev is correctly computed for n = 1", {
  n = 1
  μ = rnorm(n)
  σ = exp(rnorm(n))
  ξ = runif(n)
  p = .9
  ρ_p = function(u, p) ifelse(u >= 0, p * u, (p - 1) * u)
  F_inv = function(p) qbgev(p, μ, σ, ξ)
  quantile_integrand = function(p, y) 2 * ρ_p(y - F_inv(p), p)
  twcrps = function(y) sapply(y, function(x) integrate(quantile_integrand, p, 1, y = x)$value)
  lims = qgev(c(.00001, .99999), μ, σ, ξ)
  etwcrps1 = integrate(function(y) dbgev(y, μ, σ, ξ) * twcrps(y), lims[1], lims[2])$value
  etwcrps2 = expected_twcrps_bgev(μ, σ, ξ, p)
  expect_equal(etwcrps1, etwcrps2, tolerance = .1 * abs(etwcrps1))
})
