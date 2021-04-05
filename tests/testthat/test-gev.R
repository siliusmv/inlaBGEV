test_that("pgev is the inverse of qgev", {
  n1 = 20
  n2 = 20
  μ = rnorm(n1)
  σ = exp(rnorm(n1))
  ξ = runif(n1)
  p = matrix(runif(n1 * n2), n1, n2)
  p2 = p
  for (i in seq_len(n2)) {
    p2[, i] = qgev(p[, i], μ, σ, ξ) %>%
      pgev(μ, σ, ξ)
  }
  expect_equal(p, p2, tolerance = 1e-4)
})

test_that("dgev integrates to pgev", {
  n = 5
  μ = rnorm(n)
  σ = exp(rnorm(n))
  ξ = runif(n)
  integrals = rep(0, n)
  p_low = runif(n, 0, .1)
  p_high = runif(n, .9, 1)
  for (i in seq_len(n)) {
    lower = qgev(p_low[i], μ[i], σ[i], ξ[i])
    upper = qgev(p_high[i], μ[i], σ[i], ξ[i])
    integrals[i] = integrate(function(x) dgev(x, μ[i], σ[i], ξ[i]), lower, upper)$value
  }
  expect_equal(integrals, p_high - p_low, tolerance = 1e-3)
})


test_that("reparametrisation is invertible", {
  n = 100
  μ = rnorm(n)
  σ = exp(rnorm(n))
  ξ = runif(n)
  p1 = locscale_to_locspread(μ, σ, ξ)
  p2 = locspread_to_locscale(p1$q, p1$s, p1$ξ)
  expect_identical(p2, list(μ = μ, σ = σ, ξ = ξ), tolerance = 1e-8)
  n = 100
  q = rnorm(n)
  s = exp(rnorm(n))
  ξ = runif(n)
  p1 = locspread_to_locscale(q, s, ξ)
  p2 = locscale_to_locspread(p1$μ, p1$σ, p1$ξ)
  expect_identical(p2, list(q = q, s = s, ξ = ξ), tolerance = 1e-8)
})
