## =============================================================================
## Save data for the INLA tests
## =============================================================================
#n = 100
#μ = 1; σ = 1; ξ = .3
#set.seed(1)
#y = rgev(n, μ, σ, ξ)
#res = inla_bgev(y)
#samples = INLA::inla.posterior.sample(5, res, seed = 1)
#saveRDS(samples, file.path(here::here(), "tests", "testdata", "inla-bgev-univariate.rds"))
#n = 200
#p = 2
#set.seed(1)
#x = matrix(rnorm(n * p), n, p)
#μ = as.numeric(2 + x %*% c(-1, 1))
#σ = as.numeric(exp(2 + x %*% c(-.3, .3)))
#ξ = .2
#y = rgev(n, μ, σ, ξ)
#df = as.data.frame(cbind(x, y))
#names(df) = c("x1", "x2", "y")
#res = inla_bgev(
#  data = df,
#  response_name = "y",
#  covariate_names = list(c("x1", "x2"), c("x1", "x2"), NULL))
#samples = INLA::inla.posterior.sample(5, res, seed = 1)
#saveRDS(samples, file.path(here::here(), "tests", "testdata", "inla-bgev-with-covariates.rds"))

test_that("extracting coefficients from a bGEV model works for a univariate model", {
  samples = readRDS(file.path(here::here(), "tests", "testdata", "inla-bgev-univariate.rds"))
  coeffs = inla_bgev_coeffs(samples)
  expected_result = list(
    q = sapply(samples, function(x) x$latent[102]),
    s = sapply(samples, function(x) x$hyperpar[1]),
    ξ = sapply(samples, function(x) x$hyperpar[2]))
  expect_equal(coeffs, expected_result)
})

test_that("extracting coefficients from a bGEV model works for when we have covariates", {
  samples = readRDS(file.path(here::here(), "tests", "testdata", "inla-bgev-with-covariates.rds"))
  coeffs = inla_bgev_coeffs(samples)
  expected_result = list(
    q = sapply(samples, function(x) x$latent[401:403]),
    s = sapply(samples, function(x) x$hyperpar[c(1, 3, 4)]),
    ξ = sapply(samples, function(x) x$hyperpar[2]))
  rownames(expected_result$q) = c("intercept", "x1", "x2")
  expect_equal(coeffs, expected_result)
})

test_that("parameters are correctly estimated with a bGEV model with no covariates", {
  samples = readRDS(file.path(here::here(), "tests", "testdata", "inla-bgev-univariate.rds"))
  coeffs = inla_bgev_coeffs(samples)
  pars = inla_bgev_pars(coeffs)
  expect_equal(pars, coeffs)
})

test_that("parameters are correctly estimated with a bGEV model with covariates", {
  samples = readRDS(file.path(here::here(), "tests", "testdata", "inla-bgev-with-covariates.rds"))
  coeffs = inla_bgev_coeffs(samples)
  X = matrix(rep(0:2, 2), 3, 2)
  colnames(X) = c("x1", "x2")
  p1 = inla_bgev_pars(coeffs, X, c("x1", "x2"))
  # Expect that a change of one covariate is the same each time
  expect_equal(p1$q[3, ] - p1$q[2, ], p1$q[2, ] - p1$q[1, ])
  expect_equal(log(p1$s[3, ]) - log(p1$s[2, ]), log(p1$s[2, ]) - log(p1$s[1, ]))
  # Is q correct?
  q = cbind(1, X) %*% coeffs$q
  expect_equal(p1$q, q)
  # Is s correct?
  s_coeffs = coeffs$s
  s_coeffs[1, ] = log(s_coeffs[1, ])
  s = exp(cbind(1, X) %*% s_coeffs)
  expect_equal(p1$s, s)
  # Is ξ correct?
  expect_equal(p1$ξ, coeffs$ξ)
})

test_that("add_spde_samples works as it should", {
  stop("This is not implemented yet")
})
