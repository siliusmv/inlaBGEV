
# This script contains functions for performing univariate parameter estimation with the GEV
# and bGEV distributions


#' @export
gev_fit = function(x, method = c("pwm", "inla_bgev", "inla_gev"), ...) {
  fit = switch(method[1],
               pwm = gev_pwm_fit(x, ...),
               inla_bgev = bgev_inla_fit(x, ...),
               inla_gev = gev_inla_fit(x, ...),
               stop("This is not a valid method"))
  class(fit) = "gev-fit"
  fit$data = x
  fit
}

bgev_inla_fit = function(x, α = .5, β = .8, ...) {
  mdata = INLA::inla.mdata(x, matrix(0, length(x), 0), matrix(0, length(x), 0))
  inla_data = as.list(data.frame(x = x, intercept = 1))
  inla_data$mdata = mdata
  estimate = tryCatch({
    INLA::inla(formula = mdata ~ -1 + intercept,
         family = "bgev",
         data = inla_data,
         control.family = list(control.bgev = list(q.location = α, q.spread = β)),
         control.fixed = list(mean = list(intercept = quantile(x, α)),
                              prec = list(intercept = 20)),
         control.compute = list(config = TRUE),
         num.threads = 1,
         control.inla = list(tolerance = 1e-12, h = 1e-3),
         ...)
  },
  error = function(e) NULL)
  if (is.null(estimate)) {
    res = list(par = c(q = NA, s = NA, ξ = NA))
  } else {
    q = estimate$summary.fixed[, -c(6, 7)]
    s = estimate$summary.hyper[1, -6]
    ξ = estimate$summary.hyper[2, -6]
    res = list(
      estimate = c(q = q$mean, s = s$mean, ξ = ξ$mean),
      stats = rbind(q, s, ξ))
    res$stats$par = c("q", "s", "ξ")
    res$stats$sd = NULL
  }
  res$α = α
  res$β = β
  res
}

gev_inla_fit = function(x, ...) {
  inla_data = data.frame(x = x, intercept = 1)
  estimate = tryCatch({
    INLA::inla(formula = x ~ -1 + intercept,
         family = "gev",
         data = inla_data,
         num.threads = 1,
         control.inla = list(tolerance = 1e-12, h = 1e-3),
         ...)
  },
  error = function(e) NULL)
  if (is.null(estimate)) {
    res = list(par = c(μ = NA, σ = NA, ξ = NA))
  } else {
    μ = estimate$summary.fixed[, -7]
    σ = 1 / sqrt(estimate$summary.hyper[1, ])
    ξ = estimate$summary.hyper[2, ]
    res = list(
      estimate = c(μ = μ$mean, σ = σ$mean, ξ = ξ$mean),
      stats = rbind(μ, σ, ξ))
    res$stats$sd = NULL
    res$stats$par = c("μ", "σ", "ξ")
  }
  res
}

gev_pwm_fit = function(x, α = .5, β = .8, use_locspread = TRUE) {
  x = sort(x)
  n = length(x)
  j_vec = seq_along(x)
  θ = rep(0, 3)
  θ[1] = mean(x)
  θ[2] = 1 / (n * (n - 1)) * sum((j_vec - 1) * x)
  θ[3] = 1 / (n * (n - 1) * (n - 2)) * sum((j_vec - 1) * (j_vec - 2) * x)
  const = (3 * θ[3] - θ[1]) / (2 * θ[2] - θ[1])
  f = function(x) (1 - 3^x) / (1 - 2^x) - const
  sol = tryCatch(uniroot(f, c(-3, 1 - 1e-10)), error = function(e) NULL)
  if (is.null(sol)) {
    if (use_locspread) {
      res = list(estimate = c(q = NA, s = NA, ξ = NA))
    } else {
      res = list(estimate = c(μ = NA, σ = NA, ξ = NA))
    }
  } else {
    ξ = sol$root
    σ = ξ * (θ[1] - 2 * θ[2]) / (gamma(1 - ξ) * (1 - 2^ξ))
    μ = θ[1] + (σ / ξ) * (1 - gamma(1 - ξ))
    if (use_locspread) {
      pars = locscale_to_locspread(μ, σ, ξ, α, β)
      res = list(estimate = c(q = pars$q, s = pars$s, ξ = pars$ξ))
    } else {
      res = list(estimate = c(μ = μ, σ = σ, ξ = ξ))
    }
  }
  res
}
