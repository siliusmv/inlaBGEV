library(INLA)
library(ggplot2)
library(dplyr)
library(inlaBGEV)
library(parallel)

n_trials = 20
n_loc = 250
n = 1500
α = .5
β = .8
n_trials = 30
n_standardised_samples = 12
num_cores = 12

#inla.setOption(pardiso.license = "~/.R/licences/pardiso.lic")

tail_estimates = list()
set.seed(1, kind = "L'Ecuyer-CMRG")
for (i in 1:n_trials) {
  # Use a location parameter of 0, draw random coefficients for s and a random tail parameter
  q = rep(0, n_loc)
  n_s = sample.int(4, 1)
  s_coeffs = c(runif(1) * 2 + 1, rnorm(n_s) * .4)
  ξ = runif(1, .01, .3)

  s_cov_names = paste0("s_", 1:n_s)
  names(s_coeffs) = c("intercept", s_cov_names)

  # Simulate covariates
  X = cbind(1, matrix(rnorm(n_loc * n_s), nrow = n_loc))
  colnames(X) = c("intercept", s_cov_names)

  # Simulate the number of observations at each location
  location_indices = c(1:n_loc, sample.int(n_loc, n - n_loc, replace = TRUE))

  # Calculate the spread parameter
  s = as.numeric(X %*% s_coeffs)
  s = exp(s)

  # Simulate observations from the latent field
  locscale_pars = locspread_to_locscale(q, s, ξ, α, β)
  observations = rgev(n,
                      locscale_pars$μ[location_indices],
                      locscale_pars$σ[location_indices],
                      locscale_pars$ξ[location_indices])

  # Create a data frame with covariates and response
  data = as.data.frame(X[location_indices, ])
  data$obs = observations

  # Run R-INLA
  res1 = tryCatch(
    inla_bgev(
      data = data,
      covariate_names = list(NULL, s_cov_names, NULL),
      α = α,
      β = β,
      response_name = "obs"),
    error = function(e) NULL)
  #summary(res1)

  # Run R-inla to estimate the BGEV-parameters once for each of the SD samples
  samples = parallel::mclapply(
    X = seq_len(n_standardised_samples),
    mc.cores = num_cores,
    mc.preschedule = FALSE,
    FUN = function(i) {
      noisy_s = s[location_indices] * rnorm(n, mean = 1, sd = .15)
      res = tryCatch({
        inla_bgev(
          data = data,
          s_est = noisy_s,
          covariate_names = list(NULL, NULL, NULL),
          response_name = "obs",
          α = α,
          β = β)},
        error = function(e) NULL)
      if (is.null(res) || !res$convergence) return(NULL)
      inla.posterior.sample(100, res, seed = 1) %>%
        sapply(function(x) x$hyper[2])
    })


  tail_estimates[[i]] = list(truth = ξ, n_s = n_s)
  # res1 might some times crash form numerical errors. Sometimes it also says that
  # ξ is nondegenerate in its summary, but not when we sample from it. Therefore,
  # we might have to sample ξ from the posterior of res1. Since the crashes are not
  # 'consistent', we have to preserve the random seed here to ensure as much
  # reproducibility as possible
  seed = .Random.seed
  if (!is.null(res1)) {
    tail_estimates[[i]]$joint = summary(res1)$hyper[2, -6]
    if (any(is.na(summary(res1)$hyper[2, ]))) {
      tmp = INLA::inla.posterior.sample(100, res1) %>%
        sapply(function(x) x$hyper[2]) %>%
        data_stats(q = c(.025, .5, .975))
      names(tmp) = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
      tail_estimates[[i]]$joint = tmp
    }
  }
  .Random.seed = seed

  if (any(sapply(samples, is.null))) samples = samples[-which(sapply(samples, is.null))]
  tmp = unlist(samples) %>%
    data_stats(q = c(.025, .5, .975))
  names(tmp) = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
  tail_estimates[[i]]$twostep = tmp
  message("Done with iteration nr. ", i)
}

mydf = lapply(
  seq_along(tail_estimates),
  function(i) {
    x = tail_estimates[[i]]
    if (is.null(x$joint)) {
      x$joint = data.frame(NA, NA, NA, NA, NA)
      names(x$joint) = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
    }
    if (is.null(x$twostep)) {
      x$twostep = data.frame(NA, NA, NA, NA, NA)
      names(x$twostep) = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
    }
    res = as.data.frame(rbind(x$joint, x$twostep))
    res$tag = c("Joint model", "Two-step model")[c(!is.null(x$joint), !is.null(x$twostep))]
    res$truth = x$truth
    res$x = i
    res$n_s = x$n_s
    res
  }) %>%
  do.call(rbind, .)

plot = mydf %>%
  dplyr::mutate(n_s = factor(n_s)) %>%
  ggplot(aes(x = x)) +
  geom_errorbar(aes(ymin = `0.025quant`, ymax = `0.975quant`, col = n_s)) +
  geom_point(aes(y = truth, col = n_s)) +
  #scale_shape_manual(values = c(15, 16, 17, 18)) +
  facet_wrap(~tag, nrow = 2) +
  labs(y = "$\\xi$", x = "Simulation nr.", col = "$|\\beta_s|$") +
  theme_bw() +
  #theme_linedraw() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5),
        text = element_text(size = 20))

tikz_plot(file.path(here::here(), "inst", "extdata", "simulation-study.pdf"),
          print(plot), width = 10, height = 7, view = TRUE)
