library(INLA)

n_trials = 20
n_loc = 250
n = 1500
α = .5
β = .8
n_trials = 30
n_standardised_samples = 10

tail_estimates = list()
set.seed(123)
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
      response_name = "obs"),
    error = function(e) NULL)
  #summary(res1)

  samples = list()
  for (j in 1:n_standardised_samples) {
    # Add some noise to s and standardise the observations
    noisy_s = s[location_indices] * rnorm(n, mean = 1, sd = .1)

    res2 = tryCatch(
      inla_bgev(
        data = data,
        s_est = noisy_s,
        covariate_names = list(NULL, NULL, NULL),
        response_name = "obs"),
      error = function(e) NULL)
    #summary(res2)

    # Sample ξ from the posterior of res2
    if (!is.null(res2)) {
      samples[[length(samples) + 1]] = INLA::inla.posterior.sample(100, res2) %>%
        sapply(function(x) x$hyper[2])
    }
  }

  tail_estimates[[i]] = list(truth = ξ)

  # res1 might some times crash form numerical errors. Sometimes it also says that
  # ξ is nondegenerate in its summary, but not when we sample from it. Therefore,
  # we might have to sample ξ from the posterior of res1. Since the crashes are not
  # 'consistent', we have to preserve the random seed here to ensure as much
  # reproducibility as possible
  seed = .Random.seed
  if (!is.null(res1)) {
    tail_estimates[[i]]$with_spread = summary(res1)$hyper[2, -6]
    if (any(is.na(summary(res1)$hyper[2, ]))) {
      tmp = INLA::inla.posterior.sample(100, res1) %>%
        sapply(function(x) x$hyper[2]) %>%
        data_stats(q = c(.025, .5, .975))
      names(tmp) = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
      tail_estimates[[i]]$with_spread = tmp
    }
  }
  .Random.seed = seed

  tail_estimates[[i]]$n_s = n_s

  tmp = unlist(samples) %>%
    data_stats(q = c(.025, .5, .975))
  names(tmp) = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
  tail_estimates[[i]]$standardised = tmp
  message("Done with iteration nr. ", i)
}

mydf = lapply(
  seq_along(tail_estimates),
  function(i) {
    x = tail_estimates[[i]]
    if (is.null(x$with_spread)) {
      x$with_spread = data.frame(NA, NA, NA, NA, NA)
      names(x$with_spread) = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
    }
    if (is.null(x$standardised)) {
      x$standardised = data.frame(NA, NA, NA, NA, NA)
      names(x$standardised) = c("mean", "sd", "0.025quant", "0.5quant", "0.975quant")
    }
    res = as.data.frame(rbind(x$with_spread, x$standardised))
    res$tag = c("Modelling all parameters jointly",
                "Modelling the spread separately")[
                  c(!is.null(x$with_spread), !is.null(x$standardised))]
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
plot

#tikz_plot(file.path(here::here(), "inst", "extdata", "two-step-performance.pdf"),
#          print(plot), width = 10, height = 7, view = TRUE)
