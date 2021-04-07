library(INLA)
library(ggplot2)
library(dplyr)
library(inlaBGEV)
library(parallel)


#mle_bgev = function(x, p_a = .1, p_b = .2) {
#  f = function(par, x) {
#    res = -sum(dbgev(x, par[1], par[2], par[3], p_a, p_b, log = TRUE))
#    if (is.infinite(res)) res = 99999999 * sign(res)
#    res
#  }
#  optim(c(0, 1, .1), fn = f, lower = c(-100, .00001, .00001), upper = c(100, 100, .5),
#        method = "L-BFGS-B", x = x)
#}

# set.seed(7)
# q = 2.3
# s = 1
# ξ = .23
# α = .5
# β = .8
# locscale_pars = locspread_to_locscale(q, s, ξ, α, β)
# μ = locscale_pars$μ; σ = locscale_pars$σ
# y = rgev(10000, μ, σ, ξ)
# standardising_const = diff(quantile(y, c(.05, .95)))
# y = y / standardising_const
# mdata = inla.mdata(y, matrix(0, length(y), 0), matrix(0, length(y), 0))
# 
# r = inla(
#   formula = mdata ~ -1 + intercept,
#   data = data.frame(intercept = 1),
#   control.family = list(control.bgev = list(q.location = α, q.spread = β)),
#   family = "bgev")
# summary(r)$hyper
# truth = data.frame(q = q / standardising_const, s = s / standardising_const, ξ = ξ)
# rownames(truth) = NULL
# truth

#r2 = mle_bgev(y)
#p2 = locscale_to_locspread(r2$par[1], r2$par[2], r2$par[3], α, β)
#p2


expected_twcrps_bgev_gev = function(μ, σ, ξ, p,
                                    μ_true = μ, σ_true = σ, ξ_true = ξ,
                                    p_a = .1, p_b = .2) {
  p_min = .00001
  y_min = min(qgev(p_min, μ_true, σ_true, ξ_true))
  p_max = .99999
  y_max = max(qgev(p_max, μ_true, σ_true, ξ_true))
  if (length(c(μ_true, σ_true, ξ_true)) == 3) {
    density = function(x) dgev(x, μ_true, σ_true, ξ_true)
  } else {
    density = function(x) sapply(x, function(z) mean(dgev(z, μ_true, σ_true, ξ_true)))
  }
  integrate(function(y) density(y) * twcrps_bgev(y, μ, σ, ξ, p, p_b),
            lower = y_min, upper = y_max)$value
}


α = .5; β = .8 # Probabilities used in the location and spread parameters
n_vec = c(10, 50, 100, 500, 1000, 5000, 10000)
n_trials = 500
num_cores = 25
return_periods = c(5, 10, 25, 50, 100, 200, 500)

set.seed(1, kind = "L'Ecuyer-CMRG") # Set a seed that works for paralellisation
stats = parallel::mclapply(
  X = seq_len(n_trials),
  mc.cores = num_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {
    q = rnorm(1, 0, 5)
    s = runif(1, 0, 10)
    ξ = runif(1, 0, .4)
    locscale_pars = locspread_to_locscale(q, s, ξ, α, β)
    μ = locscale_pars$μ; σ = locscale_pars$σ
    return_levels = return_level_bgev(return_periods, μ, σ, ξ)
    stats = lapply(
      n_vec,
      function(n) {
        y = rgev(n, μ, σ, ξ)
        res = tryCatch(inla_bgev(y, α = α, β = β), error = function(e) NULL)
        if (is.null(res)) return(NULL)
        samples = INLA::inla.posterior.sample(1000, res, seed = 1)
        q_s = vapply(samples, function(x) utils::tail(x$latent, 1), 1) * res$standardising_const
        s_s = vapply(samples, function(x) x$hyperpar[1], 1) * res$standardising_const
        ξ_s = vapply(samples, function(x) x$hyperpar[2], 1)
        locscale_s = locspread_to_locscale(q_s, s_s, ξ_s, α, β)
        μ_s = locscale_s$μ; σ_s = locscale_s$σ
        est_return_levels = mapply(return_level_bgev, μ = μ_s, σ = σ_s, ξ = ξ_s,
                                   MoreArgs = list(period = return_periods))
        pars = c("q", "s", "ξ", paste(return_periods, "return period"))
        is_inside = NULL; lower = NULL; upper = NULL; mean_vals = NULL; mse = NULL
        values = c(q, s, ξ, return_levels)
        for (var in c("q", "s", "ξ")) {
          interval = quantile(get(paste0(var, "_s")), c(.025, .975))
          is_inside = c(is_inside, get(var) >= interval[1] && get(var) <= interval[2])
          lower = c(lower, interval[1]); upper = c(upper, interval[2])
          mean_vals = c(mean_vals, base::mean(get(paste0(var, "_s"))))
          mse = c(mse, base::mean((get(paste0(var, "_s")) - get(var))^2))
        }
        for (j in seq_along(return_periods)) {
          interval = quantile(est_return_levels[j, ], c(.025, .975))
          is_inside = c(is_inside, return_levels[j] >= interval[1] && return_levels[j] <= interval[2])
          lower = c(lower, interval[1]); upper = c(upper, interval[2])
          mean_vals = c(mean_vals, base::mean(est_return_levels[j, ]))
          mse = c(mse, base::mean((est_return_levels[j, ] - return_levels[j])^2))
        }
        stwcrps = base::mean(stwcrps_bgev(y, μ_s, σ_s, ξ_s, .9))
        etwcrps = expected_twcrps_bgev_gev(μ_s, σ_s, ξ_s, .9, μ_true = μ, σ_true = σ, ξ_true = ξ)
        S = abs(expected_twcrps_bgev(μ_s, σ_s, ξ_s, .9))
        estwcrps = etwcrps / S + log(S)
        #stwcrps_true = mean(stwcrps_bgev(y, μ, σ, ξ, .9))
        data.frame(
          par = pars,
          val = values,
          is_inside = is_inside,
          n = n, i = i,
          lower = lower,
          mean = mean_vals,
          mse = mse,
          stwcrps = stwcrps,
          etwcrps = etwcrps,
          estwcrps = estwcrps,
          #stwcrps_true = mean(stwcrps_true),
          upper = upper,
          CI_width = upper - lower,
          stringsAsFactors = FALSE)
      })
    stats = do.call(rbind, stats)
    message("Done with iter ", i, " out of ", n_trials)
    stats
  })

saveRDS(
  stats,
  file.path(here::here(), "results", "simulation-study-univariate.rds"))

stats %>%
  do.call(rbind, .) %>%
  dplyr::filter(n == tail(n_vec, 1), par == "s") %>%
  dplyr::mutate(below = case_when(val < lower ~ TRUE, val > upper ~ FALSE, TRUE ~ NA)) %>%
  ggplot() +
  geom_point(aes(x = i, y = val, col = below)) +
  geom_errorbar(aes(x = i, ymin = lower, ymax = upper, col = below)) +
  facet_wrap(~is_inside)

stats %>%
  do.call(rbind, .) %>%
  dplyr::filter(n == tail(n_vec, 1), par == "ξ") %>%
  dplyr::mutate(below = case_when(val < lower ~ TRUE, val > upper ~ FALSE, TRUE ~ NA)) %>%
  ggplot() +
  geom_point(aes(x = i, y = val, col = below)) +
  geom_errorbar(aes(x = i, ymin = lower, ymax = upper, col = below)) +
  facet_wrap(~is_inside)

stats %>%
  do.call(rbind, .) %>%
  dplyr::filter(n == tail(n_vec, 1), par == "s") %>%
  dplyr::mutate(below = case_when(val < lower ~ TRUE, val > upper ~ FALSE, TRUE ~ NA)) %>%
  ggplot() +
  geom_point(aes(x = i, y = val, col = below)) +
  geom_errorbar(aes(x = i, ymin = lower, ymax = upper, col = below)) +
  facet_wrap(~is_inside)

stats %>%
  do.call(rbind, .) %>%
  dplyr::filter(n == tail(n_vec, 1), par %in% c("s", "ξ"), !is_inside) %>%
  dplyr::mutate(below = case_when(val < lower ~ TRUE, val > upper ~ FALSE, TRUE ~ NA)) %>%
  ggplot() +
  geom_point(aes(x = i, y = val, col = below)) +
  geom_errorbar(aes(x = i, ymin = lower, ymax = upper, col = below)) +
  facet_wrap(~par, scales = "free_y")

stats %>%
  do.call(rbind, .) %>%
  dplyr::group_by(par, n) %>%
  dplyr::summarise(inclusion_percentage = base::mean(is_inside)) %>%
  ggplot() +
  geom_col(aes(x = as.numeric(factor(n)), y = inclusion_percentage, fill = factor(n))) +
  facet_wrap(~par) +
  geom_hline(yintercept = .95) +
  labs(fill = "n")
# We care about the return periods, and we see that they get better and better!
# We don't really care about the parameters!

stats %>%
  do.call(rbind, .) %>%
  ggplot() +
  #geom_point(aes(x = as.numeric(factor(n)), y = mse, col = factor(n))) +
  geom_jitter(aes(x = as.numeric(factor(n)), y = mse, col = factor(n)), height = 0, width = .1) +
  facet_wrap(~par, scales = "free_y") +
  labs(col = "n")
# MSE goes down with n, which is great!

stats %>%
  do.call(rbind, .) %>%
  dplyr::group_by(n) %>%
  dplyr::summarise(
    stwcrps = base::mean(stwcrps),
    estwcrps = base::mean(estwcrps),
    etwcrps = base::mean(etwcrps))
# Everything goes down, except the stwcrps. This kind of makes sense, I think.
# The reason is that we are "overfitting" to the available data or something like that...
