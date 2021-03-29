library(INLA)
library(ggplot2)
library(dplyr)
library(inlaBGEV)
library(parallel)

α = .5; β = .8 # Probabilities used in the location and spread parameters
n_vec = c(50, 100, 500, 1000, 2000)
n_trials = 2000
num_cores = 25
return_periods = c(5, 10, 25, 50)

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
    stats = list()
    for (n in n_vec) {
      y = rgev(n, μ, σ, ξ)
      res = tryCatch(inla_bgev(y, α = α, β = β), error = function(e) NULL)
      if (is.null(res)) return(NULL)
      samples = INLA::inla.posterior.sample(500, res, seed = 1)
      q_s = sapply(samples, function(x) tail(x$latent, 1)) * res$standardising_const
      s_s = sapply(samples, function(x) x$hyperpar[1]) * res$standardising_const
      ξ_s = sapply(samples, function(x) x$hyperpar[2])
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
        mean_vals = c(mean_vals, mean(get(paste0(var, "_s"))))
        mse = c(mse, mean((get(paste0(var, "_s")) - get(var))^2))
      }
      for (j in seq_along(return_periods)) {
        interval = quantile(est_return_levels[j, ], c(.025, .975))
        is_inside = c(is_inside, return_levels[j] >= interval[1] && return_levels[j] <= interval[2])
        lower = c(lower, interval[1]); upper = c(upper, interval[2])
        mean_vals = c(mean_vals, mean(est_return_levels[j, ]))
        mse = c(mse, mean((est_return_levels[j, ] - return_levels[j])^2))
      }
      #stwcrps = mean(stwcrps_bgev(y, μ_s, σ_s, ξ_s, .9))
      stwcrps = expected_twcrps_bgev(μ_s, σ_s, ξ_s, .9, μ_true = μ, σ_true = σ, ξ_true = ξ)
      #stwcrps_true = mean(stwcrps_bgev(y, μ, σ, ξ, .9))
      stats[[n]] = data.frame(
        par = pars,
        val = values,
        is_inside = is_inside,
        n = n, i = i,
        lower = lower,
        mean = mean_vals,
        mse = mse,
        stwcrps = stwcrps,
        #stwcrps_true = mean(stwcrps_true),
        upper = upper,
        CI_width = upper - lower)
    }
    stats = do.call(rbind, stats)
    message("Done with iter ", i, " out of ", n_trials)
    stats
  })

saveRDS(
  stats,
  file.path(here::here(), "inst", "extdata", "simulation-study-univariate.rds"))


stats %>%
  do.call(rbind, .) %>%
  dplyr::group_by(par, n) %>%
  dplyr::summarise(inclusion_percentage = mean(is_inside)) %>%
  ggplot() +
  geom_col(aes(x = as.numeric(factor(n)), y = inclusion_percentage, fill = factor(n))) +
  facet_wrap(~par) +
  geom_hline(yintercept = .95) +
  labs(fill = "n")

stats %>%
  do.call(rbind, .) %>%
  ggplot() +
  #geom_point(aes(x = as.numeric(factor(n)), y = mse, col = factor(n))) +
  geom_jitter(aes(x = as.numeric(factor(n)), y = mse, col = factor(n)), height = 0, width = .1) +
  facet_wrap(~par, scales = "free_y") +
  labs(col = "n")


stats %>%
  do.call(rbind, .) %>%
  #ggplot() +
  #geom_boxplot(aes(x = as.numeric(factor(n)), y = stwcrps, fill = factor(n))) +
  #facet_wrap(~par, scales = "free_y") +
  #labs(fill = "n")
  dplyr::group_by(par, n) %>%
  dplyr::summarise(mean = mean(stwcrps),
                   upper = quantile(stwcrps, .95),
                   lower = quantile(stwcrps, .05)) %>%
  ggplot() +
  geom_point(aes(x = n, y = mean)) +
  geom_ribbon(aes(x = n, ymin = lower, ymax = upper), alpha = .3) +
  facet_wrap(~par)


#stop("You should do this the other way, so we use the same q, s, ξ for all values of n. Then we can actually see that the error goes down as n increases!!!")
#inclusion_percentages = list()
#for (i in seq_along(n_vec)) {
#  n = n_vec[i]
#  set.seed(1, kind = "L'Ecuyer-CMRG") # Set a seed that works for paralellisation
#  inclusion = parallel::mclapply(
#    X = seq_len(n_trials),
#    mc.preschedule = FALSE,
#    mc.cores = num_cores,
#    FUN = function(j) {
#    })
#  inclusion_percentages[[i]] = do.call(rbind, inclusion_percentages[[i]])
#  message("Done with n = ", n, ". Number of successful runs: ", length(unique(inclusion$j)))
#  print(inclusion_percentages[[i]])
#  saveRDS(
#    inclusion_percentages,
#    file.path(here::here(), "inst", "extdata", "simulation-study-univariate.rds"))
#}
#
#saveRDS(
#  inclusion_percentages,
#  file.path(here::here(), "inst", "extdata", "simulation-study-univariate.rds"))
#
#inclusion_percentages %>%
#  do.call(rbind, .) %>%
#  dplyr::select(par, n, inclusion_percentage) %>%
#  tidyr::pivot_wider(names_from = "n", values_from = "inclusion_percentage") %>%
#  xtable::xtable()
