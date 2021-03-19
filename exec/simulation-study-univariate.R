library(INLA)
library(ggplot2)
library(dplyr)
library(inlaBGEV)
library(parallel)

α = .5; β = .8 # Probabilities used in the location and spread parameters
n_vec = c(5, 10, 25, 50, 100, 250, 500, 1000)
n_trials = 200
num_cores = 20
return_periods = c(5, 10, 25, 50)

inclusion_percentages = list()
for (i in seq_along(n_vec)) {
  n = n_vec[i]
  set.seed(1, kind = "L'Ecuyer-CMRG") # Set a seed that works for paralellisation
  inclusion = parallel::mclapply(
    X = seq_len(n_trials),
    mc.preschedule = FALSE,
    mc.cores = num_cores,
    FUN = function(j) {
      q = rnorm(1, 0, 5)
      s = runif(1, 0, 10)
      ξ = runif(1, 0, .5)
      locscale_pars = locspread_to_locscale(q, s, ξ, α, β)
      μ = locscale_pars$μ; σ = locscale_pars$σ
      return_levels = return_level_bgev(return_periods, μ, σ, ξ)
      y = rgev(n, μ, σ, ξ)
      res = tryCatch(inla_bgev(y, α = α, β = β), error = function(e) NULL)
      if (is.null(res)) return(NULL)
      samples = INLA::inla.posterior.sample(500, res, seed = 1)
      q_s = sapply(samples, function(x) tail(x$latent, 1)) * res$standardising_const
      s_s = sapply(samples, function(x) x$hyperpar[1]) * res$standardising_const
      ξ_s = sapply(samples, function(x) x$hyperpar[2])
      locscale_s = locspread_to_locscale(q_s, s_s, ξ_s)
      μ_s = locscale_s$μ; σ_s = locscale_s$σ
      est_return_levels = mapply(return_level_bgev, μ = μ_s, σ = σ_s, ξ = ξ_s,
                                 MoreArgs = list(period = return_periods))

      pars = c("q", "s", "ξ", paste(return_periods, "return period"))
      is_inside = NULL
      values = c(q, s, ξ, return_levels)
      for (var in c("q", "s", "ξ")) {
        interval = quantile(get(paste0(var, "_s")), c(.025, .975))
        is_inside = c(is_inside, get(var) >= interval[1] && get(var) <= interval[2])
      }
      for (i in seq_along(return_periods)) {
        interval = quantile(est_return_levels[i, ], c(.025, .975))
        is_inside = c(is_inside, return_levels[i] >= interval[1] && return_levels[i] <= interval[2])
      }
      data.frame(par = pars,
                 val = values,
                 inside = is_inside,
                 n = n, j = j)
    })
  inclusion = inclusion[!sapply(inclusion, is.null)]
  inclusion = do.call(rbind, inclusion)
  inclusion_percentages[[i]] = lapply(
    X = seq_along(unique(inclusion$par)),
    FUN = function(i) {
      par = unique(inclusion$par)[i]
      included = dplyr::filter(inclusion, par == !!par)$inside
      data.frame(par = par, n = inclusion$n[1], inclusion_percentage = mean(included))
    })
  inclusion_percentages[[i]] = do.call(rbind, inclusion_percentages[[i]])
  message("Done with n = ", n, ". Number of successful runs: ", length(unique(inlusion$j)))
  print(inclusion_percentages[[i]])
  saveRDS(
    inclusion_percentages,
    file.path(here::here(), "inst", "extdata", "simulation-study-univariate.rds"))
}

saveRDS(
  inclusion_percentages,
  file.path(here::here(), "inst", "extdata", "simulation-study-univariate.rds"))
