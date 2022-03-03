library(INLA)
library(ggplot2)
library(dplyr)
library(inlaBGEV)
library(parallel)
library(evd)

α = .5; β = .8 # Probabilities used in the location and spread parameters
q = 10.8
s = 1.8
ξ = .18
return_periods = c(25, 50, 100, 250, 500)
n_vec = c(25, 50, 100, 500, 1000)
n_trials = 500
num_cores = 4
locscale_pars = locspread_to_locscale(q, s, ξ, α, β)
μ = locscale_pars$μ; σ = locscale_pars$σ
return_levels = return_level_gev(return_periods, μ, σ, ξ)

plots = list()

set.seed(1, kind = "L'Ecuyer-CMRG") # Set a seed that works for paralellisation
res1 = parallel::mclapply(
  X = seq_len(n_trials),
  mc.cores = num_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {
    res = lapply(
      n_vec,
      function(n) {
        y = rgev(n, μ, σ, ξ)
        res_bgev = optim(
          par = c(μ, log(σ), log(ξ)),
          fn = function(θ, ...) {
            res = sum(dbgev(μ = rep(θ[1], n), σ = rep(exp(θ[2]), n), ξ = rep(exp(θ[3]), n), ...))
            if (any(θ[-1] < -20)) res = -1e299
            if (is.na(res)) res = -1e299
            if (is.infinite(res)) {
              res = 1e299 * sign(res)
            }
            res
          },
          x = y,
          log = TRUE,
          control = list(fnscale = -1))
        res_gev = optim(
          par = c(μ, log(σ), log(ξ)),
          fn = function(θ, ...) {
            res = sum(inlaBGEV::dgev(μ = rep(θ[1], n), σ = rep(exp(θ[2]), n), ξ = rep(exp(θ[3]), n), ...))
            if (is.infinite(res)) {
              res = 1e99 * sign(res)
            }
            res
          },
          x = y,
          log = TRUE,
          control = list(fnscale = -1))
        est_pars = exp(c(res_bgev$par, res_gev$par))
        est_pars[c(1, 4)] = log(est_pars[c(1, 4)])
        gev_rt = return_level_gev(return_periods, est_pars[4], est_pars[5], est_pars[6])
        bgev_rt = return_level_bgev(return_periods, est_pars[1], est_pars[2], est_pars[3])

        data.frame(
          value = c(gev_rt, bgev_rt),
          truth = rep(return_levels, 2),
          model = rep(c("GEV", "bGEV"), each = length(return_periods)),
          n = n,
          iter = i,
          name = paste0("R", return_periods))
      })
    message(i, " / ", n_trials)
    do.call(rbind, res)
  })
res1 = do.call(rbind, res1)

res1$tag = paste0(res1$name, ", n=", res1$n)
res1$tag = factor(res1$tag, levels = unique(res1$tag))

plots[[1]] = ggplot(res1) +
  geom_density(aes(x = value, col = model), size = 1) +
  geom_vline(aes(xintercept = truth)) +
  labs(col = "Model", y = "Density", x = "Return level") +
  facet_wrap(~tag, scales = "free", dir = "v", ncol = length(n_vec)) +
  theme(text = element_text(size = 15))

α = .5; β = .8 # Probabilities used in the location and spread parameters
σ = .5
ξ = .18
block_size = 24 * 365
block_size = 365
return_periods = c(25, 50, 100, 250, 500)
n_vec = c(25, 50, 100, 500, 1000)
n_trials = 500
num_cores = 4
return_levels = evd::qgpd(1 - 1 / (block_size * return_periods), scale = σ, shape = ξ)

set.seed(1, kind = "L'Ecuyer-CMRG") # Set a seed that works for paralellisation
res2 = parallel::mclapply(
  X = seq_len(n_trials),
  mc.cores = num_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {
    res = lapply(
      n_vec,
      function(n) {
        y = sapply(1:n, function(i) max(evd::rgpd(block_size, scale = σ, shape = ξ)))
        res_bgev = optim(
          par = c(0, 0, log(ξ)),
          fn = function(θ, ...) {
            res = sum(inlaBGEV::dbgev(μ = rep(θ[1], n), σ = rep(exp(θ[2]), n), ξ = rep(exp(θ[3]), n), ...))
            if (any(θ[-1] < -20)) res = -1e299
            if (is.na(res)) res = -1e299
            if (is.infinite(res)) {
              res = 1e299 * sign(res)
            }
            res
          },
          x = y,
          log = TRUE,
          control = list(fnscale = -1))
        res_gev = optim(
          par = c(0, 0, log(ξ)),
          fn = function(θ, ...) {
            res = sum(inlaBGEV::dgev(μ = rep(θ[1], n), σ = rep(exp(θ[2]), n), ξ = rep(exp(θ[3]), n), ...))
            if (is.infinite(res)) {
              res = 1e99 * sign(res)
            }
            res
          },
          x = y,
          log = TRUE,
          control = list(fnscale = -1))
        est_pars = exp(c(res_bgev$par, res_gev$par))
        est_pars[c(1, 4)] = log(est_pars[c(1, 4)])
        gev_rt = return_level_gev(return_periods, est_pars[4], est_pars[5], est_pars[6])
        bgev_rt = return_level_bgev(return_periods, est_pars[1], est_pars[2], est_pars[3])

        data.frame(
          value = c(gev_rt, bgev_rt),
          truth = rep(return_levels, 2),
          model = rep(c("GEV", "bGEV"), each = length(return_periods)),
          n = n,
          iter = i,
          name = paste0("R", return_periods))
      })
    message(i, " / ", n_trials)
    do.call(rbind, res)
  })
res2 = do.call(rbind, res2)

res2$tag = paste0(res2$name, ", n=", res2$n)
res2$tag = factor(res2$tag, levels = unique(res2$tag))

plots[[2]] = ggplot(res2) +
  geom_density(aes(x = value, col = model), size = 1) +
  geom_vline(aes(xintercept = truth)) +
  labs(col = "Model", y = "Density", x = "Return level") +
  facet_wrap(~tag, scales = "free", dir = "v", ncol = length(n_vec)) +
  theme(text = element_text(size = 15))


set.seed(1, kind = "L'Ecuyer-CMRG") # Set a seed that works for paralellisation
res3 = parallel::mclapply(
  X = seq_len(n_trials),
  mc.cores = num_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {
    res = lapply(
      n_vec,
      function(n) {
        y = sapply(1:n, function(i) max(evd::rgpd(block_size, scale = σ, shape = ξ)))
        res_bgev = optim(
          par = c(μ, log(σ), log(ξ)),
          fn = function(θ, ...) {
            res = sum(inlaBGEV::dbgev(μ = rep(θ[1], n), σ = rep(exp(θ[2]), n), ξ = rep(exp(θ[3]), n), ...))
            if (any(θ[-1] < -20)) res = -1e299
            if (is.na(res)) res = -1e299
            if (is.infinite(res)) {
              res = 1e299 * sign(res)
            }
            res
          },
          x = y,
          log = TRUE,
          control = list(fnscale = -1))
        res_gev = optim(
          par = c(μ, log(σ), log(ξ)),
          fn = function(θ, ...) {
            res = sum(inlaBGEV::dgev(μ = rep(θ[1], n), σ = rep(exp(θ[2]), n), ξ = rep(exp(θ[3]), n), ...))
            if (is.infinite(res)) {
              res = 1e99 * sign(res)
            }
            res
          },
          x = y,
          log = TRUE,
          control = list(fnscale = -1))
        est_pars = exp(c(res_bgev$par, res_gev$par))
        est_pars[c(1, 4)] = log(est_pars[c(1, 4)])
        gev_rt = return_level_gev(return_periods, est_pars[4], est_pars[5], est_pars[6])
        bgev_rt = return_level_bgev(return_periods, est_pars[1], est_pars[2], est_pars[3])

        data.frame(
          value = c(gev_rt, bgev_rt),
          truth = rep(return_levels, 2),
          model = rep(c("GEV", "bGEV"), each = length(return_periods)),
          n = n,
          iter = i,
          name = paste0("R", return_periods))
      })
    message(i, " / ", n_trials)
    do.call(rbind, res)
  })
res3 = do.call(rbind, res3)

res3$tag = paste0(res3$name, ", n=", res3$n)
res3$tag = factor(res3$tag, levels = unique(res3$tag))

plots[[3]] = res3 |>
  dplyr::filter(value < 50) |> # Remove the worst outliers to make visualisation easier
  ggplot() +
  geom_density(aes(x = value, col = model), size = 1) +
  geom_vline(aes(xintercept = truth)) +
  labs(col = "Model", y = "Density", x = "Return level") +
  facet_wrap(~tag, scales = "free", dir = "v", ncol = length(n_vec)) +
  theme(text = element_text(size = 15))


tikz_plot(
  file = file.path(here::here(), "results", "simulation-study-univariate-MLE.pdf"),
  plot = plots,
  width = 12, height = 12)
