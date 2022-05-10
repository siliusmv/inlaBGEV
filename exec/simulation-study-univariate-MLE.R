library(INLA)
library(ggplot2)
library(dplyr)
library(inlaBGEV)
library(parallel)
library(evd)

α = .5; β = .8 # Probabilities used in the location and spread parameters
q = exp(2.835) * .661
s = exp(2.835) * .118
ξ = .178
return_periods = c(25, 50, 100, 250, 500)
n_vec = c(25, 50, 100, 500, 1000)
n_trials = 500
num_cores = 4
locscale_pars = locspread_to_locscale(q, s, ξ, α, β)
μ = locscale_pars$μ; σ = locscale_pars$σ
return_levels = return_level_gev(return_periods, μ, σ, ξ)

plots = list()

set.seed(1, kind = "L'Ecuyer-CMRG") # Set a seed that works for paralellisation
res = parallel::mclapply(
  X = seq_len(n_trials),
  mc.cores = num_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {
    res = lapply(
      n_vec,
      function(n) {
        y = rgev(n, μ, σ, ξ)
        init = c(μ, log(σ), log(ξ))
        res_bgev = optim(
          par = init,
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
          par = init,
          fn = function(θ, ...) {
            res = sum(inlaBGEV::dgev(μ = rep(θ[1], n), σ = rep(exp(θ[2]), n), ξ = rep(exp(θ[3]), n), ...))
            if (is.infinite(res)) {
              res = 1e299 * sign(res)
            }
            res
          },
          x = y,
          log = TRUE,
          control = list(fnscale = -1))
        est_pars = exp(c(res_gev$par, res_bgev$par))
        est_pars[c(1, 4)] = log(est_pars[c(1, 4)])
        gev_rt = return_level_gev(return_periods, est_pars[1], est_pars[2], est_pars[3])
        bgev_rt = return_level_bgev(return_periods, est_pars[4], est_pars[5], est_pars[6])

        df1 = data.frame(
          value = c(gev_rt, bgev_rt),
          truth = rep(return_levels, 2),
          model = rep(c("GEV", "bGEV"), each = length(return_periods)),
          n = n,
          iter = i,
          init = "good",
          name = paste0("R", return_periods))

        df11 = data.frame(
          value = est_pars,
          truth = rep(c(μ, σ, ξ)),
          model = rep(c("GEV", "bGEV"), each = 3),
          n = n,
          iter = i,
          init = "good",
          name = rep(c("$\\mu$", "$\\sigma$", "$\\xi$"), 2))

        init = c(μ, log(.9), log(ξ))
        res_bgev = optim(
          par = init,
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
          par = init,
          fn = function(θ, ...) {
            res = sum(inlaBGEV::dgev(μ = rep(θ[1], n), σ = rep(exp(θ[2]), n), ξ = rep(exp(θ[3]), n), ...))
            if (is.infinite(res)) {
              res = 1e299 * sign(res)
            }
            res
          },
          x = y,
          log = TRUE,
          control = list(fnscale = -1))
        est_pars = exp(c(res_gev$par, res_bgev$par))
        est_pars[c(1, 4)] = log(est_pars[c(1, 4)])
        gev_rt = return_level_gev(return_periods, est_pars[1], est_pars[2], est_pars[3])
        bgev_rt = return_level_bgev(return_periods, est_pars[4], est_pars[5], est_pars[6])

        df2 = data.frame(
          value = c(gev_rt, bgev_rt),
          truth = rep(return_levels, 2),
          model = rep(c("GEV", "bGEV"), each = length(return_periods)),
          n = n,
          iter = i,
          init = "bad",
          name = paste0("R", return_periods))

        df22 = data.frame(
          value = est_pars,
          truth = rep(c(μ, σ, ξ)),
          model = rep(c("GEV", "bGEV"), each = 3),
          n = n,
          iter = i,
          init = "bad",
          name = rep(c("$\\mu$", "$\\sigma$", "$\\xi$"), 2))

        rbind(df1, df2, df11, df22)
      })
    message(i, " / ", n_trials)
    do.call(rbind, res)
  })
res = do.call(rbind, res)

return_level_index = which(substr(res$name, 1, 1) == "R")
res1 = res[return_level_index, ]
res2 = res[-return_level_index, ]

res1$tag = paste0("$T = ", substr(res1$name, 2, 999), "$, $n = ", res1$n, "$")
res1$tag = factor(res1$tag, levels = unique(res1$tag))

plots[[1]] = res1 |>
  dplyr::filter(init == "good") |>
  ggplot() +
  geom_density(aes(x = value, col = model), size = 1) +
  geom_vline(aes(xintercept = truth)) +
  facet_wrap(~tag, scales = "free", dir = "h", ncol = length(n_vec))

plots[[2]] = res1 |>
  dplyr::filter(init == "bad") |>
  ggplot() +
  geom_density(aes(x = value, col = model), size = 1) +
  geom_vline(aes(xintercept = truth)) +
  facet_wrap(~tag, scales = "free", dir = "h", ncol = length(n_vec))

res2$tag = paste0(res2$name, ", $n = ", res2$n, "$")
res2$tag = factor(res2$tag, levels = unique(res2$tag))

plots[[3]] = res2 |>
  dplyr::filter(init == "good") |>
  ggplot() +
  geom_density(aes(x = value, col = model), size = 1) +
  geom_vline(aes(xintercept = truth)) +
  facet_wrap(~tag, scales = "free", dir = "h", ncol = 3)

plots[[4]] = res2 |>
  dplyr::filter(init == "bad") |>
  ggplot() +
  geom_density(aes(x = value, col = model), size = 1) +
  geom_vline(aes(xintercept = truth)) +
  facet_wrap(~tag, scales = "free", dir = "h", ncol = 3)



for (i in 1:2) {
  plots[[i]] = plots[[i]] +
    labs(col = "Model", y = "Density", x = "Return level") +
    theme_light() +
    theme(
      text = element_text(size = 15),
      axis.text = element_text(size = 8),
      strip.text = element_text(size = 10))
}


tikz_plot(
  file = file.path(here::here(), "results", "simulation-study-univariate-MLE.pdf"),
  plot = plots,
  width = 9, height = 9)




# α = .5; β = .8 # Probabilities used in the location and spread parameters
# σ_gpd = .9
# ξ_gpd = .18
# p = .9
# block_size = 365 * 24
# return_periods = c(25, 50, 100, 250, 500)
# n_vec = c(25, 50, 100, 500, 1000)
# n_trials = 500
# num_cores = 4
# #return_levels = evd::qgpd(1 - 1 / (block_size * return_periods), scale = σ_gpd, shape = ξ_gpd)
# 
# return_levels = evd::qgpd(
#   p = (1 - 1 / (block_size * return_periods) - p) / (1 - p),
#   scale = σ_gpd,
#   shape = ξ_gpd)
# 
# r_gpd_mix = function(n, σ, ξ, p) {
#   res = rep(0, n)
#   is_positive = sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(1-p, p))
#   if (any(is_positive)) {
#     res[is_positive] = evd::rgpd(sum(is_positive), scale = σ, shape = ξ)
#   }
#   res
# }
# 
# set.seed(1, kind = "L'Ecuyer-CMRG") # Set a seed that works for paralellisation
# res2 = parallel::mclapply(
#   X = seq_len(n_trials),
#   mc.cores = num_cores,
#   mc.preschedule = FALSE,
#   FUN = function(i) {
#     res = lapply(
#       n_vec,
#       function(n) {
#         #y = sapply(1:n, function(i) max(evd::rgpd(block_size, scale = σ_gpd, shape = ξ_gpd)))
#         y = sapply(1:n, function(i) max(r_gpd_mix(block_size, σ_gpd, ξ_gpd, p)))
#         init = c(μ, log(σ), log(ξ))
#         res_bgev = optim(
#           par = init,
#           fn = function(θ, ...) {
#             res = sum(inlaBGEV::dbgev(μ = rep(θ[1], n), σ = rep(exp(θ[2]), n), ξ = rep(exp(θ[3]), n), ...))
#             if (any(θ[-1] < -20)) res = -1e299
#             if (is.na(res)) res = -1e299
#             if (is.infinite(res)) {
#               res = 1e299 * sign(res)
#             }
#             res
#           },
#           x = y,
#           log = TRUE,
#           control = list(fnscale = -1))
#         res_gev = optim(
#           par = init,
#           fn = function(θ, ...) {
#             res = sum(inlaBGEV::dgev(μ = rep(θ[1], n), σ = rep(exp(θ[2]), n), ξ = rep(exp(θ[3]), n), ...))
#             if (is.infinite(res)) {
#               res = 1e299 * sign(res)
#             }
#             res
#           },
#           x = y,
#           log = TRUE,
#           control = list(fnscale = -1))
#         est_pars = exp(c(res_bgev$par, res_gev$par))
#         est_pars[c(1, 4)] = log(est_pars[c(1, 4)])
#         gev_rt = return_level_gev(return_periods, est_pars[4], est_pars[5], est_pars[6])
#         bgev_rt = return_level_bgev(return_periods, est_pars[1], est_pars[2], est_pars[3])
# 
#         df1 = data.frame(
#           value = c(gev_rt, bgev_rt),
#           truth = rep(return_levels, 2),
#           model = rep(c("GEV", "bGEV"), each = length(return_periods)),
#           initial = "good",
#           n = n,
#           iter = i,
#           name = paste0("R", return_periods))
# 
#         init = c(μ, log(.5), log(ξ))
#         res_bgev = optim(
#           par = init,
#           fn = function(θ, ...) {
#             res = sum(inlaBGEV::dbgev(μ = rep(θ[1], n), σ = rep(exp(θ[2]), n), ξ = rep(exp(θ[3]), n), ...))
#             if (any(θ[-1] < -20)) res = -1e299
#             if (is.na(res)) res = -1e299
#             if (is.infinite(res)) {
#               res = 1e299 * sign(res)
#             }
#             res
#           },
#           x = y,
#           log = TRUE,
#           control = list(fnscale = -1))
#         res_gev = optim(
#           par = init,
#           fn = function(θ, ...) {
#             res = sum(inlaBGEV::dgev(μ = rep(θ[1], n), σ = rep(exp(θ[2]), n), ξ = rep(exp(θ[3]), n), ...))
#             if (is.infinite(res)) {
#               res = 1e299 * sign(res)
#             }
#             res
#           },
#           x = y,
#           log = TRUE,
#           control = list(fnscale = -1))
#         est_pars = exp(c(res_bgev$par, res_gev$par))
#         est_pars[c(1, 4)] = log(est_pars[c(1, 4)])
#         gev_rt = return_level_gev(return_periods, est_pars[4], est_pars[5], est_pars[6])
#         bgev_rt = return_level_bgev(return_periods, est_pars[1], est_pars[2], est_pars[3])
# 
#         df2 = data.frame(
#           value = c(gev_rt, bgev_rt),
#           truth = rep(return_levels, 2),
#           model = rep(c("GEV", "bGEV"), each = length(return_periods)),
#           initial = "bad",
#           n = n,
#           iter = i,
#           name = paste0("R", return_periods))
# 
#         rbind(df1, df2)
#       })
#     message(i, " / ", n_trials)
#     do.call(rbind, res)
#   })
# res2 = do.call(rbind, res2)
# 
# #res2$tag = paste0(res2$name, ", $n$=", res2$n)
# res2$tag = paste0("$T = ", substr(res2$name, 2, 999), "$, $n = ", res2$n, "$")
# res2$tag = factor(res2$tag, levels = unique(res2$tag))
# 
# plots[[2]] = res2 |>
#   dplyr::filter(initial == "good") |>
#   ggplot() +
#   geom_density(aes(x = value, col = model), size = 1) +
#   geom_vline(aes(xintercept = truth)) +
#   labs(col = "Model", y = "Density", x = "Return level") +
#   facet_wrap(~tag, scales = "free", dir = "h", ncol = length(n_vec)) +
#   theme_light() +
#   theme(text = element_text(size = 15),
#         axis.text = element_text(size = 11))
# 
# 
# plots[[3]] = res2 |>
#   dplyr::filter(initial == "bad") |>
#   ggplot() +
#   geom_density(aes(x = value, col = model), size = 1) +
#   geom_vline(aes(xintercept = truth)) +
#   labs(col = "Model", y = "Density", x = "Return level") +
#   facet_wrap(~tag, scales = "free", dir = "h", ncol = length(n_vec)) +
#   theme_light() +
#   theme(text = element_text(size = 15),
#         axis.text = element_text(size = 11))
# 
# for (i in 1:3) {
#   plots[[i]] = plots[[i]] + theme(axis.text = element_text(size = 8),
#                                   strip.text = element_text(size = 10))
# }
# 
# 
# tikz_plot(
#   file = file.path(here::here(), "results", "simulation-study-univariate-MLE.pdf"),
#   plot = plots,
#   width = 9, height = 9)
