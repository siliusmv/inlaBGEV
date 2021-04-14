library(INLA)
library(ggplot2)
library(dplyr)
library(inlaBGEV)
library(parallel)


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
        loglik = sum(dbgev(y, mean(μ_s), mean(σ_s), mean(ξ_s), log = TRUE))
        true_loglik = sum(dbgev(y, μ, σ, ξ, log = TRUE))
        gev_loglik = sum(dgev(y, mean(μ_s), mean(σ_s), mean(ξ_s), log = TRUE))
        gev_true_loglik = sum(dgev(y, μ, σ, ξ, log = TRUE))
        pars = c("q", "s", "ξ", "μ", "σ", paste(return_periods, "return period"))
        is_inside = NULL; lower = NULL; upper = NULL; mean_vals = NULL; mse = NULL
        values = c(q, s, ξ, μ, σ, return_levels)
        for (var in c("q", "s", "ξ", "μ", "σ")) {
          interval = quantile(get(paste0(var, "_s")), c(.025, .975))
          is_inside = c(is_inside, get(var) >= interval[1] && get(var) <= interval[2])
          lower = c(lower, interval[1]); upper = c(upper, interval[2])
          mean_vals = c(mean_vals, base::mean(get(paste0(var, "_s"))))
          mse = c(mse, base::mean((get(paste0(var, "_s")) - get(var))^2))
        }
        for (j in seq_along(return_periods)) {
          interval = quantile(est_return_levels[j, ], c(.025, .975))
          is_inside = c(is_inside,
                        return_levels[j] >= interval[1] && return_levels[j] <= interval[2])
          lower = c(lower, interval[1]); upper = c(upper, interval[2])
          mean_vals = c(mean_vals, base::mean(est_return_levels[j, ]))
          mse = c(mse, base::mean((est_return_levels[j, ] - return_levels[j])^2))
        }
        twcrps = base::mean(twcrps_bgev(y, μ_s, σ_s, ξ_s, .9))
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
          larger_than_the_truth = mean_vals > values,
          mse = mse,
          stwcrps = stwcrps,
          twcrps = twcrps,
          etwcrps = etwcrps,
          estwcrps = estwcrps,
          loglik = loglik,
          true_loglik = true_loglik,
          gev_loglik = gev_loglik,
          gev_true_loglik = gev_true_loglik,
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

# This shows the inclusion stats for different parameter/return levels
# Since we overestimate ξ and underestimate σ, that means that we overestimate the large
# return periods and underestimate the small periods. It also means that we do quite well
# on the in-between size return periods, like 10-50 years.
stats %>%
  do.call(rbind, .) %>%
  dplyr::group_by(par, n) %>%
  dplyr::summarise(inclusion_percentage = base::mean(is_inside)) %>%
  ggplot() +
  geom_col(aes(x = as.numeric(factor(n)), y = inclusion_percentage, fill = factor(n))) +
  facet_wrap(~par) +
  geom_hline(yintercept = .95) +
  labs(fill = "n") +
  coord_cartesian(y = c(.4, 1))

# Create a table of inclusion percentages
table = stats %>%
  do.call(rbind, .) %>%
  dplyr::group_by(n, par) %>%
  dplyr::summarise(inclusion_percentage = base::mean(is_inside)) %>%
  tidyr::pivot_wider(names_from = par, values_from = inclusion_percentage) %>%
  .[, c(1, 11, 13, 12, 2, 7, 4, 8)] %>%
  as.matrix()
table[, -1] = paste0("\\(", format(table[, -1] * 100, digits = 3), "\\%\\)")
table = sapply(seq_len(nrow(table)), function(i) c(paste(table[i, ], collapse = " & "), "\\\\\n"))
cat(table)


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
    twcrps = base::mean(twcrps),
    stwcrps = base::mean(stwcrps),
    estwcrps = base::mean(estwcrps),
    etwcrps = base::mean(etwcrps))
# Everything goes down, except the (s)twcrps. This kind of makes sense, I think.
# The reason is that we are "overfitting" to the available data or something like that...


# Plot MSE for all parameters/return levels
stats %>%
  do.call(rbind, .) %>%
  dplyr::group_by(n, par) %>%
  dplyr::summarise(mse = base::mean(mse)) %>%
  ggplot() +
  geom_point(aes(x = log(n), y = log(mse))) +
  facet_wrap(~par, scales = "free_y")

# Create a table of MSE values
table = stats %>%
  do.call(rbind, .) %>%
  dplyr::group_by(n, par) %>%
  dplyr::summarise(mse = base::mean(mse)) %>%
  tidyr::pivot_wider(names_from = par, values_from = mse) %>%
  .[, c(1, 11, 13, 12, 2, 7, 4, 8)] %>%
  as.matrix()
table[, -1] = paste0("\\(", format(log(table[, -1]), digits = 2, scientific = FALSE), "\\)")
table = sapply(seq_len(nrow(table)), function(i) c(paste(table[i, ], collapse = " & "), "\\\\\n"))
cat(table)




# Look at the likelihood surface for the GEV distribution and the bGEV distribution.
# We see that there is a difference between the MLE of the two distributions
μ = 0
σ = 3
ξ = .25
n = 500000
set.seed(1)
y = rgev(n, μ, σ, ξ)
ll_bgev = function(μ, σ, ξ) {
  res = pbapply::pbsapply(
    X = seq_along(μ),
    cl = 7,
    FUN = function(i) {
      sum(dbgev(y, μ[i], σ[i], ξ[i], log = TRUE))
    }
  )
  res
}
ll_gev = function(μ, σ, ξ) {
  res = pbapply::pbsapply(
    X = seq_along(μ),
    cl = 7,
    FUN = function(i) {
      sum(dgev(y, μ[i], σ[i], ξ[i], log = TRUE))
    }
  )
  res
}
df = expand.grid(μ = 0,
                 σ = seq(2.5, 3.5, length = 50),
                 ξ = seq(.2, .35, length = 50))

df$ll_bgev = ll_bgev(df$μ, df$σ, df$ξ)
df$ll_gev = ll_gev(df$μ, df$σ, df$ξ)

gg1 = ggplot(df) +
  geom_contour(aes(x = σ, y = ξ, z = log(log(-ll_bgev))), bins = 500) +
  #geom_raster(aes(x = σ, y = ξ, fill = log(log(-ll_bgev))), alpha = .5) +
  geom_point(data = data.frame(x = σ, y = ξ), aes(x = x, y = y), col = "red") +
  scale_fill_viridis_c() +
  labs(title = "Log-likelihood for the bGEV distribution")
gg2 = ggplot(df) +
  geom_contour(aes(x = σ, y = ξ, z = log(log(-ll_gev))), bins = 500) +
  #geom_raster(aes(x = σ, y = ξ, fill = log(log(-ll_gev))), alpha = .5) +
  geom_point(data = data.frame(x = σ, y = ξ), aes(x = x, y = y), col = "red") +
  scale_fill_viridis_c() +
  labs(title = "Log-likelihood for the GEV distribution")
patchwork::wrap_plots(gg1, gg2)


# This shows that we have a tendency to overestimate ξ and underestimate σ, but
# that is something we already know by now...
stats %>%
  do.call(rbind, .) %>%
  dplyr::filter(n == tail(n_vec, 1), par %in% c("σ", "ξ"), !is_inside) %>%
  dplyr::mutate(below = case_when(val < lower ~ TRUE, val > upper ~ FALSE, TRUE ~ NA)) %>%
  ggplot() +
  geom_point(aes(x = i, y = val, col = below)) +
  geom_errorbar(aes(x = i, ymin = lower, ymax = upper, col = below)) +
  facet_wrap(~par, scales = "free_y")


stats %>%
  do.call(rbind, .) %>%
  dplyr::group_by(par, n) %>%
  dplyr::summarise(larger = mean(larger_than_the_truth)) %>%
  dplyr::mutate(n = factor(n), par = factor(par)) %>%
  ggplot() +
  geom_col(aes(x = n, y = larger, fill = n)) +
  facet_wrap(~par) +
  geom_hline(yintercept = .5)

stats %>%
