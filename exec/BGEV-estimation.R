library(dplyr)
library(INLA)
library(parallel)
library(sf)
library(inlaBGEV)
library(ggplot2)
library(patchwork)


# In this script we estimate return level maps for short-term precipitation,
# using the two-step model
# Parameters are also estimated

hour_vec = c(1, 3, 6, 12, 24) # Which aggregation lengths are we examining?
α = .5; β = .8 # Probabilities used in the location and spread parameters
min_sd_years = 4L # Minimum number of years before we use the computed SD values
return_level_period = 20 # Period we are computing return levels for
n_sd_samples = 100 # Number of samples drawn from the distribution of the SD
num_cores = 25 # Number of cores used for parallel computations
num_twostep_samples = 2000 # Number of samples used in the two-step model

# A list containing covariate_names for location, spread and tail parameter
covariate_names = list(c("precipitation", "height", "x", "y", "dist_sea"),
                       c("x", "y", "dist_sea"), NULL)

stats = list()
set.seed(1, kind = "L'Ecuyer-CMRG")
for (i in seq_along(hour_vec)) {
  n_hours = hour_vec[i]

  # Filter out the data of interest and standardise the covariates in the observations data
  # and in the prediction data
  data = dplyr::left_join(observations, estimated_sd, by = c("id", "n_hours", "n_years")) %>%
    dplyr::filter(n_hours == !!n_hours)
  standardisation_stats = get_mean_and_sd_stats(data, unique(unlist(covariate_names)))
  data = standardise(data, standardisation_stats)
  prediction_data = standardise(prediction_grid, standardisation_stats)

  # Create the mesh used for modelling, and define the prior for the
  # spatial Gaussian field
  mesh = inla_mesh(data)
  spde = inla.spde2.pcmatern(
    mesh = mesh,
    prior.range = c(75, .05),
    prior.sigma = c(.5, .05))

  # Data used for modelling the SD at all observation locations
  sd_df = data %>%
    dplyr::distinct(id, .keep_all = TRUE) %>%
    dplyr::filter(n_years >= min_sd_years) %>%
    dplyr::mutate(log_sd = log(sd))

  # Estimate s^*
  sd_inla_args = inla_default_args("gaussian")
  sd_inla_args$data = sd_df
  sd_inla_args$data$intercept = 1
  sd_inla_args$formula = as.formula(paste("log_sd ~ -1 + intercept +",
                                          paste(covariate_names[[2]], collapse = " + ")))
  sd_res = do.call(inla, sd_inla_args)

  # Run the two-step model
  samples = twostep_modelling(
    data = data,
    sd_model = sd_res,
    covariate_names = covariate_names,
    response_name = "value",
    n_sd_samples = n_sd_samples,
    prediction_data = prediction_data,
    verbose = FALSE,
    num_samples = num_twostep_samples,
    spde = spde,
    num_cores = num_cores,
    α = α,
    β = β)

  # Compute parameter stats and return level stats at all locations
  stats[[i]] = inla_stats(
    sample_list = lapply(samples, `[[`, "samples"),
    data = prediction_data,
    covariate_names = list(covariate_names[[1]], NULL, NULL),
    s_list = lapply(samples, `[[`, "s_est"),
    mesh = mesh,
    family = "bgev",
    n_batches = 20,
    fun = list(
      return_level = function(pars) {
        locscale_pars = locspread_to_locscale(pars$q, pars$s, pars$ξ, α, β)
        return_level_gev(return_level_period, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ)
      },
      relative_matern_size = function(pars) {
        abs(pars$matern) / abs(pars$q)
      }))

  σ_bgev = unlist(lapply(samples, \(x) sapply(x$samples, \(y) y$hyperpar[1])))
  ρ_samples = unlist(lapply(samples, \(x) sapply(x$samples, \(y) y$hyperpar[3])))
  σ_samples = unlist(lapply(samples, \(x) sapply(x$samples, \(y) y$hyperpar[4])))
  stats[[i]]$ρ = data_stats(ρ_samples)
  stats[[i]]$σ = data_stats(σ_samples)
  stats[[i]]$σ_bgev = data_stats(σ_bgev)
  stats[[i]]$sd_fixed = sd_res$summary.fixed
  stats[[i]]$sd_hyper = sd_res$summary.hyperpar
  stats[[i]]$sd_fixed_with_standardising_const = sapply(
    samples,
    function(x) x$s_coeffs + c(log(x$standardising_const), rep(0, length(x$s_coeffs) - 1)))
  stats[[i]]$sd_fixed_with_standardising_const = data_stats(stats[[i]]$sd_fixed_with_standardising_const)

  beta_vals = lapply(
    seq_along(samples),
    function(j) {
      contents = attr(samples[[j]]$samples, ".contents")
      indices = which(contents$tag %in% c("intercept", covariate_names[[1]]))
      res = vapply(
        samples[[j]]$samples,
        function(z) z$latent[contents$start[indices], ],
        numeric(length(indices)))
      rownames(res) = vapply(rownames(res), function(n) strsplit(n, ":")[[1]][1], character(1))
      res
    }) %>%
    do.call(cbind, .)
  percentage_positive = apply(beta_vals, 1, function(x) mean(x > 0))
  stats[[i]]$beta_q = data_stats(beta_vals)
  stats[[i]]$beta_q$percentage_positive = percentage_positive

  message("Done with ", n_hours, " hours")
  message("Number of succesful runs: ", length(samples), " of ", n_sd_samples)
}
saveRDS(stats, file.path(here::here(), "results", "return-level-stats.rds"))

# Plot the results ===================================================
stats = readRDS(file.path(here::here(), "results", "return-level-stats.rds"))


# Filter out the data of interest and standardise the covariates in the observations data
# and in the prediction data
data = dplyr::left_join(observations, estimated_sd, by = c("id", "n_hours", "n_years")) %>%
  dplyr::filter(n_hours == 1)
standardisation_stats = get_mean_and_sd_stats(data, unique(unlist(covariate_names)))
data = standardise(data, standardisation_stats)
prediction_data = standardise(prediction_grid, standardisation_stats)


# Return levels =======================
my_breaks = seq(17, by = 2, length = 6)
CI_breaks = seq(2, by = 1, length = 6)
p1 = stats[[1]]$return_level %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(breaks = my_breaks, CI_breaks = CI_breaks, use_tex = TRUE,
             size = .3, axis_text = FALSE, add_stations = FALSE)
p1[[1]] = p1[[1]] + labs(title = "1 hour precipitation")
p1 = patchwork::wrap_plots(p1[[1]], p1[[2]], nrow = 2)

my_breaks = seq(27, by = 4, length = 6)
CI_breaks = seq(3.5, by = 1, length = 6)
p2 = stats[[2]]$return_level %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(breaks = my_breaks, CI_breaks = CI_breaks, use_tex = TRUE,
             size = .3, axis_text = FALSE, add_stations = FALSE)
p2[[1]] = p2[[1]] + labs(title = "3 hour precipitation")
p2 = patchwork::wrap_plots(p2[[1]], p2[[2]], nrow = 2)

my_breaks = seq(35, by = 6, length = 6)
CI_breaks = seq(4, by = 1.5, length = 6)
p3 = stats[[3]]$return_level %>%
  cbind(st_geometry(prediction_data)) %>%
  st_as_sf() %>%
  plot_stats(breaks = my_breaks, CI_breaks = CI_breaks, use_tex = TRUE,
             size = .3, axis_text = FALSE, add_stations = FALSE)
p3[[1]] = p3[[1]] + labs(title = "6 hour precipitation")
p3 = patchwork::wrap_plots(p3[[1]], p3[[2]], nrow = 2)

text_size = 12
myplot = patchwork::wrap_plots(
  p1 * theme(text = element_text(size = text_size)),
  p2 * theme(text = element_text(size = text_size)),
  p3 * theme(text = element_text(size = text_size)),
  nrow = 1)

tikz_plot(file.path(here::here(), "results", "return-level-maps.pdf"),
          myplot, width = 10, height = 6)

# Examine the matern field for 1/3/6 hour precipitation ================
plots = list()
for (i in 1:3) {
  plots[[i]] = stats[[i]]$matern %>%
    cbind(st_geometry(prediction_data)) %>%
    st_as_sf() %>%
    plot_stats(axis_text = FALSE, use_tex = TRUE, size = .3)
  plots[[i]] = plots[[i]][[1]] +
    theme(text = element_text(size = 16)) +
    labs(title = paste(c(1, 3, 6)[i], "hour precipitation"))# +
    #scico::scale_fill_scico(palette = "oslo", direction = -1, end = .9, begin = .1, limits = c(-3.7, 3.7))
}
plot = patchwork::wrap_plots(plots, nrow = 1)

tikz_plot(file.path(here::here(), "results", "matern-plot-1,3,6_hours.pdf"), plot, width = 10, height = 5)


# Tail parameter summary =========================

table = NULL
for (i in seq_along(hour_vec)) {
  message(paste(hour_vec[i], "hour ξ:"))
  print(stats[[i]]$ξ)
  table[i] = stats[[i]]$ξ[c(1, 3, 5, 7)] %>%
    round(digits = 3) %>%
    format() %>%
    paste0("\\(", ., "\\)", collapse = " & ")
  table[i] = paste(hour_vec[i], ifelse(hour_vec[i] == 1, "hour", "hours"), "&", table[i], "\\\\\n")
}
cat(unlist(table))


# Range summary =========================
table = NULL
for (i in seq_along(hour_vec)) {
  message(paste(hour_vec[i], "hour ρ:"))
  print(stats[[i]]$ρ)
  table[i] = stats[[i]]$ρ[c(1, 3, 5, 7)] %>%
    round(digits = 0) %>%
    format() %>%
    paste0("\\(", ., "\\)", collapse = " & ")
  table[i] = paste(hour_vec[i], ifelse(hour_vec[i] == 1, "hour", "hours"), "&", table[i], "\\\\\n")
}
cat(unlist(table))

# Regression coefficient summary

table = list()
for (i in seq_along(hour_vec)) {
  message(paste(hour_vec[i], "hour beta_q:"))
  print(stats[[i]]$beta_q)
  cat("\n")
  message(paste(hour_vec[i], "hour beta_s:"))
  print(stats[[i]]$sd_fixed_with_standardising_const)
  cat("\n")
  beta_q = stats[[i]]$beta_q[, c(1, 2, 3, 5, 7)] %>%
    round(digits = 3) %>%
    format()
  beta_s = stats[[i]]$sd_fixed_with_standardising_const[, c(1, 2, 3, 4, 5)] %>%
    round(digits = 3) %>%
    format()
  rownames(beta_s) = c("Intercept", covariate_names[[2]])
  sigmas = matrix(nrow = 2, ncol = 5)
  sigmas[1, ] = as.numeric(stats[[i]]$σ[, c(1, 2, 3, 5, 7)])
  sigmas[2, ] = as.numeric(stats[[i]]$σ_bgev[c(1, 2, 3, 5, 7)])
  sigmas = sigmas %>%
    round(digits = 3) %>%
    format()
  rownames(sigmas) = c("sigma_q", "s_beta")
  tmp_q = NULL
  for (j in seq_len(nrow(beta_q))) {
    tmp_q[j] = paste("&", rownames(beta_q)[j], "&",
                   paste0("\\(", beta_q[j, ], "\\)", collapse = " & "), "\\\\\n")
    if (j == 1) tmp_q[j] = paste("\\(\\bm{\\beta}_\\mu\\)", tmp_q[j], collapse = "")
  }
  tmp_s = NULL
  for (j in seq_len(nrow(beta_s))) {
    tmp_s[j] = paste("&", rownames(beta_s)[j], "&",
                   paste0("\\(", beta_s[j, ], "\\)", collapse = " & "), "\\\\\n")
    if (j == 1) tmp_s[j] = paste("\\(\\bm{\\beta}_\\sigma\\)", tmp_s[j], collapse = "")
  }
  tmp_σ = NULL
  for (j in seq_len(nrow(sigmas))) {
    tmp_σ[j] = paste(rownames(sigmas)[j], "& &",
                   paste0("\\(", sigmas[j, ], "\\)", collapse = " & "), "\\\\\n")
  }
  table[[i]] = c(tmp_q, tmp_σ, tmp_s)
}
final_table = list()
for (i in seq_along(table)) {
  final_table[[i]] = paste("&", table[[i]])
  final_table[[i]] = gsub("& \\\\midrule", "\\\\midrule", final_table[[i]])
  final_table[[i]][1] = paste(hour_vec[i],
                              ifelse(hour_vec[i] == 1, "hour", "hours"),
                              final_table[[i]][1])
  final_table[[i]] = paste(final_table[[i]], collapse = "")
}
final_table = unlist(final_table)
selected_tables = c(1, 2, 3)
final_table = paste(final_table[selected_tables], collapse = "\\midrule\n")
final_table = gsub("dist_sea", "Distance to the open sea", final_table)
final_table = gsub("height", "Altitude", final_table)
final_table = gsub("intercept", "Intercept", final_table)
final_table = gsub("x", "Easting", final_table)
final_table = gsub("y", "Northing", final_table)
final_table = gsub("precipitation", "Mean annual precipitation", final_table)
final_table = gsub("sigma_q", "\\\\(\\\\text{SD}(u_\\\\mu)\\\\)", final_table)
final_table = gsub("s_beta", "\\\\(\\\\sigma^*_\\\\beta\\\\)", final_table)

cat(final_table)
