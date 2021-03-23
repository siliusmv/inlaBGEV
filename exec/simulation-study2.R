library(INLA)
library(ggplot2)
library(dplyr)
library(inlaBGEV)
library(parallel)
library(sf)
library(mvtnorm)


get_return_level_function = function(period) {
  function(pars) {
    locscale_pars = locspread_to_locscale(pars$q, pars$s, pars$ξ, α, β)
    return_level_bgev(period, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ)
  }
}


n = 2000 # Number of samples
n_loc = 250 # Number of "locations" that the data are sampled from
n_leave_out_loc = 50
α = .5; β = .8 # Probabilities used in the location and spread parameters
n_trials = 200
block_size = 24 * 365
num_cores = 25
verbose = FALSE

#μ = rnorm(1)
#σ = runif(1)
#ξ = runif(1)
#k = 50
#μ_k = μ + σ / ξ * (k^ξ - 1)
#σ_k = σ * k^ξ
#x = seq(qgev(.1, μ_k, σ_k, ξ), qgev(.9, μ_k, σ_k, ξ), length = 10)
#y1 = pgev(x, μ, σ, ξ)^k
#y2 = pbgev(x, μ_k, σ_k, ξ)
#y1 - y2

set.seed(1, kind = "L'Ecuyer-CMRG")
inclusion = parallel::mclapply(
  X = seq_len(n_trials),
  mc.cores = num_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {

    inclusion = list()
    μ = rep(0, n_loc)
    n_σ = sample.int(4, 1)
    σ_coeffs = c(runif(1, 1, 3), rnorm(n_σ, 0, .2)) / 10
    ξ = rep(runif(1, .01, .4), n_loc)
    covariate_names = list(NULL, paste0("σ_", seq_len(n_σ)), NULL)

    X = cbind(1, matrix(rnorm(n_σ * n_loc), nrow = n_loc))
    colnames(X) = c("intercept", covariate_names[[2]])

    # Simulate coordinates
    #proj = "+proj=utm +zone=33 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=km +no_defs"
    #coords = data.frame(x = runif(n_loc) + 58, y = runif(n_loc) + 60) %>%
    #  st_as_sf(coords = c("x", "y"), crs = st_crs(4326)) %>%
    #  st_transform(proj)

    ## Create the covariance matrix of a spatial Gaussian field
    #distances = as.matrix(st_distance(coords))
    #units(distances) = NULL
    #matern_correlation = function(distances, ρ, ν = 1) {
    #  κ = sqrt(8 * ν) / ρ
    #  Σ = besselK(distances * κ, ν) * (κ * distances) ^ ν * 2 ^ (1 - ν) / gamma(ν)
    #  diag(Σ) = 1 + 1e-10 # Ensure numerical stability
    #  Σ
    #}
    #ρ = 30
    #matern_σ = 1
    #Σ = matern_σ^2 * matern_correlation(distances, ρ)

    σ = as.numeric(exp(X %*% σ_coeffs))
    # We want the block location par to be zero
    μ = μ - σ / ξ * (block_size^ξ * (-log(α))^-ξ - 1)

    location_indices = c(seq_len(n_loc), sample.int(n_loc, n - n_loc, replace = TRUE))

    x = matrix(nrow = block_size, ncol = n)
    for (j in seq_len(n)) {
      x[, j] = rgev(
        block_size, μ[location_indices[j]], σ[location_indices[j]], ξ[location_indices[j]])
    }
    y = apply(x, 2, max)

    df = data.frame(X[location_indices, ])
    df$y = y
    in_sample_df = df[which(location_indices > n_leave_out_loc), ]
    leave_out_df = df[which(location_indices <= n_leave_out_loc), ]

    μ_k = μ + σ / ξ * (block_size^ξ - 1)
    σ_k = σ * block_size^ξ
    q_k = locscale_to_locspread(μ_k, σ_k, ξ, α, β)$q
    s_k = locscale_to_locspread(μ_k, σ_k, ξ, α, β)$s
    r10 = return_level_bgev(10, μ_k, σ_k, ξ)
    r25 = return_level_bgev(25, μ_k, σ_k, ξ)
    r50 = return_level_bgev(50, μ_k, σ_k, ξ)

    # Run R-INLA with the joint model
    joint_res = tryCatch(
      inla_bgev(
        data = in_sample_df,
        response_name = "y",
        diagonal = .1,
        covariate_names = covariate_names,
        α = α,
        β = β),
      error = function(e) NULL)

    if (!is.null(joint_res)) {
      joint_samples = INLA::inla.posterior.sample(1000, joint_res, seed = 1)
      joint_stats = inla_stats(
        sample_list = list(joint_samples),
        data = dplyr::distinct(df, σ_1, .keep_all = TRUE),
        covariate_names = covariate_names,
        s_list = list(rep(joint_res$standardising_const, n_loc)),
        verbose = verbose,
        fun = list(
          r10 = get_return_level_function(10),
          r25 = get_return_level_function(25),
          r50 = get_return_level_function(50)))

      joint_pars = inla_bgev_pars(
        samples = joint_samples,
        data = df,
        covariate_names = covariate_names,
        s_est = rep(joint_res$standardising_const, n_loc))
      joint_score = vector("numeric", n)
      for (j in seq_len(n_loc)) {
        obs = y[which(location_indices == j)]
        par = locspread_to_locscale(joint_pars$q[j, ], joint_pars$s[j, ],
                                    joint_pars$ξ[j, ], α, β)
        if (verbose) message(j)
        joint_score[which(location_indices == j)] = stwcrps_bgev(obs, par$μ, par$σ, par$ξ, .9)
      }

      inclusion$joint = lapply(
        c("q", "s", "ξ", "r10", "r25", "r50"),
        function(name) {
          if (name %in% c("q", "s")) {
            value_name = paste0(name, "_k")
          } else {
            value_name = name
          }
          res = data.frame(
            name = name,
            value = get(value_name),
            lower = joint_stats[[name]]$`2.5%`,
            upper = joint_stats[[name]]$`97.5%`,
            mean = joint_stats[[name]]$mean,
            score = joint_score,
            in_sample = location_indices > n_leave_out_loc,
            n_σ = n_σ,
            model = "joint")
          res$err = res$value - res$mean
          res$included = get(value_name) > res$lower & get(value_name) < res$upper
          res
        }) %>%
        do.call(rbind, .)
    }

    # Run R-inla with the two-step model
    s_est = sapply(
      seq_len(n_loc),
      function(i) {
        obs = as.numeric(x[, which(location_indices == i)])
        sd(obs[obs >= quantile(obs, .8)])
      }
    )

    # Data used for modelling the SD at all observation locations
    sd_df = data.frame(X[-seq_len(n_leave_out_loc), ]) %>%
      dplyr::mutate(log_sd = log(s_est[-seq_len(n_leave_out_loc)]))

    # Estimate s^*
    sd_inla_args = inla_default_args("gaussian")
    sd_inla_args$formula = as.formula(
      paste("log_sd ~ -1 + intercept +", paste(covariate_names[[2]], collapse = " + ")))
    sd_inla_args$data = sd_df
    sd_res = do.call(inla, sd_inla_args)

    log_sd_samples = INLA::inla.posterior.sample(1000, sd_res, seed = 1)
    log_sd_stats = inla_stats(
      sample_list = list(log_sd_samples),
      data = as.data.frame(X),
      covariate_names = covariate_names[[2]],
      verbose = verbose,
      fun = function(pars) exp(rnorm(length(pars$μ), pars$μ, 1 / sqrt(pars$τ))),
      family = "gaussian")
    sd_samples = log_sd_stats$fun$mean

    twostep_res = tryCatch(
      inla_bgev(
        data = in_sample_df,
        response_name = "y",
        s_est = sd_samples[location_indices][location_indices > n_leave_out_loc],
        diagonal = .1,
        covariate_names =  covariate_names,
        α = α,
        β = β),
      error = function(e) NULL)

    if (!is.null(twostep_res)) {
      twostep_samples = INLA::inla.posterior.sample(1000, twostep_res, seed = 1)
      twostep_stats = inla_stats(
        sample_list = list(twostep_samples),
        data = dplyr::distinct(df, σ_1, .keep_all = TRUE),
        covariate_names = covariate_names,
        s_list = list(sd_samples * twostep_res$standardising_const),
        verbose = verbose,
        fun = list(
          r10 = get_return_level_function(10),
          r25 = get_return_level_function(25),
          r50 = get_return_level_function(50)))

      twostep_pars = inla_bgev_pars(
        samples = twostep_samples,
        data = df,
        covariate_names = covariate_names,
        s_est = sd_samples * twostep_res$standardising_const)
      twostep_score = vector("numeric", n)
      for (j in seq_len(n_loc)) {
        obs = y[which(location_indices == j)]
        par = locspread_to_locscale(twostep_pars$q[j, ], twostep_pars$s[j, ],
                                    twostep_pars$ξ[j, ], α, β)
        if (verbose) message(j)
        twostep_score[which(location_indices == j)] = stwcrps_bgev(obs, par$μ, par$σ, par$ξ, .9)
      }

      inclusion$twostep = lapply(
        c("q", "s", "ξ", "r10", "r25", "r50"),
        function(name) {
          if (name %in% c("q", "s")) {
            value_name = paste0(name, "_k")
          } else {
            value_name = name
          }
          res = data.frame(
            name = name,
            value = get(value_name),
            lower = twostep_stats[[name]]$`2.5%`,
            upper = twostep_stats[[name]]$`97.5%`,
            mean = twostep_stats[[name]]$mean,
            in_sample = location_indices > n_leave_out_loc,
            score = twostep_score,
            n_σ = n_σ,
            model = "twostep")
          res$err = res$value - res$mean
          res$included = get(value_name) > res$lower & get(value_name) < res$upper
          res
        }) %>%
        do.call(rbind, .)
    }

    message("Done with iter nr. ", i)
    #message("twostep stats:")
    #print(apply(inclusion$twostep[, 1:6], 2, mean))
    #message("joint stats:")
    #print(apply(inclusion$joint[, 1:6], 2, mean))

    inclusion = rbind(inclusion[[1]], inclusion[[2]])
    inclusion$i = i

    inclusion
  })
inclusion = do.call(rbind, inclusion)

saveRDS(inclusion, file.path(here::here(), "inst", "extdata", "simulation2.rds"))

score_stats = inclusion %>%
  dplyr::group_by(model, i, in_sample) %>%
  dplyr::summarise(score = mean(score), n_σ = unique(n_σ)) %>%
  tidyr::pivot_wider(values_from = score, names_from = model) %>%
  dplyr::mutate(diff = joint - twostep, n_σ = factor(n_σ))

score_stats$diff %>% summary()
score_stats$twostep %>% summary()
score_stats$joint %>% summary()

ggplot(score_stats) +
  geom_point(aes(x = i, y = diff, col = n_σ))

message("Joint StwCRPS in-sample:")
print(summary(dplyr::filter(inclusion, model == "joint", in_sample)$score))
message("Twostep StwCRPS in-sample:")
print(summary(dplyr::filter(inclusion, model == "twostep", in_sample)$score))
message("Joint StwCRPS out-of-sample:")
print(summary(dplyr::filter(inclusion, model == "joint", !in_sample)$score))
message("Twostep StwCRPS out-of-sample:")
print(summary(dplyr::filter(inclusion, model == "twostep", !in_sample)$score))

percentages = inclusion %>%
  dplyr::group_by(name, n_σ, model, in_sample) %>%
  dplyr::mutate(percentage = mean(included)) %>%
  dplyr::select(name, n_σ, model, percentage, in_sample) %>%
  dplyr::slice(1)

ggplot(percentages) +
  geom_col(aes(x = name, y = percentage, fill = model), position = "dodge") +
  facet_wrap(~n_σ + in_sample) +
  geom_hline(yintercept = .95)

message("joint inclusion stats")
dplyr::filter(inclusion, model == "joint") %>%
  dplyr::group_by(in_sample, name) %>%
  dplyr::summarise(percentage = mean(included)) %>%
  print()
message("twostep inclusion stats")
dplyr::filter(inclusion, model == "twostep") %>%
  dplyr::group_by(in_sample, name) %>%
  dplyr::summarise(percentage = mean(included)) %>%
  print()
