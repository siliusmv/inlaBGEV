library(INLA)
library(ggplot2)
library(dplyr)
library(inlaBGEV)
library(parallel)
library(sf)
library(mvtnorm)

# Function for computing the expected twCRPS of a forecast when
# the forecast is a mixture of bGEV distributions with parameters (μ, σ, ξ)
# and the true distribution is a GEV distribution with parameters (μ_true, σ_true, ξ_true)
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

# Function factory for creating functions that compute return levels
get_return_level_function = function(period) {
  function(pars) {
    locscale_pars = locspread_to_locscale(pars$q, pars$s, pars$ξ, α, β)
    return_level_bgev(period, locscale_pars$μ, locscale_pars$σ, locscale_pars$ξ)
  }
}


n = 1500 # Number of samples
n_loc = 250 # Number of "locations" that the data are sampled from
n_leave_out_loc = 50 # Number of locations that are only used for testing
α = .5; β = .8 # Probabilities used in the location and spread parameters
n_trials = 200 # Number of times to perform the simulation study
block_size = 24 * 365 # The size of a block
num_cores = 25 # Number of cores for parallelisation
B = 100 # The number of bootstraps used in the two-step model for estimating s^*
verbose = FALSE # Print a lot of progress messages?

# Set a seed that works when doing parallelisation
set.seed(1, kind = "L'Ecuyer-CMRG")
res = parallel::mclapply(
  X = seq_len(n_trials),
  mc.cores = num_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {

    # Prepare the list of results
    res = list(
      inclusion = list(),
      score = list(),
      time = list())

    # Sample ξ and the regression coefficients for σ for every single observation
    n_σ = sample.int(4, 1)
    σ_coeffs = c(runif(1, 1, 3), rnorm(n_σ, 0, .2)) / 10
    ξ = rep(runif(1, .01, .4), n_loc)
    covariate_names = list(NULL, paste0("σ_", seq_len(n_σ)), NULL)

    # Create the design matrix
    X = cbind(1, matrix(rnorm(n_σ * n_loc), nrow = n_loc))
    colnames(X) = c("intercept", covariate_names[[2]])

    # Compute σ, and choose μ such that the distribution of the block maxima have a location
    # parameter that is equal to zero
    σ = as.numeric(exp(X %*% σ_coeffs))
    μ = - σ / ξ * (block_size^ξ * (-log(α))^-ξ - 1)

    # How many block maxima will there be at the different locations?
    location_indices = c(seq_len(n_loc), sample.int(n_loc, n - n_loc, replace = TRUE))

    # Sample all observations at all locations (x) and compute the block maxima (y)
    x = matrix(nrow = block_size, ncol = n)
    for (j in seq_len(n)) {
      x[, j] = rgev(
        block_size, μ[location_indices[j]], σ[location_indices[j]], ξ[location_indices[j]])
    }
    y = apply(x, 2, max)

    # Create a data frame with the covariates and block maxima at all locations
    df = data.frame(X[location_indices, ])
    df$y = y
    # Create a data frame for all training locations and all testing locations
    in_sample_df = df[base::which(location_indices > n_leave_out_loc), ]
    leave_out_df = df[base::which(location_indices <= n_leave_out_loc), ]

    # Create a list containing the true block maxima parameters and return level values
    # at all locations
    truth = list(μ = 0, σ = σ * block_size^ξ, ξ = ξ)
    truth$q = locscale_to_locspread(truth$μ, truth$σ, ξ, α, β)$q
    truth$s = locscale_to_locspread(truth$μ, truth$σ, ξ, α, β)$s
    truth$r10 = return_level_bgev(10, truth$μ, truth$σ, ξ)
    truth$r25 = return_level_bgev(25, truth$μ, truth$σ, ξ)
    truth$r50 = return_level_bgev(50, truth$μ, truth$σ, ξ)

    # ==========================================================================
    # Evaluate the performance of the joint model
    # ==========================================================================

    # Perform inference and report the time it takes
    joint_time = proc.time()
    joint_res = tryCatch(
      inla_bgev(
        data = in_sample_df,
        response_name = "y",
        diagonal = .1,
        covariate_names = covariate_names,
        α = α,
        β = β),
      error = function(e) NULL)
    joint_time = proc.time() - joint_time
    joint_time = sum(joint_time[-3]) # That is just how it works...

    # Estimate parameters, return levels and compute scoring rules, if INLA did not crash
    if (!is.null(joint_res)) {
      # Sample from the estimated posterior
      joint_samples = INLA::inla.posterior.sample(1000, joint_res, seed = 1)
      # Use the samples to compute estimated parameters and return levels at all locations
      joint_pars = inla_bgev_pars(
        samples = joint_samples,
        data = df,
        covariate_names = covariate_names,
        s_est = rep(joint_res$standardising_const, n_loc),
        fun = list(
          r10 = get_return_level_function(10),
          r25 = get_return_level_function(25),
          r50 = get_return_level_function(50)))
      # Use the samples to compute estimated parameters and return levels statistics at all locations
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
      # Compute scoring rules at all testing locations
      joint_score = list()
      joint_etwcrps = vector("numeric", n_leave_out_loc)
      joint_estwcrps = vector("numeric", n_leave_out_loc)
      for (j in seq_len(n_leave_out_loc)) {
        obs = y[base::which(location_indices == j)]
        par = locspread_to_locscale(joint_pars$q[j, ], joint_pars$s[j, ],
                                    joint_pars$ξ[j, ], α, β)
        if (verbose) message(j)
        joint_score[[j]] = stwcrps_bgev(obs, par$μ, par$σ, par$ξ, .9)
        etwcrps = expected_twcrps_bgev_gev(par$μ, par$σ, par$ξ, .9,
                                       μ_true = truth$μ, σ_true = truth$σ, ξ_true = truth$ξ)
        joint_etwcrps[j] = etwcrps
        S = expected_twcrps_bgev(par$μ, par$σ, par$ξ, .9)
        joint_estwcrps[j] = etwcrps / S + log(S)
      }
      joint_score = unlist(joint_score)
      # Test how good the interval estimates are
      res$inclusion$joint = lapply(
        c("q", "s", "ξ", "r10", "r25", "r50"),
        function(name) {
          res = data.frame(
            name = name,
            value = truth[[name]],
            lower = joint_stats[[name]]$`2.5%`,
            upper = joint_stats[[name]]$`97.5%`,
            mean = joint_stats[[name]]$mean,
            in_sample = location_indices > n_leave_out_loc,
            n_σ = n_σ,
            model = "joint",
            stringsAsFactors = FALSE)
          res$err = res$value - res$mean
          res$included = truth[[name]] > res$lower & truth[[name]] < res$upper
          res$mse = base::mean((truth[[name]] - joint_pars[[name]])^2)
          res
        }) %>%
        do.call(rbind, .)
      # Return the mean of all the computed scores
      res$score$joint = data.frame(
        score = base::mean(joint_score),
        etwcrps = base::mean(joint_etwcrps),
        estwcrps = base::mean(joint_estwcrps),
        model = "joint",
        n_σ = n_σ,
        stringsAsFactors = FALSE)
      # Return the time it took to run the inference
      res$time$joint = data.frame(
        time = joint_time,
        model = "joint",
        n_σ = n_σ,
        stringsAsFactors = FALSE)
    }


    # ==========================================================================
    # Evaluate the performance of the two-step model, with bootstrapping
    # ==========================================================================

    # Compute the standard deviation of large observations
    s_est = vapply(
      seq_len(n_loc),
      function(i) {
        obs = as.numeric(x[, base::which(location_indices == i)])
        sd(obs[obs >= quantile(obs, .8)])
      },
      numeric(1))

    # Create a data frame with covariates and SD at all locations
    sd_df = as.data.frame(X) %>%
      dplyr::mutate(log_sd = log(s_est))
    # Set the SD to NA at all testing locations
    sd_df$log_sd[seq_len(n_leave_out_loc)] = NA

    # Prepare for estimating the posterior of s^* at all locations with R-INLA
    sd_inla_args = inla_default_args("gaussian")
    sd_inla_args$formula = as.formula(
      paste("log_sd ~ -1 + intercept +", paste(covariate_names[[2]], collapse = " + ")))
    sd_inla_args$data = sd_df

    # Estimate the posterior of s^* and draw B realisation from it.
    # Use each of these realisations for standardising the response and performing inference
    # with INLA
    twostep_time = proc.time()
    num_samples = 1000 / B
    samples = list()
    s_vals = list()
    sd_res = do.call(inla, sd_inla_args)
    sd_samples = INLA::inla.posterior.sample(B, sd_res, seed = 1)
    sd_pars = inla_gaussian_pars(
      samples = sd_samples,
      data = dplyr::distinct(df, σ_1, .keep_all = TRUE),
      covariate_names = covariate_names[[2]])
    for (b in seq_len(B)) {
      s_est2 = sd_pars$μ[, b]
      twostep_res = tryCatch(
        inla_bgev(
          data = in_sample_df,
          response_name = "y",
          s_est = s_est2[location_indices][location_indices > n_leave_out_loc],
          diagonal = .1,
          covariate_names =  covariate_names,
          α = α,
          β = β),
        error = function(e) NULL)
      # Sample from the posterior of the standardised model
      if (!is.null(twostep_res)) {
        samples[[b]] = INLA::inla.posterior.sample(num_samples, twostep_res, seed = 1)
        s_vals[[b]] = s_est2 * twostep_res$standardising_const
      } else {
        samples[[b]] = NULL
        s_vals[[b]] = NULL
      }
      if (verbose) message("sd_samle nr. ", b)
    }
    twostep_time = proc.time() - twostep_time
    twostep_time = sum(twostep_time[-3]) # That is just how it works...

    # Use the samples to compute estimated parameters and return levels at all locations
    twostep_pars = inla_multiple_sample_pars(
      sample_list = samples,
      data = dplyr::distinct(df, σ_1, .keep_all = TRUE),
      covariate_names = covariate_names,
      s_list = s_vals,
      fun = list(
        r10 = get_return_level_function(10),
        r25 = get_return_level_function(25),
        r50 = get_return_level_function(50)))
    # Use the samples to compute estimated parameters and return levels statistics at all locations
    twostep_stats = inla_stats(
      sample_list = samples,
      data = dplyr::distinct(df, σ_1, .keep_all = TRUE),
      covariate_names = covariate_names,
      s_list = s_vals,
      verbose = verbose,
      fun = list(
        r10 = get_return_level_function(10),
        r25 = get_return_level_function(25),
        r50 = get_return_level_function(50)))
    # Compute scoring rules at all testing locations
    twostep_score = list()
    twostep_etwcrps = vector("numeric", n_leave_out_loc)
    twostep_estwcrps = vector("numeric", n_leave_out_loc)
    for (j in seq_len(n_leave_out_loc)) {
      obs = y[base::which(location_indices == j)]
      par = locspread_to_locscale(twostep_pars$q[j, ], twostep_pars$s[j, ],
                                  twostep_pars$ξ[j, ], α, β)
      if (verbose) message(j)
      twostep_score[[j]] = stwcrps_bgev(obs, par$μ, par$σ, par$ξ, .9)
      etwcrps = expected_twcrps_bgev_gev(par$μ, par$σ, par$ξ, .9,
                                     μ_true = truth$μ, σ_true = truth$σ, ξ_true = truth$ξ)
      twostep_etwcrps[j] = etwcrps
      S = expected_twcrps_bgev(par$μ, par$σ, par$ξ, .9)
      twostep_estwcrps[j] = etwcrps / S + log(S)
    }
    twostep_score = unlist(twostep_score)
    # Test how good the interval estimates are
    res$inclusion$twostep = lapply(
      c("q", "s", "ξ", "r10", "r25", "r50"),
      function(name) {
        res = data.frame(
          name = name,
          value = truth[[name]],
          lower = twostep_stats[[name]]$`2.5%`,
          upper = twostep_stats[[name]]$`97.5%`,
          mean = twostep_stats[[name]]$mean,
          in_sample = location_indices > n_leave_out_loc,
          n_σ = n_σ,
          model = "twostep",
          stringsAsFactors = FALSE)
        res$err = res$value - res$mean
        res$included = truth[[name]] > res$lower & truth[[name]] < res$upper
        res$mse = base::mean((truth[[name]] - twostep_pars[[name]])^2)
        res
      }) %>%
      do.call(rbind, .)
    # Return the mean of all the computed scores
    res$score$twostep = data.frame(
      score = base::mean(twostep_score),
      etwcrps = base::mean(twostep_etwcrps),
      estwcrps = base::mean(twostep_estwcrps),
      model = "twostep",
      n_σ = n_σ,
      stringsAsFactors = FALSE)
    # Return the time it took to run the inference
    res$time$twostep = data.frame(
      time = twostep_time,
      model = "twostep",
      n_σ = n_σ,
      stringsAsFactors = FALSE)

    # ==========================================================================
    # Evaluate the performance of the two-step model, without bootstrapping
    # ==========================================================================

    # Performe inference and report the time it takes
    twostep_one_time = proc.time()
    sd_res = do.call(inla, sd_inla_args)
    sd_samples = exp(sd_res$summary.linear.predictor$mean)
    twostep_one_res = tryCatch(
      inla_bgev(
        data = in_sample_df,
        response_name = "y",
        s_est = sd_samples[location_indices][location_indices > n_leave_out_loc],
        diagonal = .1,
        covariate_names =  covariate_names,
        α = α,
        β = β),
      error = function(e) NULL)
    twostep_one_time = proc.time() - twostep_one_time
    twostep_one_time = sum(twostep_one_time[-3]) # That is just how it works...

    # Estimate parameters, return levels and compute scoring rules, if INLA did not crash
    if (!is.null(twostep_one_res)) {
      # Sample from the estimated posterior
      twostep_one_samples = INLA::inla.posterior.sample(1000, twostep_one_res, seed = 1)
      # Use the samples to compute estimated parameters and return levels at all locations
      twostep_one_stats = inla_stats(
        sample_list = list(twostep_one_samples),
        data = dplyr::distinct(df, σ_1, .keep_all = TRUE),
        covariate_names = covariate_names,
        s_list = list(twostep_one_res$standardising_const * sd_samples),
        verbose = verbose,
        fun = list(
          r10 = get_return_level_function(10),
          r25 = get_return_level_function(25),
          r50 = get_return_level_function(50)))
      # Use the samples to compute estimated parameters and return levels statistics at all locations
      twostep_one_pars = inla_bgev_pars(
        samples = twostep_one_samples,
        data = df,
        covariate_names = covariate_names,
        s_est = twostep_one_res$standardising_const * sd_samples,
        fun = list(
          r10 = get_return_level_function(10),
          r25 = get_return_level_function(25),
          r50 = get_return_level_function(50)))
      # Compute scoring rules at all locations
      twostep_one_score = list()
      twostep_one_etwcrps = vector("numeric", n_leave_out_loc)
      twostep_one_estwcrps = vector("numeric", n_leave_out_loc)
      for (j in seq_len(n_leave_out_loc)) {
        obs = y[base::which(location_indices == j)]
        par = locspread_to_locscale(twostep_one_pars$q[j, ], twostep_one_pars$s[j, ],
                                    twostep_one_pars$ξ[j, ], α, β)
        if (verbose) message(j)
        twostep_one_score[[j]] = stwcrps_bgev(obs, par$μ, par$σ, par$ξ, .9)
        etwcrps = expected_twcrps_bgev_gev(par$μ, par$σ, par$ξ, .9,
                                       μ_true = truth$μ, σ_true = truth$σ, ξ_true = truth$ξ)
        twostep_one_etwcrps[j] = etwcrps
        S = expected_twcrps_bgev(par$μ, par$σ, par$ξ, .9)
        twostep_one_estwcrps[j] = etwcrps / S + log(S)
      }
      twostep_one_score = unlist(twostep_one_score)
      # Test how good the interval estimates are
      res$inclusion$twostep_one = lapply(
        c("q", "s", "ξ", "r10", "r25", "r50"),
        function(name) {
          res = data.frame(
            name = name,
            value = truth[[name]],
            lower = twostep_one_stats[[name]]$`2.5%`,
            upper = twostep_one_stats[[name]]$`97.5%`,
            mean = twostep_one_stats[[name]]$mean,
            in_sample = location_indices > n_leave_out_loc,
            n_σ = n_σ,
            model = "twostep_one",
            stringsAsFactors = FALSE)
          res$err = res$value - res$mean
          res$included = truth[[name]] > res$lower & truth[[name]] < res$upper
          res$mse = base::mean((truth[[name]] - twostep_pars[[name]])^2)
          res
        }) %>%
        do.call(rbind, .)
      # Return the mean of all the computed scores
      res$score$twostep_one = data.frame(
        score = base::mean(twostep_one_score),
        etwcrps = base::mean(twostep_one_etwcrps),
        estwcrps = base::mean(twostep_one_estwcrps),
        model = "twostep_one",
        n_σ = n_σ,
        stringsAsFactors = FALSE)
      # Return the time it took to run the inference
      res$time$twostep_one = data.frame(
        time = twostep_one_time,
        model = "twostep_one",
        n_σ = n_σ,
        stringsAsFactors = FALSE)
    }

    message("Done with iter nr. ", i)

    # Tidy the results
    for (n in names(res)) {
      res[[n]] = do.call(rbind, res[[n]])
      if (!is.null(res[[n]])) res[[n]]$i = i
      res[[n]]$i = i
    }

    res
  })

saveRDS(res, file.path(here::here(), "results", "simulation-out-of-sample.rds"))

tmp = list()
for (n in names(res[[1]])) {
  good_runs = vapply(res, function(x) !is.null(x[[n]]), logical(1))
  tmp[[n]] = do.call(rbind, lapply(res[good_runs], function(x) x[[n]]))
}
res = tmp

res$score %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(score = base::mean(score),
                   etwcrps = base::mean(etwcrps),
                   estwcrps = base::mean(estwcrps))
# This shows that we do better on the observations and on the estwcrps.

# Add a permutation test to be certain
joint_estwcrps = dplyr::filter(res$score, model == "joint")$estwcrps
twostep_estwcrps = dplyr::filter(res$score, model == "twostep")$estwcrps
diff_estwcrps = joint_estwcrps - twostep_estwcrps
summary(diff_estwcrps)
diff_bootstraps = vapply(
  1:1000,
  function(i) base::mean(sample(diff_estwcrps, length(diff_estwcrps), replace = TRUE)),
  numeric(1))
summary(diff_bootstraps)
# Permutation test shows that we have a significant difference!!!!!!!

res$score %>%
  tidyr::pivot_longer(c("score", "etwcrps", "estwcrps")) %>%
  ggplot() +
  geom_point(aes(x = i, y = value, col = model)) +
  facet_wrap(~name, scales = "free_y")

res$score %>%
  tidyr::pivot_longer(c("score", "etwcrps", "estwcrps")) %>%
  tidyr::pivot_wider(names_from = model) %>%
  dplyr::mutate(diff = joint - twostep) %>%
  ggplot() +
  geom_point(aes(x = i, y = diff)) +
  facet_wrap(~name)

res$score %>%
  tidyr::pivot_longer(c("score", "etwcrps", "estwcrps")) %>%
  tidyr::pivot_wider(names_from = model) %>%
  dplyr::mutate(diff = joint - twostep) %>%
  dplyr::group_by(name) %>%
  dplyr::summarise(summ = list(summary(diff))) %>%
  dplyr::pull(summ)

time_stats = res$time %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(mean_time = base::mean(time),
                   lower_time = quantile(time, .1),
                   upper_time = quantile(time, .9))
base::print(time_stats)

percentages = res$inclusion %>%
  dplyr::group_by(name, n_σ, model, in_sample) %>%
  dplyr::mutate(percentage = base::mean(included)) %>%
  dplyr::select(name, n_σ, model, percentage, in_sample) %>%
  dplyr::slice(1)

ggplot(percentages) +
  geom_col(aes(x = name, y = percentage, fill = model), position = "dodge") +
  facet_wrap(~n_σ + in_sample) +
  geom_hline(yintercept = .95)

message("joint inclusion stats")
dplyr::filter(res$inclusion, model == "joint") %>%
  dplyr::group_by(in_sample, name) %>%
  dplyr::summarise(percentage = base::mean(included)) %>%
  base::print()
message("twostep inclusion stats")
dplyr::filter(res$inclusion, model == "twostep") %>%
  dplyr::group_by(in_sample, name) %>%
  dplyr::summarise(percentage = base::mean(included)) %>%
  base::print()
