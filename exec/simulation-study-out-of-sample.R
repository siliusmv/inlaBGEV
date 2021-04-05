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


n = 1500 # Number of samples
n_loc = 250 # Number of "locations" that the data are sampled from
n_leave_out_loc = 50
α = .5; β = .8 # Probabilities used in the location and spread parameters
n_trials = 200
block_size = 24 * 365
num_cores = 25
verbose = FALSE
B = 100

#verbose = TRUE
#num_cores = 1
#n_trials = 1
#n_leave_out_loc = 10
#B = 100

set.seed(1, kind = "L'Ecuyer-CMRG")
res = parallel::mclapply(
  X = seq_len(n_trials),
  mc.cores = num_cores,
  mc.preschedule = FALSE,
  FUN = function(i) {

    res = list(
      inclusion = list(),
      score = list(),
      time = list())

    μ = rep(0, n_loc)
    n_σ = sample.int(4, 1)
    σ_coeffs = c(runif(1, 1, 3), rnorm(n_σ, 0, .2)) / 10
    ξ = rep(runif(1, .01, .4), n_loc)
    covariate_names = list(NULL, paste0("σ_", seq_len(n_σ)), NULL)

    X = cbind(1, matrix(rnorm(n_σ * n_loc), nrow = n_loc))
    colnames(X) = c("intercept", covariate_names[[2]])

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
      joint_score = list()
      joint_etwcrps = vector("numeric", n_leave_out_loc)
      joint_estwcrps = vector("numeric", n_leave_out_loc)
      for (j in seq_len(n_leave_out_loc)) {
        obs = y[which(location_indices == j)]
        par = locspread_to_locscale(joint_pars$q[j, ], joint_pars$s[j, ],
                                    joint_pars$ξ[j, ], α, β)
        if (verbose) message(j)
        joint_score[[j]] = stwcrps_bgev(obs, par$μ, par$σ, par$ξ, .9)
        etwcrps = expected_twcrps_bgev(par$μ, par$σ, par$ξ, .9,
                                       μ_true = μ_k, σ_true = σ_k, ξ_true = ξ)
        joint_etwcrps[j] = etwcrps
        S = expected_twcrps_bgev(par$μ, par$σ, par$ξ, .9)
        joint_estwcrps[j] = etwcrps / S + log(S)
      }
      joint_score = unlist(joint_score)

      res$inclusion$joint = lapply(
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
            in_sample = location_indices > n_leave_out_loc,
            n_σ = n_σ,
            model = "joint")
          res$err = res$value - res$mean
          res$included = get(value_name) > res$lower & get(value_name) < res$upper
          res
        }) %>%
        do.call(rbind, .)

      res$score$joint = data.frame(
        score = mean(joint_score),
        etwcrps = mean(joint_etwcrps),
        estwcrps = mean(joint_estwcrps),
        model = "joint",
        n_σ = n_σ)

      res$time$joint = data.frame(
        time = joint_time,
        model = "joint",
        n_σ = n_σ)
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
    sd_df = as.data.frame(X) %>%
      dplyr::mutate(log_sd = log(s_est))
    sd_df$log_sd[seq_len(n_leave_out_loc)] = NA

    # Estimate s^*
    sd_inla_args = inla_default_args("gaussian")
    sd_inla_args$formula = as.formula(
      paste("log_sd ~ -1 + intercept +", paste(covariate_names[[2]], collapse = " + ")))
    sd_inla_args$data = sd_df

    twostep_time = proc.time()
    #sd_res = do.call(inla, sd_inla_args)
    #sd_samples = exp(sd_res$summary.linear.predictor$mean)
    #sd_samples = s_est
    #twostep_res = tryCatch(
    #  inla_bgev(
    #    data = in_sample_df,
    #    response_name = "y",
    #    s_est = sd_samples[location_indices][location_indices > n_leave_out_loc],
    #    diagonal = .1,
    #    covariate_names =  covariate_names,
    #    α = α,
    #    β = β),
    #  error = function(e) NULL)
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
      #s_est2 = sapply(
      #  seq_len(n_loc),
      #  function(i) {
      #    obs = as.numeric(x[, which(location_indices == i)])
      #    obs = sample(obs, length(obs), replace = TRUE)
      #    sd(obs[obs >= quantile(obs, .8)])
      #  })
      s_est2 = sd_pars$μ[, b]
      twostep_res2 = tryCatch(
        inla_bgev(
          data = df,
          response_name = "y",
          s_est = s_est2[location_indices],
          diagonal = .1,
          covariate_names =  covariate_names,
          α = α,
          β = β),
        error = function(e) NULL)
      if (!is.null(twostep_res2)) {
        #samples[[b]] = INLA::inla.posterior.sample(10, twostep_res2, seed = 1)
        #s_vals[[b]] = s_est2 * twostep_res2$standardising_const
        samples[[b]] = INLA::inla.posterior.sample(num_samples, twostep_res2, seed = 1)
        s_vals[[b]] = s_est2 * twostep_res2$standardising_const
      } else {
        samples[[b]] = NULL
        s_vals[[b]] = NULL
      }
      if (verbose) message("sd_samle nr. ", b)
    }
    twostep_time = proc.time() - twostep_time
    twostep_time = sum(twostep_time[-3]) # That is just how it works...

    twostep_pars = inla_multiple_sample_pars(
      sample_list = samples,
      data = dplyr::distinct(df, σ_1, .keep_all = TRUE),
      covariate_names = covariate_names,
      s_list = s_vals,
      fun = list(
        r10 = get_return_level_function(10),
        r25 = get_return_level_function(25),
        r50 = get_return_level_function(50)))

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

    twostep_score = list()
    twostep_etwcrps = vector("numeric", n_leave_out_loc)
    twostep_estwcrps = vector("numeric", n_leave_out_loc)
    for (j in seq_len(n_leave_out_loc)) {
      obs = y[which(location_indices == j)]
      par = locspread_to_locscale(twostep_pars$q[j, ], twostep_pars$s[j, ],
                                  twostep_pars$ξ[j, ], α, β)
      if (verbose) message(j)
      twostep_score[[j]] = stwcrps_bgev(obs, par$μ, par$σ, par$ξ, .9)
      etwcrps = expected_twcrps_bgev(par$μ, par$σ, par$ξ, .9,
                                     μ_true = μ_k, σ_true = σ_k, ξ_true = ξ)
      twostep_etwcrps[j] = etwcrps
      S = expected_twcrps_bgev(par$μ, par$σ, par$ξ, .9)
      twostep_estwcrps[j] = etwcrps / S + log(S)
    }
    twostep_score = unlist(twostep_score)

    res$inclusion$twostep = lapply(
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
          n_σ = n_σ,
          model = "twostep")
        res$err = res$value - res$mean
        res$included = get(value_name) > res$lower & get(value_name) < res$upper
        res
      }) %>%
      do.call(rbind, .)

    res$score$twostep = data.frame(
      score = mean(twostep_score),
      etwcrps = mean(twostep_etwcrps),
      estwcrps = mean(twostep_estwcrps),
      model = "twostep",
      n_σ = n_σ)

    res$time$twostep = data.frame(
      time = twostep_time,
      model = "twostep",
      n_σ = n_σ)

    message("Done with iter nr. ", i)
    #message("twostep stats:")
    #print(apply(inclusion$twostep[, 1:6], 2, mean))
    #message("joint stats:")
    #print(apply(inclusion$joint[, 1:6], 2, mean))

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
  good_runs = sapply(res, function(x) !is.null(x[[n]]))
  tmp[[n]] = do.call(rbind, lapply(res[good_runs], function(x) x[[n]]))
}
res = tmp

res$score %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(score = mean(score), etwcrps = mean(etwcrps), estwcrps = mean(estwcrps))


score_stats = res$score %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(score = mean(score), n_σ = unique(n_σ)) %>%
  tidyr::pivot_wider(values_from = score, names_from = model) %>%
  dplyr::mutate(diff = joint - twostep, n_σ = factor(n_σ))

message("diff")
score_stats$diff %>% summary()
message("twostep StwCRPS")
score_stats$twostep %>% summary()
message("joint StwCRPS")
score_stats$joint %>% summary()

ggplot(score_stats) +
  geom_point(aes(x = i, y = diff, col = n_σ))

time_stats = res$time %>%
  dplyr::group_by(model) %>%
  dplyr::summarise(mean_time = mean(time),
                   lower_time = quantile(time, .1),
                   upper_time = quantile(time, .9))
print(time_stats)

message("Joint StwCRPS :")
print(summary(dplyr::filter(res$score, model == "joint")$score))
message("Twostep StwCRPS :")
print(summary(dplyr::filter(res$score, model == "twostep")$score))

percentages = res$inclusion %>%
  dplyr::group_by(name, n_σ, model, in_sample) %>%
  dplyr::mutate(percentage = mean(included)) %>%
  dplyr::select(name, n_σ, model, percentage, in_sample) %>%
  dplyr::slice(1)

ggplot(percentages) +
  geom_col(aes(x = name, y = percentage, fill = model), position = "dodge") +
  facet_wrap(~n_σ + in_sample) +
  geom_hline(yintercept = .95)

message("joint inclusion stats")
dplyr::filter(res$inclusion, model == "joint") %>%
  dplyr::group_by(in_sample, name) %>%
  dplyr::summarise(percentage = mean(included)) %>%
  print()
message("twostep inclusion stats")
dplyr::filter(res$inclusion, model == "twostep") %>%
  dplyr::group_by(in_sample, name) %>%
  dplyr::summarise(percentage = mean(included)) %>%
  print()


res$inclusion %>%
  dplyr::group_by(model, name, in_sample) %>%
  dplyr::summarise(mean(included))
