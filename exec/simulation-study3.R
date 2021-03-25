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
n_trials = 100
block_size = 24 * 365
num_cores = 25
verbose = FALSE

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

    n_q = sample.int(5, 1)
    n_s = sample.int(4, 1)
    s_coeffs = c(runif(1, 1, 3), rnorm(n_s, 0, .2))
    q_coeffs = c(rnorm(1, 5, 1), rnorm(n_q, 0, 1))
    ξ = runif(1, .01, .4)
    covariate_names = list(paste0("q_", seq_len(n_q)), paste0("s_", seq_len(n_s)), NULL)

    X = cbind(1, matrix(rnorm((n_s + n_q) * n_loc), nrow = n_loc))
    colnames(X) = c("intercept", unlist(covariate_names))

    # Simulate coordinates
    proj = inlaBGEV:::get_proj_xy()
    coords = data.frame(x = runif(n_loc) + 58, y = runif(n_loc) + 60) %>%
      st_as_sf(coords = c("x", "y"), crs = st_crs(4326)) %>%
      st_transform(proj)

    # Create the covariance matrix of a spatial Gaussian field
    distances = as.matrix(st_distance(coords))
    units(distances) = NULL
    matern_correlation = function(distances, ρ, ν = 1) {
      κ = sqrt(8 * ν) / ρ
      Σ = besselK(distances * κ, ν) * (κ * distances) ^ ν * 2 ^ (1 - ν) / gamma(ν)
      diag(Σ) = 1 + 1e-10 # Ensure numerical stability
      Σ
    }
    ρ = 30
    matern_σ = 1
    Σ = matern_σ^2 * matern_correlation(distances, ρ)
    u = as.numeric(mvtnorm::rmvnorm(1, rep(0, n_loc), Σ))

    # Create the mesh used for modelling, and define the prior for the
    # spatial Gaussian field
    boundary = inla.nonconvex.hull(points = st_coordinates(coords), convex = -.15)
    mesh_crs = sp::CRS(as.character(proj)[1])
    mesh = inla.mesh.2d(st_coordinates(coords), boundary = boundary,
                        max.edge = c(5, 15), cutoff = 2, offset = c(15, 50),
                        crs = mesh_crs)
    #plot(mesh)
    spde = inla.spde2.pcmatern(
      mesh = mesh,
      prior.range = c(15, .05),
      prior.sigma = c(.5, .05)) # This is because we standardise stuff

    s = as.numeric(exp(X[, c("intercept", covariate_names[[2]])] %*% s_coeffs))
    q = as.numeric(X[, c("intercept", covariate_names[[1]])] %*% q_coeffs) + u

    location_indices = c(seq_len(n_loc), sample.int(n_loc, n - n_loc, replace = TRUE))

    x = matrix(nrow = block_size, ncol = n)
    for (j in seq_len(n)) {
      # This results in block maxima with distribution y ~ GEV(μ = 0, σ = 1, ξ)
      x[, j] = rgev(block_size, (block_size^-ξ - 1) / ξ, block_size^-ξ, ξ)
    }
    # Because of INLA we prefer to have q = 0, s = 1 over μ = 0, σ = 1
    x = x * ξ / ((-log(1 - β / 2))^-ξ - (-log(β / 2))^-ξ)
    x = x - ((-log(α))^-ξ - 1) / ((-log(1 - β / 2))^-ξ - (-log(β / 2))^-ξ)
    # Locate the block maxima
    y = apply(x, 2, max)

    df = data.frame(X[location_indices, ])
    df$y = (y * s[location_indices]) + q[location_indices]
    df = st_as_sf(cbind(df, coords[location_indices, ]))
    df$s = s[location_indices]
    df$q = q[location_indices]
    df$ξ = ξ
    locscale_pars = locspread_to_locscale(q[location_indices], s[location_indices], ξ, α, β)
    df$r10 = return_level_bgev(10, locscale_pars$μ, locscale_pars$σ, ξ)
    df$r25 = return_level_bgev(25, locscale_pars$μ, locscale_pars$σ, ξ)
    df$r50 = return_level_bgev(50, locscale_pars$μ, locscale_pars$σ, ξ)
    in_sample_df = df[which(location_indices > n_leave_out_loc), ]
    leave_out_df = df[which(location_indices <= n_leave_out_loc), ]

    # Run R-INLA with the joint model
    joint_time = proc.time()
    joint_res = tryCatch(
      inla_bgev(
        data = in_sample_df,
        response_name = "y",
        diagonal = .1,
        spde = spde,
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
        data = dplyr::distinct(df, s_1, .keep_all = TRUE),
        covariate_names = covariate_names,
        mesh = mesh,
        s_list = list(rep(joint_res$standardising_const, n_loc)),
        verbose = verbose,
        fun = list(
          r10 = get_return_level_function(10),
          r25 = get_return_level_function(25),
          r50 = get_return_level_function(50)))

      joint_pars = inla_bgev_pars(
        samples = joint_samples,
        data = df,
        coords = coords,
        mesh = mesh,
        covariate_names = covariate_names,
        s_est = rep(joint_res$standardising_const, n_loc))
      joint_score = list()
      for (j in seq_len(n_leave_out_loc)) {
        obs = df$y[which(location_indices == j)]
        par = locspread_to_locscale(joint_pars$q[j, ], joint_pars$s[j, ],
                                    joint_pars$ξ[j, ], α, β)
        if (verbose) message(j)
        joint_score[[j]] = stwcrps_bgev(obs, par$μ, par$σ, par$ξ, .9)
      }
      joint_score = unlist(joint_score)

      res$inclusion$joint = lapply(
        c("q", "s", "ξ", "r10", "r25", "r50"),
        function(name) {
          res = data.frame(
            name = name,
            value = df[[name]],
            lower = joint_stats[[name]]$`2.5%`,
            upper = joint_stats[[name]]$`97.5%`,
            mean = joint_stats[[name]]$mean,
            in_sample = location_indices > n_leave_out_loc,
            n_s = n_s,
            n_q = n_q,
            model = "joint")
          res$err = res$value - res$mean
          res$included = df[[name]] > res$lower & df[[name]] < res$upper
          res
        }) %>%
        do.call(rbind, .)

      res$score$joint = data.frame(
        score = joint_score,
        model = "joint",
        n_s = n_s, n_q = n_q)

      res$time$joint = data.frame(
        time = joint_time,
        model = "joint",
        n_s = n_s, n_q = n_q)
    }

    # Run R-inla with the two-step model
    s_est = sapply(
      seq_len(n_loc),
      function(i) {
        obs = as.numeric(x[, which(location_indices == i)])
        obs = (obs + q[i]) * s[i]
        sd(obs[obs >= quantile(obs, .9)])
      }
    )

    # Data used for modelling the SD at all observation locations
    sd_df = data.frame(X) %>%
      dplyr::mutate(log_sd = log(s_est))
    sd_df$log_sd[seq_len(n_leave_out_loc)] = NA

    # Estimate s^*
    sd_inla_args = inla_default_args("gaussian")
    sd_inla_args$formula = as.formula(
      paste("log_sd ~ -1 + intercept +", paste(covariate_names[[2]], collapse = " + ")))
    sd_inla_args$data = sd_df

    df = data.frame(X[location_indices, ])
    df$y = (y + q[location_indices]) * s[location_indices]
    df = st_as_sf(cbind(df, coords[location_indices, ]))
    df$s = s[location_indices]
    df$q = q[location_indices] * df$s
    df$ξ = ξ
    locscale_pars = locspread_to_locscale(df$q, df$s, ξ, α, β)
    df$r10 = return_level_bgev(10, locscale_pars$μ, locscale_pars$σ, ξ)
    df$r25 = return_level_bgev(25, locscale_pars$μ, locscale_pars$σ, ξ)
    df$r50 = return_level_bgev(50, locscale_pars$μ, locscale_pars$σ, ξ)
    in_sample_df = df[which(location_indices > n_leave_out_loc), ]
    leave_out_df = df[which(location_indices <= n_leave_out_loc), ]

    twostep_time = proc.time()
    sd_res = do.call(inla, sd_inla_args)
    sd_samples = exp(sd_res$summary.linear.predictor$mean)
    twostep_res = tryCatch(
      inla_bgev(
        data = in_sample_df,
        response_name = "y",
        s_est = sd_samples[location_indices][location_indices > n_leave_out_loc],
        diagonal = .1,
        spde = spde,
        covariate_names = list(covariate_names[[1]], NULL, NULL),
        α = α,
        β = β),
      error = function(e) NULL)
    twostep_time = proc.time() - twostep_time
    twostep_time = sum(twostep_time[-3]) # That is just how it works...

    if (!is.null(twostep_res)) {
      twostep_samples = INLA::inla.posterior.sample(1000, twostep_res, seed = 1)
      twostep_stats = inla_stats(
        sample_list = list(twostep_samples),
        data = dplyr::distinct(df, s_1, .keep_all = TRUE),
        covariate_names = covariate_names,
        mesh = mesh,
        s_list = list(sd_samples * twostep_res$standardising_const),
        verbose = verbose,
        fun = list(
          r10 = get_return_level_function(10),
          r25 = get_return_level_function(25),
          r50 = get_return_level_function(50)))

      twostep_pars = inla_bgev_pars(
        samples = twostep_samples,
        data = df,
        coords = coords,
        mesh = mesh,
        covariate_names = covariate_names,
        s_est = sd_samples * twostep_res$standardising_const)
      twostep_score = list()
      for (j in seq_len(n_leave_out_loc)) {
        obs = df$y[which(location_indices == j)]
        par = locspread_to_locscale(twostep_pars$q[j, ], twostep_pars$s[j, ],
                                    twostep_pars$ξ[j, ], α, β)
        if (verbose) message(j)
        twostep_score[[j]] = stwcrps_bgev(obs, par$μ, par$σ, par$ξ, .9)
      }
      twostep_score = unlist(twostep_score)

      res$inclusion$twostep = lapply(
        c("q", "s", "ξ", "r10", "r25", "r50"),
        function(name) {
          res = data.frame(
            name = name,
            value = df[[name]],
            lower = twostep_stats[[name]]$`2.5%`,
            upper = twostep_stats[[name]]$`97.5%`,
            mean = twostep_stats[[name]]$mean,
            in_sample = location_indices > n_leave_out_loc,
            n_s = n_s, n_q = n_q,
            model = "twostep")
          res$err = res$value - res$mean
          res$included = df[[name]] > res$lower & df[[name]] < res$upper
          res
        }) %>%
        do.call(rbind, .)

      res$score$twostep = data.frame(
        score = twostep_score,
        model = "twostep",
        n_s = n_s, n_q = n_q)

      res$time$twostep = data.frame(
        time = twostep_time,
        model = "twostep",
        n_s = n_s, n_q = n_q)
    }

    message("Done with iter nr. ", i)
    #message("twostep stats:")
    #print(apply(inclusion$twostep[, 1:6], 2, mean))
    #message("joint stats:")
    #print(apply(inclusion$joint[, 1:6], 2, mean))

    for (n in names(res)) {
      res[[n]] = do.call(rbind, res[[n]])
      if (!is.null(res[[n]])) res[[n]]$i = i
    }

    res
  })

saveRDS(res, file.path(here::here(), "inst", "extdata", "simulation3.rds"))


tmp = list()
for (n in names(res[[1]])) {
  good_runs = sapply(res, function(x) !is.null(x[[n]]))
  tmp[[n]] = do.call(rbind, lapply(res[good_runs], function(x) x[[n]]))
}
res = tmp

score_stats = res$score %>%
  dplyr::group_by(model, i) %>%
  dplyr::summarise(score = mean(score), n_s = unique(n_s), n_q = unique(n_q)) %>%
  tidyr::pivot_wider(values_from = score, names_from = model) %>%
  dplyr::mutate(diff = joint - twostep, n_par = factor(n_q + n_s),
                n_q = factor(n_q), n_s = factor(n_s))

message("diff")
score_stats$diff %>% summary()
message("twostep")
score_stats$twostep %>% summary()
message("joint")
score_stats$joint %>% summary()

ggplot(score_stats) +
  geom_point(aes(x = i, y = diff, col = n_par))

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
  dplyr::mutate(n_par = n_q + n_s) %>%
  dplyr::group_by(name, n_par, model, in_sample) %>%
  dplyr::mutate(percentage = mean(included)) %>%
  dplyr::select(name, n_par, model, percentage, in_sample) %>%
  dplyr::slice(1)

ggplot(percentages) +
  geom_col(aes(x = name, y = percentage, fill = model), position = "dodge") +
  facet_wrap(~n_par + in_sample) +
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
