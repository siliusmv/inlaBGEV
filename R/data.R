
data_stats = function(x, q = c(.025, .25, .5, .75, .975)) {
  if (all(is.null(dim(x)))) {
    quantiles = quantile(x, probs = q, na.rm = TRUE)
    mean = mean(x, na.rm = TRUE)
    sd = sd(x, na.rm = TRUE)
    res = c(mean = mean, sd = sd, quantiles)
    res = as.data.frame(t(res))
  } else {
    quantiles = matrixStats::rowQuantiles(x, probs = q, na.rm = TRUE)
    mean = matrixStats::rowMeans2(x, na.rm = TRUE)
    sd = matrixStats::rowSds(x, na.rm = TRUE)
    res = cbind(mean = mean, sd = sd, quantiles)
    res = as.data.frame(res)
  }
  res
}

standardise = function(df, stats) {
  for (name in colnames(stats)) {
    df[[name]] = (df[[name]] - stats["mean", name]) / stats["sd", name]
  }
  df
}

un_standardise = function(df, stats) {
  for (name in colnames(stats)) {
    df[[name]] = df[[name]] * stats["sd", name] + stats["mean", name]
  }
  df
}

get_stats_for_standardisation = function(df, covariate_names) {
  statistics = matrix(NA, 2, length(covariate_names))
  rownames(statistics) = c("mean", "sd")
  colnames(statistics) = covariate_names
  for (name in covariate_names) {
    statistics["mean", name] = mean(df[[name]], na.rm = TRUE)
    statistics["sd", name] = sd(df[[name]], na.rm = TRUE)
  }
  statistics
}
