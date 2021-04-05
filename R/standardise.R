
#' @export
standardise = function(x, stats = NULL, names = NULL) {
  if (is.null(dim(x))) {
    standardise_vec(x, stats)
  } else {
    standardise_df(x, stats, names)
  }
}

#' @export
get_mean_and_sd_stats = function(x, names = NULL) {
  if (is.null(names)) {
    numeric_cols = sapply(x, is.numeric)
    cols = seq_len(ncol(x))[numeric_cols]
    names = names(x)[cols]
  } else {
    cols = which(names(x) %in% names)
  }
  res = matrix(NA, 2, length(cols))
  rownames(res) = c("mean", "sd")
  colnames(res) = names
  for (i in seq_len(ncol(res))) {
    res[1, i] = mean(x[[cols[i]]], na.rm = TRUE)
    res[2, i] = sd(x[[cols[i]]], na.rm = TRUE)
  }
  res
}

standardise_df = function(x, stats = NULL, names = NULL) {
  names = c(colnames(stats), names)
  if (is.null(names)) {
    numeric_cols = sapply(x, is.numeric)
    cols = seq_len(ncol(x))[numeric_cols]
    names = names(x)[cols]
  } else {
    cols = which(names(x) %in% names)
  }
  for (i in seq_along(cols)) {
    if (is.null(stats)) {
      x[[cols[i]]] = standardise_vec(x[[cols[i]]])
    } else {
      x[[cols[i]]] = standardise_vec(x[[cols[i]]], stats[, colnames(stats) == names[i]])
    }
  }
  x
}

standardise_vec = function(x, stats = NULL) {
  if (length(x) == 0) stop("length of x is zero")
  if (is.logical(x)) {
    warning("x is logical. Are you sure you want to standardise this?")
  } else if (!is.numeric(x)) {
    stop("x is not numeric")
  }
  if (is.null(stats)) {
    μ = mean(x, na.rm = TRUE)
    σ = sd(x, na.rm = TRUE)
  } else {
    if (!is.numeric(stats)) stop("stats are not numerical")
    if (length(stats) != 2) stop("stats does not have length 2")
    if (any(is.na(stats))) stop("stats contains NA or NaN")
    μ = stats[1]
    σ = stats[2]
  }
  if (σ == 0) warning("x has a standard deviation of 0. Returned NaN")
  (x - μ) / σ
}
