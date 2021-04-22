
# This script provides the necessary functions for working with the GEV distribution


#' @export
pgev = function(x, μ, σ, ξ) {
  fix_lengths(x, μ, σ, ξ)
  ifelse(ξ == 0,
         exp(-exp(- (x - μ) / σ)),
         exp(-pmax(0, 1 + ξ * (x - μ) / σ) ^ (-1 / ξ)))
}

#' @export
qgev = function(p, μ, σ, ξ) {
  fix_lengths(p, μ, σ, ξ)
  ifelse(ξ == 0,
         μ - σ * log(-log(p)),
         μ - σ * (1 / ξ) * (1 - (- log(p)) ^ (-ξ)))
}

#' @export
rgev = function(n, μ, σ, ξ) {
  lengths = sapply(list(μ, σ, ξ), length)
  if (any(lengths > n)) stop("Bad input lengths")
  qgev(runif(n), μ, σ, ξ)
}

#' @export
dgev = function(x, μ, σ, ξ, log = FALSE) {
  fix_lengths(x, μ, σ, ξ)
  res = ifelse(ξ == 0,
               -exp(- (x - μ) / σ),
               -pmax(0, 1 + ξ * (x - μ) / σ) ^ (-1 / ξ))
  res = res - log(σ) +
    ifelse(ξ == 0,
           - (x - μ) / σ,
           ifelse(1 + ξ * (x - μ) / σ > 0,
                  - (1 / ξ + 1) * log(1 + ξ * (x - μ) / σ),
                  -Inf))
  if (!log) res = exp(res)
  res
}

#' @export
return_level_gev = function(period, μ, σ, ξ) {
  if (any(period <= 1)) warning("invalid period")
  p = ifelse(period > 1, 1 - 1 / period, NA)
  qgev(p, μ, σ, ξ)
}

#' @export
locscale_to_locspread = function(μ, σ, ξ, α = .5, β = .8) {
  fix_lengths(μ, σ, ξ)
  s = σ * (ℓ(1 - β / 2, ξ) - ℓ(β / 2, ξ)) / ifelse(ξ == 0, 1, ξ)
  q = μ + s * ℓ(α, ξ) / (ℓ(1 - β / 2, ξ) - ℓ(β / 2, ξ))
  list(q = q, s = s, ξ = ξ)
}

#' @export
locspread_to_locscale = function(q, s, ξ, α = .5, β = .8) {
  fix_lengths(q, s, ξ)
  μ = q - s * ℓ(α, ξ) / (ℓ(1 - β / 2, ξ) - ℓ(β / 2, ξ))
  σ = s / (ℓ(1 - β / 2, ξ) - ℓ(β / 2, ξ)) * ifelse(ξ == 0, 1, ξ)
  list(μ = μ, σ = σ, ξ = ξ)
}

ℓ = function(x, ξ) {
  ifelse(ξ == 0, -log(-log(x)), (-log(x))^-ξ - 1)
}
