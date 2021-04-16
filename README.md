
The goal of this package is to perform spatial Bayesian hierarchical modelling of extreme short-term
precipitation in Norway using the bGEV distribution

## Installation

You can install the package directly from github

``` r
devtools::install_github("siliusmv/inlaBGEV")
```

## Scripts

There are a number of scripts available in the `exec/` folder. These are described below.

### Univariate data exploration

In the script `local-gev-estiation.R`, the method of probability weighted moments is used for
univariate modelling of yearly maximum precipitation at all the available weather
stations. Exploratory analysis is performed to look for trends in the data.

### Examination of priors

In the script `pc-priors.R` we compute PC-priors for the tail parameter of the GEV, bGEV and GP
distributions and compare them for different penalty parameters.

### Cross validation

The script `cross-validation.R` performs k-fold cross-validation using the scaled threshold
weighted CRPS.

### Estimation of return levels

The script `BGEV-estimation.R` is used for modelling short-term precipitation maxima and
estimating return levels in a map.

### Simulation studies

The scripts `simulation-study-univariate.R`, `simulation-study-in-sample.R` and
`simulation-study-out-of-sample.R` contains simulation studies that evaluate the performance of the
bGEV distribution and the two-step model.

### Visualisation of all weather stations

The script `visualise-data.R` simply plots the location of all weather stations on a map.

## Data

This package contains three data-sets, saved under the `data/` folder. Those are `observations`,
`estimated_sd` and `prediction_grid`. The `observations` is an `sf`-object containing all yearly
precipitation maxima for periods of 1, 3, 6, 12 and 24 hours. The object also contains coordinates
for all stations and the values of the available explanatory variables. The `estimated_sd` is a
`data.frame` that contains the standard deviation of all observations larger than the 99% quantile
from all available weather stations. The `prediction_grid` object is an `sf`-object containing
coordinates and explanatory variables for all locations in a 1x1 km^2 grid over the south of
Norway.

## Numerical dificulties

If you are experiencing numerical difficulties, you might have to stop using inla.pardiso. The inla.pardiso library speeds up computations, but it is also more numerically unstable.
``` r
INLA::inla.setOption(pardiso.license = NULL)
```
