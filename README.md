
This repository contains the code used for creating all the results in the paper [Modelling Sub-daily Precipitation Extremes with the Blended Generalised Extreme Value Distribution](https://doi.org/10.1007/s13253-022-00500-7),
which performs spatial Bayesian hierarchical modelling of extreme
sub-daily precipitation in Norway using the bGEV distribution. 

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

The script `simulation-study-univariate-MLE.R` contains simulation studies for comparing the
performance of the bGEV with that of the GEV.

### Visualisation of all weather stations

The script `visualise-data.R` simply plots the location of all weather stations on a map.

## Data

This package contains three data-sets, saved under the `data/` folder. Those are `observations`,
`estimated_sd` and `prediction_grid`. `observations` is an `sf`-object containing all yearly
precipitation maxima for periods of 1, 3, 6, 12 and 24 hours. The object also contains coordinates
for all stations and the values of the available explanatory variables. `estimated_sd` is a
`data.frame` that contains the standard deviation of all observations larger than the 99% quantile
from all available weather stations. The `prediction_grid` object is an `sf`-object containing
coordinates and explanatory variables for all locations in a 1x1 km^2 grid over the south of
Norway.

## R package versions
The version of all R packages used for running this code can be found in the file `renv.lock`
