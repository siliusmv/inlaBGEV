
# inlaBGEV

The goal of inlaBGEV is to ...

## Installation

You can install the released version of inlaBGEV from github

``` r
devtools::install_github("siliusmv/inlaBGEV")
```

## Scripts

### Simulation study

The script `exec/simulation-study.R` contains a simulation study.

### Return level maps

### PC priors


## Numerical dificulties

If you are experiencing numerical difficulties, you might have to stop using inla.pardiso. The inla.pardiso library speeds up computations, but it is also more numerically unstable.
``` r
INLA::inla.setOption(pardiso.license = NULL)
```