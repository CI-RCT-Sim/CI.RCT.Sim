# Create Analyse Functions for

Create Analyse Functions for

## Usage

``` r
analyse_oncology_rpsftm(recensor = TRUE, alpha = 0.05, B = 100)
```

## Arguments

- recensor:

  logical, whether to recensor the data after adjusting for treatment
  switching. Default is TRUE.

- alpha:

  numeric, significance level for confidence intervals. Default is 0.05.

- B:

  integer, number of bootstrap samples for estimating standard errors.
  Default is 100.

## Value

an analyse function that can be used in runSimulation

## Examples

``` r
setting <- oncology_scenario()[1, ]
dat <- generate_oncology(setting)
analyse_oncology_rpsftm()(setting, dat)
#> $HR
#> [1] 0.3323076
#> 
#> $SElogHR
#> [1] 1.966864
#> 
#> $low
#> [1] 0.006708503
#> 
#> $up
#> [1] 16.46095
#> 
#> $p
#> [1] 0.5766579
#> 
#> $N_pat
#> [1] 199
#> 
#> $N_evt
#> [1] 66
#> 
```
