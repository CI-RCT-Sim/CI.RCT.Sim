# Create Analyse Functions for

Create Analyse Functions for

## Usage

``` r
analyse_oncology_TSE(recensor = TRUE, alpha = 0.05, B = 100)
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
analyse_oncology_TSE()(setting, dat)
#> $HR
#> [1] 0.687268
#> 
#> $SElogHR
#> [1] 0.3923182
#> 
#> $low
#> [1] 0.3155382
#> 
#> $up
#> [1] 1.496926
#> 
#> $p
#> [1] 0.3414333
#> 
#> $N_pat
#> [1] 211
#> 
#> $N_evt
#> [1] 66
#> 
```
