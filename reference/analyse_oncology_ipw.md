# Create Analyse Functions for IPW

Create Analyse Functions for IPW

## Usage

``` r
analyse_oncology_ipw(X)
```

## Arguments

- X:

  input can be used to pass parameters to the analyse function

## Value

an analyse function that can be used in runSimulation

## Examples

``` r
setting <- oncology_scenario()[1, ]
dat <- generate_oncology(setting)
analyse_oncology_ipw()(setting, dat)
#> $HR
#>       trt 
#> 0.6402642 
#> 
#> $SElogHR
#> [1] 0.2800608
#> 
#> $low
#>     2.5 % 
#> 0.3698035 
#> 
#> $up
#>  97.5 % 
#> 1.10853 
#> 
#> $p
#> [1] 0.1113705
#> 
#> $N_pat
#> [1] 209
#> 
#> $N_evt
#> [1] 66
#> 
```
