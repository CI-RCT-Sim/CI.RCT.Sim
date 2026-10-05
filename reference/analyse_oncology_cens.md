# Create Analyse Functions for

Create Analyse Functions for

## Usage

``` r
analyse_oncology_cens(X)
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

analyse_oncology_cens()(setting, dat)
#> $HR
#>       trt 
#> 0.4219662 
#> 
#> $SElogHR
#> [1] 0.2762358
#> 
#> $low
#>    2.5 % 
#> 0.245553 
#> 
#> $up
#>    97.5 % 
#> 0.7251201 
#> 
#> $p
#> [1] 0.00178697
#> 
#> $N_pat
#> [1] 203
#> 
#> $N_evt
#> [1] 66
#> 
```
