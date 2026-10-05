# Create Analyse Functions for ITT analysis

Create Analyse Functions for ITT analysis

## Usage

``` r
analyse_oncology_itt(X)
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
analyse_oncology_itt()(setting, dat)
#> $HR
#>       trt 
#> 0.7861532 
#> 
#> $SElogHR
#> [1] 0.2524533
#> 
#> $low
#>     2.5 % 
#> 0.4793123 
#> 
#> $up
#>   97.5 % 
#> 1.289424 
#> 
#> $p
#> [1] 0.3405588
#> 
#> $N_pat
#> [1] 215
#> 
#> $N_evt
#> [1] 66
#> 
```
