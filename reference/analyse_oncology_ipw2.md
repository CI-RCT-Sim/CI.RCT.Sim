# Create Analyse Functions for

Create Analyse Functions for

## Usage

``` r
analyse_oncology_ipw2(
  use_censoring_IPW = FALSE,
  trunc_weights = 5,
  requ_n_cens = 5
)
```

## Arguments

- use_censoring_IPW:

  logical, default is FALSE

- trunc_weights:

  numeric

- requ_n_cens:

  numeric

## Value

an analyse function that can be used in runSimulation

## Examples

``` r
setting <- oncology_scenario()[1, ]
dat <- generate_oncology(setting)
analyse_oncology_ipw2()(setting, dat)
#> $HR
#> [1] 0.6896052
#> 
#> $SElogHR
#> [1] 0.2518295
#> 
#> $low
#> [1] 0.420962
#> 
#> $up
#> [1] 1.129687
#> 
#> $p
#> [1] 0.1400126
#> 
```
