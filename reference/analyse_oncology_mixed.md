# Create Analyse Functions for

Create Analyse Functions for

## Usage

``` r
analyse_oncology_mixed(
  method = c("RPSFTM", "TSE")[1],
  recensor = TRUE,
  alpha = 0.05,
  B = 100,
  trunc_weights = 5,
  show_progress = FALSE,
  use_censoring_IPW = FALSE,
  requ_n_cens = 5
)
```

## Arguments

- method:

  character, method to adjust for treatment switching, either "RPSFTM"
  or "TSE". Default is "RPSFTM".

- recensor:

  logical, whether to recensor the data after adjusting for treatment
  switching. Default is TRUE.

- alpha:

  numeric, significance level for confidence intervals. Default is 0.05.

- B:

  integer, number of bootstrap samples for estimating standard errors.
  Default is 100.

- trunc_weights:

  numeric, maximum value for inverse probability weights to avoid
  extreme weights. Default is 5.

- show_progress:

  logical, whether to print progress during bootstrapping. Default is
  FALSE.

- use_censoring_IPW:

  logical, whether to use inverse probability weighting to adjust for
  censoring in the final Cox model. Default is FALSE.

- requ_n_cens:

  numeric, minimum number of censoring events required to use inverse
  probability weighting for censoring. Default is 5.

## Value

an analyse function that can be used in runSimulation

## Examples

``` r
setting <- oncology_scenario()[1, ]
dat <- generate_oncology(setting)
analyse_oncology_mixed(method = "RPSFTM")(setting, dat)
#> $HR
#> [1] 0.4267059
#> 
#> $SElogHR
#> [1] 0.5091401
#> 
#> $low
#> [1] 0.1573069
#> 
#> $up
#> [1] 1.157469
#> 
#> $p
#> [1] 0.09437798
#> 
analyse_oncology_mixed(method = "TSE")(setting, dat)
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> $HR
#> [1] 0.4136363
#> 
#> $SElogHR
#> [1] 2.601461
#> 
#> $low
#> [1] 0.002524942
#> 
#> $up
#> [1] 67.76197
#> 
#> $p
#> [1] 0.734357
#> 
```
