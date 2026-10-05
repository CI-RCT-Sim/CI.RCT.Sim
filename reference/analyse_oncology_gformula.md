# Create Analyse Functions for g formula

Create Analyse Functions for g formula

## Usage

``` r
analyse_oncology_gformula(
  B = 20,
  reps = 1,
  n_ev_cutoff_no_bootstrap = 100,
  use_censoring_IPW = FALSE,
  requ_n_cens = 5,
  trunc_weights = 5
)
```

## Arguments

- B:

  Number of bootstrap samples

- reps:

  Number of repetitions for the g-formula simulation (per bootstrap
  sample)

- n_ev_cutoff_no_bootstrap:

  If the observed number of events exceeds this cutoff, no bootstrap
  will be used. Only a point estimate for the hazard ratio is then
  calculated, but no further inferential statistics. The default value
  is 100.

- use_censoring_IPW:

  Logical, indicating whether to use inverse probability of censoring
  weighting when calculating the models. Default is FALSE.

- requ_n_cens:

  The minimum number of random censoring events required to use inverse
  probability of censoring weighting when calculating the models.
  Default is 5.

- trunc_weights:

  Weights exceeding this value will be truncated to this value. Default
  is 5.

## Value

an analyse function that can be used in runSimulation

## Examples

``` r
# \donttest{
setting <- oncology_scenario()[1, ]

dat <- generate_oncology(setting)

analyse_oncology_gformula()(setting, dat)
#> $HR
#> [1] 0.4768824
#> 
#> $SElogHR
#> [1] 0.4208236
#> 
#> $low
#> [1] 0.2090286
#> 
#> $up
#> [1] 1.08797
#> 
#> $p
#> [1] 0.07847398
#> 
#> $N_pat
#> [1] 180
#> 
#> $N_evt
#> [1] 66
#> 
# }
```
