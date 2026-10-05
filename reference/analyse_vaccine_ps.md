# Analyse dataset from vaccine scenario using principal score weighting

Analyse dataset from vaccine scenario using principal score weighting

## Usage

``` r
analyse_vaccine_ps(
  ci_level = 0.95,
  VE_margin = 0.3,
  covariates_in_outcomes_model = TRUE,
  V_unobserved = FALSE,
  W_unobserved = FALSE
)
```

## Arguments

- ci_level:

  the confidence level for the CIs (defaults to 0.95)

- VE_margin:

  vaccine efficacy margin for the super-superiority test

- covariates_in_outcomes_model:

  should the covariates to estimate the principal score also be included
  in the outcomes model

- V_unobserved:

  consider covariate V unobserved (don't use it in analysis)

- W_unobserved:

  consider covariate W unobserved (don't use it in analysis)

## Value

an analyse function that returns a list with the elements

- `p` the p-value of the super-superiority test

- `VE` the point estimate for the vaccine efficacy

- `VE_lower` the lower CI limit for vaccine efficacy

- `VE_upper` the upper CI limit for vaccine efficacy

- `OR` the point estimate for the odds-ratio for infection

- `OR_lower` the lower CI limit for the odds-ratio for infection

- `OR_upper` the upper CI limit for the odds-ratio for infection

## Examples

``` r
Design <- vaccine_scenario() |>
  vaccine_scenario_set_beta_A1_relative() |>
  vaccine_scenario_set_gamma_0() |>
  vaccine_scenario_set_true_eff() |>
  vaccine_scenario_set_samplesize()

dat <- generate_vaccine(Design[3,])
my_analyse <- analyse_vaccine_ps(ci_level=0.95)
my_analyse(Design[3, ], dat)
#> $p
#> [1] 7.604086e-05
#> 
#> $VE
#> [1] 0.8338516
#> 
#> $VE_lower
#> [1] 0.6897293
#> 
#> $VE_upper
#> [1] 1
#> 
#> $OR
#> [1] 0.1652874
#> 
#> $OR_lower
#> [1] 0.078466
#> 
#> $OR_upper
#> [1] 0.3481755
#> 
```
