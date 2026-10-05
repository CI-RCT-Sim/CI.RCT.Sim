# Analyse dataset from vaccine scenario using instrumental variable regression

Analyse dataset from vaccine scenario using instrumental variable
regression

## Usage

``` r
analyse_vaccine_ivreg(
  ci_level = 0.95,
  VE_margin = 0.3,
  V_unobserved = FALSE,
  W_unobserved = FALSE
)
```

## Arguments

- ci_level:

  the confidence level for the CIs (defaults to 0.95)

- VE_margin:

  vaccine efficacy margin for the super-superiority test

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

- `RD` the point estimate for risk difference

- `RD_lower` the lower CI limit for risk difference

- `RD_upper` the upper CI limit for risk difference

## Examples

``` r
Design <- vaccine_scenario() |>
  vaccine_scenario_set_beta_A1_relative() |>
  vaccine_scenario_set_gamma_0() |>
  vaccine_scenario_set_true_eff() |>
  vaccine_scenario_set_samplesize()

dat <- generate_vaccine(Design[3,])
my_analyse <- analyse_vaccine_ivreg(ci_level=0.95)
my_analyse(Design[3, ], dat)
#> $p
#> [1] 0.0116831
#> 
#> $VE
#> [1] 0.6404559
#> 
#> $VE_lower
#> [1] 0.4170186
#> 
#> $VE_upper
#> [1] 1
#> 
```
