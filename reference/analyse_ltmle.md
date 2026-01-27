# Analyse Dataset with Longitudinal Targeted Minimum Loss/(Maximum Likelihood) Estimation

Analyse Dataset with Longitudinal Targeted Minimum Loss/(Maximum
Likelihood) Estimation

## Usage

``` r
analyse_ltmle(level = 0.95, alternative = "two.sided")
```

## Arguments

- level:

  confidence level for CI computation

- alternative:

  alternative hypothesis for the tests "two.sided" or "one.sided"

## Value

an analyse function that returns a list with the elements

- `p` p value of the score test (two.sided) or the Wald test (one.sided)

- `alternative` the alternative used

- `coef` coefficient for `trt`

- `hr` hazard ratio for `trt`

- `hr_lower` lower 95% confidence intervall boundary for the hazard
  ratio for `trt`

- `hr_upper`lower 95% confidence intervall boundary for the hazard ratio
  for `trt`

- `CI_level` the CI level used

- `N_pat` number of patients

- `N_evt` number of events

## Details

`alternative` can be "two.sided" for a two sided test of equality of the
summary statistic or "one.sided" for a one sided test testing H0:
treatment has equal or shorter survival than control vs. H1 treatment
has longer survival than control.

## Examples

``` r
condition <- merge(
  assumptions_delayed_effect(),
  design_fixed_followup(),
  by = NULL
) |>
  head(1)
#> Error in assumptions_delayed_effect(): could not find function "assumptions_delayed_effect"
dat <- generate_delayed_effect(condition)
#> Error in generate_delayed_effect(condition): could not find function "generate_delayed_effect"
analyse_ltmle()(condition, dat)
#> Error: object 'dat' not found
analyse_ltmle()(assumptions_diabetes_rescue, generate_diabetes_rescue(assumptions_diabetes_rescue()))
#> Error in dplyr::mutate(dat, rescue_start = ifelse(is.na(rescue_start),     condition$k[1] + 2, rescue_start)): ℹ In argument: `rescue_start = ifelse(is.na(rescue_start),
#>   condition$k[1] + 2, rescue_start)`.
#> Caused by error in `condition$k`:
#> ! object of type 'closure' is not subsettable
```
