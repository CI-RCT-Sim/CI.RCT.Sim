# Analyse data set with G-computation

Analyse data set with G-computation

## Usage

``` r
analyse_diabetes_gcomputation()
```

## Value

A function that, when called with `condition` and `dat`, returns a list
with:

- `coef` estimated difference in mean change in HbA1c between treatment
  groups

- `se` standard error for coef

- `p` p-value for coef

- `ci_lower` lower bound of 95% confidence interval for coef

- `ci_upper` upper bound of 95% confidence interval for coef

## Details

This function implements G-computation to estimate the average treatment
effect (ATE) on the **change in HbA1c** under a hypothetical scenario
where rescue medication is not available, regardless of treatment
discontinuation.

Two scenarios are considered depending on `setup`:

- `setup = 1`: the effect of rescue goes on top of the effect of the
  assigned treatment

- `setup = 0`: once rescue starts, treatment is deterministically
  recoded to control

The estimand is:

- **Endpoint**: Change in HbA1c at final visit k

- **Treatments**: Experimental vs. control

- **Summary measure**: Difference in means

- **Population**: People with type 2 diabetes

- **Intercurrent events**: Use of rescue medication and treatment
  discontinuation

- **Intercurrent event strategy**:

  - *Hypothetical strategy*: Assume rescue medication was not available
    (i.e., rescue = 0 always)

  - *Treatment policy strategy*: Ignore treatment discontinuation (i.e.,
    continue treatment in counterfactual)

The null hypothesis of no treatment effect is rejected if the 95% CI
does not include 0.

## Examples

``` r
# \donttest{
Design <- diabetes_scenario()[1, ] |>
  diabetes_scenario_set_truevalues()

dat <- generate_diabetes(Design)

analyse_diabetes_gcomputation()(Design, dat)
#> $coef
#> [1] -0.6337911
#> 
#> $se
#> [1] 0.1994487
#> 
#> $p
#> [1] 0.004975124
#> 
#> $ci_lower
#> [1] -0.966468
#> 
#> $ci_upper
#> [1] -0.261006
#> 
# }
```
