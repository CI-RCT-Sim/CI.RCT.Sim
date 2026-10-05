# Analyse data set with De-mediation

Analyse data set with De-mediation

## Usage

``` r
analyse_diabetes_demediation(separate = TRUE)
```

## Arguments

- separate:

  Default is TRUE which means that the imputation will be performed
  separately for the treatment and control groups. If FALSE, the
  imputation will be performed on the combined data set, with treatment
  as a predictor variable.

## Value

A function that, when called with `condition` and `dat`, returns a list
with:

- `coef` estimated difference in mean change in HbA1c between treatment
  groups

- `se` standard error for coef

- `ci_lower` lower bound of 95% confidence interval for coef

- `ci_upper` upper bound of 95% confidence interval for coef

## Details

This function implements the de-mediation approach described in the
paper by Olarte Parra et al (2025) for handling rescue medication in
longitudinal diabetes trials. The function first imputes the missing
rescue medication data using the `mice` package, and then performs a
de-mediation analysis by fitting a series of models to estimate the
effect of treatment on the outcome while accounting for the potential
mediation through rescue medication use. The final estimate of the
treatment effect is obtained by pooling the results from the multiple
imputations using Rubin's rules.

## Examples

``` r
# \donttest{
Design <- diabetes_scenario()[1, ] |>
  diabetes_scenario_set_truevalues()

dat <- generate_diabetes(Design)

analyse_diabetes_demediation(separate = TRUE)(Design, dat)
#> $coef
#> [1] -0.6115408
#> 
#> $p
#> [1] 0.0007644099
#> 
#> $se
#> [1] 0.1928813
#> 
#> $ci_lower
#> [1] -0.9896543
#> 
#> $ci_upper
#> [1] -0.2334273
#> 
#> $two_sided
#> [1] 0.00152882
#> 
# }
```
