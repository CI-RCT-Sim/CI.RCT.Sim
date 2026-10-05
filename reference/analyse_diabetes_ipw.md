# Analyse data set with the Inverse probability weighting

Analyse data set with the Inverse probability weighting

## Usage

``` r
analyse_diabetes_ipw(strategy = "hypothetical")
```

## Arguments

- strategy:

  the estimand targeted by the CI method, default: "hypothetical".
  Allowed also: "treatment_policy".

## Value

A function that, when called with `condition` and `dat`, returns a list
with:

- `coef` coefficient for `trt`

- `se` standard error for coef

- `p` p-value for coef

- `ci_lower` lower bound of 95% confidence interval for coef

- `ci_upper` upper bound of 95% confidence interval for coef

## Details

This is a function to implement the inverse probability weighting for
estimation of the treatment effect in a diabetes trial. The function
first reshapes the data to a long format, where each row corresponds to
a visit for a patient. Then, depending on the chosen strategy, it
creates an indicator variable for missing outcomes due to treatment
discontinuation (treatment policy) or due to both treatment
discontinuation and rescue medication use (hypothetical). It fits a
logistic regression model to estimate the probability of missingness and
computes inverse probability weights. Finally, it fits a weighted linear
regression model to estimate the treatment effect on the change in HbA1c
at the final visit, using robust standard errors to account for the
weighting.

## Examples

``` r
Design <- diabetes_scenario()[1, ] |>
  diabetes_scenario_set_truevalues()

dat <- generate_diabetes(Design)

analyse_diabetes_ipw(strategy = "treatment_policy")(Design, dat)
#> $coef
#> [1] -0.245836
#> 
#> $se
#> [1] 0.1953588
#> 
#> $p
#> [1] 0.1057911
#> 
#> $ci_lower
#> [1] -0.6340707
#> 
#> $ci_upper
#> [1] 0.1423986
#> 
#> $n
#> [1] 0
#> 
#> $m
#> [1] 1.633881
#> 
#> $k
#> [1] 1094.964
#> 
#> $s
#> [1] 0
#> 
#> $j
#> [1] 0
#> 
analyse_diabetes_ipw(strategy = "hypothetical")(Design, dat)
#> $coef
#> [1] -0.3173733
#> 
#> $se
#> [1] 0.2102685
#> 
#> $p
#> [1] 0.06767515
#> 
#> $ci_lower
#> [1] -0.7361594
#> 
#> $ci_upper
#> [1] 0.1014128
#> 
#> $n
#> [1] 0
#> 
#> $m
#> [1] 1.621233
#> 
#> $k
#> [1] 181.5266
#> 
#> $s
#> [1] 0
#> 
#> $j
#> [1] 0
#> 
```
