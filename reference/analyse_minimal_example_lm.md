# Analyse dataset from minimal example scenario with linear regression

Analyse dataset from minimal example scenario with linear regression

## Usage

``` r
analyse_minimal_example_lm(ci_level = 0.95)
```

## Arguments

- ci_level:

  the confidence level for the CIs (defaults to 0.95)

## Value

an analyse function that returns a list with the elements

- `p` the p-value of the F-test

- `coef` the point estimate

- `ci_lower` the lower CI limit

- `ci_upper` the upper CI limit

## Examples

``` r
Design <- assumptions_minimal_example() |>
  true_summary_statistics_minimal_example()

condition <- Design[1, ]

dat <- generate_minimal_example(condition)

my_analyse_lm <- analyse_minimal_example_lm(ci_level=0.9)
my_analyse_lm(condition, dat)
#> $p
#> [1] 0.0001640878
#> 
#> $coef
#>    group 
#> 0.734723 
#> 
#> $ci_lower
#> [1] 0.4234935
#> 
#> $ci_upper
#> [1] 1.045952
#> 
```
