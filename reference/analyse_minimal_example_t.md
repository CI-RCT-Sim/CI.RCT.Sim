# Analyse dataset from minimal example scenario with t-test

Analyse dataset from minimal example scenario with t-test

## Usage

``` r
analyse_minimal_example_t()
```

## Value

an analyse function that returns a list with the elements

- `p` the p-value

## Examples

``` r
Design <- assumptions_minimal_example() |>
  true_summary_statistics_minimal_example()

condition <- Design[1, ]

dat <- generate_minimal_example(condition)

my_analyse_t  <- analyse_minimal_example_t()
my_analyse_t(condition, dat)
#> $p
#> [1] 7.212631e-05
#> 
```
