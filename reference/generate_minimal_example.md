# Create an empty assumtions tibble for generate_minimal_example

Create an empty assumtions tibble for generate_minimal_example

Calculate true summary statistics value for generate_minimal_example

Generate dataset for minimal example

## Usage

``` r
assumptions_minimal_example(print = interactive())

true_summary_statistics_minimal_example(Design)

generate_minimal_example(condition, fixed_objects = NULL)
```

## Arguments

- print:

  print code to generate parameter set?

- Design:

- condition:

  row of Design dataset

- fixed_objects:

  list of parameters that are fixed across simulations

## Value

For assumptions_mimimal_example: a design tibble with default values
invisibly

For true_summary_statistics_minimal_example: the Design tibble with
added column eff_size

a simulated dataset

## Functions

- `assumptions_minimal_example()`: generate default assumptions `tibble`

- `true_summary_statistics_minimal_example()`: calculate true effect
  size

## Examples

``` r
Design <- assumptions_minimal_example() |>
true_summary_statistics_minimal_example()

dat <- generate_minimal_example(Design[1,])
head(dat)
#>   group          y
#> 1     1  0.3601949
#> 2     1  2.5666906
#> 3     1 -0.4491557
#> 4     1  0.2084960
#> 5     1  0.4955193
#> 6     1  1.4018267
```
