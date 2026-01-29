# Generate dataset for minimal example

Generate dataset for minimal example

Calculate true summary statistics value for generate_minimal_example

Create an empty assumtions tibble for generate_minimal_example

## Usage

``` r
generate_minimal_example(condition, fixed_objects = NULL)

true_summary_statistics_minimal_example(Design)

assumptions_minimal_example(print = interactive())
```

## Arguments

- condition:

  row of Design dataset

- fixed_objects:

  list of parameters that are fixed across simulations

- Design:

  data.frame with parameter values

- print:

  print code to generate parameter set?

## Value

a simulated dataset

For true_summary_statistics_minimal_example: the Design tibble with
added column eff_size

For assumptions_mimimal_example: a design tibble with default values
invisibly

## Functions

- `true_summary_statistics_minimal_example()`: calculate true effect
  size

- `assumptions_minimal_example()`: generate default assumptions `tibble`

## Examples

``` r
Design <- assumptions_minimal_example() |>
true_summary_statistics_minimal_example()

dat <- generate_minimal_example(Design[1,])
head(dat)
#>   group          y
#> 1     1  1.4864143
#> 2     1  0.1427555
#> 3     1  0.7911104
#> 4     1  0.9781716
#> 5     1  1.3717289
#> 6     1 -1.3006149
```
