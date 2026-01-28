# Generate Grid of Parameter Values

Generate Grid of Parameter Values

## Usage

``` r
params_scenarios_grid(...)
```

## Arguments

- ...:

  vectors of parameter values

## Value

Returns a tibble of parameter values for the simulations.

## Details

The first value of each argument is taken to be the reference value. In
addition to the reference scenario the output contains one line per
additional value in each of the input arguments In each of those lines
only one parameter is varied. Each line differs from the first line in
exactly one column.

List columns can be used as they can be used in tibbles.

## Examples

``` r
params_scenarios_grid(x=1:3, y=11:12, z=21)
#> # A tibble: 4 × 3
#>       x     y     z
#>   <int> <int> <dbl>
#> 1     1    11    21
#> 2     2    11    21
#> 3     3    11    21
#> 4     1    12    21
params_scenarios_grid(x=1:2, y=list(c(1), c(1, 2)))
#> # A tibble: 3 × 2
#>       x y        
#>   <int> <list>   
#> 1     1 <dbl [1]>
#> 2     2 <dbl [1]>
#> 3     1 <dbl [2]>
```
