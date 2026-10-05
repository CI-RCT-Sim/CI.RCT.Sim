# Analysis function with parameter extraction : cox model without covariates

Analysis function with parameter extraction : cox model without
covariates

## Usage

``` r
analyse_mace_cox_nocov(condition, dat, fixed_objects = NULL)
```

## Arguments

- condition:

  SimDesign object

- dat:

  SimDesign object

- fixed_objects:

  SimDesign object

## Value

named list of results coef = estimated treatment effect se = standard
error of the coef HR = hazard ratio CI.l CI.u = lower and upper CI of HR
