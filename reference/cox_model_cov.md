# Cox model with covariates

Cox model with covariates

## Usage

``` r
cox_model_cov(data)
```

## Arguments

- data:

  data in wide format, output of censor_buffer window

## Value

a coxph fit

## Examples

``` r
Design <- mace_scenario() |>
  mace_scenario_set_sample_size()
toy_data <- generate_mace(Design[1,])
buffer <- as.numeric( Design[1,"true_window"])
toy_data_trunc <- censor_buffer_window(toy_data,buffer)
cox_model_cov(toy_data_trunc)
#> Call:
#> coxph(formula = Surv(t_mace, event_mace) ~ A + X + Z, data = data, 
#>     id = ID)
#> 
#>       coef exp(coef) se(coef)      z        p
#> A -0.61884   0.53857  0.16080 -3.849 0.000119
#> X  0.72272   2.06003  0.08568  8.436  < 2e-16
#> Z  0.55798   1.74715  0.08334  6.695 2.16e-11
#> 
#> Likelihood ratio test=111.8  on 3 df, p=< 2.2e-16
#> n= 1910, unique id= 1910, number of events= 162 
```
