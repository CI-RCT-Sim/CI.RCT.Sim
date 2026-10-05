# Cox model without covariates

Cox model without covariates

## Usage

``` r
cox_model_nocov(data)
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
cox_model_nocov(toy_data_trunc)
#> Call:
#> coxph(formula = Surv(t_mace, event_mace) ~ A, data = data, id = ID)
#> 
#>      coef exp(coef) se(coef)      z     p
#> A -0.2180    0.8042   0.1519 -1.435 0.151
#> 
#> Likelihood ratio test=2.05  on 1 df, p=0.1523
#> n= 1910, unique id= 1910, number of events= 174 
```
