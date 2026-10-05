# Cox model with IPW weights with covariates

Cox model with IPW weights with covariates

## Usage

``` r
cox_model_ipw(data)
```

## Arguments

- data:

  simulated data with long format, output of data_process_start_stop
  function

## Value

a coxph fit

## Examples

``` r
Design <- mace_scenario() |>
  mace_scenario_set_sample_size()
toy_data <- generate_mace(Design[1,])
buffer <- as.numeric( Design[1,"true_window"])
toy_data_trunc <- censor_buffer_window(toy_data,buffer)
toy_data_trunc_long <- data_process_start_stop(toy_data_trunc)
cox_model_ipw(toy_data_trunc_long)
#> Call:
#> coxph(formula = Surv(t_mace_start, t_mace_stop, event_mace) ~ 
#>     A + X + Z, data = data, weights = IPW, robust = T, id = ID)
#> 
#>       coef exp(coef) se(coef) robust se      z        p
#> A -0.46596   0.62753  0.15001   0.19708 -2.364   0.0181
#> X  0.72358   2.06179  0.07900   0.09180  7.882 3.21e-15
#> Z  0.79856   2.22235  0.08372   0.12345  6.469 9.89e-11
#> 
#> Likelihood ratio test=163.2  on 3 df, p=< 2.2e-16
#> n= 2275, unique id= 1910, number of events= 160 
```
