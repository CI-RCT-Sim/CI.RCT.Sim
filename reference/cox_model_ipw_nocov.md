# Cox model with IPW weights without covariates

Cox model with IPW weights without covariates

## Usage

``` r
cox_model_ipw_nocov(data)
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
buffer <-as.numeric( Design[1,"true_window"])
toy_data_trunc <- censor_buffer_window(toy_data,buffer)
toy_data_trunc_long <- data_process_start_stop(toy_data_trunc)
cox_model_ipw_nocov(toy_data_trunc_long)
#> Call:
#> coxph(formula = Surv(t_mace_start, t_mace_stop, event_mace) ~ 
#>     A, data = data, weights = IPW, robust = T, id = ID)
#> 
#>      coef exp(coef) se(coef) robust se    z       p
#> A -0.5024    0.6051   0.1485    0.1570 -3.2 0.00137
#> 
#> Likelihood ratio test=11.65  on 1 df, p=0.000642
#> n= 2288, unique id= 1910, number of events= 177 
```
