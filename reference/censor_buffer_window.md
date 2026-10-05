# Function that censor MACE events after discontinuation + buffer window t_mace = T_disc + buffer is censored at one year

Function that censor MACE events after discontinuation + buffer window
t_mace = T_disc + buffer is censored at one year

## Usage

``` r
censor_buffer_window(data, buffer)
```

## Arguments

- data:

  data in wide format, output of generate_mace function

- buffer:

  real non-negative value, duration of the buffer window

## Value

data in wide format, with MACE events censored after the buffer window

## Examples

``` r
Design <- mace_scenario() |>
  mace_scenario_set_sample_size()
toy_data <- generate_mace(Design[1,])
buffer <- as.numeric( Design[1,"true_window"])
toy_data_trunc <- censor_buffer_window(toy_data,buffer)
```
