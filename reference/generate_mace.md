# Generate Data for MACE Scenario

Generate Data for MACE Scenario

## Usage

``` r
mace_scenario()

mace_scenario_set_sample_size(Design)

generate_mace(condition, fixed_objects = c())

true_trt_mace(Design)
```

## Arguments

- Design:

  Design dataset

- condition:

  row of Design dataset

- fixed_objects:

  list of parameters that are fixed across simulations

## Value

A data.frame with parameter Values for scenarios

for `mace_scenario_set_sample_size`: The Design data set with additional
column sample_size

a simulated dataset

Design data set with additional column `true_trt`

## Functions

- `mace_scenario()`: generate default assumptions `tibble`

- `mace_scenario_set_sample_size()`: Calculate sample size for
  logHRassumed

- `true_trt_mace()`: Calculate true treatment effect for MACE scenario

## Examples

``` r
Design <- mace_scenario() |>
  mace_scenario_set_sample_size()
dat <- generate_mace(Design[1,])
head(dat)
#>   ID          X          Z A      T_disc    t_mace event_mace censor withdraw
#> 1  1 -0.3934121  2.6773716 1  1.70039070 0.1647049          1      0        0
#> 2  2 -2.0561015 -1.5951391 1 59.97530674 1.0000000          0      1        0
#> 3  3 -0.9202663 -0.2330685 1  0.03169908 1.0000000          0      1        0
#> 4  4 -0.5391314  0.3695043 0  1.95904969 0.4659114          1      0        0
#> 5  5 -0.7038993  0.4208998 1  2.66334688 1.0000000          0      1        0
#> 6  6 -1.9404775 -0.3272956 0  6.06857554 1.0000000          0      1        0
#>       t_disc event_disc      T_mace
#> 1 0.16470488          0   0.1647049
#> 2 1.00000000          0 191.0983871
#> 3 0.03169908          1   3.4827215
#> 4 0.46591141          0   0.4659114
#> 5 1.00000000          0  16.2862496
#> 6 1.00000000          0  71.6110630
```
