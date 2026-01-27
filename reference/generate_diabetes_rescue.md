# Generate Dataset ...

Generate Dataset ...

Create an empty assumptions data.frame for generate_diabetes_rescue

Calculate true summary statistics for scenarios with delayed treatment
effect

## Usage

``` r
generate_diabetes_rescue(condition, fixed_objects = NULL)

assumptions_diabetes_rescue(print = interactive())

true_summary_statistics_diabetes_rescue(
  Design,
  cutoff_stats = 10,
  fixed_objects = NULL
)
```

## Arguments

- condition:

  condition row of Design dataset

- fixed_objects:

  fixed objects not used for now

- print:

  print code to generate parameter set?

- Design:

  Design data.frame for x

- cutoff_stats:

  Cutoff time for rmst and average hazard ratios

## Value

For generate_diabetes_rescue: A data set with the columns id, trt
(1=treatment, 0=control), evt (event, currently TRUE for all
observations)

For assumptions_diabetes_rescue: a design tibble with default values
invisibly

For true_summary_statistics_x: the design data.frame passed as argument
with the additional columns:

- `rmst_trt` rmst in the treatment group

- `median_surv_trt` median survival in the treatment group

- `rmst_ctrl` rmst in the control group

- `median_surv_ctrl` median survial in the control group

- `gAHR` geometric average hazard ratio

- `AHR` average hazard ratio

## Details

Condition has to contain the following columns:

- eff the effect size

- 

- ...

assumptions_diabetes_rescue generates a default design `data.frame` for
use with generate_diabetes_rescue If print is `TRUE` code to produce the
template is also printed for copying, pasting and editing by the user.
(This is the default when run in an interactive session.)

## Functions

- `generate_diabetes_rescue()`: simulates a dataset with ...

- `assumptions_diabetes_rescue()`: generate default design tibble

- `true_summary_statistics_diabetes_rescue()`: calculate true summary
  statistics for ...

## Examples

``` r
Design <- assumptions_diabetes_rescue()
Design
#>    eff rescue_effect k mean_bl mean_age sd_age      b_age pr_rescue      h_y
#> 1    0             0 5       8       60     10 0.06931472      0.05 1.098612
#> 2    1             0 5       8       60     10 0.06931472      0.05 1.098612
#> 3    2             0 5       8       60     10 0.06931472      0.05 1.098612
#> 4    3             0 5       8       60     10 0.06931472      0.05 1.098612
#> 5    4             0 5       8       60     10 0.06931472      0.05 1.098612
#> 6    5             0 5       8       60     10 0.06931472      0.05 1.098612
#> 7    0            -2 5       8       60     10 0.06931472      0.05 1.098612
#> 8    1            -2 5       8       60     10 0.06931472      0.05 1.098612
#> 9    2            -2 5       8       60     10 0.06931472      0.05 1.098612
#> 10   3            -2 5       8       60     10 0.06931472      0.05 1.098612
#> 11   4            -2 5       8       60     10 0.06931472      0.05 1.098612
#> 12   5            -2 5       8       60     10 0.06931472      0.05 1.098612
#> 13   0            -4 5       8       60     10 0.06931472      0.05 1.098612
#> 14   1            -4 5       8       60     10 0.06931472      0.05 1.098612
#> 15   2            -4 5       8       60     10 0.06931472      0.05 1.098612
#> 16   3            -4 5       8       60     10 0.06931472      0.05 1.098612
#> 17   4            -4 5       8       60     10 0.06931472      0.05 1.098612
#> 18   5            -4 5       8       60     10 0.06931472      0.05 1.098612
#> 19   0            -6 5       8       60     10 0.06931472      0.05 1.098612
#> 20   1            -6 5       8       60     10 0.06931472      0.05 1.098612
#> 21   2            -6 5       8       60     10 0.06931472      0.05 1.098612
#> 22   3            -6 5       8       60     10 0.06931472      0.05 1.098612
#> 23   4            -6 5       8       60     10 0.06931472      0.05 1.098612
#> 24   5            -6 5       8       60     10 0.06931472      0.05 1.098612
#> 25   0            -8 5       8       60     10 0.06931472      0.05 1.098612
#> 26   1            -8 5       8       60     10 0.06931472      0.05 1.098612
#> 27   2            -8 5       8       60     10 0.06931472      0.05 1.098612
#> 28   3            -8 5       8       60     10 0.06931472      0.05 1.098612
#> 29   4            -8 5       8       60     10 0.06931472      0.05 1.098612
#> 30   5            -8 5       8       60     10 0.06931472      0.05 1.098612
#> 31   0           -10 5       8       60     10 0.06931472      0.05 1.098612
#> 32   1           -10 5       8       60     10 0.06931472      0.05 1.098612
#> 33   2           -10 5       8       60     10 0.06931472      0.05 1.098612
#> 34   3           -10 5       8       60     10 0.06931472      0.05 1.098612
#> 35   4           -10 5       8       60     10 0.06931472      0.05 1.098612
#> 36   5           -10 5       8       60     10 0.06931472      0.05 1.098612
#>           h_age pr_missing       g_y      g_age  g_rescue
#> 1  -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 2  -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 3  -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 4  -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 5  -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 6  -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 7  -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 8  -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 9  -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 10 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 11 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 12 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 13 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 14 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 15 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 16 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 17 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 18 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 19 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 20 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 21 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 22 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 23 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 24 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 25 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 26 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 27 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 28 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 29 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 30 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 31 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 32 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 33 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 34 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 35 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
#> 36 -0.009950331       0.02 0.4054651 0.01980263 0.4054651
```
