# Analyse Dataset with the Inverse probability weighting

Analyse Dataset with the Inverse probability weighting

## Usage

``` r
analyse_ipw(estimand = "tp", level = 0.95, alternative = "two.sided")
```

## Arguments

- estimand:

  the estimand targeted by the CI method "tp" for treatment policy

- level:

  confidence level for CI computation

- alternative:

  alternative hypothesis for the tests "two.sided" or "one.sieded"

## Value

an analyse function that returns a list with the elements

- `p` p value of the score test (two.sided) or the Wald test (one.sided)

- `alternative` the alternative used

- `coef` coefficient for `trt`

- `hr` hazard ratio for `trt`

- `hr_lower` lower 95% confidence intervall boundary for the hazard
  ratio for `trt`

- `hr_upper`lower 95% confidence intervall boundary for the hazard ratio
  for `trt`

- `CI_level` the CI level used

- `N_pat` number of patients

- `N_evt` number of events

## Details

`alternative` can be "two.sided" for a two sided test of equality of the
summary statistic or "one.sided" for a one sided test testing H0:
treatment has equal or shorter survival than control vs. H1 treatment
has longer survival than control.
