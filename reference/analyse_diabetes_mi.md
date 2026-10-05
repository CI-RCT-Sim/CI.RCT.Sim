# Analyse data set with Multiple Imputation

Creates an analysis function for the diabetes rescue simulation.
Depending on the chosen estimand, either a hypothetical strategy with
rescue-based censoring or a treatment policy strategy using observed
rescue information is applied.

## Usage

``` r
analyse_diabetes_mi(
  strategy = c("hypothetical", "treatment_policy"),
  m = 10,
  maxit = 10,
  ci_level = 0.95,
  seed = 123
)
```

## Arguments

- strategy:

  Either `"hypothetical"` or `"treatment_policy"`.

- m:

  Number of imputations.

- maxit:

  Maximum number of MICE iterations.

- ci_level:

  Confidence level for the two-sided confidence interval.

- seed:

  Random seed for imputation (used locally via withr).

## Value

A function that, when called with `condition` and `dat`, returns a list
with:

- `coef` estimated treatment effect for `trt`

- `p` one-sided p-value for the hypothesis \\\beta\_{trt} \< 0\\

- `ci_lower` lower bound of the two-sided confidence interval

- `ci_upper` upper bound of the two-sided confidence interval

## Details

### Multiple Imputation

Missing HbA1c outcomes are imputed using the mice package.

Imputation is performed separately within each treatment group,
consistent with the trial protocol and avoiding cross-treatment
borrowing of information.

### Estimand-specific handling

#### Treatment policy strategy

Under the treatment policy estimand, observed post-rescue outcomes are
retained.

Rescue indicators may be included as predictors in the imputation model
to preserve associations between rescue use and future HbA1c
trajectories.

#### Hypothetical strategy

Under the hypothetical estimand, the analysis targets the scenario in
which rescue medication had not been available.

For subjects initiating rescue medication at visit \\s\\, all
post-rescue outcomes are set to missing:

\$\$ y\_{i,j} = \mathrm{NA} \quad \text{for all } j \> s \$\$

The outcome at the rescue visit itself (\\j=s\\) is retained.

Rescue indicators are excluded from the HbA1c imputation models under
the hypothetical strategy to avoid conditioning on post-rescue
information.

### Imputation model

HbA1c values are imputed using predictive mean matching (`"pmm"`).

The imputation model includes:

- baseline HbA1c (\\y_0\\),

- age at baseline,

- adjacent HbA1c visit history,

- optionally rescue indicators under the treatment policy strategy.

Rescue indicators are imputed using logistic regression.

A total of `m` imputations are generated.

### Analysis model

For each imputed dataset, an ANCOVA model is fitted:

\$\$ \Delta y_i = \beta_0 + \beta_1 \cdot \mathrm{trt}\_i + \beta_2
\cdot \mathrm{age}\_i + \beta_3 \cdot y\_{0,i} + \varepsilon_i \$\$

where:

- \\\Delta y_i\\ is the change from baseline in HbA1c,

- \\y\_{0,i}\\ is baseline HbA1c.

### Hypothesis testing

One-sided p-values are reported for the treatment effect:

\$\$ H_0: \beta\_{trt} \ge 0 \$\$

versus

\$\$ H_A: \beta\_{trt} \< 0 \$\$

since negative changes in HbA1c are considered favorable.

Two-sided confidence intervals are retained.

### Pooling

Estimates from the imputed datasets are pooled using Rubin's rules via
[`mice::pool()`](https://amices.org/mice/reference/pool.html).

### Reproducibility

Randomness from the multiple imputation procedure is controlled using
[`withr::with_seed()`](https://withr.r-lib.org/reference/with_seed.html),
ensuring reproducibility without affecting the global random number
generator state.

## Examples

``` r
# \donttest{

Design <- diabetes_scenario()[1, ] |>
  diabetes_scenario_set_truevalues()

dat <- generate_diabetes(Design)

## ----------------------------
## Treatment policy estimand
## ----------------------------
res_tp <- analyse_diabetes_mi(
  strategy = "treatment_policy"
)(Design, dat)
#> Warning: Number of logged events: 6
#> Warning: Number of logged events: 7

res_tp
#> $coef
#> [1] -0.3104466
#> 
#> $p
#> [1] 0.06138293
#> 
#> $se
#> [1] 0.1992206
#> 
#> $ci_lower
#> [1] -0.70638
#> 
#> $ci_upper
#> [1] 0.08548682
#> 

## ----------------------------
## Hypothetical estimand
## (censor after rescue + MI)
## ----------------------------
res_hyp <- analyse_diabetes_mi(
  strategy = "hypothetical"
)(Design, dat)

res_hyp
#> $coef
#> [1] -0.3240827
#> 
#> $p
#> [1] 0.06167555
#> 
#> $se
#> [1] 0.2081304
#> 
#> $ci_lower
#> [1] -0.7382111
#> 
#> $ci_upper
#> [1] 0.09004575
#> 

## Compare estimated effects
c(
  treatment_policy = res_tp$coef,
  hypothetical     = res_hyp$coef
)
#> treatment_policy     hypothetical 
#>       -0.3104466       -0.3240827 
# }
```
