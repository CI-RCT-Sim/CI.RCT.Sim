# Analyse diabetes endpoint using a Mixed Model for Repeated Measures (MMRM)

This function fits a Mixed Model for Repeated Measures (MMRM) to
estimate the treatment effect at the final visit in a longitudinal
clinical trial setting.

## Usage

``` r
analyse_diabetes_mmrm(
  ci_level = 0.95,
  strategy = c("treatment_policy", "hypothetical")
)
```

## Arguments

- ci_level:

  Confidence level for the two-sided confidence interval (default 0.95)

- strategy:

  Strategy for handling rescue medication:

  - `"treatment_policy"`

  - `"hypothetical"`

## Value

A function that, when called with `condition` and `dat`, returns a list
with:

- `coef` estimated treatment effect at final visit

- `p` one-sided p-value for \\\beta\_{trt} \< 0\\

- `ci_lower` lower bound of the two-sided confidence interval

- `ci_upper` upper bound of the two-sided confidence interval

- `converged` logical indicating whether the model converged

- `covariance` covariance structure used (`us`, `cs`, or `diag`)

- `fallback` logical indicating whether fallback covariance was used

## Details

### Model specification

The following linear mixed model is fitted:

\$\$ y\_{ij} = \beta_0 + \beta_1 \cdot \text{trt}\_i + \sum\_{j}
\gamma_j \cdot \mathbb{1}(\text{visit}\_j) + \sum\_{j} \delta_j \cdot
\left(\text{trt}\_i \times \mathbb{1}(\text{visit}\_j)\right) +
\sum\_{j} \alpha_j \cdot y\_{0,i} \cdot \mathbb{1}(\text{visit}\_j) +
\sum\_{j} \eta_j \cdot \text{age}\_i \cdot \mathbb{1}(\text{visit}\_j) +
\varepsilon\_{ij} \$\$

where:

- \\y\_{ij}\\ is the outcome for subject \\i\\ at visit \\j\\

- \\\text{trt}\_i\\ is the treatment indicator

- \\\text{visit}\_j\\ is a categorical visit factor

- \\y\_{0,i}\\ is baseline value

- \\\text{age}\_i\\ is a baseline covariate

The within-subject covariance is modeled using an unstructured
covariance:

\$\$ \varepsilon_i \sim \mathcal{N}(0, \Sigma) \$\$

where \\\Sigma\\ is an unstructured covariance matrix across visits.

### Estimand and interpretation

The visit factor is re-leveled such that the final visit (visit = k) is
the reference category. Therefore, the main treatment coefficient
corresponds directly to the treatment effect at the final visit.

### Handling of intercurrent events (rescue medication)

Two strategies are supported:

#### Treatment policy

All observed post-baseline data are used regardless of rescue
medication.

#### Hypothetical strategy

For subjects who initiate rescue medication at visit \\s\\, all outcomes
strictly after rescue are set to missing:

\$\$ y\_{i,j} = \mathrm{NA} \quad \text{for all } j \> s \$\$

Subjects without rescue (or with rescue after the final visit) remain
unchanged.

This targets a hypothetical estimand corresponding to outcomes that
would have been observed had rescue medication not been initiated.

### Inference

Treatment effects are estimated using restricted maximum likelihood
(REML).

Degrees of freedom are computed using Satterthwaite approximation.

#### Hypothesis testing

A one-sided test is used for the treatment effect:

\$\$ H_0: \beta\_{trt} \ge 0 \quad \text{vs} \quad H_A: \beta\_{trt} \<
0 \$\$

since a reduction in HbA1c (negative change) is considered beneficial.

Two-sided confidence intervals are retained and reported.

### Covariance structures

The model attempts the following covariance structures in order:
unstructured (`us`), compound symmetry (`cs`), and diagonal (`diag`),
falling back if convergence fails.

## Examples

``` r
# \donttest{
Design <- diabetes_scenario()[1, ] |>
  diabetes_scenario_set_truevalues()

dat <- generate_diabetes(Design)

analyse_diabetes_mmrm(strategy = "treatment_policy")(Design, dat)
#> $p
#>         trt 
#> 0.001130536 
#> 
#> $coef
#>        trt 
#> -0.6927603 
#> 
#> $se
#> [1] 0.2058947
#> 
#> $ci_lower
#>       trt 
#> -1.101842 
#> 
#> $ci_upper
#>        trt 
#> -0.2836783 
#> 
#> $converged
#> [1] TRUE
#> 
#> $covariance
#> [1] "us"
#> 
#> $fallback
#> [1] FALSE
#> 
analyse_diabetes_mmrm(strategy = "hypothetical")(Design, dat)
#> $p
#>       trt 
#> 0.0013923 
#> 
#> $coef
#>        trt 
#> -0.6995528 
#> 
#> $se
#> [1] 0.2117798
#> 
#> $ci_lower
#>       trt 
#> -1.120535 
#> 
#> $ci_upper
#>        trt 
#> -0.2785711 
#> 
#> $converged
#> [1] TRUE
#> 
#> $covariance
#> [1] "us"
#> 
#> $fallback
#> [1] FALSE
#> 
# }
```
