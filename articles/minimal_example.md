# Minimal Example

## Intro

This vignette is a short explanation of the structure of the
`CI.RCT.Sim` package. And gives an explanation on how to contribute code
for the data-generation of scenarios as well as for analysis functions.

The example uses a simple example of a test and an estimator for the
difference in means of two normal distributions.

For to run the simulations the package `SimDesign` is used. Some
functions for sumarisation of simulation results are re-used from the
`SimNPH` package written for a previous simulation study.

## Data Generation

Let’s look at the code for the data generation. The simulation study
uses the `SimDesign` pacakge. To run simulations with SimDesign, one has
to define three functions for *Data Generation*, *Analysis* and
*Summarisation*. Each of those functions has to use fixed parameters.
The function signatures can be displayed with
[`SimDesign::SimFunctions()`](http://philchalmers.github.io/SimDesign/reference/SimFunctions.md).
This is the first function `generate_minimal_example` displayed in the
code below.

In `CI.RCT.Sim` there’s two additional functions for each
data-generating function. `true_summary_statistics_minimal_example`,
calculates the true (or assymptotic) values of summary statistics from
the parameters given in a “Design” dataset. This is specific to the data
generating model and the estimand.

The third function is a function that prints and outputs an example of a
Design dataset that can be used with the data generating functions. This
function gives a dataset that is directly useable with the data
generating function for testing, documentation of the parameter and for
users to have boilerplate code to be able to cerate custom parameter
ranges.

`R/generate_minimal_example.R`

``` r
#' Generate dataset for minimal example
#'
#' @param condition row of Design dataset
#' @param fixed_objects list of parameters that are fixed across simulations
#'
#' @returns a simulated dataset
#' @export
#'
#' @examples
#' Design <- assumptions_minimal_example() |>
#' true_summary_statistics_minimal_example()
#'
#' dat <- generate_minimal_example(Design[1,])
#' head(dat)
generate_minimal_example <- function(condition, fixed_objects=NULL){
  data.frame(
    group = c(rep(1, condition$n), rep(0, condition$n)),
    y = c(rnorm(condition$n, condition$mean1), rnorm(condition$n, condition$mean0))
  )
}

#' Calculate true summary statistics value for generate_minimal_example
#'
#' @param Design
#'
#' @returns For true_summary_statistics_minimal_example: the Design tibble with added column eff_size
#' @export
#' @describeIn generate_minimal_example calculate true effect size
true_summary_statistics_minimal_example <- function(Design){
  Design$eff_size <- Design$mean1 - Design$mean0

  Design
}

#' Create an empty assumtions tibble for generate_minimal_example
#'
#' @param print print code to generate parameter set?
#'
#' @returns For assumptions_mimimal_example: a design tibble with default values invisibly
#' @export
#' @describeIn generate_minimal_example generate default assumptions `tibble`
assumptions_minimal_example <- function(print=interactive()){
  skel <- "params_scenarios_grid(
    n = c(50, 100), # n, pre group
    mean1 = c(1,0), # mean, group 1
    mean0 = 0       # mean, group 0
  )"

  if(print){
    cat(skel)
  }

  invisible(
    skel |>
      str2expression() |>
      eval()
  )
}
```

The true summary statistics is really simple here, it is just the
difference of group means.

The utility function `params_scenario_grid` expands a dataset in the
following way: the first row is the reference scenario with the first
value of each defined parameter. The subsequent rows contain scenarios
in which each parameter is varied across all given values and all other
parameters are kept at the reference scenario value.

## Analysis

The analysis functions used in `SimDesign` take the arguments
`condition` and `dat`, additional arguments can be passed in the
optional argument `fixed_objects`. When running the simulations the
current line of the Design dataset, containing the parameters for the
current simulation scenario is passed as the `condition` argument. The
simulated dataset is passed as the `dat` argument. Each analysis
function is called once for each simulation replication for each
scenario.

To be able to use hyper-parameters for the analysis methods, CI.RCT.Sim
implements the analysis functions such that they are functions with the
hyper-parameters as arguments that return a function that can be used in
`SimDesign::RunSimulation`. This way one analysis method can easily be
used multiple times with different hyper-parameters in one simulation.

For the minimal example two analysis functions are defined,
`analyse_minimal_example_lm` uses a linear model and
`analyse_minimal_example_t` uses a t-test. The return value for each
method should be a named list to enable summarisation later.

`R/analyse_minimal_example.R`

``` r
#' Analyse dataset from minimal example scenario with linear regression
#'
#' @param ci_level the confidence level for the CIs (defaults to 0.95)
#'
#' @returns an analyse function that returns a list with the elements
#'  * `p` the p-value of the F-test
#'  * `coef` the point estimate
#'  * `ci_lower` the lower CI limit
#'  * `ci_upper` the upper CI limit
#' @export
#'
#' @examples
#' Design <- assumptions_minimal_example() |>
#'   true_summary_statistics_minimal_example()
#'
#' condition <- Design[1, ]
#'
#' dat <- generate_minimal_example(condition)
#'
#' my_analyse_lm <- analyse_minimal_example_lm(ci_level=0.9)
#' my_analyse_lm(condition, dat)
analyse_minimal_example_lm <- function(ci_level=0.95){
  function(condition, dat, fixed_objects = NULL){
    model <- lm(y~group, dat)
    CI <- confint(model, level=ci_level)

    list(
      p = anova(model)$`Pr(>F)`[1],
      coef = coefficients(model)["group"],
      ci_lower = CI["group", 1],
      ci_upper = CI["group", 2]
    )
  }
}



#' Analyse dataset from minimal example scenario with t-test
#'
#' @returns an analyse function that returns a list with the elements
#'  * `p` the p-value
#' @export
#'
#' @examples
#' Design <- assumptions_minimal_example() |>
#'   true_summary_statistics_minimal_example()
#'
#' condition <- Design[1, ]
#'
#' dat <- generate_minimal_example(condition)
#'
#' my_analyse_t  <- analyse_minimal_example_t()
#' my_analyse_t(condition, dat)
analyse_minimal_example_t <- function(){
  function(condition, dat, fixed_objects = NULL){
    ttest <- t.test(dat$y[dat$group==1], dat$y[dat$group==0])

    list(
      p = ttest$p.value
    )
  }
}
```

Some names of the names of the analysis results are directly used by
analysis functions, others can be choosen freely. “reserved” names are:

- `p` the p-value of tests
- `N_pat` number of recruited patients, relevant for sequential or
  adaptive designs
- `N_evt` number of observed events, relevant for time-to event outcomes
  under censoring

## Running Simulations

To actually run the simulations, first a dataset of simulation
parameters is created and true values of summary statistics are
calculated if necessary. A list of analysis functions is created, here
hyper-parameters are set, in this example the coverage of the confidence
intervals.

A list of summary functions is defined, those analysis functions are the
applied to the list of outputs of all replications of each scenario.
Which summary functions are applied to the output of which analysis
functions is matched by name by the function
`SimNPH::create_summary_function`, names can be repeated in the list of
summary functions. The functions
[`SimNPH::summarise_estimator`](https://simnph.github.io/SimNPH/reference/summarise_estimator.html)
and
[`SimNPH::summarise_test`](https://simnph.github.io/SimNPH/reference/summarise_test.html)
provide common summaries for estimators and tests, see the documentation
in the SimNPH package for details.

Finally
[`SimDesign::runSimulation`](http://philchalmers.github.io/SimDesign/reference/runSimulation.md)
is used to conduct the actual simulations. This function takes care of
parallelisation, initialisation of random seeds etc.

`scripts/minimal_example.R`

``` r
library("CI.RCT.Sim")
#> Loading required package: SimDesign

# Define the parameter valued, calculate derived quantities from the parameters
sim_parameters <- assumptions_minimal_example() |>
  true_summary_statistics_minimal_example()

# constants for simulation and aggregation
N_sim <- 1000
alpha <- 0.1

# list of analysis functions
my_analyse <- list(
  lm    = analyse_minimal_example_lm(ci_level=1-alpha),
  ttest = analyse_minimal_example_t()
)

# list of summarisation functions
# those are applied to the results of the analysis functions by the same name,
# repeating one name multiple times is allowed
# summarise_estimator and summarise_test are generic summarisation function
# from SimNPH, but own functions can be defined.
my_summarise <- create_summarise_function(
  lm    = summarise_estimator(est=coef, real=eff_size, lower=ci_lower, upper=ci_upper, null=0),
  ttest = summarise_test(alpha),
  lm    = function(condition, results, fixed_objects = NULL){
    data.frame(mean_ci_width=mean(results$ci_upper - results$ci_lower))
  }
)

# run the simulations using SimDesign::runSimulation
results <- runSimulation(
  sim_parameters,
  replications = N_sim,
  generate = generate_minimal_example,
  analyse = my_analyse,
  summarise = my_summarise
)

# inspect the results
results |>
  subset(select=c(names(sim_parameters), "lm.bias", "lm.sd_est", "ttest.rejection_0.1", "lm.1.mean_ci_width"))
#> # A tibble: 3 × 8
#>       n mean1 mean0 eff_size    lm.bias lm.sd_est ttest.rejection_0.1
#>   <dbl> <dbl> <dbl>    <dbl>      <dbl>     <dbl>               <dbl>
#> 1    50     1     0        1 -0.0022700   0.19598               1    
#> 2   100     1     0        1 -0.0038767   0.14157               1    
#> 3    50     0     0        0  0.010151    0.20418               0.104
#> # ℹ 1 more variable: lm.1.mean_ci_width <dbl>
```
