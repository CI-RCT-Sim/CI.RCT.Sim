#' Analyse Dataset with G-estimation
#'
#' @param level confidence level for CI computation
#' @param alternative alternative hypothesis for the tests "two.sided" or "one.sided"
#'
#' @return an analyse function that returns a list with the elements
#' * `p` p value of the score test (two.sided) or the Wald test (one.sided)
#' * `alternative` the alternative used
#' * `coef` coefficient for `trt`
#' * `hr` hazard ratio for `trt`
#' * `hr_lower` lower 95% confidence interval boundary for the hazard ratio for `trt`
#' * `hr_upper`lower 95% confidence interval boundary for the hazard ratio for `trt`
#' * `CI_level` the CI level used
#' * `N_pat` number of patients
#' * `N_evt` number of events
#'
#' @export
#'
#' @importFrom gfoRmula gformula_continuous_eof
#' @details
#' `alternative` can be "two.sided" for a two sided test of equality of the
#' summary statistic or "one.sided" for a one sided test testing H0: treatment
#' has equal or shorter survival than control vs. H1 treatment has longer
#' survival than control.
#'
#' @examples
#' Design <- assumptions_diabetes_rescue() |>
#'   true_summary_statistics_diabetes_rescue()
#'
#' condition <- Design[1, ]
#'
#' dat <- generate_diabetes_rescue(condition)
#'
#' analyse_gestimation  <- analyse_gestimation()
#' analyse_gestimation (condition, dat)


analyse_gestimation <- function(level = 0.95, alternative = "two.sided") {
  stopifnot(alternative %in% c("two.sided", "one.sided"))

  # What still remains to be done is
  # write function
  # documentation
  # determine if the "alternative" statement can be incorporated

  function(condition, dat, fixed_objects = NULL) {

  #  model the full joint distribution of age at baseline, HbA1c trajectories
    # and administration of rescue over time
    # by fitting respective models for data at each visit
    # conditional on data at the previous visits.

    # model 1: y ∼ lag1 y + age and
    # model 2: R ∼ y + age
    # A restriction will be implemented
    # such that the rescue indicator Rij equals 1 if Rij−1 = 1
    # to take into account that rescue medication is taken continuously once started.

    #The fitted models are used to sample for each treatment group
    # a number of trajectories equal to the number of patients originally included,
    # starting from baseline covariate values and treatment group,
    # albeit under the assumption that rescue medication is not possible.
    # The mean difference of the outcome between groups is estimated from the sampled data.

    # bootstrap samples 500 (may be reduced to 200):
    # nsamples=500
    # Bootstrap 95% confidence intervals for the mean difference
    # will be calculated using the percentile method:
    #ci_method="percentile"
      CI <- confint(model, level=level)

      list(
        p = "NA", # No formal p-value will be calculated for this method
        coef = summary(model)$coefficients["trt", "Estimate"],
        ci_lower = CI["trt", 1],
        ci_upper = CI["trt", 2]
      )
}}


