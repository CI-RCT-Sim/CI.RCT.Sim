#' Analyse dataset from diabetes rescue scenario using MMRM
#'
#' @param ci_level the confidence level for the CIs (defaults to 0.95)
#'
#' @returns an analyse function that returns a list with the elements
#'  * `p` p-value for the treatment effect at the final visit
#'  * `coef` estimated treatment effect at the final visit
#'  * `ci_lower` the lower CI limit
#'  * `ci_upper` the upper CI limit
#'
#' @export
#'
#' @importFrom mmrm mmrm
#' @importFrom tidyselect all_of
#'
#' @examples
#' Design <- assumptions_diabetes_rescue() |>
#'   true_summary_statistics_diabetes_rescue()
#'
#' condition <- Design[1, ]
#'
#' dat <- generate_diabetes_rescue(condition)
#'
#' analyse_mmrm <- analyse_diabetes_rescue_mmrm(ci_level = 0.95)
#' analyse_mmrm(condition, dat)
analyse_diabetes_rescue_mmrm <- function(ci_level = 0.95) {
  function(condition, dat, fixed_objects = NULL) {

    # reshape post-baseline data to long format
    visit_vars <- paste0("y", 1:condition$k)

    long <- tidyr::pivot_longer(
      dat,
      cols = all_of(visit_vars),
      names_to = "visit",
      values_to = "y"
    )

    # ensure factors for MMRM
    long$id <- factor(long$id)
    long$visit <- factor(
      as.integer(sub("y", "", long$visit)),
      levels = 1:condition$k
    )

    # fit MMRM with baseline covariates and interactions
    fit <- mmrm::mmrm(
      y ~ trt * visit +
        y0 * visit +
        age * visit +
        us(visit | id),
      data = long
    )

    # extract treatment effect at final visit
    term <- paste0("trt:visit", condition$k)

    est <- coef(fit)[[term]]
    se  <- sqrt(vcov(fit)[term, term])

    z <- qnorm(1 - (1 - ci_level) / 2)

    list(
      p = 2 * (1 - pnorm(abs(est / se))),
      coef = est,
      ci_lower = est - z * se,
      ci_upper = est + z * se
    )
  }
}
