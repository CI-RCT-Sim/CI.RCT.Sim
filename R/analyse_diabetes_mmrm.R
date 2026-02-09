#' Analyse dataset from diabetes rescue scenario using MMRM
#'
#' @param ci_level confidence level for the CI (default 0.95)
#'
#' @returns an analyse function that returns a list with
#' * `p` p-value for treatment effect at final visit
#' * `coef` estimated treatment effect at final visit
#' * `ci_lower` lower confidence limit
#' * `ci_upper` upper confidence limit
#'
#' @export
#' @importFrom mmrm mmrm
#' @importFrom stats coef vcov pnorm
#'
#' @examples
#' Design <- assumptions_diabetes_rescue()
#' dat <- generate_diabetes_rescue(Design[1, ])
#' analyse_mmrm <- analyse_diabetes_rescue_mmrm()
#' analyse_mmrm(Design[1, ], dat)
analyse_diabetes_rescue_mmrm <- function(ci_level = 0.95) {
  function(condition, dat, fixed_objects = NULL) {

    # keep baseline
    baseline <- dat$y0

    # reshape post-baseline only
    visit_vars <- paste0("y", 1:condition$k)

    long <- tidyr::pivot_longer(
      dat,
      cols = tidyr::all_of(visit_vars),
      names_to = "visit",
      values_to = "y"
    )

    # convert id and visit to factors for MMRM
    long$id <- factor(long$id)
    long$visit <- factor(as.integer(sub("y", "", long$visit)))
    long$y0 <- baseline[match(long$id, dat$id)]

    # fit MMRM — covariance structure in formula
    fit <- mmrm::mmrm(
      formula = y ~ trt * visit + y0 + us(visit | id),
      data = long
    )

    # treatment effect at final visit
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

