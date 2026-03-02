#' Analyse diabetes endpoint using MMRM
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
analyse_diabetes_rescue_mmrm <- function(
  ci_level = 0.95,
  strategy = c("treatment_policy", "hypothetical")
) {
  strategy <- match.arg(strategy)

  function(condition, dat, fixed_objects = NULL) {
    # keep baseline
    baseline <- dat$y0
    dat_work <- dat
    baseline <- dat$y0

    # --- Hypothetical strategy: remove post-rescue values ---
    if (strategy == "hypothetical") {
      for (i in seq_len(nrow(dat_work))) {
        start <- dat_work$rescue_start[i]
        if (!is.na(start)) {
          dat_work[i, paste0("y", start:condition$k)] <- NA
        }
      }
    }

    # --- Reshape to long ---
    visit_vars <- paste0("y", 1:condition$k)

    long <- tidyr::pivot_longer(
      dat_work,
      cols = tidyselect::all_of(visit_vars),
      names_to = "visit",
      values_to = "y"
    )

    # --- Prepare variables ---
    long$id <- factor(long$id)
    long$visit <- factor(
      as.integer(sub("y", "", long$visit)),
      levels = 1:condition$k
    )
    long$y0 <- baseline[match(long$id, dat$id)]

    # --- Fit MMRM ---
    fit <- mmrm::mmrm(
      formula = y ~ trt * visit + y0 + us(visit | id),
      data = long
    )

    # --- Extract treatment effect at final visit ---
    term <- paste0("trt:visit", condition$k)

    est <- coef(fit)[[term]]
    se <- sqrt(vcov(fit)[term, term])

    z <- qnorm(1 - (1 - ci_level) / 2)

    list(
      p = 2 * (1 - pnorm(abs(est / se))),
      coef = est,
      ci_lower = est - z * se,
      ci_upper = est + z * se
    )
  }
}
