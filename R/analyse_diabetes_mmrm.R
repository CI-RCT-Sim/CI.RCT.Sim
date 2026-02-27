#' Analyse dataset from diabetes rescue scenario using MMRM
#'
#' Supports two intercurrent event strategies:
#' * "treatment_policy" (default) uses all observed data
#' * "hypothetical" sets post-rescue data to missing
#'
#' @param ci_level Confidence level for CIs (default 0.95)
#' @param strategy "treatment_policy" or "hypothetical"
#'
#' @return A function that takes (condition, dat) and returns a list:
#' * `p` p-value for treatment effect at final visit
#' * `coef` estimated treatment effect at final visit
#' * `ci_lower` lower CI
#' * `ci_upper` upper CI
#'
#' @export
#' @importFrom mmrm mmrm
analyse_diabetes_rescue_mmrm <- function(ci_level = 0.95, strategy = c("treatment_policy", "hypothetical")) {

  strategy <- match.arg(strategy)

  function(condition, dat, fixed_objects = NULL) {

    dat_work <- dat

    if(strategy == "hypothetical") {
      # set post-rescue measurements to NA
      visit_vars <- paste0("y", 1:condition$k)
      for(i in seq_len(nrow(dat_work))) {
        start <- dat_work$rescue_start[i]
        if(!is.na(start)) {
          # set all visits after rescue_start to NA
          dat_work[i, paste0("y", start:condition$k)] <- NA
        }
      }
    }

    # reshape to long format
    visit_vars <- paste0("y", 1:condition$k)
    long <- tidyr::pivot_longer(
      dat_work,
      cols = all_of(visit_vars),
      names_to = "visit",
      values_to = "y"
    )

    # ensure factors for MMRM
    long$id <- factor(long$id)
    long$visit <- factor(as.integer(sub("y", "", long$visit)), levels = 1:condition$k)

    # fit MMRM
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

#How to use
# Design <- assumptions_diabetes_rescue() |>
#   true_summary_statistics_diabetes_rescue()
#
# condition <- Design[1, ]
# dat <- generate_diabetes_rescue(condition)
#
# # Treatment policy strategy (all observed data)
# analyse_mmrm_tp <- analyse_diabetes_rescue_mmrm(strategy = "treatment_policy")
# analyse_mmrm_tp(condition, dat)
#
# # Hypothetical strategy (set post-rescue measurements to NA)
# analyse_mmrm_hyp <- analyse_diabetes_rescue_mmrm(strategy = "hypothetical")
# analyse_mmrm_hyp(condition, dat)
