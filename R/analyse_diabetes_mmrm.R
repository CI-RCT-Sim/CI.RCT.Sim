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
#' @importFrom stats vcov pnorm pt qt
analyse_diabetes_rescue_mmrm <- function(
  ci_level = 0.95,
  strategy = c("treatment_policy", "hypothetical")
) {
  strategy <- match.arg(strategy)

  function(condition, dat, fixed_objects = NULL) {
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

    # Set FINAL visit as reference so trt coefficient = effect at final visit
    long$visit <- stats::relevel(long$visit, ref = as.character(condition$k))

    long$y0 <- baseline[match(long$id, dat$id)]

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
    term <- "trt" # Main effect of treatment
    int_term <- paste0("trt:visit", condition$k) # Effect of the interaction of trt and the final visit

    est <- coef(fit)[term] + coef(fit)[int_term]
    se <- sqrt(vcov(fit)[term, term] + vcov(fit)[int_term, int_term] + 2 * vcov(fit)[term, int_term])

    # Satterthwaite df from mmrm
    # df <- fit$beta_vcov_denom_df[term]
    df <- (summary(fit)$coefficients[term, "df"] + summary(fit)$coefficients[int_term, "df"]) / 2

    tcrit <- qt(1 - (1 - ci_level) / 2, df)

    z <- qnorm(1 - (1 - ci_level) / 2)

    list(
      p = 2 * (1 - pt(abs(est / se), df)),
      coef = est,
      ci_lower = est - tcrit * se,
      ci_upper = est + tcrit * se
    )
  }
}
