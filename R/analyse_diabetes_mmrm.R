#' Analyse diabetes endpoint using MMRM
#'
#' @param ci_level confidence level for the CI (default 0.95)
#' @param strategy strategy to handle rescue medication: either "treatment_policy" (use all observed data) or "hypothetical" (set post-rescue values to missing and ignore them in the analysis)
#'
#' @returns an analyse function that returns a list with
#' * `p` p-value for treatment effect at final visit
#' * `coef` estimated treatment effect at final visit
#' * `ci_lower` lower confidence limit
#' * `ci_upper` upper confidence limit
#'
#' @export
#' @importFrom mmrm mmrm
#' @importFrom stats vcov pnorm pt qt relevel
#'
#' @examples
#' setting <- diabetes_scenario()[1, ] |> diabetes_scenario_set_truevalues()
#' dat <- generate_diabetes(setting)
#' analyse_diabetes_mmrm()(setting, dat)
#'
analyse_diabetes_mmrm <- function(
    ci_level = 0.95,
    strategy = c("treatment_policy", "hypothetical")
) {
  strategy <- match.arg(strategy)

  function(condition, dat, fixed_objects = NULL) {
    dat_work <- dat
    baseline <- dat$y0

    # --- Hypothetical strategy: remove post-rescue values ---
    if (strategy == "hypothetical") {
      dat_work$rescue_start[is.na(dat_work$rescue_start)] <- condition$k+2
      for (i in seq_len(nrow(dat_work))) {
        start <- dat_work$rescue_start[i]
        if (start < (condition$k+2)) {
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
    long$visit <- relevel(long$visit, ref = as.character(condition$k))

    long$y0 <- baseline[match(long$id, dat$id)]

    # --- Fit MMRM ---
    fit <- mmrm(
      y ~ trt * visit +
        y0 * visit +
        age * visit +
        us(visit | id),
      data = long
    )

    # --- Extract treatment effect at final visit ---
    term <- "trt" # Main effect of treatment

    est <- coef(fit)[term]
    se <- sqrt(vcov(fit)[term, term])

    # Satterthwaite df from mmrm
    df <- summary(fit)$coefficients[term, "df"]

    tcrit <- qt(1 - (1 - ci_level) / 2, df)

    list(
      p = 2 * (1 - pt(abs(est / se), df)),
      coef = est,
      ci_lower = est - tcrit * se,
      ci_upper = est + tcrit * se
    )
  }
}
