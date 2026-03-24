#' Analyse diabetes endpoint using a Mixed Model for Repeated Measures (MMRM)
#'
#' This function constructs an analysis function that fits a Mixed Model for
#' Repeated Measures (MMRM) to longitudinal continuous outcomes in a diabetes
#' trial setting. The model estimates the treatment effect at the **final visit**
#' while adjusting for baseline and covariates.
#'
#' ## Model specification
#'
#' The fitted model has the following form:
#'
#' \deqn{
#'   y_{ij} = \mu + \beta_{\text{trt}} \cdot \text{trt}_i
#'          + \gamma_j \cdot \text{visit}*j
#'          + (\beta\gamma)*{\text{trt},j} \cdot (\text{trt}_i \times \text{visit}_j)
#'          + \delta_j \cdot (y0_i \times \text{visit}_j)
#'          + \alpha_j \cdot (\text{age}_i \times \text{visit}*j)
#'          + \varepsilon*{ij}
#' }
#'
#' where:
#' * y_ij is the outcome for subject i at visit j,
#' * y0_i is the baseline value,
#' * trt_i is the treatment indicator,
#' * visit_j is treated as a categorical factor,
#' * errors follow a multivariate normal distribution with an
#'   unstructured covariance matrix within subject.
#'
#' In implementation:
#'   y ~ trt * visit + y0 * visit + age * visit + us(visit | id)
#'
#' ## Treatment effect at final visit
#'
#' The factor `visit` is re-leveled so that the **final visit is the reference**.
#' With this parametrization, the coefficient for `trt` directly represents the
#' treatment effect at the final visit.
#'
#' ## Handling of rescue medication
#'
#' strategy = "treatment_policy":
#'   All observed data are used.
#'
#' strategy = "hypothetical":
#'   For a subject with rescue_start = r, all outcomes with visit >= r are set
#'   to NA and ignored by the MMRM (MAR assumption).
#'
#' @param ci_level Confidence level (default 0.95)
#' @param strategy "treatment_policy" or "hypothetical"
#'
#' @returns A function returning a list with:
#'   p, coef, ci_lower, ci_upper, converged
#'
#' @export
#'
#' @examples
#' setting <- diabetes_scenario()[1, ] |>
#'   diabetes_scenario_set_truevalues()
#'
#' dat <- generate_diabetes(setting)
#'
#' # Treatment policy estimand
#' analyse_diabetes_mmrm()(setting, dat)
#'
#' # Hypothetical estimand (censor after rescue)
#' analyse_diabetes_mmrm(strategy = "hypothetical")(setting, dat)
analyse_diabetes_mmrm <- function(
    ci_level = 0.95,
    strategy = c("treatment_policy", "hypothetical")
) {
  strategy <- match.arg(strategy)

  function(condition, dat, fixed_objects = NULL) {

    term <- "trt"

    safe_result <- list(
      p = NA_real_,
      coef = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_,
      converged = FALSE
    )

    dat_work <- dat
    baseline <- dat$y0

    # --- Hypothetical strategy ---
    if (strategy == "hypothetical") {

      has_rescue <- !is.na(dat_work$rescue_start) &
        dat_work$rescue_start <= condition$k

      for (v in seq_len(condition$k)) {
        var <- paste0("y", v)

        dat_work[[var]] <- ifelse(
          has_rescue & dat_work$rescue_start <= v,
          NA,
          dat_work[[var]]
        )
      }
    }

    # --- Check required columns ---
    visit_vars <- paste0("y", seq_len(condition$k))
    if (!all(visit_vars %in% names(dat_work))) {
      return(safe_result)
    }

    # --- Reshape ---
    long <- tryCatch(
      tidyr::pivot_longer(
        dat_work,
        cols = tidyselect::all_of(visit_vars),
        names_to = "visit",
        values_to = "y"
      ),
      error = function(e) NULL
    )

    if (is.null(long)) return(safe_result)

    # --- Prepare variables ---
    long$id <- factor(long$id)

    long$visit <- factor(
      as.integer(sub("y", "", long$visit)),
      levels = seq_len(condition$k)
    )

    # Final visit as reference
    long$visit <- stats::relevel(long$visit, ref = as.character(condition$k))

    # Baseline mapping
    long$y0 <- baseline[match(long$id, dat$id)]

    # --- Fit model ---
    fit <- tryCatch(
      mmrm::mmrm(
        y ~ trt * visit +
          y0 * visit +
          age * visit +
          mmrm::us(visit | id),
        data = long
      ),
      error = function(e) NULL
    )

    if (is.null(fit)) return(safe_result)

    # --- Extract components safely ---
    coefs <- tryCatch(stats::coef(fit), error = function(e) NULL)
    vc <- tryCatch(stats::vcov(fit), error = function(e) NULL)
    summ <- tryCatch(summary(fit), error = function(e) NULL)

    if (
      is.null(coefs) ||
      is.null(vc) ||
      is.null(summ) ||
      !(term %in% names(coefs)) ||
      !(term %in% rownames(vc)) ||
      !(term %in% rownames(summ$coefficients))
    ) {
      return(safe_result)
    }

    est <- coefs[term]

    se <- tryCatch(sqrt(vc[term, term]), error = function(e) NA_real_)
    df <- summ$coefficients[term, "df"]

    if (is.na(est) || is.na(se) || is.na(df) || se <= 0 || df <= 0) {
      return(safe_result)
    }

    tcrit <- stats::qt(1 - (1 - ci_level) / 2, df)

    list(
      p = 2 * (1 - stats::pt(abs(est / se), df)),
      coef = est,
      ci_lower = est - tcrit * se,
      ci_upper = est + tcrit * se,
      converged = TRUE
    )

  }
}
