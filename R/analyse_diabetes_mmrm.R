#' Analyse diabetes endpoint using a Mixed Model for Repeated Measures (MMRM)
#'
#' This function fits a Mixed Model for Repeated Measures (MMRM) to estimate
#' the treatment effect at the final visit in a longitudinal clinical trial setting.
#'
#' ## Model specification
#'
#' The following linear model is fitted:
#'
#' \deqn{
#' y_{ij} = \beta_0
#' + \beta_1 \cdot \text{trt}_i
#' + \gamma_j \cdot \text{visit}_j
#' + \delta_j \cdot (\text{trt}i \times \text{visit}j)
#' + \alpha_j \cdot y{0,i}
#' + \eta_j \cdot \text{age}i
#' + \varepsilon{ij}
#' }
#'
#' where:
#' * \eqn{y{ij}} is the outcome for subject \eqn{i} at visit \eqn{j}
#' * \eqn{\text{trt}_i} is the treatment indicator
#' * \eqn{\text{visit}j} is a categorical visit effect
#' * \eqn{y{0,i}} is the baseline value
#' * \eqn{\text{age}i} is a baseline covariate
#'
#' The within-subject covariance is modeled using an unstructured covariance matrix:
#'
#' \deqn{
#' \varepsilon_i \sim \mathcal{N}(0, \Sigma)
#' }
#'
#' with \eqn{\Sigma} fully unstructured across visits.
#'
#' ## Interpretation of treatment effect
#'
#' The visit factor is re-leveled such that the final visit (visit = k) is the
#' reference category. As a result, the main effect coefficient for trt
#' corresponds directly to the treatment effect at the final visit.
#'
#' Without this releveling, the treatment effect would correspond to the reference
#' visit (typically the first visit), and additional contrasts would be required
#' to obtain the effect at the final visit.
#'
#' ## Handling of intercurrent events (rescue medication)
#'
#' Two strategies are supported:
#'
#' * "treatment_policy":
#' All observed post-baseline data are used regardless of rescue medication.
#'
#' * "hypothetical":
#' For subjects who initiate rescue medication during the study
#' (i.e. rescue_start <= k), all outcomes from the rescue visit onward are
#' set to missing. This results in a monotone missingness pattern, where:
#'
#' \deqn{
#' y{i, j} = \text{NA for all } j \geq \text{rescue_start}_i
#' }
#'
#' Subjects without rescue (i.e. rescue_start > k or NA) remain unchanged.
#'
#' This approach targets a hypothetical estimand corresponding to outcomes
#' that would have been observed had rescue medication not been initiated.
#'
#' ## Inference
#'
#' Treatment effects are estimated using restricted maximum likelihood (REML),
#' and inference is based on Satterthwaite approximations for the degrees of freedom.
#'
#' Confidence intervals are constructed using a t-distribution with the estimated
#' degrees of freedom.
#'
#' @param ci_level Confidence level for the confidence interval (default 0.95)
#' @param strategy Strategy for handling rescue medication:
#' "treatment_policy" or "hypothetical"
#'
#' @returns A function returning a list with:
#' * p: p-value for the treatment effect at the final visit
#' * coef: estimated treatment effect at the final visit
#' * ci_lower: lower confidence limit
#' * ci_upper: upper confidence limit
#'
#' @examples
# setting <- diabetes_scenario()[1, ] |>
# diabetes_scenario_set_truevalues()
#
# dat <- generate_diabetes(setting)
#
#  #Treatment policy estimand
# analyse_diabetes_mmrm()(setting, dat)
#
#  #Hypothetical estimand (censor after rescue)
# analyse_diabetes_mmrm(strategy = "hypothetical")(setting, dat)
#' @export
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
      converged = FALSE,
      covariance = NA_character_,
      fallback = FALSE
    )

    dat_work <- dat
    baseline <- dat$y0

    # ============================================================
    # HYPOTHETICAL STRATEGY (HARMONIZED)
    # ============================================================
    if (strategy == "hypothetical") {

      # encode no-rescue consistently with generator
      dat_work$rescue_start[is.na(dat_work$rescue_start)] <- condition$k + 2

      for (i in seq_len(nrow(dat_work))) {

        start <- dat_work$rescue_start[i]

        # CENSOR ONLY AFTER rescue visit
        if (start < condition$k) {

          post_visits <- (start + 1):condition$k

          dat_work[i, paste0("y", post_visits)] <- NA

        }
      }
    }

    # ============================================================
    # RESHAPE TO LONG
    # ============================================================
    visit_vars <- paste0("y", seq_len(condition$k))

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

    long$id <- factor(long$id)

    long$visit <- factor(
      as.integer(sub("y", "", long$visit)),
      levels = seq_len(condition$k)
    )

    # Reference = final visit
    long$visit <- stats::relevel(long$visit, ref = as.character(condition$k))

    long$y0 <- baseline[match(long$id, dat$id)]

    # ============================================================
    # FIT WITH COVARIANCE FALLBACK
    # ============================================================
    fit_mmrm <- function(cov_type) {

      formula_str <- switch(
        cov_type,
        "us" = y ~ trt * visit + y0 * visit + age * visit + us(visit | id),
        "cs" = y ~ trt * visit + y0 * visit + age * visit + cs(visit | id),
        "diag" = y ~ trt * visit + y0 * visit + age * visit + diag(visit | id)
      )

      tryCatch(
        mmrm::mmrm(formula_str, data = long),
        error = function(e) NULL
      )
    }

    fit <- fit_mmrm("us")
    covariance_used <- "us"

    if (is.null(fit)) {
      fit <- fit_mmrm("cs")
      covariance_used <- "cs"
      safe_result$fallback <- TRUE
    }

    if (is.null(fit)) {
      fit <- fit_mmrm("diag")
      covariance_used <- "diag"
      safe_result$fallback <- TRUE
    }

    if (is.null(fit)) return(safe_result)

    safe_result$covariance <- covariance_used

    # ============================================================
    # SAFE EXTRACTION
    # ============================================================
    coefs <- tryCatch(coef(fit), error = function(e) NULL)
    vc <- tryCatch(vcov(fit), error = function(e) NULL)
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

    se <- tryCatch(
      sqrt(vc[term, term]),
      error = function(e) NA_real_
    )

    df <- summ$coefficients[term, "df"]

    if (
      is.na(est) || is.na(se) || is.na(df) ||
      se <= 0 || df <= 0 ||
      !is.finite(se) || !is.finite(df)
    ) {
      return(safe_result)
    }

    tcrit <- stats::qt(1 - (1 - ci_level) / 2, df)

    list(
      p = 2 * (1 - stats::pt(abs(est / se), df)),
      coef = est,
      ci_lower = est - tcrit * se,
      ci_upper = est + tcrit * se,
      converged = TRUE,
      covariance = covariance_used,
      fallback = safe_result$fallback
    )
  }
}
