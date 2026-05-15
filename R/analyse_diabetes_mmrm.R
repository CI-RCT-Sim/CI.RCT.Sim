#' Analyse diabetes endpoint using a Mixed Model for Repeated Measures (MMRM)
#'
#' This function fits a Mixed Model for Repeated Measures (MMRM) to estimate
#' the treatment effect at the final visit in a longitudinal clinical trial setting.
#'
#' @details
#' ## Model specification
#'
#' The following linear mixed model is fitted:
#'
#' \deqn{
#' y_{ij} =
#' \beta_0
#' + \beta_1 \cdot \text{trt}_i
#' + \sum_{j} \gamma_j \cdot \mathbb{1}(\text{visit}_j)
#' + \sum_{j} \delta_j \cdot \left(\text{trt}_i \times \mathbb{1}(\text{visit}_j)\right)
#' + \sum_{j} \alpha_j \cdot y_{0,i} \cdot \mathbb{1}(\text{visit}_j)
#' + \sum_{j} \eta_j \cdot \text{age}_i \cdot \mathbb{1}(\text{visit}_j)
#' + \varepsilon_{ij}
#' }
#'
#' where:
#' * \eqn{y_{ij}} is the outcome for subject \eqn{i} at visit \eqn{j}
#' * \eqn{\text{trt}_i} is the treatment indicator
#' * \eqn{\text{visit}_j} is a categorical visit factor
#' * \eqn{y_{0,i}} is baseline value
#' * \eqn{\text{age}_i} is a baseline covariate
#'
#' The within-subject covariance is modeled using an unstructured covariance:
#'
#' \deqn{
#' \varepsilon_i \sim \mathcal{N}(0, \Sigma)
#' }
#'
#' where \eqn{\Sigma} is an unstructured covariance matrix across visits.
#'
#' ## Estimand and interpretation
#'
#' The visit factor is re-leveled such that the final visit (visit = k)
#' is the reference category. Therefore, the main treatment coefficient
#' corresponds directly to the treatment effect at the final visit.
#'
#' ## Handling of intercurrent events (rescue medication)
#'
#' Two strategies are supported:
#'
#' ### Treatment policy
#'
#' All observed post-baseline data are used regardless of rescue medication.
#'
#' ### Hypothetical strategy
#'
#' For subjects who initiate rescue medication at visit \eqn{s}, all outcomes
#' strictly after rescue are set to missing:
#'
#' \deqn{
#' y_{i,j} = \mathrm{NA} \quad \text{for all } j > s
#' }
#'
#' Subjects without rescue (or with rescue after the final visit)
#' remain unchanged.
#'
#' This targets a hypothetical estimand corresponding to outcomes that
#' would have been observed had rescue medication not been initiated.
#'
#' ## Inference
#'
#' Treatment effects are estimated using restricted maximum likelihood (REML).
#'
#' Degrees of freedom are computed using Satterthwaite approximation.
#'
#' ### Hypothesis testing
#'
#' A one-sided test is used for the treatment effect:
#'
#' \deqn{
#' H_0: \beta_{trt} \ge 0
#' \quad \text{vs} \quad
#' H_A: \beta_{trt} < 0
#' }
#'
#' since a reduction in HbA1c (negative change) is considered beneficial.
#'
#' Two-sided confidence intervals are retained and reported.
#'
#' ## Covariance structures
#'
#' The model attempts the following covariance structures in order:
#' unstructured (`us`), compound symmetry (`cs`), and diagonal (`diag`),
#' falling back if convergence fails.
#'
#' @importFrom mmrm mmrm
#' @importFrom stats vcov pt qt relevel
#'
#' @param ci_level Confidence level for the two-sided confidence interval (default 0.95)
#' @param strategy Strategy for handling rescue medication:
#'   * `"treatment_policy"`
#'   * `"hypothetical"`
#'
#' @returns A function that, when called with `condition` and `dat`, returns a list with:
#' * `coef` estimated treatment effect at final visit
#' * `p` one-sided p-value for \eqn{\beta_{trt} < 0}
#' * `ci_lower` lower bound of the two-sided confidence interval
#' * `ci_upper` upper bound of the two-sided confidence interval
#' * `converged` logical indicating whether the model converged
#' * `covariance` covariance structure used (`us`, `cs`, or `diag`)
#' * `fallback` logical indicating whether fallback covariance was used
#'
#' @examples
#' \donttest{
#' Design <- diabetes_scenario()[1, ] |>
#'   diabetes_scenario_set_truevalues()
#'
#' dat <- generate_diabetes(Design)
#'
#' analyse_diabetes_mmrm(strategy = "treatment_policy")(Design, dat)
#' analyse_diabetes_mmrm(strategy = "hypothetical")(Design, dat)
#' }
#'
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

    if (is.null(long)) {
      return(safe_result)
    }

    long$id <- factor(long$id)

    long$visit <- factor(
      as.integer(sub("y", "", long$visit)),
      levels = seq_len(condition$k)
    )

    # Reference = final visit
    long$visit <- relevel(long$visit, ref = as.character(condition$k))

    long$y0 <- baseline[match(long$id, dat$id)]

    # ============================================================
    # FIT WITH COVARIANCE FALLBACK
    # ============================================================
    fit_mmrm <- function(cov_type) {
      formula_str <- switch(cov_type,
        "us" = y ~ trt * visit + y0 * visit + age * visit + us(visit | id),
        "cs" = y ~ trt * visit + y0 * visit + age * visit + cs(visit | id),
        "diag" = y ~ trt * visit + y0 * visit + age * visit + diag(visit | id)
      )

      tryCatch(
        mmrm(formula_str, data = long),
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

    if (is.null(fit)) {
      return(safe_result)
    }

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

    tcrit <- qt(1 - (1 - ci_level) / 2, df)

    list(
      p = 2 * (1 - pt(abs(est / se), df)),
      coef = est,
      ci_lower = est - tcrit * se,
      ci_upper = est + tcrit * se,
      converged = TRUE,
      covariance = covariance_used,
      fallback = safe_result$fallback
    )
  }
}
