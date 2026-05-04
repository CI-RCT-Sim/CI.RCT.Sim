#' Analyse data set with Multiple imputation
#'
#' Creates an analysis function for the diabetes rescue simulation.
#' Depending on the chosen estimand, either a hypothetical strategy
#' with multiple imputation after rescue initiation is applied,
#' or a treatment policy strategy using observed data only.
#'
#' ## Multiple Imputation (Hypothetical Strategy)
#'
#' Multiple imputation is used as a comparator method for the hypothetical
#' estimand, targeting the scenario in which rescue medication had not been
#' available.
#'
#' ### Handling of rescue medication
#'
#' For subjects initiating rescue medication at visit \eqn{s}, all post-rescue
#' outcomes are set to missing:
#'
#' \deqn{
#' y_{i,j} = \text{NA} \quad \text{for all } j > s
#' }
#'
#' The outcome at the rescue visit itself (\eqn{j = s}) is retained.
#'
#' This ensures consistency with the hypothetical estimand and aligns with the
#' data-generating mechanism.
#'
#' ### Imputation model (Protocol-aligned)
#'
#' The imputation is performed using the \pkg{mice} package.
#'
#' * HbA1c at each visit is imputed if missing.
#' * The imputation model includes:
#'   - all available HbA1c measurements,
#'   - baseline HbA1c (\eqn{y_0}),
#'   - age at baseline.
#'
#' The imputation is performed **separately within each treatment group**,
#' as specified in the protocol.
#'
#' A total of 10 imputations are generated.
#'
#' ### Analysis model
#'
#' For each imputed dataset, an ANCOVA model is fitted:
#'
#' \deqn{
#' \Delta y_i = \beta_0 + \beta_1 \cdot \text{trt}_i + \beta_2 \cdot \text{age}_i + \beta_3 \cdot y_{0,i} + \varepsilon_i
#' }
#'
#' where:
#' * \eqn{\Delta y_i} is the change from baseline in HbA1c
#' * \eqn{y_{0,i}} is baseline HbA1c
#'
#' ### Pooling
#'
#' Estimates from the imputed datasets are pooled using Rubin’s rules via
#' \code{mice::pool}.
#'
#' ### Reproducibility
#'
#' Randomness from the multiple imputation procedure is controlled using
#' \code{withr::with_seed()}, ensuring reproducibility without affecting the
#' global random number generator state. This avoids unintended interference
#' with simulation data generation.
#'
#' @param strategy Either `"hypothetical"` or `"treatment_policy"`.
#' @param m Number of imputations (used for hypothetical only).
#' @param maxit Maximum number of MICE iterations.
#' @param ci_level confidence level for the CI (default 0.95)
#' @param seed Random seed for imputation (used locally via withr).
#'
#' @return A function that, when called with `condition` and `dat`, returns a list with:
#' * `coef` coefficient for `trt`
#' * `p` p-value for coef
#' * `ci_lower` lower bound of `ci-level`% confidence interval for coef
#' * `ci_upper` upper bound of `ci-level`% confidence interval for coef
#'
#' @importFrom mice mice make.method make.predictorMatrix complete as.mira pool
#' @importFrom dplyr filter select bind_rows all_of
#' @importFrom stats lm confint coef
#' @importFrom withr with_seed
#' @export
#'
#' @examples
#' \donttest{
#' Design <- diabetes_scenario()[1, ] |>
#'   diabetes_scenario_set_truevalues()
#'
#' dat <- generate_diabetes(Design)
#'
#' analyse_diabetes_mi(strategy = "treatment_policy")(Design, dat)
#' analyse_diabetes_mi(strategy = "hypothetical")(Design, dat)
#'
#' dat <- generate_diabetes(Design)
#'
#' ## ----------------------------
#' ## Treatment policy estimand
#' ## ----------------------------
#' res_tp <- analyse_diabetes_mi(
#'   strategy = "treatment_policy"
#' )(Design, dat)
#'
#' res_tp
#'
#' ## ----------------------------
#' ## Hypothetical estimand
#' ## (censor after rescue + MI)
#' ## ----------------------------
#' res_hyp <- analyse_diabetes_mi(
#'   strategy = "hypothetical"
#' )(Design, dat)
#'
#' res_hyp
#'
#' ## Compare results
#' c(
#'   treatment_policy = res_tp$coef,
#'   hypothetical     = res_hyp$coef
#' )
#' }
analyse_diabetes_mi <- function(
    strategy = c("hypothetical", "treatment_policy"),
    m = 10,
    maxit = 10,
    ci_level = 0.95,
    seed = 123
) {
  strategy <- match.arg(strategy)

  function(condition, dat, fixed_objects = NULL) {

    k <- condition$k

    ############################################################
    # TREATMENT POLICY (unchanged, protocol OK)
    ############################################################
    if (strategy == "treatment_policy") {

      dat$chg <- dat[[paste0("y", k)]] - dat$y0

      fit <- lm(chg ~ trt + age + y0, data = dat)

      ci <- confint(fit, level = ci_level)["trt", ]

      return(list(
        coef     = coef(fit)["trt"],
        p        = summary(fit)$coefficients["trt", "Pr(>|t|)"],
        ci_lower = ci[1],
        ci_upper = ci[2]
      ))
    }

    ############################################################
    # HYPOTHETICAL ESTIMAND (PROTOCOL-COMPLIANT MI)
    ############################################################

    dat_hyp <- dat

    ## 1. Censor post-rescue outcomes
    for (i in seq_len(nrow(dat_hyp))) {
      rs <- dat_hyp$rescue_start[i]
      if (!is.na(rs) && rs < k) {
        dat_hyp[i, paste0("y", (rs + 1):k)] <- NA
      }
    }

    ############################################################
    # VARIABLES
    ############################################################
    vars_y <- paste0("y", 0:k)
    vars_R <- if (k > 1) paste0("R", 1:(k - 1)) else character(0)
    vars_all <- c(vars_y, vars_R, "age", "trt")

    dat_mi <- dat_hyp[, vars_all]

    ############################################################
    # Ensure correct types
    ############################################################
    if (length(vars_R) > 0) {
      dat_mi[vars_R] <- lapply(dat_mi[vars_R], function(x)
        factor(x, levels = c(0, 1))
      )
    }

    ############################################################
    # MICE METHOD SPECIFICATION (PROTOCOL-COMPLIANT)
    ############################################################
    meth <- mice::make.method(dat_mi)

    ## HbA1c: predictive mean matching
    meth[vars_y] <- "pmm"

    ## Rescue indicators: logistic regression
    if (length(vars_R) > 0) {
      meth[vars_R] <- "logreg"
    }

    ## No imputation for fully observed covariates
    meth[c("age", "trt")] <- ""

    ############################################################
    # PASSIVE CONSTRAINT: RESCUE ONCE STARTED STAYS ON
    ############################################################
    if (length(vars_R) > 1) {
      for (j in 2:length(vars_R)) {
        r_prev <- vars_R[j - 1]
        r_curr <- vars_R[j]

        meth[r_curr] <- paste0("~ I(pmax(", r_prev, ", ", r_curr, "))")
      }
    }

    ############################################################
    # PREDICTOR MATRIX (protocol-aligned)
    ############################################################
    pred <- mice::make.predictorMatrix(dat_mi)
    pred[,] <- 0

    ## HbA1c depends on:
    pred[vars_y, c(vars_y, vars_R, "age")] <- 1
    diag(pred[vars_y, vars_y]) <- 0

    ## Rescue depends on:
    if (length(vars_R) > 0) {
      pred[vars_R, c(vars_y, vars_R, "age")] <- 1
      diag(pred[vars_R, vars_R]) <- 0
    }

    ## treatment never used in imputation model
    pred[, "trt"] <- 0

    ############################################################
    # IMPUTATION (TRUE mice, stable)
    ############################################################
    imp <- withr::with_seed(
      seed,
      mice::mice(
        dat_mi,
        m = m,
        maxit = maxit,
        method = meth,
        predictorMatrix = pred,
        printFlag = FALSE,
        ridge = 1e-6
      )
    )

    ############################################################
    # ANALYSIS (ANCOVA per imputation)
    ############################################################
    fit_models <- lapply(1:m, function(i) {

      d <- mice::complete(imp, i)

      d$chg <- d[[paste0("y", k)]] - d$y0

      lm(chg ~ trt + age + y0, data = d)
    })

    imp_obj <- mice::as.mira(fit_models)
    pooled <- mice::pool(imp_obj)

    sum_pooled <- summary(pooled, conf.int = TRUE, conf.level = ci_level)

    trt_row <- sum_pooled[sum_pooled$term == "trt", ]

    list(
      coef     = trt_row$estimate,
      p        = trt_row$p.value,
      ci_lower = trt_row$`2.5 %`,
      ci_upper = trt_row$`97.5 %`
    )
  }
}
