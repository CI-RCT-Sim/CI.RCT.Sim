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
# Design <- diabetes_scenario()[1, ] |>
#   diabetes_scenario_set_truevalues()
#
# dat <- generate_diabetes(Design)
#
# analyse_diabetes_mi(strategy = "treatment_policy")(Design, dat)
# analyse_diabetes_mi(strategy = "hypothetical")(Design, dat)
#
# dat <- generate_diabetes(Design)
#
# ## ----------------------------
# ## Treatment policy estimand
# ## ----------------------------
# res_tp <- analyse_diabetes_mi(
#   strategy = "treatment_policy"
# )(Design, dat)
#
# res_tp
#
# ## ----------------------------
# ## Hypothetical estimand
# ## (censor after rescue + MI)
# ## ----------------------------
# res_hyp <- analyse_diabetes_mi(
#   strategy = "hypothetical"
# )(Design, dat)
#
# res_hyp
#
# ## Compare results
# c(
#   treatment_policy = res_tp$coef,
#   hypothetical     = res_hyp$coef
# )
#' }
analyse_diabetes_mi <- function(
    strategy = c("hypothetical", "treatment_policy"),
    m = 10,
    maxit = 20,
    ci_level = 0.95,
    seed = 123
) {

  strategy <- match.arg(strategy)

  function(condition, dat, fixed_objects = NULL) {

    k <- condition$k

    vars_y <- paste0("y", 0:k)
    vars_R <- if (k > 1) paste0("R", 1:(k - 1)) else character(0)

    dat_hyp <- dat

    ############################################################
    # OPTIONAL: scaling for numerical stability
    ############################################################
    dat_hyp$age <- as.numeric(scale(dat_hyp$age))
    dat_hyp$y0  <- as.numeric(scale(dat_hyp$y0))

    ############################################################
    # HYPOTHETICAL STRATEGY
    ############################################################
    if (strategy == "hypothetical") {

      for (i in seq_len(nrow(dat_hyp))) {

        rs <- dat_hyp$rescue_start[i]

        if (!is.na(rs) && rs < k) {
          dat_hyp[i, paste0("y", (rs + 1):k)] <- NA
        }
      }
    }

    ############################################################
    # VARIABLES
    ############################################################
    vars_imp <- c(vars_y, vars_R, "age", "trt")

    ############################################################
    # METHODS
    ############################################################
    meth <- mice::make.method(dat_hyp[vars_imp])
    meth[vars_y] <- "pmm"

    if (length(vars_R) > 0) {
      meth[vars_R] <- ""
    }

    meth[c("age", "trt")] <- ""

    ############################################################
    # PREDICTOR MATRIX (ENHANCED LAG STRUCTURE)
    ############################################################
    pred <- mice::make.predictorMatrix(dat_hyp[vars_imp])
    pred[,] <- 0

    # baseline + covariates
    pred[vars_y, c("y0", "age")] <- 1

    # 🔧 STRONGER LAG STRUCTURE
    if (k >= 1) {
      pred["y1", c("y0", "age")] <- 1
    }

    if (k > 1) {
      for (j in 2:k) {

        lag1 <- paste0("y", j - 1)
        lag2 <- if (j > 2) paste0("y", j - 2) else NULL

        preds <- c(lag1, "y0", "age")

        if (!is.null(lag2)) {
          preds <- c(preds, lag2)
        }

        pred[paste0("y", j), preds] <- 1
      }
    }

    # rescue as predictor (treatment policy consistency)
    if (length(vars_R) > 0) {
      pred[vars_y, vars_R] <- 1
    }

    pred[, "trt"] <- 0

    ############################################################
    # IMPUTATION BY TREATMENT ARM
    ############################################################
    imp_list <- vector("list", 2)

    for (g in 0:1) {

      dat_g <- dat_hyp |>
        dplyr::filter(trt == g) |>
        dplyr::select(dplyr::all_of(vars_imp))

      imp_list[[g + 1]] <- withr::with_seed(
        seed + g,
        mice::mice(
          dat_g,
          m = m,
          method = meth,
          predictorMatrix = pred,
          maxit = maxit,
          ridge = 1e-4,
          visitSequence = sort(vars_y),

          # 🔧 KEY STABILITY FIX
          donors = 10,

          remove.collinear = FALSE,
          printFlag = FALSE
        )
      )
    }

    ############################################################
    # COMBINE IMPUTATIONS
    ############################################################
    imp_full <- lapply(seq_len(m), function(i) {

      d <- dplyr::bind_rows(
        mice::complete(imp_list[[1]], i),
        mice::complete(imp_list[[2]], i)
      )

      d
    })

    ############################################################
    # ANALYSIS (UNCHANGED)
    ############################################################
    fits <- lapply(imp_full, function(d) {

      d$chg <- d[[paste0("y", k)]] - d$y0

      lm(chg ~ trt + age + y0, data = d)
    })

    imp_obj <- mice::as.mira(fits)
    pooled <- mice::pool(imp_obj)

    sm <- summary(pooled, conf.int = TRUE, conf.level = ci_level)

    trt_row <- sm[sm$term == "trt", ]

    list(
      coef = trt_row$estimate,
      p = trt_row$p.value,
      ci_lower = trt_row$`2.5 %`,
      ci_upper = trt_row$`97.5 %`
    )
  }
}
