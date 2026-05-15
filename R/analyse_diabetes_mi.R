#' Analyse data set with Multiple Imputation
#'
#' Creates an analysis function for the diabetes rescue simulation.
#' Depending on the chosen estimand, either a hypothetical strategy
#' with rescue-based censoring or a treatment policy strategy using
#' observed rescue information is applied.
#'
#' ## Multiple Imputation
#'
#' Missing HbA1c outcomes are imputed using the \pkg{mice} package.
#'
#' Imputation is performed separately within each treatment group,
#' consistent with the trial protocol and avoiding cross-treatment
#' borrowing of information.
#'
#' ## Estimand-specific handling
#'
#' ### Treatment policy strategy
#'
#' Under the treatment policy estimand, observed post-rescue outcomes
#' are retained.
#'
#' Rescue indicators may be included as predictors in the imputation
#' model to preserve associations between rescue use and future HbA1c
#' trajectories.
#'
#' ### Hypothetical strategy
#'
#' Under the hypothetical estimand, the analysis targets the scenario
#' in which rescue medication had not been available.
#'
#' For subjects initiating rescue medication at visit \eqn{s},
#' all post-rescue outcomes are set to missing:
#'
#' \deqn{
#' y_{i,j} = \mathrm{NA}
#' \quad \text{for all } j > s
#' }
#'
#' The outcome at the rescue visit itself (\eqn{j=s}) is retained.
#'
#' Rescue indicators are excluded from the HbA1c imputation models
#' under the hypothetical strategy to avoid conditioning on post-rescue
#' information.
#'
#' ## Imputation model
#'
#' HbA1c values are imputed using predictive mean matching (`"pmm"`).
#'
#' The imputation model includes:
#'
#' * baseline HbA1c (\eqn{y_0}),
#' * age at baseline,
#' * adjacent HbA1c visit history,
#' * optionally rescue indicators under the treatment policy strategy.
#'
#' Rescue indicators are imputed using logistic regression.
#'
#' A total of `m` imputations are generated.
#'
#' ## Analysis model
#'
#' For each imputed dataset, an ANCOVA model is fitted:
#'
#' \deqn{
#' \Delta y_i =
#' \beta_0 +
#' \beta_1 \cdot \mathrm{trt}_i +
#' \beta_2 \cdot \mathrm{age}_i +
#' \beta_3 \cdot y_{0,i} +
#' \varepsilon_i
#' }
#'
#' where:
#'
#' * \eqn{\Delta y_i} is the change from baseline in HbA1c,
#' * \eqn{y_{0,i}} is baseline HbA1c.
#'
#' ## Hypothesis testing
#'
#' One-sided p-values are reported for the treatment effect:
#'
#' \deqn{
#' H_0: \beta_{trt} \ge 0
#' }
#'
#' versus
#'
#' \deqn{
#' H_A: \beta_{trt} < 0
#' }
#'
#' since negative changes in HbA1c are considered favorable.
#'
#' Two-sided confidence intervals are retained.
#'
#' ## Pooling
#'
#' Estimates from the imputed datasets are pooled using Rubin's rules
#' via \code{mice::pool()}.
#'
#' ## Reproducibility
#'
#' Randomness from the multiple imputation procedure is controlled using
#' \code{withr::with_seed()}, ensuring reproducibility without affecting
#' the global random number generator state.
#'
#' @param strategy Either `"hypothetical"` or `"treatment_policy"`.
#' @param m Number of imputations.
#' @param maxit Maximum number of MICE iterations.
#' @param ci_level Confidence level for the two-sided confidence interval.
#' @param seed Random seed for imputation (used locally via withr).
#'
#' @return A function that, when called with `condition` and `dat`,
#' returns a list with:
#'
#' * `coef` estimated treatment effect for `trt`
#' * `p` one-sided p-value for the hypothesis \eqn{\beta_{trt} < 0}
#' * `ci_lower` lower bound of the two-sided confidence interval
#' * `ci_upper` upper bound of the two-sided confidence interval
#'
#' @importFrom mice mice make.method make.predictorMatrix
#'   complete as.mira pool
#' @importFrom dplyr filter select bind_rows all_of
#' @importFrom stats lm pt
#' @importFrom withr with_seed
#'
#' @export
#'
#' @examples
#' \donttest{
#'
#' Design <- diabetes_scenario()[1, ] |>
#'   diabetes_scenario_set_truevalues()
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
#' ## Compare estimated effects
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

    vars_y <- paste0("y", 0:k)
    vars_R <- if (k > 1) paste0("R", 1:(k - 1)) else character(0)

    vars_imp <- c(vars_y, vars_R, "age", "trt")

    dat_hyp <- dat

    ############################################################
    # Hypothetical censoring (rescue-based)
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
    # METHOD SPECIFICATION
    ############################################################
    meth <- mice::make.method(dat_hyp[vars_imp])
    meth[vars_y] <- "pmm"
    meth[c("age", "trt")] <- ""

    if (length(vars_R) > 0) {
      meth[vars_R] <- "logreg"
    }

    ############################################################
    # PREDICTOR MATRIX
    ############################################################
    pred <- mice::make.predictorMatrix(dat_hyp[vars_imp])
    pred[,] <- 0

    pred[vars_y, c("y0", "age")] <- 1

    if (strategy == "treatment_policy" && length(vars_R) > 0) {
      pred[vars_y, vars_R] <- 1
    }

    if (strategy == "hypothetical" && length(vars_R) > 0) {
      pred[vars_y, vars_R] <- 0
    }

    if (k > 1) {
      for (j in 2:k) {
        pred[paste0("y", j), paste0("y", j - 1)] <- 1
      }
    }

    pred[, "trt"] <- 0

    ############################################################
    # IMPUTATION PER TREATMENT ARM
    ############################################################
    imp_list <- vector("list", 2)

    for (g in 0:1) {

      dat_g <- dat_hyp |>
        dplyr::filter(trt == g) |>
        dplyr::select(dplyr::all_of(vars_imp))

      ##########################################################
      # ✅ Minimal robust fix: keep ALL R as factors
      ##########################################################
      if (length(vars_R) > 0) {
        dat_g[vars_R] <- lapply(dat_g[vars_R], function(x)
          factor(as.integer(x), levels = c(0, 1))
        )
      }

      ##########################################################
      # ✅ Disable imputation for degenerate R variables
      ##########################################################
      meth_g <- meth

      if (length(vars_R) > 0) {
        for (r in vars_R) {
          vals <- unique(stats::na.omit(dat_g[[r]]))
          if (length(vals) < 2) {
            meth_g[r] <- ""   # no variation → skip imputation
          }
        }
      }

      ##########################################################
      # MICE call
      ##########################################################
      imp_list[[g + 1]] <- withr::with_seed(
        seed + g,
        mice::mice(
          dat_g,
          m = m,
          method = meth_g,
          predictorMatrix = pred,
          maxit = maxit,
          ridge = 5e-5,
          donors = 5,
          pmm.k = 5,
          visitSequence = sort(vars_y),
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

      if (length(vars_R) > 0) {
        d[vars_R] <- lapply(d[vars_R], function(x)
          factor(as.integer(x), levels = c(0, 1))
        )
      }

      d
    })

    ############################################################
    # ANALYSIS
    ############################################################
    fits <- lapply(imp_full, function(d) {

      d$chg <- d[[paste0("y", k)]] - d$y0

      lm(chg ~ trt + age + y0, data = d)
    })

    imp_obj <- mice::as.mira(fits)
    pooled <- mice::pool(imp_obj)

    sm <- summary(pooled, conf.int = TRUE, conf.level = ci_level)

    trt_row <- sm[sm$term == "trt", ]

    if (nrow(trt_row) == 0) {
      return(list(
        coef = NA_real_,
        p = NA_real_,
        ci_lower = NA_real_,
        ci_upper = NA_real_
      ))
    }

    ############################################################
    # ONE-SIDED P-VALUE (H1: trt < 0)
    ############################################################
    t_stat <- trt_row$statistic
    df <- trt_row$df

    p_one_sided <- stats::pt(t_stat, df = df)

    list(
      coef = trt_row$estimate,
      p = p_one_sided,
      ci_lower = trt_row[["2.5 %"]],
      ci_upper = trt_row[["97.5 %"]]
    )
  }
}
