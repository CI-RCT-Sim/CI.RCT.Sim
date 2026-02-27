#' Analysis for hypothetical and treatment policy estimands
#'
#' Creates an analysis function for the diabetes rescue simulation.
#' Depending on the chosen estimand, either a hypothetical strategy
#' with multiple imputation after rescue initiation is applied,
#' or a treatment policy strategy using observed data only.
#'
#' @param estimand Either `"hypothetical"` or `"treatment_policy"`.
#' @param m Number of imputations (used for hypothetical only).
#' @param maxit Maximum number of MICE iterations.
#' @param ci_level Confidence level.
#' @param seed Random seed for imputation.
#'
#' @return A function with arguments `(condition, dat, fixed_objects = NULL)`
#'   returning a list with elements `coef`, `p`, `ci_lower`, `ci_upper`.
#'
#' @importFrom mice mice make.method make.predictorMatrix complete as.mira pool
#' @importFrom dplyr filter select bind_rows all_of
#' @importFrom stats lm
#' @export
analyse_diabetes_rescue <- function(
    estimand = c("hypothetical", "treatment_policy"),
    m = 10,
    maxit = 20,
    ci_level = 0.95,
    seed = 123) {

  estimand <- match.arg(estimand)

  function(condition, dat, fixed_objects = NULL) {

    k <- condition$k

    ############################################################
    # TREATMENT POLICY ESTIMAND
    ############################################################
    if (estimand == "treatment_policy") {

      dat$chg <- dat[[paste0("y", k)]] - dat$y0

      fit <- stats::lm(chg ~ trt + age + y0, data = dat)

      sum_fit <- summary(fit)
      ci <- confint(fit, level = ci_level)["trt", ]

      return(list(
        coef     = coef(fit)["trt"],
        p        = sum_fit$coefficients["trt", "Pr(>|t|)"],
        ci_lower = ci[1],
        ci_upper = ci[2]
      ))
    }

    ############################################################
    # HYPOTHETICAL ESTIMAND
    # Set post-rescue data to missing, then MI
    ############################################################

    dat_hyp <- dat

    # Set all post-rescue Y values to NA
    for (i in seq_len(nrow(dat_hyp))) {

      rs <- dat_hyp$rescue_start[i]

      if (!is.na(rs) && rs < k) {

        miss_visits <- (rs + 1):k

        for (v in miss_visits) {
          dat_hyp[i, paste0("y", v)] <- NA
        }
      }
    }

    vars_y <- paste0("y", 0:k)
    vars_R <- if (k > 1) paste0("R", 1:(k - 1)) else character(0)
    vars_imp <- c(vars_y, vars_R, "age", "trt")

    meth <- mice::make.method(dat_hyp[vars_imp])
    meth[vars_y] <- "pmm"

    if (k >= 2) {
      meth["R1"] <- "logreg"
    }

    if (k > 2) {
      for (j in 2:(k - 1)) {
        meth[paste0("R", j)] <-
          paste0("~ I(pmax(R", j - 1, ", R", j, "))")
      }
    }

    meth[c("age", "trt")] <- ""

    pred <- mice::make.predictorMatrix(dat_hyp[vars_imp])
    diag(pred) <- 0
    pred[, "trt"] <- 0
    pred[, "y0"] <- 1

    imp_list <- vector("list", 2)

    for (g in 0:1) {

      dat_g <- dat_hyp |>
        dplyr::filter(trt == g) |>
        dplyr::select(dplyr::all_of(vars_imp))

      if (length(vars_R) > 0) {
        dat_g[vars_R] <- lapply(
          dat_g[vars_R],
          function(x) factor(x, levels = c(0, 1))
        )
      }

      imp_list[[g + 1]] <- mice::mice(
        dat_g,
        m = m,
        method = meth,
        predictorMatrix = pred,
        maxit = maxit,
        seed = seed + g,
        printFlag = FALSE
      )
    }

    imp_full <- lapply(seq_len(m), function(i) {
      dplyr::bind_rows(
        mice::complete(imp_list[[1]], i),
        mice::complete(imp_list[[2]], i)
      )
    })

    fit_models <- lapply(imp_full, function(d) {
      d$chg <- d[[paste0("y", k)]] - d$y0
      stats::lm(chg ~ trt + age + y0, data = d)
    })

    imp_obj <- mice::as.mira(fit_models)
    pooled <- mice::pool(imp_obj)

    sum_pooled <- summary(
      pooled,
      conf.int = TRUE,
      conf.level = ci_level
    )

    trt_row <- sum_pooled[sum_pooled$term == "trt", ]

    list(
      coef     = trt_row$estimate,
      p        = trt_row$p.value,
      ci_lower = trt_row[[grep("%", names(trt_row))[1]]],
      ci_upper = trt_row[[grep("%", names(trt_row))[2]]]
    )
  }
}

# Design <- assumptions_diabetes_rescue()
# dat <- generate_diabetes_rescue(Design[1, ])
#
# analyse_mi <- analyse_diabetes_rescue(estimand = "treatment_policy",m = 5)
# analyse_mi(Design[1, ], dat)
