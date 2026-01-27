#' Analyse Dataset with Longitudinal Targeted Minimum Loss/(Maximum Likelihood) Estimation
#'
#' @param level confidence level for CI computation
#' @param alternative alternative hypothesis for the tests "two.sided" or "one.sided"
#'
#' @return an analyse function that returns a list with the elements
#' * `p` p value of the score test (two.sided) or the Wald test (one.sided)
#' * `alternative` the alternative used
#' * `coef` coefficient for `trt`
#' * `hr` hazard ratio for `trt`
#' * `hr_lower` lower 95% confidence intervall boundary for the hazard ratio for `trt`
#' * `hr_upper`lower 95% confidence intervall boundary for the hazard ratio for `trt`
#' * `CI_level` the CI level used
#' * `N_pat` number of patients
#' * `N_evt` number of events
#'
#' @export
#'
#' @details
#'
#' `alternative` can be "two.sided" for a two sided test of equality of the
#' summary statistic or "one.sided" for a one sided test testing H0: treatment
#' has equal or shorter survival than control vs. H1 treatment has longer
#' survival than control.
#'
#' @examples
#' condition <- merge(
#'   assumptions_delayed_effect(),
#'   design_fixed_followup(),
#'   by = NULL
#' ) |>
#'   head(1)
#' dat <- generate_delayed_effect(condition)
#' analyse_ltmle()(condition, dat)
#' analyse_ltmle()(assumptions_diabetes_rescue, generate_diabetes_rescue(assumptions_diabetes_rescue()))
analyse_ltmle <- function(level = 0.95, alternative = "two.sided") {
  stopifnot(alternative %in% c("two.sided", "one.sided"))

  # What still remains to be done is to
    # documentation
    # determine if the "alternative" statement can be incorporated

  function(condition, dat, fixed_objects = NULL) {
    dat <- dat |> dplyr::mutate(rescue_start = ifelse(is.na(rescue_start), condition$k[1]+2, rescue_start))
    pred <- mice::make.predictorMatrix(dat)
    pred[upper.tri(pred)] <- 0
    dats <- mice::mice(dat,
                 m = 3, print = FALSE, seed = 2026,
                 formulas = make.formulas(dat, blocks = mice::make.blocks(dat), predictorMatrix = pred)
    )
    t <- lapply(1:dats$m, function(i) {
      dat <- mice::complete(dats, i)
      dat_comp <- dat |> mutate(
        across(
          matches("^y[0-9]+$") & !y0,
          ~ .x - y0,
          .names = "yc{gsub('y', '', .col)}"
        ))
      for (j in 1:condition$k[1]) {
        dat_comp[[paste0("R", j)]] <- ifelse(dat_comp$rescue_start <= j, 1L, 0L)
      }

      model <- ltmle::ltmle(
        data = dat_comp |>
          dplyr::select(
            age,
            y0,
            trt,
            unlist(
              lapply(
                seq_len(condition$k[1]),
                \(i) c(paste0("R", i), paste0("yc", i))
              )
            )),
        Anodes <- c(
          "trt",
          paste0("R", seq_len(condition$k[1]))
        ),
        Ynodes = c(paste0("yc", seq_len(condition$k[1]))),
        SL.library = "glm",
        abar = list(
          c(1, rep(0, condition$k[1])),
          rep(0, condition$k[1] + 1)
        )
      )
      summary(model)$effect.measures$ATE$estimate
      out <- data.frame(
        ATE = summary(model)$effect.measures$ATE$estimate,
        CI_lower = summary(model)$effect.measures$ATE$CI[1],
        CI_upper = summary(model)$effect.measures$ATE$CI[2],
        se = summary(model)$effect.measures$ATE$std.dev,
        p_value = summary(model)$effect.measures$ATE$pvalue,
        iteration = i
      )
      out
    })
    t
  }
}

<<<<<<< HEAD
=======
# analyse_ltmle()(assumptions_diabetes_rescue(print = FALSE),
#                 generate_diabetes_rescue(assumptions_diabetes_rescue(print = FALSE)))

>>>>>>> tmp_fix_TF
