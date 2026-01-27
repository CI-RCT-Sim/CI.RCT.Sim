#' Analyse Dataset with the Inverse probability weighting
#'
#' @param level confidence level for CI computation
#' @param alternative alternative hypothesis for the tests "two.sided" or "one.sieded"
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
#' analyse_ipw()(condition, dat)
analyse_ipw <- function(estimand = "tp", level = 0.95, alternative = "two.sided") {
  stopifnot(alternative %in% c("two.sided", "one.sided"))

  # What still remains to be done is to
  # documentation
  # Statements to be able to switch between tp and hyp

  function(condition, dat, fixed_objects = NULL) {
    if (estimand == "tp") {
      temp <- ipw::ipwpoint(
        exposure = trt,
        family = "binomial",
        link = "logit",
        numerator = ~1,
        denominator = ~age,
        data = dat
      )
      # browser()
      model1 <- stats::lm(as.formula(paste0("y", condition$k[1], " ~ trt + y0")), data = dat)
      model <- stats::lm(as.formula(paste0("y", condition$k[1], " ~ trt + y0")),
        weights = temp$ipw.weights, data = dat
      )
      rbind(
        summary(model1)$coefficients["trt", ],
        summary(model)$coefficients["trt", ]
      )
    } else {
      rweight <- ipw::ipwtm(
        exposure = R,
        family = "binomial",
        link = "logit", # other options: "logit", "probit", "cauchit", "log" and "cloglog"
        denominator = ~ A + Y + X2 + X3,
        id = id,
        timevar = Visit,
        type = "first",
        data = dat
      )
      ipw.mod <- lm(Y ~ A,
        data = d %>% filter(Visit == 3, R == 0),
        weights = rweight$ipw.weights[d$Visit == 3 & d$R == 0]
      )

      estimate <- ipw.mod$coefficients[["A"]]
      estimate
    }
  }
}

analyse_ipw(estimand = "tp")(assumptions_diabetes_rescue(print = FALSE),
  generate_diabetes_rescue(assumptions_diabetes_rescue(print = FALSE)))
