#' Analyse Dataset with the Inverse probability weighting
#'
#' @param estimand the estimand targeted by the CI method "tp" for treatment policy, otherwise hypothetical
#' @param level confidence level for CI computation
#' @param alternative alternative hypothesis for the tests "two.sided" or "one.sided"
#'
#' @return an analyse function that returns a list with the elements
#' * `p` p value of the score test (two.sided) or the Wald test (one.sided)
#' * `alternative` the alternative used
#' * `coef` coefficient for `trt`
#' * `hr` hazard ratio for `trt`
#' * `hr_lower` lower 95% confidence interval boundary for the hazard ratio for `trt`
#' * `hr_upper`lower 95% confidence interval boundary for the hazard ratio for `trt`
#' * `CI_level` the CI level used
#' * `N_pat` number of patients
#' * `N_evt` number of events
#'
#' @export
#'
#' @importFrom stats as.formula
#' @importFrom estimatr lm_robust
#' @importFrom dplyr filter lag
#' @importFrom tidyr pivot_longer
#' @importFrom ipw ipwtm
#' @details
#' `alternative` can be "two.sided" for a two sided test of equality of the
#' summary statistic or "one.sided" for a one sided test testing H0: treatment
#' has equal or shorter survival than control vs. H1 treatment has longer
#' survival than control.
#'
#' @examples
#' Design <- assumptions_diabetes_rescue() #|>
#'   #true_summary_statistics_diabetes_rescue()
#'
#' condition <- Design[1, ]
#'
#' dat <- generate_diabetes_rescue(condition)
#'
#' analyse_ipw_tp  <- analyse_ipw(estimand="tp")
#' analyse_ipw_tp (condition, dat)
#'
#' analyse_ipw_hyp  <- analyse_ipw(estimand="hyp")
#' analyse_ipw_hyp (condition, dat)

analyse_ipw <- function(estimand = "tp", level = 0.95, alternative = "two.sided") {
  stopifnot(alternative %in% c("two.sided", "one.sided"))

  # documentation

  function(condition, dat, fixed_objects = NULL) {
    dat_long <- tidyr::pivot_longer(dat, y0:y12, names_to = "visit", values_to = "y")
    dat_long<- dat_long %>%
      mutate(visit = as.numeric(sub('y', '', visit))) %>%
      mutate(r = ifelse(!is.na(rescue_start)&rescue_start<=visit, 1, 0)) %>% #rescue
      mutate(r_lag = dplyr::lag(r, default = NA)) %>% #rescue j-1
      arrange(id, visit) %>% group_by(id) %>%
      mutate(y_lag = dplyr::lag(y, default = NA)) %>% #HbA1 j-1
      mutate(y0 = y[visit == 0]) %>% # y0
      filter(visit!=0) #do not include baseline visits

    if (estimand == "tp") {
      dat_long$exposure<- ifelse(is.na(dat_long$y), 1, 0) #mark missing outcomes (somehow this only works if variable named exactly "exposure")

      temp <- ipw::ipwtm(
        exposure = exposure, # indicator for missing data at visit j
        family = "binomial",
        link = "logit",
        denominator = ~ trt + age + y_lag + r_lag,
        id = id,
        timevar = visit,
        type = "first",
        data = dat_long
      )
      model <- estimatr::lm_robust( # OLS with HC2 variance estimator
        as.formula(paste0("y ~ trt + y0 + age")),
        weights = temp$ipw.weights[dat_long$visit==12 & dat_long$exposure==0], # maybe use condition$k[1] instead of "12" if we want to vary the follow-up length
        data = dat_long[dat_long$visit==12 & dat_long$exposure==0,]
        )
      CI <- confint(model, level=level)

      list(
        p = summary(model)$coefficients["trt", "Pr(>|t|)"],
        coef = summary(model)$coefficients["trt", "Estimate"],
        ci_lower = CI["trt", 1],
        ci_upper = CI["trt", 2]
      )
    } else {
      dat_long$exposure <- ifelse(is.na(dat_long$y) | dat_long$r == 1, 1, 0) #indicator for missing outcomes and rescue medication

      temp <- ipw::ipwtm(
        exposure = exposure, #indicator for missingness or rescue
        family = "binomial",
        link = "logit",
        denominator = ~ trt + age + y_lag,
        id = id,
        timevar = visit,
        type = "first", #  models are fitted only on observations up to and including the first time point where y is missing, afterwards weights will be constant
        data = dat_long
      )
      model <- estimatr::lm_robust( # OLS with HC2 variance estimator
        as.formula(paste0("y ~ trt + y0 + age")),
        weights = temp$ipw.weights[dat_long$visit==12 & dat_long$exposure==0],
        data = dat_long[dat_long$visit==12 & dat_long$exposure==0,]
      )
      CI <- confint(model, level=level)

      list(
        p = summary(model)$coefficients["trt", "Pr(>|t|)"],
        coef = summary(model)$coefficients["trt", "Estimate"],
        ci_lower = CI["trt", 1],
        ci_upper = CI["trt", 2]
      )
    }
  }
}

# Testing
#analyse_ipw(estimand = "tp")(assumptions_diabetes_rescue(print = FALSE),
#   generate_diabetes_rescue(assumptions_diabetes_rescue(print = FALSE)))
#analyse_ipw(estimand = "hyp")(assumptions_diabetes_rescue(print = FALSE),
#   generate_diabetes_rescue(assumptions_diabetes_rescue(print = FALSE)))
