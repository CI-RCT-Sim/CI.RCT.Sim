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
#' @importFrom dplyr filter lag arrange group_by
#' @importFrom tidyr pivot_longer
#' @importFrom ipw ipwtm
#' @importFrom magrittr `%>%`
#' @details
#' `alternative` can be "two.sided" for a two sided test of equality of the
#' summary statistic or "one.sided" for a one sided test testing H0: treatment
#' has equal or shorter survival than control vs. H1 treatment has longer
#' survival than control.
#'
#' @examples
#' Design <- assumptions_diabetes_rescue() |>
#'   true_summary_statistics_diabetes_rescue()
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

  # What still remains to be done is
  # documentation
  # determine if the "alternative" statement can be incorporated

  function(condition, dat, fixed_objects = NULL) {
    k<-condition$k[1] #number of last visit

    #reformate dat to long format with one outcome column 'y' and a time variable 'visit'
    dat_long <- pivot_longer(dat, y0:paste0("y",k), names_to = "visit", values_to = "y")

    dat_long<- dat_long %>%
      mutate(visit = as.numeric(sub('y', '', visit))) %>%
      mutate(r = ifelse(!is.na(rescue_start)&rescue_start<=visit, 1, 0)) %>% #new variable for rescue at visit j
      mutate(r_lag = lag(r, default = NA)) %>% #rescue at visit j-1
      arrange(id, visit) %>% group_by(id) %>%
      mutate(y_lag = lag(y, default = NA)) %>% #HbA1 at visit j-1
      mutate(y0 = y[visit == 0]) %>% # HbA1 at baseline
      filter(visit!=0) #do not include baseline visits

    if (estimand == "tp") {
      dat_long$exposure<- ifelse(is.na(dat_long$y), 1L, 0L) #indicator for missing outcomes (somehow this only works if variable named exactly "exposure")

      temp <- ipwtm(
        exposure = exposure, # indicator for missing data at visit j
        family = "binomial",
        link = "logit",
        denominator = ~ trt + age + y_lag + r_lag,
        id = id,
        timevar = visit,
        type = "first",
        data = dat_long
      )
      model <- lm_robust( # OLS with HC2 variance estimator
        as.formula(paste0("y ~ trt + y0 + age")),
        weights = temp$ipw.weights[dat_long$visit==k & dat_long$exposure==0],
        data = dat_long[dat_long$visit==k & dat_long$exposure==0,]
        )
      CI <- confint(model, level=level)

      list(
        p = summary(model)$coefficients["trt", "Pr(>|t|)"],
        coef = summary(model)$coefficients["trt", "Estimate"],
        ci_lower = CI["trt", 1],
        ci_upper = CI["trt", 2]
      )
    } else {
      dat_long$exposure <- ifelse(is.na(dat_long$y) | dat_long$r == 1, 1L, 0L) #indicator for missing outcomes and rescue medication

      temp <- ipwtm(
        exposure = exposure, #indicator for missingness or rescue
        family = "binomial",
        link = "logit",
        denominator = ~ trt + age + y_lag,
        id = id,
        timevar = visit,
        type = "first", #  models are fitted only on observations up to and including the first time point where y is missing, afterwards weights will be constant
        data = dat_long
      )
      model <- lm_robust( # OLS with HC2 variance estimator
        as.formula(paste0("y ~ trt + y0 + age")),
        weights = temp$ipw.weights[dat_long$visit==k & dat_long$exposure==0],
        data = dat_long[dat_long$visit==k & dat_long$exposure==0,]
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


