#' Analyse Dataset with the Inverse probability weighting
#'
#' @param estimand the estimand targeted by the CI method, default: "hyp" (hypothetical). Allowed also: "tp" (treatment policy).
#'
#' @return an analyse function that can be used in runSimulation that returns a list with the elements
#' * `coef` coefficient for `trt`
#' * `sd` standard deviation for coef
#'
#' @export
#'
#' @importFrom stats as.formula
#' @importFrom estimatr lm_robust
#' @importFrom dplyr filter lag arrange group_by mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ipw ipwtm
#' @importFrom magrittr `%>%`
#'
#' @examples
#' Design <- assumptions_diabetes_rescue() |>
#'   true_summary_statistics_diabetes_rescue()
#'
#' condition <- Design[1, ]
#'
#' dat <- generate_diabetes_rescue(condition)
#'
#' analyse_diabetes_ipw_tp <- analyse_diabetes_ipw(estimand = "tp")
#' analyse_diabetes_ipw_tp(condition, dat)
#'
#' analyse_diabetes_ipw_hyp <- analyse_diabetes_ipw(estimand = "hyp")
#' analyse_diabetes_ipw_hyp(condition, dat)
analyse_diabetes_ipw <- function(estimand = "hyp") {
  # What still remains to be done is
  # documentation

  function(condition, dat, fixed_objects = NULL) {
    stopifnot(estimand %in% c("hyp", "tp"))

    k <- condition$k # number of last visit

    # reformate dat to long format with one Hba1c column, one outcome column 'y' (change of Hba1c) and a time variable 'visit'
    dat_long <- tidyr::pivot_longer(dat,
                                    cols = matches("^[yR]\\d+$"),
                                    names_to = c(".value", "visit"),
                                    names_pattern = "([yR])(\\d+)"
    ) %>%
      mutate(
        hba1c = y,
        visit = as.numeric(sub("y", "", visit)),
        R = ifelse(visit == 0, 0, ifelse(visit == condition$k, dplyr::lag(R), R)), # set rescue at baseline to 0
        R_lag = dplyr::lag(R, default = NA)
      ) %>%
      arrange(id, visit) %>%
      group_by(id) %>% # make sure table is grouped by id and ordered by visit
      mutate(
        hba1c_0 = hba1c[visit == 0], # HbA1 at baseline
        hba1c_lag = dplyr::lag(hba1c, default = NA), # HbA1 at visit j-1
        y = hba1c - hba1c_0
      ) %>% # HbA1c change
      filter(visit != 0) # do not include baseline visits

    if (estimand == "tp") {
      dat_long$exposure <- ifelse(is.na(dat_long$y), 1L, 0L) # indicator for missing outcomes (somehow this only works if variable named exactly "exposure")

      temp <- ipw::ipwtm(
        exposure = exposure, # indicator for missing data at visit j
        family = "binomial",
        link = "logit",
        denominator = ~ trt + age + hba1c_lag + R_lag,
        id = id,
        timevar = visit,
        type = "first",
        data = dat_long
      )
      model <- lm_robust( # OLS with HC2 variance estimator
        as.formula(paste0("y ~ trt + hba1c_0 + age")),
        weights = temp$ipw.weights[dat_long$visit == k & dat_long$exposure == 0],
        data = dat_long[dat_long$visit == k & dat_long$exposure == 0,]
      )

      list(
        coef = summary(model)$coefficients["trt", "Estimate"],
        sd = summary(model)$coefficients["trt", "Std. Error"],
        p = summary(model)$coefficients["trt", "Pr(>|t|)"],
        ci_lower = summary(model)$coefficients["trt", "CI Lower"],
        ci_upper = summary(model)$coefficients["trt", "CI Upper"]
      )
    } else if (estimand == "hyp") {
      dat_long$exposure <- ifelse(is.na(dat_long$y) | dat_long$R == 1, 1L, 0L) # indicator for missing outcomes and rescue medication

      temp <- ipw::ipwtm(
        exposure = exposure, # indicator for missingness or rescue
        family = "binomial",
        link = "logit",
        denominator = ~ trt + age + hba1c_lag,
        id = id,
        timevar = visit,
        type = "first", #  models are fitted only on observations up to and including the first time point where y is missing, afterwards weights will be constant
        data = dat_long
      )
      model <- estimatr::lm_robust( # OLS with HC2 variance estimator
        as.formula(paste0("y ~ trt + hba1c_0 + age")),
        weights = temp$ipw.weights[dat_long$visit == k & dat_long$exposure == 0],
        data = dat_long[dat_long$visit == k & dat_long$exposure == 0,]
      )

      list(
        coef = summary(model)$coefficients["trt", "Estimate"],
        sd = summary(model)$coefficients["trt", "Std. Error"],
        p = summary(model)$coefficients["trt", "Pr(>|t|)"],
        ci_lower = summary(model)$coefficients["trt", "CI Lower"],
        ci_upper = summary(model)$coefficients["trt", "CI Upper"]
      )
    }
  }
}
