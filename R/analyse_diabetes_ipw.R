#' Analyse data set with the Inverse probability weighting
#'
#' @param strategy the estimand targeted by the CI method, default: "hypothetical". Allowed also: "treatment_policy".
#'
#' @return A function that, when called with `condition` and `dat`, returns a list with:
#' * `coef` coefficient for `trt`
#' * `se` standard error for coef
#' * `p` p-value for coef
#' * `ci_lower` lower bound of 95% confidence interval for coef
#' * `ci_upper` upper bound of 95% confidence interval for coef
#'
#' @details
#' This is a function to implement the inverse probability weighting for estimation of the treatment effect in a diabetes trial.
#' The function first reshapes the data to a long format, where each row corresponds to a visit for a patient.
#' Then, depending on the chosen strategy, it creates an indicator variable for missing outcomes due to treatment discontinuation (treatment policy) or due to both treatment discontinuation and rescue medication use (hypothetical).
#' It fits a logistic regression model to estimate the probability of missingness and computes inverse probability weights.
#' Finally, it fits a weighted linear regression model to estimate the treatment effect on the change in HbA1c at the final visit, using robust standard errors to account for the weighting.
#'
#'
#'
#'
#' @export
#'
#' @importFrom stats as.formula lm confint
#' @importFrom lmtest coeftest
#' @importFrom sandwich vcovHC
#' @importFrom dplyr filter lag arrange group_by mutate
#' @importFrom tidyr pivot_longer
#' @importFrom ipw ipwtm
#' @importFrom magrittr `%>%`
#'
#' @examples
#' Design <- diabetes_scenario()[1, ] |>
#'   diabetes_scenario_set_truevalues()
#'
#' dat <- generate_diabetes(Design)
#'
#' analyse_diabetes_ipw(strategy = "treatment_policy")(Design, dat)
#' analyse_diabetes_ipw(strategy = "hypothetical")(Design, dat)
#'
analyse_diabetes_ipw <- function(strategy = "hypothetical") {
  function(condition, dat, fixed_objects = NULL) {
    stopifnot(strategy %in% c("hypothetical", "treatment_policy"))

    k <- condition$k # number of last visit

    # reformat dat to long format with one HbA1c column, one outcome column 'y' (change of HbA1c) and a time variable 'visit'
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
        hba1c_0 = hba1c[visit == 0], # HbA1c at baseline
        hba1c_lag = dplyr::lag(hba1c, default = NA), # HbA1c at visit j-1
        y = hba1c - hba1c_0
      ) %>% # HbA1c change
      filter(visit != 0) # do not include baseline visits

    if (strategy == "treatment_policy") {
      dat_long$exposure <- ifelse(is.na(dat_long$y), 1L, 0L) # indicator for missing outcomes (somehow this only works if variable named exactly "exposure")
    } else if (strategy == "hypothetical") {
      dat_long$exposure <- ifelse(is.na(dat_long$y) | dat_long$R == 1, 1L, 0L) # indicator for missing outcomes and rescue medication
    }

    if (nrow(dat_long[dat_long$visit == k & dat_long$exposure == 1, ]) > 1) { # there need to be more than one missing value due to discontinuation (both estimands) or rescue (in case of hypothetical estimand only)
    } else if (strategy == "hypothetical") {
      dat_long$exposure <- ifelse(is.na(dat_long$y) | dat_long$R == 1, 1L, 0L) # indicator for missing outcomes and rescue medication
    }

    if (nrow(dat_long[dat_long$visit == k & dat_long$exposure == 1, ]) > 1) { # there need to be more than one missing value due to discontinuation (both estimands) or rescue (in case of hypothetical estimand only)

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
      fit <- lm( # OLS with HC2 variance estimator
        as.formula(paste0("y ~ trt + hba1c_0 + age")),
        weights = temp$ipw.weights[dat_long$visit == k & dat_long$exposure == 0],
        data = dat_long[dat_long$visit == k & dat_long$exposure == 0, ]
      )
    } else {
      fit <- lm( # OLS with HC2 variance estimator
        as.formula(paste0("y ~ trt + hba1c_0 + age")),
        data = dat_long[dat_long$visit == k & dat_long$exposure == 0, ]
      )
    }

    model <- lmtest::coeftest(fit, vcov = sandwich::vcovHC(fit, type = "HC2"))
    ci <- stats::confint(model)

    list(
      coef = model["trt", "Estimate"],
      se = model["trt", "Std. Error"],
      p = model["trt", "Pr(>|t|)"],
      ci_lower = ci[2, 1],
      ci_upper = ci[2, 2]
    )
  }
}
