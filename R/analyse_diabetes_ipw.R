#' Analyse Dataset with the Inverse probability weighting
#'
#' @param strategy the estimand targeted by the CI method, default: "hypothetical". Allowed also: "treatment_policy".
#'
#' @return an analyse function that can be used in runSimulation that returns a list with the elements
#' * `coef` coefficient for `trt`
#' * `sd` standard deviation for coef
#' * `p` p-value for coef
#' * `ci_lower` lower bound of confidence interval for coef
#' * `ci_upper` upper bound of confidence interval for coef
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
#' Design <- diabetes_scenario() |>
#'   diabetes_scenario_set_truevalues()
#'
#' condition <- Design[1, ]
#'
#' dat <- generate_diabetes(condition)
#'
#' analyse_diabetes_ipw(strategy = "treatment_policy")(condition, dat)
#' analyse_diabetes_ipw(estimand = "hypothetical")(condition, dat)
#'
analyse_diabetes_ipw <- function(strategy = "hypothetical") {
  function(condition, dat, fixed_objects = NULL) {
    stopifnot(strategy %in% c("hypothetical", "treatment_policy"))

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

    if (strategy == "treatment_policy") {
      dat_long$exposure <- ifelse(is.na(dat_long$y), 1L, 0L) # indicator for missing outcomes (somehow this only works if variable named exactly "exposure")
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
      sd = model["trt", "Std. Error"],
      p = model["trt", "Pr(>|t|)"],
      ci_lower = ci[2, 1],
      ci_upper = ci[2, 2]
    )
  }
}
