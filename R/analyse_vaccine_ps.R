#' Analyse dataset from vaccine scenario using principal score weighting
#'
#' @param ci_level the confidence level for the CIs (defaults to 0.95)
#' @param VE_margin vaccine efficacy margin for the super-superiority test
#'
#' @returns an analyse function that returns a list with the elements
#'  * `p` the p-value of the super-superiority test
#'  * `VE` the point estimate for the vaccine efficacy
#'  * `VE_lower` the lower CI limit for vaccine efficacy
#'  * `VE_upper` the upper CI limit for vaccine efficacy
#'  * `OR` the point estimate for the odds-ratio for infection
#'  * `OR_lower` the lower CI limit for the odds-ratio for infection
#'  * `OR_upper` the upper CI limit for the odds-ratio for infection
#' @export
#'
#' @importFrom dplyr case_when mutate
#' @importFrom emmeans emmeans regrid contrast
#' @importFrom graphics pairs
#'
#' @examples
#' Design <- vaccine_scenario() |>
#'   vaccine_scenario_set_gamma_0() |>
#'   vaccine_scenario_set_true_eff() |>
#'   vaccine_scenario_set_samplesize()
#'
#' dat <- generate_vaccine(Design[3,])
#' my_analyse <- analyse_vaccine_ps(ci_level=0.95)
#' my_analyse(Design[3, ], dat)
analyse_vaccine_ps <- function(ci_level=0.95, VE_margin=0.3){
  function(condition, dat, fixed_objects = NULL){

    dat1 <- dat |>
      within({
        C <- factor(C, levels=c("0", "1"))
        trt <- factor(trt, levels=c("1", "0"))
      })
    mod_ps <- glm(C~V+W, subset = (trt==1), data=dat1, family=binomial())
    odds <- predict(mod_ps, newdata = dat1, type="response") |>
      exp()

    dat1 <- dat1 |>
      mutate(
        weight = case_when(
          ((trt==1) & (C=="1")) ~ 1,
          ((trt==0) & (C=="0")) ~ 0,
          ((trt==1) & (C=="0")) ~ 0,
          ((trt==0) & (C=="1")) ~ odds
        )
      )

    # suppressWarnings is used to ignore warnings about non-integer outcomes
    # this is expected due to using non-integer weights
    outcome_mod <- suppressWarnings({
      glm(evt ~ trt + V + W, weights=weight, family=binomial(), data = dat1)
    })

    emm <- emmeans(outcome_mod, ~ trt)
    # results on the log-odds scale (odds-ratio)
    ci_or <- emm |>
      pairs(type="response") |>
      confint(level=ci_level)
    # results on the log scale (risk-ratio)
    lemm <- regrid(emm, "log")
    ci_rr <- lemm |>
      contrast(method="pairwise", type="response") |>
      summary(null=log(1-VE_margin), side="<", infer=TRUE, level=ci_level)

    list(
      p = ci_rr$p.value,
      VE = 1-ci_rr$ratio,
      VE_lower = 1-ci_rr$asymp.UCL,
      VE_upper = 1-ci_rr$asymp.LCL,
      OR = ci_or$odds.ratio,
      OR_lower = ci_or$asymp.LCL,
      OR_upper = ci_or$asymp.UCL
    )
  }
}
