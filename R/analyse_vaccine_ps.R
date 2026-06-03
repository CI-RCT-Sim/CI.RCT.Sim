#' Analyse dataset from vaccine scenario using principal score weighting
#'
#' @param ci_level the confidence level for the CIs (defaults to 0.95)
#' @param VE_margin vaccine efficacy margin for the super-superiority test
#' @param covariates_in_outcomes_model should the covariates to estimate the principal score also be included in the outcomes model
#' @param V_unobserved consider covariate V unobserved (don't use it in analysis)
#' @param W_unobserved consider covariate W unobserved (don't use it in analysis)
#' @param W_interaction add treatment x W interaction in outcomes model
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
#'   vaccine_scenario_set_beta_A1_relative() |>
#'   vaccine_scenario_set_gamma_0() |>
#'   vaccine_scenario_set_true_eff() |>
#'   vaccine_scenario_set_samplesize()
#'
#' dat <- generate_vaccine(Design[3,])
#' my_analyse <- analyse_vaccine_ps(ci_level=0.95)
#' my_analyse(Design[3, ], dat)
analyse_vaccine_ps <- function(ci_level=0.95, VE_margin=0.3, covariates_in_outcomes_model=TRUE, V_unobserved=FALSE, W_unobserved=FALSE, W_interaction=FALSE){
  if(W_unobserved && W_interaction){
    stop("Cannot include treatement x W interaction if W is not observed.")
  }

  function(condition, dat, fixed_objects = NULL){

    formula_ps <- C ~ 1
    if(!V_unobserved){
      formula_ps <- update.formula(formula_ps, .~.+V)
    }

    if(!W_unobserved){
      formula_ps <- update.formula(formula_ps, .~.+W)
    }

    dat1 <- dat |>
      within({
        C <- factor(C, levels=c("0", "1"))
        trt <- factor(trt, levels=c("1", "0"))
      })
    mod_ps <- glm(formula_ps, subset = (trt==1), data=dat1, family=binomial())
    odds <- predict(mod_ps, newdata = dat1, type="link") |>
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

    formula_outcome <- evt ~ trt
    if(covariates_in_outcomes_model & (!V_unobserved)){
      formula_outcome <- update.formula(formula_outcome, .~.+V)
    }

    if(covariates_in_outcomes_model & (!W_unobserved)){
      if(W_interaction){
        formula_outcome <- update.formula(formula_outcome, .~.+W*trt)
      } else {
        formula_outcome <- update.formula(formula_outcome, .~.+W)
      }
    }

    outcome_mod <- suppressWarnings({
        glm(formula_outcome, weights=weight, family=binomial(), data = dat1)
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
      summary(infer=TRUE, level=ci_level)

    test_rr <- lemm |>
      contrast(method="pairwise", type="response") |>
      summary(null=log(1-VE_margin), side="<")

    emm_sandwich <- emmeans(outcome_mod, ~ trt, vcov. = sandwich::vcovHAC)
    # results on the log-odds scale (odds-ratio)
    ci_or_sandwich <- emm_sandwich |>
      pairs(type="response") |>
      confint(level=ci_level)
    # results on the log scale (risk-ratio)
    lemm_sandwich <- regrid(emm_sandwich, "log")

    ci_rr_sandwich <- lemm_sandwich |>
      contrast(method="pairwise", type="response") |>
      summary(infer=TRUE, level=ci_level)

    test_rr_sandwich <- lemm_sandwich |>
      contrast(method="pairwise", type="response") |>
      summary(null=log(1-VE_margin), side="<", vcov. = sandwich::vcovHAC)

    list(
      p_model = test_rr$p.value,
      VE = 1-test_rr$ratio,
      VE_lower = 1-ci_rr$asymp.UCL,
      VE_upper = 1-ci_rr$asymp.LCL,
      OR = ci_or$odds.ratio,
      OR_lower = ci_or$asymp.LCL,
      OR_upper = ci_or$asymp.UCL,
      OR_sandwich = ci_or_sandwich$odds.ratio,
      OR_lower_sandwich = ci_or_sandwich$asymp.LCL,
      OR_upper_sandwich = ci_or_sandwich$asymp.UCL,
      p = test_rr_sandwich$p.value,
      VE_sandwich = 1-test_rr_sandwich$ratio,
      VE_lower_sandwich = 1-ci_rr_sandwich$asymp.UCL,
      VE_upper_sandwich = 1-ci_rr_sandwich$asymp.LCL
    )
  }
}
