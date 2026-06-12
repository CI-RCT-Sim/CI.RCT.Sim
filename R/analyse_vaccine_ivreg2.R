#' Analyse dataset from vaccine scenario using instrumental variable regression
#'
#' @param ci_level the confidence level for the CIs (defaults to 0.95)
#' @param VE_margin vaccine efficacy margin for the super-superiority test
#' @param V_unobserved consider covariate V unobserved (don't use it in analysis)
#' @param W_unobserved consider covariate W unobserved (don't use it in analysis)
#'
#' @returns an analyse function that returns a list with the elements
#'  * `p` the p-value of the super-superiority test
#'  * `VE` the point estimate for the vaccine efficacy
#'  * `VE_lower` the lower CI limit for vaccine efficacy
#'  * `VE_upper` the upper CI limit for vaccine efficacy
#'  * `RD` the point estimate for risk difference
#'  * `RD_lower` the lower CI limit for risk difference
#'  * `RD_upper` the upper CI limit for risk difference
#' @export
#'
#' @details
#' This implementation uses a different first stage than `analyse_vaccine_ivreg`.
#'
#'
#' @importFrom emmeans emmeans regrid contrast
#' @importFrom graphics pairs
#' @importFrom stats lm.fit model.matrix update
#'
#' @examples
#' Design <- vaccine_scenario()
#'
#' condition <- Design[3, ]
#' condition$p_W <- 0.3
#'
#' condition <- condition  |>
#'   vaccine_scenario_set_beta_A1_relative() |>
#'   vaccine_scenario_set_gamma_0() |>
#'   vaccine_scenario_set_true_eff() |>
#'   vaccine_scenario_set_samplesize()
#'
#' dat <- generate_vaccine(condition)
#' my_analyse  <- analyse_vaccine_ivreg (ci_level=0.95)
#' my_analyse2 <- analyse_vaccine_ivreg2(ci_level=0.95)
#'
#' res1 <- my_analyse(conditon, dat)
#' res2 <- my_analyse2(conditon, dat)
#'
#' res1
#' res2
#'
#' lapply(names(res1), \(i){res1[[i]]-res2[[i]]})
analyse_vaccine_ivreg2 <- function(ci_level=0.95, VE_margin=0.3, V_unobserved=FALSE, W_unobserved=FALSE){
  function(condition, dat, fixed_objects = NULL){

    dat1 <- dat |>
      within({
        T <- ifelse(((trt==1) & (C==1)), 0L, 1L)
        V <- factor(V)
        W <- factor(W)
        trt <- factor(trt)
        evt <- as.integer(evt)
      })

    formula_stage1 <- T ~ trt

    # lm drops (multi-)colinear values but but we need to know which columns are
    # included in the second stage of the estimation therefore we build our
    # formulas and variable names here
    formula_stage2 <- evt ~ stage1_T

    if((length(unique(dat$V)) > 1) & (!V_unobserved)){
      formula_stage1 <- update(formula_stage1, . ~ . + V)
      formula_stage2 <- update(formula_stage2, . ~ . + V)
    }

    if((length(unique(dat$W)) > 1) & (!W_unobserved)){
      formula_stage1 <- update(formula_stage1, . ~ . + W)
      formula_stage2 <- update(formula_stage2, . ~ . + W)
    }

    lm_stage1 <- lm(formula_stage1, dat1)
    dat1$stage1_T <- predict(lm_stage1)

    stage2 <- glm(formula_stage2, data=dat1, family=poisson(link="log"))

    emm <- emmeans(stage2, ~ stage1_T, at=list(stage1_T=c(0,1)))
    # results on the log scale (risk-ratio)
    res <- emm |>
      pairs(type="response") |>
      summary(null=log(1-VE_margin), side="<")

    res_ci <- emm |>
      pairs(type="response") |>
      summary(infer=TRUE, level=ci_level)

    emm_sandwich <- emmeans(stage2, ~ stage1_T, at=list(stage1_T=c(0,1)), vcov. = sandwich::vcovHAC)
    res_sandwich <- emm_sandwich |>
      pairs(type="response") |>
      summary(null=log(1-VE_margin), side="<")

    res_ci_sandwich <- emm_sandwich |>
      pairs(type="response") |>
      summary(infer=TRUE, level=ci_level)

    # lower and upper exchanged because VE = 1-RR
    list(
      p_model = res$p.value,
      VE = 1-res$ratio,
      VE_lower = 1-res_ci$asymp.UCL,
      VE_upper = 1-res_ci$asymp.LCL,
      p = res_sandwich$p.value,
      VE_sandwich = 1-res_sandwich$ratio,
      VE_lower_sandwich = 1-res_ci_sandwich$asymp.UCL,
      VE_upper_sandwich = 1-res_ci_sandwich$asymp.LCL,
      N_pat = nrow(dat1),
      N_evt = sum(dat1$evt)
    )
  }
}
