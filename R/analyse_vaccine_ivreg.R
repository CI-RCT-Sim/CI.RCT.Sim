#' Analyse dataset from vaccine scenario using instrumental variable regression
#'
#' @param ci_level the confidence level for the CIs (defaults to 0.95)
#' @param VE_margin vaccine efficacy margin for the super-superiority test
#'
#' @returns an analyse function that returns a list with the elements
#'  * `p` the p-value of the exact binomial test
#'  * `VE` the point estimate for the vaccine efficacy
#'  * `VE_lower` the lower CI limit
#'  * `VE_upper` the upper CI limit
#' @export
#'
#' @importFrom ivreg ivreg
#' @importFrom emmeans emmeans regrid
#'
#' @examples
#' Design <- vaccine_scenario() |>
#'   vaccine_scenario_set_gamma_0() |>
#'   vaccine_scenario_set_true_eff() |>
#'   vaccine_scenario_set_samplesize()
#'
#' dat <- generate_vaccine(Design[1,])
#' my_analyse <- analyse_vaccine_ivreg(ci_level=0.95)
#' my_analyse(Design[1, ], dat)
analyse_vaccine_ivreg <- function(ci_level=0.95, VE_margin=0.3){
  function(condition, dat, fixed_objects = NULL){

    dat1 <- dat |>
      within({
        T <- factor(ifelse(((trt==1) & (C==1)), 1L, 0L), levels=c("1", "0"))
        V <- factor(V)
        trt <- factor(trt)
        evt <- as.integer(evt)
      })

    # lm drops (multi-)colinear values but but we need to know which columns are
    # included in the second stage of the estimation therefore we build our
    # formulas and variable names here
    formula_stage1x <- ~ trt
    formula_stage1y <- ~ T
    vars_stage1 <- c("T0")
    vars_stage2 <- c("stage1_T")
    formula_stage2 <- evt ~ stage1_1 + stage1_T

    if(length(unique(dat$V)) > 1){
      formula_stage1x <- update(formula_stage1x, ~ . + V)
      formula_stage1y <- update(formula_stage1y, ~ . + V)
      formula_stage2  <- update(formula_stage2 , ~ . + stage1_V)
      vars_stage1 <- c(vars_stage1, "V1")
      vars_stage2 <- c(vars_stage2, "stage1_V")
    }

    if(length(unique(dat$W)) > 1){
      formula_stage1x <- update(formula_stage1x, ~ . + W)
      formula_stage1y <- update(formula_stage1y, ~ . + W)
      formula_stage2  <- update(formula_stage2 , ~ . + stage1_W)
      vars_stage1 <- c(vars_stage1, "V1")
      vars_stage2 <- c(vars_stage2, "stage1_W")
    }

    vars_stage1 <- c(vars_stage1, "(Intercept)")
    vars_stage2 <- c(vars_stage2, "stage1_1")
    formula_stage2  <- update(formula_stage2 , ~ . -1)

    x_stage1 <- model.matrix(formula_stage1x, data=dat1)
    y_stage1 <- model.matrix(formula_stage1y, data=dat1)

    lm_stage1 <- lm.fit(x=x_stage1, y=y_stage1)
    dat1[, vars_stage2] <- lm_stage1$fitted.values[, vars_stage1]

    stage2 <- lm(formula_stage2, data=dat1)

    # marginalize over other variables
    emm <- emmeans(stage2, ~ stage1_T, at=list(stage1_T=c(0,1)))
    lemm <- regrid(emm, "log")
    res <- pairs(lemm, type = "response")
    ci <- confint(res)

    list(
      p = summary(res)$p.value,
      VE = 1-ci$ratio,
      VE_lower = 1-ci$upper.CL,
      VE_upper = 1-ci$lower.CL
    )
  }
}
