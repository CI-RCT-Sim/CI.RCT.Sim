#' Analyse dataset from vaccine scenario, per-protocol analysis
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
#' @importFrom stats glm poisson
#' @importFrom multcomp glht
#'
#' @examples
#' Design <- vaccine_scenario() |>
#'   vaccine_scenario_set_gamma_0() |>
#'   vaccine_scenario_set_true_eff() |>
#'   vaccine_scenario_set_samplesize()
#' Design
#'
#' dat <- generate_vaccine(Design[1,])
#' my_analyse <- analyse_vaccine_pp(ci_level=0.95)
#' my_analyse(Design[1, ], dat)
analyse_vaccine_pp <- function(ci_level=0.95, VE_margin=0.3){
  function(condition, dat, fixed_objects = NULL){

    # glm drops (multi-)colinear values but glht throws an error, when columns
    # are missing, therefore construct formula and contrasts from only the
    # columsn with more than one value
    my_formula <- evt ~ trt
    my_contrast <- c(0,1)

    if(length(unique(dat$V)) > 1){
      my_formula <- update(my_formula, ~ . + V)
      my_contrast <- c(my_contrast, 0)
    }

    if(length(unique(dat$W)) > 1){
      my_formula <- update(my_formula, ~ . + W)
      my_contrast <- c(my_contrast, 0)
    }

    # filter compliant participants and calculate risk-ratio
    mod_ve <- dat |>
      subset(C==1) |>
      glm(my_formula, data = _, family=poisson(link="log"))

    p_val <- mod_ve |>
      glht(linfct = matrix(my_contrast, nrow=1), rhs=log(1-VE_margin), alternative="less") |>
      summary() |>
      _$test |>
      _$pvalues

    attributes(p_val) <- NULL

    ci <- mod_ve |>
      glht(linfct = matrix(my_contrast, nrow=1)) |>
      confint(level=ci_level) |>
      _$confint |>
      exp() |>
      (\(x) 1-x)()

    # lower and upper exchanged because VE = 1-RR
    list(
      p = p_val,
      VE = ci[1,"Estimate"],
      VE_lower = ci[1,"upr"],
      VE_upper = ci[1,"lwr"]
    )
  }
}
