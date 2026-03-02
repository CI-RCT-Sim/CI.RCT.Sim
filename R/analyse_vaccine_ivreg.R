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
#'
#' @examples
#' Design <- vaccine_scenario() |>
#'   vaccine_scenario_set_gamma_0() |>
#'   vaccine_scenario_set_true_eff() |>
#'   vaccine_scenario_set_samplesize()
#' Design
#'
#' dat <- generate_vaccine(Design[1,])
#' my_analyse <- analyse_vaccine_ivreg(ci_level=0.95)
#' my_analyse(Design[1, ], dat)
analyse_vaccine_ivreg <- function(ci_level=0.95, VE_margin=0.3){
  function(condition, dat, fixed_objects = NULL){

    dat <- dat |>
      within({
        T <- ifelse(((trt==1) & (C==1)), 1L, 0L)
      })

    ivreg::ivreg(evt ~ T + V | trt, data=dat)

    list(
      p = p_val,
      VE = ci[1,"Estimate"],
      VE_lower = ci[1,"upr"],
      VE_upper = ci[1,"lwr"]
    )
  }
}
