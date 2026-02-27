#' Analyse dataset from vaccine scenario, per-protocol analysis
#'
#' @param ci_level the confidence level for the CIs (defaults to 0.95)
#'
#' @returns an analyse function that returns a list with the elements
#'  * `p` the p-value of the exact binomial test
#'  * `coef` the point estimate
#'  * `ci_lower` the lower CI limit
#'  * `ci_upper` the upper CI limit
#' @export
#'
#' @examples
#' Design <- vaccine_scenario() |>
#'   vaccine_scenario_set_gamma_0() |>
#'   vaccine_scenario_set_true_eff() |>
#'   vaccine_scenario_set_samplesize()
#' Design
#'
#' Design$beta_A1 <- 0.1
#'
#' dat <- generate_vaccine(Design[1,])
#' my_analyse <- analyse_vaccine_pp(ci_level=0.95)
#' my_analyse(Design[1, ], dat)
analyse_vaccine_pp <- function(ci_level=0.95){
  function(condition, dat, fixed_objects = NULL){

    list(
      p = NA_real_,
      coef = NA_real_,
      ci_lower = NA_real_,
      ci_upper = NA_real_
    )
  }
}
