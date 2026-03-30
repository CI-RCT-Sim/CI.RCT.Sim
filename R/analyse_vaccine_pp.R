#' Analyse dataset from vaccine scenario, per-protocol analysis
#'
#' Estimates and tests are based on marginal estimates from a Poisson regression
#' adjusted for `V` and `W`.
#'
#' @param ci_level the confidence level for the CIs (defaults to 0.95)
#' @param VE_margin vaccine efficacy margin for the super-superiority test
#'
#' @returns an analyse function that returns a list with the elements
#'  * `p` the p-value of the super-superiority test
#'  * `VE` the point estimate for the vaccine efficacy
#'  * `VE_lower` the lower CI limit
#'  * `VE_upper` the upper CI limit
#' @export
#'
#' @importFrom stats glm poisson
#' @importFrom emmeans emmeans
#' @importFrom graphics pairs
#'
#' @examples
#' Design <- vaccine_scenario() |>
#'   vaccine_scenario_set_beta_A1_relative() |>
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

    dat1 <- dat |>
      within({
        trt <- factor(trt, levels = c("1", "0"))
      })

    # filter compliant participants and calculate risk-ratio
    mod_ve <- dat1 |>
      subset(C==1) |>
      glm(evt ~ trt + V + W, data = _, family=poisson(link="log"))

    emm <- emmeans(mod_ve, ~ trt)
    # results on the log scale (risk-ratio)
    res <- emm |>
      pairs(type="response") |>
      summary(null=log(1-VE_margin), side="<", infer=TRUE, level=ci_level)

    # lower and upper exchanged because VE = 1-RR
    list(
      p = res$p.value,
      VE = 1-res$ratio,
      VE_lower = 1-res$asymp.UCL,
      VE_upper = 1-res$asymp.LCL
    )
  }
}
