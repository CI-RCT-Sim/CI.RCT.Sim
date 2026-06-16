#' Analyse dataset from vaccine scenario, per-protocol analysis
#'
#' Estimates and tests are based on an exact conditional test. As described in
#' Nauta et. al.
#'
#' @param ci_level the confidence level for the CIs (defaults to 0.95)
#' @param VE_margin vaccine efficacy margin for the super-superiority test
#'
#' @returns an analyse function that returns a list with the elements
#'  * `VE` the point estimate for the vaccine efficacy
#'  * `VE_lower` the lower CI limit for vaccine efficacy
#'  * `VE_upper` the upper CI limit for vaccine efficacy
#'  * `p` the p-value of the exact super-superiority test
#'  * `N_pat` number of patients in the dataset
#'  * `N_evt` number of observed events in the dataset
#'
#' @export
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
#' my_analyse <- analyse_vaccine_pp2(ci_level=0.95)
#' my_analyse(Design[1, ], dat)
analyse_vaccine_pp2 <- function(ci_level=0.95, VE_margin=0.3){
  # cases0, cases1: cases control, trt
  # t0, t1: participants/person-time control, trt
  ci_ve <- function(cases0, cases1, t0, t1, conf.level=0.95){
    test <- binom.test(x=cases1, n=(cases0+cases1), conf.level = conf.level)
    ci <- unname(c(test$conf.int[2], test$estimate, test$conf.int[1]))
    r <- t1/t0

    VE <- 1-(ci/(r*(1-ci)))
    VE
  }

  # cases0, cases1: cases control, trt
  # t0, t1: participants/person-time control, trt
  test_ve <- function(cases0, cases1, t0, t1, VE_margin){
    r <- t1/t0
    # inversion of the formula to convert p to VE
    p_margin <- ((1-VE_margin)*r) / ((1-VE_margin)*r + 1)
    test <- binom.test(x=cases1, n=(cases0+cases1), p = p_margin, alternative="less")
    test$p.value
  }

  function(condition, dat, fixed_objects = NULL){

    # filter compliant participants
    dat1 <- dat |>
      within({
        trt <- factor(trt, levels = c("1", "0"))
      }) |>
      subset(C==1)

    cases <- tapply(dat$evt, dat$trt, sum)
    times <- tapply(dat$t,   dat$trt  , sum)
    ci   <-   ci_ve(cases["0"], cases["1"], times["0"], times["1"], conf.level = ci_level)
    test <- test_ve(cases["0"], cases["1"], times["0"], times["1"], VE_margin = VE_margin)


    # lower and upper exchanged because VE = 1-RR
    list(
      VE = ci[2],
      VE_lower = ci[1],
      VE_upper = ci[3],
      p = test,
      N_pat = nrow(dat1),
      N_evt = sum(dat1$evt)
    )
  }
}
