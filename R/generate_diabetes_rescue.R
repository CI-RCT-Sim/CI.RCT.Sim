#' Generate Dataset ...
#'
#' @param condition condition row of Design dataset
#' @param fixed_objects fixed objects of Design dataset
#'
#' @details
#' Condition has to contain the following columns:
#'
#'   * eff the effect size
#'   *
#'   * ...
#'
#' @return
#' For generate_diabetes_rescue: A data set with the columns id, trt
#' (1=treatment, 0=control), evt (event, currently TRUE for all observations)
#'
#' @export
#' @describeIn generate_diabetes_rescue simulates a dataset with ...
#'
#' @examples
generate_diabetes_rescue <- function(condition, fixed_objects=NULL){

  n <- 100 # Should be edited according to proposed power and effect size
  # n <- ((qnorm(1 - alpha / 2) + qnorm(power))^2)*sigma^2/(eff^2)

  visit <- 0:condition$k[1]
  # browser()
  id=1:n
  trt=rbinom(n, 1, 0.5)
  age = rnorm(n, mean = condition$mean_age, sd = condition$sd_age)
  age_slope = 2 * exp(-condition$b_age[1] * (age - 30))
  response_trt_i = runif(n)
  mu_ctr <- matrix(NA, nrow = n, ncol = length(visit))
  for (i in 1:length(visit)){
    mu_ctr[,i] <- condition$mean_bl[1] + visit[i] / condition$k[1] * age_slope + condition$eff[i] * response_trt_i * trt
  }

  resid <- mvtnorm::rmvnorm(n, rep(0, length(visit)), diag(length(visit)))

  Y0 = rnorm(n, mean = condition$mean_bl)
  data.frame(id, trt, age, age_slope, mu_ctr, Y0, resid)
}

#' Create an empty assumptions data.frame for generate_diabetes_rescue
#'
#' @param print print code to generate parameter set?
#'
#' @return For assumptions_diabetes_rescue: a design tibble with default values invisibly
#'
#' @details assumptions_diabetes_rescue generates a default design `data.frame`
#'   for use with generate_diabetes_rescue If print is `TRUE` code to produce
#'   the template is also printed for copying, pasting and editing by the user.
#'   (This is the default when run in an interactive session.)
#'
#' @export
#' @describeIn generate_diabetes_rescue generate default design tibble
#'
#' @examples
#' Design <- assumptions_diabetes_rescue()
#' Design
assumptions_diabetes_rescue <- function(print=interactive()){
  skel <- "expand.grid(
  eff = c(0,5,10),             # treatment effect
  k = 2,                       # Number of visits post baseline
  mean_bl=8,                   # mean hbalc value at baseline
  mean_age=60,                 # mean of the variable age
  sd_age=10,                   # standard deviation of the variable age
  b_age=log(2)/10,             # coefficient
  delay=m2d(seq(0, 10, by=2)), # delay of 0, 1, ..., 10 months
  hazard_ctrl=m2r(24),         # median survival control of 24 months
  hazard_trt=m2r(36),          # median survival treatment of 36 months
  random_withdrawal=m2r(120)   # median time to random withdrawal 10 years
)
"

if(print){
  cat(skel)
}

invisible(
  skel |>
    str2expression() |>
    eval()
)
}

#' Calculate true summary statistics for scenarios with delayed treatment effect
#'
#' @param Design Design data.frame for x
#' @param cutoff_stats Cutoff time for rmst and average hazard ratios
#' @param fixed_objects fixed objects not used for now
#'
#' @return For true_summary_statistics_x: the design data.frame
#'   passed as argument with the additional columns:
#' * `rmst_trt` rmst in the treatment group
#' * `median_surv_trt` median survival in the treatment group
#' * `rmst_ctrl` rmst in the control group
#' * `median_surv_ctrl` median survial in the control group
#' * `gAHR` geometric average hazard ratio
#' * `AHR` average hazard ratio
#'
#' @export
#'
#' @describeIn generate_diabetes_rescue  calculate true summary statistics for ...
#'
#' @examples
true_summary_statistics_diabetes_rescue <- function(Design, cutoff_stats=10, fixed_objects=NULL){

  true_summary_statistics_diabetes_rescue_rowwise <- function(condition, cutoff_stats){
    res <- data.frame(
      rmst_trt = NA_real_,
      medial_surv_trt = NA_real_,
      rmst_ctrl = NA_real_,
      median_surv_ctrl = NA_real_,
      gAHR = NA_real_,
      AHR = NA_real_
    )

    res
  }

  Design <- Design |>
    split(1:nrow(Design)) |>
    mapply(FUN=true_summary_statistics_x_rowwise, cutoff_stats = cutoff_stats, SIMPLIFY = FALSE)

  Design <- do.call(rbind, Design)

  Design
}
