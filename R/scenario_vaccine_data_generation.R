#' Create an empty assumtions data.frame for generate_vaccine
#'
#' @param print print code to generate parameter set?
#'
#' @return For vaccine_scenario: a design tibble with default values invisibly
#'
#' @details vaccine_scenario generates a default design `data.frame`
#'   for use with `generate_vaccine`. If print is `TRUE` code to produce
#'   the template is also printed for copying, pasting and editing by the user.
#'   (This is the default when run in an interactive session.)
#'
#' @export
#' @describeIn vaccine_scenario generate default assumptions `data.frame`
#'
#' @examples
#' Design <- vaccine_scenario() |>
#'   vaccine_callibrate_gamma0() |>
#'   vaccine_callibrate_beta_A2()
#' Design
vaccine_scenario <- function(print=interactive()){
  skel <- "expand.grid(
  p_V = c(0, 0.1, 0.3), # probability for binary covariate prognostic for ICE and infection risk
  p_W = c(0, 0.1, 0.3), # probability for binary covariate prognostic for ICE and infection risk and modifying treatment effect
  lambda_post = -log(1-(1/c(500, 1000, 2000)))/365.25, # force of infection (baseline infection hazard) after 14 days, yearly incidence of 1/500, 1/1000, 1/2000
  # gamma_0  = c(), # callibrated from other parameters to give typical compliance
  gamma_W  = c(-0.8, 0), # regression parameters for compliance
  gamma_V  = c(0.5, 0),
  gamma_A  = c(-0.357, 0),
  gamma_AW = c(-0.3, 0),
  beta_V  = c(log(1.5), 0), # regression parameters for time to infection
  beta_W  = c(log(1.2), 0),
  # beta_A1  = c(), # choosen according to beta_A2
  beta_A2  = log(1-c(0, 0.5, 0.8, 0.9)), #
  beta_AW = c(log(0.8), 0),
  n_trt     = c(NA_real_), # study design parameters
  n_ctrl    = c(NA_real_),
  follow_up = c(NA_real_)
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

#' @param design previously generated desing data.frame
#' @param dose_2_completion target percentage of completion of dose 2
#'
#' @returns design with added column gamma_0
#' @export
#' @describeIn vaccine_scenario callibrate gamma_0
vaccine_callibrate_gamma0 <- function(design, dose_2_completion = 0.95){
  # only calculate for unique values of used variables
  tmp <- merge(
    design[, c("p_W", "p_V", "gamma_W", "gamma_V", "gamma_A", "gamma_AW")] |>
      unique(),
    data.frame(tmp_dose_2_completion=dose_2_completion),
    by=NULL, all=TRUE
  )
  tmp <- within(tmp, {
    gamma_0 = binomial()$linkfun(tmp_dose_2_completion) - gamma_W*p_W - gamma_V*p_V - gamma_A - gamma_AW
  })
  tmp$tmp_dose_2_completion <- NULL
  # merge calculated values
  design <- merge(design, tmp, by=c("p_W", "p_V", "gamma_W", "gamma_V", "gamma_A", "gamma_AW"), all.x=TRUE)
  design
}

#' @param design previously generated desing data.frame
#' @param VE target VE of never compliers
#'
#' @returns design with added column beta_A1
#' @export
#' @describeIn vaccine_scenario callibrate beta_A2
vaccine_callibrate_beta_A2 <- function(design, VE=c(0, 0.3, 0.5, 0.6)){
  # only calculate for unique values of used variables
  tmp <- expand.grid(beta_A2=unique(design$beta_A2), tmp_hr=1-VE)
  tmp$beta_A1 = log(-exp(tmp$beta_A2)+exp(tmp$tmp_hr))
  tmp$beta_A1[tmp$beta_A1==-Inf] <- 0
  tmp$tmp_hr <- NULL
  # merge calculated values
  design <- merge(design, tmp, by="beta_A2", all.x=TRUE)
  design
}
