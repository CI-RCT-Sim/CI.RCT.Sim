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
#'
#' generate_vaccine(Design[1,])
vaccine_scenario <- function(print=interactive()){
  skel <- "params_scenarios_grid(
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
  n_trt     = c(50), # study design parameters
  n_ctrl    = c(50),
  follow_up = c(365)
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
vaccine_callibrate_gamma0 <- function(design, dose_2_completion = 0.95){ # TODO: check calculations
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
vaccine_callibrate_beta_A2 <- function(design, VE=c(0, 0.3, 0.5, 0.6)){ # TODO: check calculations
  # only calculate for unique values of used variables
  tmp <- expand.grid(beta_A2=unique(design$beta_A2), tmp_hr=1-VE)
  tmp$beta_A1 = log(-exp(tmp$beta_A2)+exp(tmp$tmp_hr))
  tmp$beta_A1[tmp$beta_A1==-Inf] <- 0
  tmp$tmp_hr <- NULL
  # merge calculated values
  design <- merge(design, tmp, by="beta_A2", all.x=TRUE)
  design
}




#' @param condition Row of parameters dataset
#' @param fixed_objects other objects passed to simulation runs
#'
#' @describeIn vaccine_scenario simulate a dataset from the vaccine scenario
#' @importFrom stats binomial runif rbinom
#' @importFrom miniPCH rpch
#'
#' @returns a simulated dataset
#' @export
generate_vaccine <- function(condition, fixed_objects = NULL){
  N <- condition$n_ctrl + condition$n_trt
  V <- rbinom(N, 1, condition$p_V)
  W <- rbinom(N, 1, condition$p_W)
  A <- c(rep(1, condition$n_trt), rep(0, condition$n_ctrl))

  # regression model for compliance propensity
  p <- function(A_){
    binomial()$linkinv(
      condition$gamma_0  +
        condition$gamma_W * W +
        condition$gamma_V * V +
        condition$gamma_A * A_ +
        condition$gamma_AW * A_ * W
    )
  }

  # latent variable
  U <- runif(N)
  # counterfactual compliances
  C_1 <- as.numeric(U < p(1))
  C_0 <- as.numeric(U < p(0))
  # observed compliance
  C <- as.numeric(U < p(A))

  # TODO: check if beta_AW is the same for d=1,2 or if there should be one more parameter
  theta_1 <- condition$beta_A1 + condition$beta_AW * W
  theta_2 <- condition$beta_A2 + condition$beta_AW * W
  theta_early <- theta_1 * W
  theta_late  <- theta_1 * W + C * (theta_2*W - theta_1*W)


  t_ <- c(0, 14)
  lambda_0 <- diag(c(0, condition$lambda_post))

  # under observed treamtent allocation
  # time to event hazard
  lambda <- cbind(
    exp(condition$beta_V * V + condition$beta_W * W + A * theta_early),
    exp(condition$beta_V * V + condition$beta_W * W + A * theta_late)
  ) %*% lambda_0
  # time to event outcome
  T_ <- sapply(1:N, \(i){
    rpch(1, t=t_, lambda = lambda[i,])
  })
  # binary outcome
  Y <- T_ < condition$follow_up
  T_ <- pmin(T_, condition$follow_up)

  # under potential treatment allocation to control
  # time to event hazard
  lambda_a0 <- cbind(
    exp(condition$beta_V * V + condition$beta_W * W + 0 * theta_early),
    exp(condition$beta_V * V + condition$beta_W * W + 0 * theta_late)
  ) %*% lambda_0
  # time to event outcome
  T_a0 <- sapply(1:N, \(i){
    rpch(1, t=t_, lambda = lambda_a0[i,])
  })
  # binary outcome
  Y_a0 <- T_a0 < condition$follow_up
  T_a0 <- pmin(T_a0, condition$follow_up)

  # under potential treatment allocation to active
  # time to event hazard
  lambda_a1 <- cbind(
    exp(condition$beta_V * V + condition$beta_W * W + 1 * theta_early),
    exp(condition$beta_V * V + condition$beta_W * W + 1 * theta_late)
  ) %*% lambda_0
  # time to event outcome
  T_a1 <- sapply(1:N, \(i){
    rpch(1, t=t_, lambda = lambda_a1[i,])
  })
  # binary outcome
  Y_a1 <- T_a1 < condition$follow_up
  T_a1 <- pmin(T_a1, condition$follow_up)

  # TODO: check which counterfactuals are needed in the returned dataset
  data.frame(
    V=V,
    W=W,
    trt=A,
    C=C,
    t=T_,
    evt=Y
  )
}
