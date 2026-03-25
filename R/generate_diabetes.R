#' Generate Dataset that simulates a diabetes trial in which rescue medication
#' can be introduced
#'
#' @param condition condition row of Design dataset
#' @param fixed_objects fixed objects of Design dataset
#'
#' @details
#' Condition has to contain the following columns:
#'
#'   * k: the number of visits post baseline (visit 0 is baseline, visit 1 to k are post-baseline visits)
#'   * mean_age: mean age of the patients in the trial
#'   * sd_age: standard deviation of age in the trial
#'   * b_age: coefficient for the age effect on the outcome
#'   * mean_bl: mean HbA1c value at baseline
#'   * sd_bl: standard deviation of HbA1c at baseline
#'   * rho: correlation between repeated HbA1c measurements
#'   * delta: maximal treatment effect
#'   * lambda: rate of increasing treatment effect
#'   * delta_resc: maximal effect of rescue medication
#'   * lambda_resc: rate of increasing effect of rescue medication
#'   * resc_0: intercept for the probability of receiving rescue medication
#'   * resc_y: coefficient for the effect of the current HbA1c value on the probability of receiving rescue medication
#'   * resc_age: coefficient for the effect of age on the probability of receiving rescue medication
#'   * setup: determines whether rescue medication is switched to (setup = 0) or put on top of active treatment (setup = 1)
#'   * miss: a list where each element is a vector of length 4 containing the coefficients for the probability of dropout (intercept, coefficient for the effect of the current HbA1c value, coefficient for the effect of age, coefficient for the effect of rescue medication)
#'
#' The data generation process is as follows:
#' 1. For each patient, the age is generated from a normal distribution with mean `mean_age` and standard deviation `sd_age`.
#' 2. The baseline HbA1c value (y0) is generated from a normal distribution with mean `mean_bl` and standard deviation `sd_bl`.
#' 3. The treatment effect at each visit is calculated based on the parameters `delta and `lambda`, and the effect of rescue medication is calculated based on `delta_resc` and `lambda_resc`.
#' 4. The repeated measurements of the outcome (y1, ..., yk) are generated based on the baseline value, the age effect, the treatment effect (if the patient is in the treatment group), and a residual error term that follows a multivariate normal distribution with mean 0 and a covariance matrix determined by `sd_bl` and `rho`.
#' 5. The probability of receiving rescue medication at each visit is calculated based on the current HbA1c value, age, and the parameters `resc_0`, `resc_y`, and `resc_age`. Rescue medication is then assigned based on this probability, and if a patient receives rescue medication, their subsequent outcome measurements are adjusted based on the effect of rescue medication.
#' 6. The visit at which rescue medication is started is recorded in the `rescue_start` column. If a patient never receives rescue medication, `rescue_start` is set to k + 2.
#' 7. The probability of dropout at each visit is calculated based on the current HbA1c value, age, rescue medication status, and the parameters in the `miss` list. Dropout is then assigned based on this probability, and if a patient drops out, their subsequent outcome measurements and rescue medication status are set to NA.
#'
#'
#' @return
#' A data set with n rows and the columns id, trt (1 = treatment, 0 = control), age, y0, y1, ..., yk (the repeated measurements of the outcome), R0, R1, ..., Rk (indicators for rescue medication at each visit), and rescue_start (the visit at which rescue medication was started).
#'
#' @importFrom stats rbinom rnorm runif qnorm qlogis plogis
#'
#' @export
#' @describeIn generate_diabetes simulates a data set with n rows.
#'
#' @examples
#' # first create a design data.frame with diabetes_scenario()
#' # and calculate true values with
#' # diabetes_scenario_set_truevalues()
#' Design <- diabetes_scenario() |>
#'   diabetes_scenario_set_truevalues()
#' # second call generate_diabetes() with a condition row
#' # of the design data.frame to create the data set
#' generate_diabetes(Design[1, ])
generate_diabetes <- function(condition, fixed_objects = NULL) {
  # sequence with the visits
  visit <- 0:condition$k

  # calculation of visit specific effects and the true effect
  eff <- condition$delta * (1 - exp(-condition$lambda * visit))
  rescue_effect <- condition$delta_resc * (1 - exp(-condition$lambda_resc * (visit)))
  delta_true <- condition$delta / 2 * (1 - exp(-condition$lambda * condition$k))

  # specifying parameters for sample size calculation
  alpha <- 0.05
  power <- 0.8
  rho <- condition$rho
  sigma <- condition$sd_bl
  sigma_adj <- sqrt(sigma^2 * (1 - rho^2))

  ifelse(is.null(condition$nfix),
    n <- 2 * round(((qnorm(1 - alpha / 2) + qnorm(power))^2) * sigma_adj^2 * 2 / (delta_true^2)),
    n <- condition$nfix
  )
  if (condition$hyp == 0) {
    eff <- rep(0, condition$k + 1)
  }

  # data generation
  id <- 1:n
  trt <- rbinom(n, 1, 0.5)
  age <- rnorm(n, mean = condition$mean_age, sd = condition$sd_age)
  age_slope <- 2 * exp(-condition$b_age * (age - 30))
  response_trt <- runif(n)
  # Implement the correlation structure of the repeated measurements using mvtnorm
  mu_resid <- rep(0, length(visit))
  sigma_resid <- diag(x = sigma^2, length(visit))
  sigma_resid[upper.tri(sigma_resid)] <- rho * sigma^2
  sigma_resid[lower.tri(sigma_resid)] <- rho * sigma^2

  resid <- mvtnorm::rmvnorm(n, mu_resid, sigma_resid)
  if (condition$setup == 0) {
    mu <- matrix(NA, nrow = n, ncol = length(visit))
    Y <- matrix(NA, nrow = n, ncol = length(visit))
    for (i in 1:length(visit)) {
      mu[, i] <- condition$mean_bl +
        visit[i] / condition$k * age_slope
      Y[, i] <- mu[, i] + eff[i] * response_trt * trt + resid[, i]
    }
  } else {
    mu <- matrix(NA, nrow = n, ncol = length(visit))
    for (i in 1:length(visit)) {
      mu[, i] <- condition$mean_bl +
        visit[i] / condition$k * age_slope +
        eff[i] * response_trt * trt
    }
    Y <- mu + resid
  }

  # Implement rescue medication and its effect
  p_rescue <- plogis(condition$resc_0 + (Y - 10) * condition$resc_y + (age - condition$mean_age) * condition$resc_age)
  p_rescue[, 1] <- 0
  p_rescue[, condition$k + 1] <- 0

  resc <- matrix(rbinom((condition$k + 1) * n, size = 1, prob = p_rescue), nrow = n)
  rescue <- t(apply(resc, 1, cumsum)) > 0
  rescue_start <- rowSums(!rescue)
  k_rescue <- rowSums(rescue)
  response_rescue <- runif(n)
  any_rescue <- c()

  # Based on the VISIT at which rescue is given, the rescue medication is
  # assumed to have an effect on all subsequent visits, implying a shift of one
  # visit. But, since the Y matrix contains visit 0 to visit k, we have to shift
  # once more. Starting rescue medication at visit 8, for example, would lead to
  # an effect on visits 9 to k, which corresponds to the indices 10 to k+1 in
  # the Y matrix.
  for (i in 1:n) {
    if (k_rescue[i] > 0) {
      rescue_set <- (rescue_start[i] + 2):(condition$k + 1) # equivalent would be rescue_set <- (rescue_start[i] + 1):condition$k + 1
      Y[i, rescue_set] <- mu[i, rescue_set] +
        response_rescue[i] * rescue_effect[rescue_set - rescue_start[i] + 1] +
        resid[i, rescue_set]
      any_rescue[i] <- TRUE
    } else {
      rescue_start[i] <- condition$k + 2
      any_rescue[i] <- FALSE
    }
  }

  # Implement dropout
  p_miss <- plogis(condition$miss[[1]][1] + (Y - 10) * condition$miss[[1]][2] +
    (age - condition$mean_age) * condition$miss[[1]][3] +
    rescue * condition$miss[[1]][4]) # actual prob. to drop out
  p_miss[, 1] <- 0 # we assume complete data at baseline

  wd <- matrix(rbinom((condition$k + 1) * n, size = 1, prob = p_miss), nrow = n)
  wd1 <- t(apply(wd, 1, cumsum)) > 0

  for (i in 1:n) {
    miss_start <- sum(!wd1[i, ]) + 1
    if (miss_start <= (condition$k + 1)) {
      Y[i, miss_start:(condition$k + 1)] <- NA
      rescue[i, miss_start:(condition$k + 1)] <- NA
      if ((miss_start - 1) <= rescue_start[i]) {
        rescue_start[i] <- NA
      }
    }
  }
  m_start <- rowSums(!wd1)
  out <- data.frame(id, trt, age, Y, rescue_start, rescue[, 1:(condition$k + 1)] * 1, m_start)
  names(out) <- c("id", "trt", "age", paste("y", visit, sep = ""), "rescue_start", paste("R", visit[1:(condition$k + 1)], sep = ""), "m_start")
  out
}

#' Create an empty assumptions data.frame for generate_diabetes
#'
#' @param print print code to generate parameter set?
#'
#' @return For diabetes_scenario: a design tibble with default values invisibly
#'
#' @details diabetes_scenario generates a default design `data.frame`
#'   for use with generate_diabetes If print is `TRUE` code to produce
#'   the template is also printed for copying, pasting and editing by the user.
#'   (This is the default when run in an interactive session.)
#'
#' @export
#' @describeIn generate_diabetes generate default design tibble
#'
diabetes_scenario <- function(print = interactive()) {
  skel <- "params_scenarios_grid(
  k           = 12,                   # Number of visits post baseline
  mean_age    = 60,                       # mean of the variable age
  sd_age      = 10,                       # standard deviation of the variable age
  b_age       = log(2)/10,                # age coefficient
  mean_bl     = 8,                        # mean HbA1c value at baseline
  sd_bl       = 1,                        # standard deviation of HbA1c at baseline
  rho         = c(0.5,0),                 # Correlation between repeated HbA1c measurements
  delta       = -c(1,0.5),                # Maximal treatment effect
  lambda      = log(2)/2,                 # Rate of increasing treatment effect
  delta_resc  = -0.75,                    # Maximal effect of rescue medication
  lambda_resc = log(2),                        # Rate of increasing effect of rescue medication
  resc_0      = qlogis(c(0.05,0.02)),     # probability for rescue medication
  resc_y      = log(c(3,150)),            # strong effect due to high HbA1c
  resc_age    = -log(1.01),               # weaker age effect than for dropout
  setup       = c(0,1), # determines whether rescue medication is switched to (setup = 0) or put on top of active treatment (setup = 1)
  miss        = list(
  c(qlogis(0.02), log(3),log(1.02),log(1.5)), # probability for missing data in the core scenario
  c(-100000,0,0,0),                                # probability for missing data in the scenario with no dropout
  c(qlogis(0.04),log(150),log(1.02),log(1.5))# probability for missing data in the scenario with stronger dropout
  )) |>
  merge(data.frame(hyp=c(1,0)), by=NULL)
"


  if (print) {
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
#' @param Design Design data.frame for diabetes scenarios, e.g. created with diabetes_scenario()
#'
#' @return For diabetes_scenario_set_truevalues: a design tibble with true values for the treatment effect at final visit and the sample size needed to achieve 80% power for a two-sided test at alpha = 0.05
#'
#' @export
#'
#' @describeIn generate_diabetes calculate true treatment effect and sample size
#'
diabetes_scenario_set_truevalues <- function(Design) {
  Design$eff_true <- Design$delta / 2 * (1 - exp(-Design$lambda * Design$k))

  # specifying parameters for sample size calculation
  alpha <- 0.05
  power <- 0.8

  if (!"nfix" %in% colnames(Design)) {
    Design$n <- 2 * round(((qnorm(1 - alpha / 2) + qnorm(power))^2) *
      Design$sd_bl^2 * (1 - Design$rho^2) * 2 /
      (Design$eff_true^2))
  } else {
    Design$n <- Design$nfix
  }

  Design$eff_true <- ifelse(Design$hyp == 1, Design$eff_true, 0)
  Design
}
