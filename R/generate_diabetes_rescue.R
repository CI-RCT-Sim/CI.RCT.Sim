#' Generate Dataset that simulates a diabetes trial in which rescue medication
#' can be introduced
#'
#' @param condition condition row of Design dataset
#' @param fixed_objects fixed objects of Design dataset
#'
#' @details
#' Condition has to contain the following columns:
#'
#'   * k
#'   * mean_age
#'   * sd_age
#'   * b_age
#'   * mean_bl
#'   * sd_bl
#'   * rho
#'   * delta
#'   * lambda
#'   * delta_resc
#'   * lambda_resc
#'   * resc_0
#'   * resc_y
#'   * resc_age
#'   * miss_0
#'   * miss_y
#'   * miss_age
#'   * miss_resc
#'
#' @return
#' For generate_diabetes_rescue: A data set with n rows and the columns id, trt
#' (1=treatment, 0=control), age, y0, y1, ..., yk
#' (the repeated measurements of the outcome), R1, ..., Rk-1 (indicators for
#' rescue medication at each visit except baseline and the last visit),
#' and rescue_start (the visit at which rescue medication was started,
#' NA if no rescue medication was received).
#'
#' @importFrom stats rbinom rnorm runif qnorm qlogis plogis
#'
#' @export
#' @describeIn generate_diabetes_rescue simulates a data set with n rows.
#'
#' @examples
#' Design <- assumptions_diabetes_rescue()
#' generate_diabetes_rescue(Design[1, ])
generate_diabetes_rescue <- function(condition, fixed_objects = NULL) {
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
    for (i in 1:length(visit)) {
      mu[, i] <- condition$mean_bl +
        visit[i] / condition$k * age_slope
    }
    Y <- mu + eff[i] * response_trt * trt + resid
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

  for (i in 1:n) {
    if (k_rescue[i] > 0) {
      rescue_set <- (rescue_start[i] + 2):(condition$k + 1)
      Y[i, rescue_set] <- mu[i, rescue_set] +
        response_rescue[i] * rescue_effect[rescue_set - rescue_start[i] + 1] +
        resid[rescue_set]
      any_rescue[i] <- TRUE
    } else {
      rescue_start[i] <- NA
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
    }
  }

  out <- data.frame(id, trt, age, Y, rescue_start, rescue[, 1:condition$k + 1] * 1)
  names(out) <- c("id", "trt", "age", paste("y", visit, sep = ""), "rescue_start", paste("R", visit[1:condition$k + 1], sep = ""))
  out
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
#' @describeIn assumptions_diabetes_rescue generate default design tibble
#'
#' @examples
#' Design <- assumptions_diabetes_rescue()[1, ]
#' Design
assumptions_diabetes_rescue <- function(print = interactive()) {
  skel <- "params_scenarios_grid(
  k           = 12,                   # Number of visits post baseline
  mean_age    = 60,                       # mean of the variable age
  sd_age      = 10,                       # standard deviation of the variable age
  b_age       = log(2)/10,                # age coefficient
  mean_bl     = 8,                        # mean hbalc value at baseline
  sd_bl       = 1,                        # standard deviation of hba1c at baseline
  rho         = c(0.5,0),                 # Correlation between repeated HbA1c measurements
  delta       = -c(1,0.5),                # Maximal treatment effect
  lambda      = log(2)/2,                 # Rate of increasing treatment effect
  delta_resc  = -0.75,                    # Maximal effect of rescue medication
  lambda_resc = log(2),                        # Rate of increasing effect of rescue medication
  resc_0      = qlogis(c(0.05,0.02)),     # probability for rescue medication
  resc_y      = log(c(3,150)),            # strong effect due to high hba1c
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
#' true_summary_statistics_diabetes_rescue(assumptions_diabetes_rescue())
true_summary_statistics_diabetes_rescue <- function(Design, cutoff_stats = 10, fixed_objects = NULL) {
  # true_summary_statistics_diabetes_rescue_rowwise <- function(condition, cutoff_stats) {
  #   res <- data.frame(
  #     eff_true <- condition$delta / 2 * (1 - exp(-condition$lambda * condition$k))
  #   )
  #   res
  # }
  #
  # Design <- Design |>
  #   split(1:nrow(Design)) |>
  #   mapply(FUN = true_summary_statistics_diabetes_rescue_rowwise, cutoff_stats = cutoff_stats, SIMPLIFY = FALSE)
  #
  # Design <- do.call(rbind, Design)
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
